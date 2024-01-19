#include <lightwave.hpp>

namespace lightwave
{
    class PathTracer : public SamplingIntegrator
    {
      private:
        /**
         * @brief Compute the contribution of a camera-sampled ray using iterative path tracing.
         */
        Color iterativePathTracing(const Ray &ray, Sampler &rng)
        {
            Ray iterRay = ray;
            Color L_di = Color(0.0f);
            Color throughput(1.0f);
            for (int depth = 0; depth < maxDepth; depth++)
            {
                Intersection its = m_scene->intersect(iterRay, rng);
                if (!its)
                {
                    // The ray misses and hits the background
                    L_di += m_scene->evaluateBackground(iterRay.direction).value * throughput;
                    break;
                }

                // Evaluate direct emission from the hit point with the following rules:
                // 1. If NEE is not active, evaluate the emission from the intersected surface regardless of what it is
                // 2. If NEE is active and the current intersected surface is not an area light, also evaluate the
                //    intersection emission. If it is an area light, skip emission from light sources
                if (its.instance->light() == nullptr)
                {
                    L_di += its.evaluateEmission() * throughput;
                }
                else
                {
                    if (its.instance->light() != nullptr)
                    {
                        if (!nee || depth == 0)
                        {
                            L_di += its.evaluateEmission() * throughput;
                        }
                    }
                }
                
                // Use Next-Event Estimation to sample a light if NEE is active
                if (nee)
                {
                    if (iterRay.depth < maxDepth - 1)
                    {
                        if (m_scene->hasLights())
                        {
                            LightSample lightSample = m_scene->sampleLight(rng);
                            DirectLightSample directLightSample =
                                lightSample.light->sampleDirect(its.position, rng, its);
                            Ray shadowRay(its.position, directLightSample.wi);
                            if (!m_scene->intersect(shadowRay, directLightSample.distance, rng))
                            {
                                // Evaluate the BSDF at the hit point for the light direction
                                Color bsdfVal = its.evaluateBsdf(directLightSample.wi).value;

                                // Modulate the light's contribution
                                Color lightContribution = bsdfVal * directLightSample.weight;

                                // Final color contribution multiplying by the throughput and dividing by the light
                                // probability
                                L_di += lightContribution / lightSample.probability * throughput;
                            }
                        }
                    }
                    else
                        continue;
                }

                // Sample the BSDF to get the new direction and the weight
                BsdfSample bsdfSample = its.sampleBsdf(rng);
                if (bsdfSample.isInvalid())
                    break;

                throughput *= bsdfSample.weight;

                // c) Trace a secondary ray in the direction determined by the BSDF sample.
                iterRay = Ray(its.position, bsdfSample.wi.normalized(), iterRay.depth + 1);
            }

            return L_di;
        }

      public:
        int maxDepth; // Maximum depth to explore with pathtracing
        bool nee;     // Use Next-Event Estimation
        bool mis;     // Use Multiple Importance Sampling

        PathTracer(const Properties &properties) : SamplingIntegrator(properties)
        {
            maxDepth = properties.get<int>("depth", 2);
            nee = properties.get<bool>("nee", false);
            mis = properties.get<bool>("mis", false);
        }

        /**
         * @brief Compute the contribution of a camera-sampled ray using path tracing.
         */
        Color Li(const Ray &ray, Sampler &rng) override
        {
            return iterativePathTracing(ray, rng);
        }

        /// @brief An optional textual representation of this class, useful for debugging.
        std::string toString() const override
        {
            return tfm::format("PathTracer[\n"
                               "  sampler = %s,\n"
                               "  image = %s,\n"
                               "]",
                               indent(m_sampler), indent(m_image));
        }
    };
} // namespace lightwave

// Register the PathTracer class for the "pathtracer" integrator type.
REGISTER_INTEGRATOR(PathTracer, "pathtracer")