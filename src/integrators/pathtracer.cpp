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

                // Evaluate direct emission from the hit point
                L_di += its.evaluateEmission() * throughput;

                // Compute the direct lighting using next-event estimation
                if (iterRay.depth < maxDepth - 1)
                {
                    if (m_scene->hasLights())
                    {
                        LightSample lightSample = m_scene->sampleLight(rng);
                        if (!lightSample.light->canBeIntersected())
                        {
                            DirectLightSample directLightSample = lightSample.light->sampleDirect(its.position, rng);
                            Ray shadowRay(its.position, directLightSample.wi);
                            if (!m_scene->intersect(shadowRay, directLightSample.distance, rng))
                            {
                                Color bsdfVal = its.evaluateBsdf(directLightSample.wi).value;
                                float lightSampleProb = lightSample.probability;
                                L_di += bsdfVal * directLightSample.weight / lightSampleProb * throughput;
                            }
                        }
                    }
                }
                else
                    continue;

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

        /**
         * @brief Compute the contribution of a camera-sampled ray using recursive path tracing.
         */
        Color recursivePathTracing(const Ray &ray, Sampler &rng, Color throughput)
        {
            if (ray.depth >= maxDepth)
            {
                // If maximum depth is reached, terminate the path
                return Color(0.0f);
            }

            Intersection its = m_scene->intersect(ray, rng);
            if (!its)
            {
                // If the ray misses the scene, return the background color
                return m_scene->evaluateBackground(ray.direction).value * throughput;
            }

            Color L_direct = Color(0.0f);

            // Evaluate direct emission from the hit point
            L_direct += its.evaluateEmission() * throughput;
            Color L_indirect = Color(0.0f);

            if (ray.depth < maxDepth - 1)
            {
                // Compute the direct lighting using next-event estimation
                if (m_scene->hasLights())
                {
                    LightSample lightSample = m_scene->sampleLight(rng);
                    if (!lightSample.light->canBeIntersected())
                    {
                        DirectLightSample directLightSample = lightSample.light->sampleDirect(its.position, rng);
                        Ray shadowRay(its.position, directLightSample.wi);
                        if (!m_scene->intersect(shadowRay, directLightSample.distance, rng))
                        {
                            Color bsdfVal = its.evaluateBsdf(directLightSample.wi).value;
                            L_direct += bsdfVal * directLightSample.weight / lightSample.probability * throughput;
                        }
                    }
                }
            }
            else
            {
                return L_direct;
            }

            // Sample the BSDF to get the new direction and the weight
            BsdfSample bsdfSample = its.sampleBsdf(rng);
            if (bsdfSample.isInvalid())
            {
                // If the BSDF sampling is invalid, return the direct emission
                return L_direct + L_indirect;
            }

            throughput *= bsdfSample.weight;

            // Trace a secondary ray in the direction determined by the BSDF sample
            Ray secondaryRay(its.position, bsdfSample.wi.normalized(), ray.depth + 1);

            L_indirect += recursivePathTracing(secondaryRay, rng, throughput);

            // Accumulate the direct and indirect lighting
            return L_direct + L_indirect;
        }

      public:
        int maxDepth;

        PathTracer(const Properties &properties) : SamplingIntegrator(properties)
        {
            maxDepth = properties.get<int>("depth", 2);
        }

        /**
         * @brief Compute the contribution of a camera-sampled ray using path tracing.
         */
        Color Li(const Ray &ray, Sampler &rng) override
        {
            if (UseRecursivePathTracer)
                return recursivePathTracing(ray, rng, Color(1.0f));
            else
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