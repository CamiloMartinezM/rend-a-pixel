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
            Intersection prevIts;
            Color L(0.0f), throughput(1.0f);
            DirectLightSample directLightSample;
            LightSample lightSample;
            BsdfEval bsdfVal;
            float p_bsdf = 1.0f;
            for (int depth = 0; depth < maxDepth; depth++)
            {
                Intersection its = m_scene->intersect(iterRay, rng);
                if (!its)
                {
                    // The ray misses and hits the background
                    Color Le = m_scene->evaluateBackground(iterRay.direction).value;

                    // TODO: Check if background light pdf is 1.0f
                    float misWeight = (mis && depth > 0) ? powerHeuristic(p_bsdf, 1.0f) : 1.0f;

                    // Incoporate the background contribution with MIS if depth > 0
                    L += Le * misWeight * throughput;
                    break;
                }

                // Evaluate direct emission from the hit point with the following rules:
                // 1. If NEE is not active or if it is the initial iteration (depth == 0), evaluate the emission from
                // the intersected surface regardless of what the intersected surface is;
                // 2. Regardless of if NEE is active or not, always account for the intersection emission, if it is not
                // an area light (that is, skip emission from light sources).
                // if ((its.instance->light() == nullptr) && (!nee || depth == 0))
                if (its.instance->light() == nullptr)
                {
                    L += its.evaluateEmission() * throughput;
                }
                else
                {
                    if (!nee || depth == 0)
                    {
                        L += its.evaluateEmission() * throughput;
                    }
                    else if (mis) // Compute MIS weight for area light if the intersected instance is an area light
                    {
                        // Compute the PDF given by the instance associated with the area light of having sampled the
                        // ray direction that hit the area light.
                        float p_light =
                            its.instance->light()->sampledDirectionPdf(iterRay.direction) * lightSample.probability;
                        float misWeight = powerHeuristic(p_bsdf, p_light);
                        L += its.evaluateEmission() * misWeight * throughput;
                    }
                }

                // Use Next-Event Estimation to sample a light if:
                // 1. NEE is active;
                // 2. The loop has not reached the final depth (we don't do NEE in the last iteration);
                // 3. The scene has lights.
                if (nee && iterRay.depth < maxDepth - 1 && m_scene->hasLights())
                {
                    lightSample = m_scene->sampleLight(rng);
                    directLightSample = lightSample.light->sampleDirect(its.position, rng, its);
                    Ray shadowRay(its.position, directLightSample.wi);
                    if (!m_scene->intersect(shadowRay, directLightSample.distance, rng))
                    {
                        // Evaluate the BSDF at the hit point for the light direction
                        bsdfVal = its.evaluateBsdf(directLightSample.wi);

                        // Modulate the light's contribution
                        Color lightContribution = bsdfVal.value * directLightSample.weight;

                        // Calculate the MIS weight using the light's PDF and the BSDF's PDF if MIS is active
                        float misWeight =
                            (mis) ? powerHeuristic(directLightSample.pdf * lightSample.probability, bsdfVal.pdf) : 1.0f;

                        // Final color contribution calculated by multiplying by the throughput, the MIS weight and
                        // dividing by light choosing probability
                        L += lightContribution * misWeight * throughput / lightSample.probability;
                    }
                }

                // Sample the BSDF to get the new direction and the weight
                BsdfSample bsdfSample = its.sampleBsdf(rng);
                if (bsdfSample.isInvalid())
                    break;

                // Save and update variables for next iteration
                p_bsdf = bsdfSample.pdf;
                prevIts = its;
                throughput *= bsdfSample.weight;

                // c) Trace a secondary ray in the direction determined by the BSDF sample.
                iterRay = Ray(its.position, bsdfSample.wi.normalized(), iterRay.depth + 1);
            }

            return L;
        }

        /**
         * @brief Compute the Multiple Importance Sampling weight using the power heuristic function (exponent = 2).
         * @param pdf_a Probability density function value of technique A.
         * @param pdf_b Probability density function value of technique B.
         * @return MIS weight for technique A.
         */
        float powerHeuristic(float pdf_a, float pdf_b)
        {
            float weight_a = sqr(pdf_a);
            float weight_b = sqr(pdf_b);
            return weight_a / (weight_a + weight_b);
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