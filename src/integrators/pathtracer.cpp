#include <lightwave.hpp>

namespace lightwave
{
    class PathTracer : public SamplingIntegrator
    {
      private:
        /// @brief Compute the contribution of a camera-sampled ray using iterative path tracing.
        Color iterativePathTracing(const Ray &ray, Sampler &rng)
        {
            Ray iterRay = ray;
            Color L(0.0f), throughput(1.0f);
            float p_bsdf = Infinity, p_light = 0.0f;
            float lightSelectionProb = m_scene->lightSelectionProbability(nullptr);
            for (int depth = 0; depth < maxDepth; depth++)
            {
                Intersection its = m_scene->intersect(iterRay, rng);
                if (!its)
                {
                    // The ray misses and hits the background:
                    // 1. Add background contribution for first bounce always.
                    // 2. If UseImprovedEnvSampling is not active, always account for background contribution.
                    if (!UseImprovedEnvSampling || depth == 0)
                    {
                        BackgroundLightEval bEval = m_scene->evaluateBackground(iterRay.direction);
                        L += bEval.value * throughput;
                    }
                    break;
                }

                // Evaluate direct emission from the hit point with the following rules:
                // 1. If it is the first intersection (depth == 0), evaluate it regardless;
                // 2. If NEE nor MIS are active (only do BSDF sampling), evaluate it regardless;
                // 3. If NEE is active and MIS is not, evaluate the emission on the surface if it is not an area
                // light;
                // 4. If MIS is active and it is not an area light, weigh the emission contrib. by the Power
                // Heuristic.
                if (depth == 0 || !nee || its.instance->light() == nullptr)
                    L += its.evaluateEmission() * throughput;
                else if (mis)
                {
                    // Do MIS weighing of the BSDF sample contribution based on its own PDF and the PDF of having
                    // sampled this direction with NEE (its.pdf)
                    p_light = pdfToSolidAngleMeasure(its.pdf, its.t, its.frame.normal, its.wo) * lightSelectionProb;
                    float misWeight = powerHeuristic(p_bsdf, p_light);
                    L += its.evaluateEmission() * misWeight * throughput;
                }

                if (iterRay.depth >= maxDepth - 1)
                    break;

                // Use Next-Event Estimation to sample a light if:
                // 1. NEE is active or MIS is active;
                // 2. The loop has not reached the final depth (we don't do NEE in the last iteration);
                // 3. The scene has lights.
                if (nee && m_scene->hasLights())
                {
                    LightSample lightSample = m_scene->sampleLight(rng);
                    if (mis || (nee && (!lightSample.light->canBeIntersected() ||
                                        lightSample.light->getLightType() == LightType::AreaLight)))
                    {
                        DirectLightSample directLightSample = lightSample.light->sampleDirect(its.position, rng, its);
                        Ray shadowRay(its.position, directLightSample.wi);

                        // Check for occlusion only if the sampled light is not an environment map
                        if (!m_scene->intersect(shadowRay, directLightSample.distance, rng))
                        {
                            // Evaluate the BSDF at the hit point for the light direction
                            BsdfEval bsdfVal = its.evaluateBsdf(directLightSample.wi);

                            // Modulate the light's contribution
                            Color lightContribution = bsdfVal.value * directLightSample.weight;

                            // Do MIS weighing of the NEE sample contribution based on its own PDF and the PDF of
                            // having sampled this direction with BSDF (bsdfVal.pdf)
                            float misWeight = 1.0f;
                            if (mis)
                            {
                                p_light = directLightSample.pdf * lightSelectionProb;
                                misWeight = powerHeuristic(p_light, bsdfVal.pdf);
                            }

                            L += lightContribution * misWeight * throughput / lightSelectionProb;
                        }
                    }
                }

                // Sample the BSDF to get the new direction and the weight
                BsdfSample bsdfSample = its.sampleBsdf(rng);
                if (bsdfSample.isInvalid())
                    break;

                // Save and update variables for next iteration
                p_bsdf = bsdfSample.pdf;
                throughput *= bsdfSample.weight;

                // Russian-Roulette to terminate small-contributing paths
                if (UseRussianRoulette && depth > RussianRouletteMinDepth && throughput.maxComponent() < 1)
                {
                    // Calculate the termination probability as 1 - the maximum component of the throughput
                    float q = max(0.0f, 1.0f - throughput.maxComponent());

                    // Randomly decide whether to terminate the path
                    if (rng.next() < q)
                        break;

                    // Prevent division by zero or very small numbers
                    if (q > 1.0f - Epsilon)
                        q = 1.0f - Epsilon;

                    // If the path continues, scale the throughput to account for the survival probability
                    throughput /= (1 - q);
                }

                // c) Trace a secondary ray in the direction determined by the BSDF sample.
                iterRay = Ray(its.position, bsdfSample.wi, iterRay.depth + 1);
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
            // Clamp the Pdfs to avoid floating point innacuracies
            pdf_a = clamp(pdf_a, MachineEpsilon, Infinity);
            pdf_b = clamp(pdf_b, MachineEpsilon, Infinity);

            // The sample from technique A should be trusted the most
            if (pdf_a == Infinity)
                return 1.0f;
            // The sample from technique B should be trusted, therefore it is 0.0f for technique A
            else if (pdf_b == Infinity)
                return 0.0f;
            else
            {
                float weight_a = sqr(pdf_a);
                float weight_b = sqr(pdf_b);
                return weight_a / (weight_a + weight_b);
            }
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