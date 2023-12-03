#include <lightwave.hpp>

namespace lightwave {

/**
 * @brief Direct lighting integrator simulating a single bounce of light in the scene.
 */
    class DirectIntegrator : public SamplingIntegrator {
        public:
        DirectIntegrator(const Properties& properties)
            : SamplingIntegrator(properties) {}

        /**
         * @brief Compute the contribution of a camera-sampled ray using direct lighting.
         */
        Color Li(const Ray& ray, Sampler& rng) override {
            // Initialize the accumulated weight of the ray.
            Color accumulatedWeight = Color(1.0f);

            // Step a: Determine if the ray intersects any surfaces in the scene.
            Intersection its;
            m_scene->intersect(ray, its, rng);

            // Ray doesn't hit any surface, return background color
            if (!its) {
                return accumulatedWeight * m_scene->evaluateBackground(ray.direction).value;
            }

            // Step b: If a surface intersection occurs, generate a new direction based on the BSDF
            BsdfSample bsdfSample = its.sampleBsdf(rng);

            // Check for invalid sample
            if (bsdfSample.isInvalid()) {
                return Color(0.0f);
            }

            // Update the ray's weight by multiplying it by the sampled BSDF weight.
            accumulatedWeight *= bsdfSample.weight;

            // Step c: Trace a secondary ray in the direction determined by the BSDF sample.
            Ray secondaryRay(its.position, bsdfSample.wi.normalized());

            // Step d: If the secondary ray escapes the scene, return its contribution.
            Intersection secondaryIts;
            m_scene->intersect(secondaryRay, secondaryIts, rng);

            if (!secondaryIts) {
                return accumulatedWeight * m_scene->evaluateBackground(secondaryRay.direction).value;
            }

            // If the secondary ray hits another surface, continue the recursion
            return accumulatedWeight;
        }

        /// @brief An optional textual representation of this class, useful for debugging.
        std::string toString() const override {
            return tfm::format(
                "DirectLightingIntegrator[\n"
                "  sampler = %s,\n"
                "  image = %s,\n"
                "]",
                indent(m_sampler),
                indent(m_image)
            );
        }
    };
} // namespace lightwave

// Register the DirectLightingIntegrator class for the "direct" integrator type.
REGISTER_INTEGRATOR(DirectIntegrator, "direct")