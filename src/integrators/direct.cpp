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
            Color accumulatedWeight = Color(1.0f);

            Intersection its = m_scene->intersect(ray, rng);
            if (!its) {
                // The ray misses all objects and hits the background.
                // return accumulatedWeight * m_scene->evaluateBackground(ray.direction).value;
                return accumulatedWeight * m_scene->evaluateBackground(ray.direction).value;
            }
            
            // Sample the BSDF to get the new direction and the weight.
            BsdfSample bsdfSample = its.sampleBsdf(rng);
            if (bsdfSample.isInvalid()) {
                return Color(0.0f); // Invalid BSDF sample, return black.
            }

            accumulatedWeight *= bsdfSample.weight;

            // Construct the secondary ray from the intersection point along the sampled direction.
            Ray secondaryRay(its.position, bsdfSample.wi.normalized());
            Intersection secondaryIts = m_scene->intersect(secondaryRay, rng);

            if (!secondaryIts) {
                // The secondary ray misses all objects and hits the background.
                return accumulatedWeight * m_scene->evaluateBackground(secondaryRay.direction).value;
            }

            return accumulatedWeight;
        }

        /// @brief An optional textual representation of this class, useful for debugging.
        std::string toString() const override {
            return tfm::format(
                "DirectIntegrator[\n"
                "  sampler = %s,\n"
                "  image = %s,\n"
                "]",
                indent(m_sampler),
                indent(m_image)
            );
        }
    };
} // namespace lightwave

// Register the DirectIntegrator class for the "direct" integrator type.
REGISTER_INTEGRATOR(DirectIntegrator, "direct")