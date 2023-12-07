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
                return m_scene->evaluateBackground(ray.direction).value;
            }

            Color itsEmission = its.evaluateEmission();

            // b) Sample the BSDF to get the new direction and the weight.
            BsdfSample bsdfSample = its.sampleBsdf(rng);
            if (bsdfSample.isInvalid()) {
                return itsEmission;
            }

            Color sampledColor = bsdfSample.weight;
            accumulatedWeight *= sampledColor;
            // accumulatedWeight += itsEmission;

            // c) Trace a secondary ray in the direction determined by the BSDF sample.
            Ray secondaryRay(its.position, bsdfSample.wi.normalized(), ray.depth + 1);

            // d) If this secondary ray escapes the scene (i.e., it doesn’t hit any other surfaces),
            Intersection secondaryIts = m_scene->intersect(secondaryRay, rng);
            if (!secondaryIts) {
                // multiply its weight with the background’s emission and return the ray’s contribution.
                return accumulatedWeight * m_scene->evaluateBackground(secondaryRay.direction).value;
            }
            
            // Since there's no further bounce, we return the accumulated color which includes the BSDF weight
            // and the secondary ray contribution
            Color secondaryItsEmission = secondaryIts.evaluateEmission();
            return (itsEmission + secondaryItsEmission) * accumulatedWeight;
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