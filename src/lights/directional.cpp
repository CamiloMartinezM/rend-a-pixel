#include <lightwave.hpp>

namespace lightwave {

    class DirectionalLight final : public Light {

        Vector m_direction; // Direction from which the light comes
        Color m_intensity;  // Intensity of the light

        public:

        DirectionalLight(const Properties& properties) {
            m_direction = properties.get<Vector>("direction").normalized();
            m_intensity = properties.get<Color>("intensity");
        }

        DirectLightSample sampleDirect(const Point& origin, Sampler& rng) const override {
            return DirectLightSample{
                .wi = m_direction,  // Reverse the direction for incoming light
                .weight = m_intensity, // Intensity of the light remains constant
                .distance = Infinity // Distance is infinity for a directional light
            };
        }

        bool canBeIntersected() const override { return false; }

        std::string toString() const override {
            return tfm::format(
                "DirectionalLight[\n"
                "  direction = %s,\n"
                "  intensity = %s\n"
                "]",
                indent(m_direction), indent(m_intensity)
            );
        }
    };
} // namespace lightwave

REGISTER_LIGHT(DirectionalLight, "directional")