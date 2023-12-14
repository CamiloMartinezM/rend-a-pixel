#include <lightwave.hpp>

namespace lightwave {

    class PointLight final : public Light {

        Point m_position; // Position of the point light in world space
        Color m_power;    // Power (flux) of the point light

        public:

        PointLight(const Properties& properties) {
            m_position = properties.get<Point>("position");
            m_power = properties.get<Color>("power");
        }

        DirectLightSample sampleDirect(const Point& origin,
                                       Sampler& rng) const override {
            // Direction from the point on the surface to the light source
            Vector direction = m_position - origin;

            // Squared distance from the point on the surface to the light source
            float distanceSquared = direction.lengthSquared();
            direction = direction.normalized();

            // Convert power to intensity (power per unit area)
            // Assuming isotropic point light, intensity falls off as the square of the distance
            // I = P / (4 * pi * r^2)
            Color intensity = m_power / (4 * Pi * distanceSquared);

            return DirectLightSample{
                .wi = direction,        
                .weight = intensity,
                .distance = std::sqrt(distanceSquared)
            };
        }

        bool canBeIntersected() const override { return false; }

        std::string toString() const override {
            return tfm::format("PointLight[\n"
                               "]");
        }
    };
} // namespace lightwave

REGISTER_LIGHT(PointLight, "point")
