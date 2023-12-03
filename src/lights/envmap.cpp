#include <lightwave.hpp>

namespace lightwave {

    class EnvironmentMap final : public BackgroundLight {
        /// @brief The texture to use as background
        ref<Texture> m_texture;

        /// @brief An optional transform from local-to-world space
        ref<Transform> m_transform;

        public:
        EnvironmentMap(const Properties& properties) {
            m_texture = properties.getChild<Texture>();
            m_transform = properties.getOptionalChild<Transform>();
        }

        BackgroundLightEval evaluate(const Vector& direction) const override {
            // hints:
            // * if (m_transform) { transform direction vector from world to local
            // coordinates }
            // * find the corresponding pixel coordinate for the given local
            // direction
            // Apply the environment map transform to the input world direction
            Vector transformedDirection = direction;
            if (m_transform)
                transformedDirection = m_transform->inverse(direction);
            transformedDirection = transformedDirection.normalized();
        
            // Convert 3D Cartesian coordinates to spherical coordinates (θ, φ)
            float theta = std::acos(transformedDirection.y());                          // θ: polar angle
            float phi = std::atan2(transformedDirection.z(), transformedDirection.x()); // φ: azimuthal angle

            // Normalize θ and φ to [0, 1] range for texture coordinates
            float u = 0.5f - phi / (2 * Pi);
            float v = 1.0f - theta / Pi;

            // Return the evaluated texture value at the mapped coordinates
            Vector2 warped(u, v);
            return {
                .value = m_texture->evaluate(warped),
            };
        }

        DirectLightSample sampleDirect(const Point& origin,
                                       Sampler& rng) const override {
            Vector direction = squareToUniformSphere(rng.next2D());
            auto E = evaluate(direction);

            // implement better importance sampling here, if you ever need it
            // (useful for environment maps with bright tiny light sources, like the
            // sun for example)

            return {
                .wi = direction,
                .weight = E.value / Inv4Pi,
                .distance = Infinity,
            };
        }

        std::string toString() const override {
            return tfm::format("EnvironmentMap[\n"
                               "  texture = %s,\n"
                               "  transform = %s\n"
                               "]",
                               indent(m_texture), indent(m_transform));
        }
    };

} // namespace lightwave

REGISTER_LIGHT(EnvironmentMap, "envmap")