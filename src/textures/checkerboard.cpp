#include <lightwave.hpp>

namespace lightwave {

    class CheckerboardTexture : public Texture {
        Color m_color0;
        Color m_color1;
        Vector2 m_scale;

        public:
        CheckerboardTexture(const Properties& properties) {
            m_color0 = properties.get<Color>("color0", Color(0));
            m_color1 = properties.get<Color>("color1", Color(1));
            m_scale = properties.get<Vector2>("scale", Vector2(1, 1));
        }

        Color evaluate(const Point2& uv) const override {
            // Rescale UV coordinates
            Point2 scaledUV = Point2(uv.x() * m_scale.x(), uv.y() * m_scale.y());

            // Use floor and modulo to create checkerboard pattern
            int checkerboard = static_cast<int>(std::floor(scaledUV.x()) + std::floor(scaledUV.y())) % 2;

            // Return color based on checkerboard pattern
            return (checkerboard == 0) ? m_color0 : m_color1;
        }

        std::string toString() const override {
            return tfm::format("CheckerboardTexture[\n"
                               "  color0 = %s\n"
                               "  color1 = %s\n"
                               "  scale = %s\n"
                               "]",
                               indent(m_color0), indent(m_color1), indent(m_scale));
        }
    };
} // namespace lightwave

REGISTER_TEXTURE(CheckerboardTexture, "checkerboard")