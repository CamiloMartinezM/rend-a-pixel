#include <lightwave.hpp>

namespace lightwave {

    class ImageTexture : public Texture {
        enum class BorderMode {
            Clamp,
            Repeat,
        };

        enum class FilterMode {
            Nearest,
            Bilinear,
        };

        ref<Image> m_image;
        float m_exposure;
        BorderMode m_border;
        FilterMode m_filter;

        public:
        ImageTexture(const Properties& properties) {
            if (properties.has("filename")) {
                m_image = std::make_shared<Image>(properties);
            }
            else {
                m_image = properties.getChild<Image>();
            }
            m_exposure = properties.get<float>("exposure", 1);

            m_border =
                properties.getEnum<BorderMode>("border", BorderMode::Repeat,
                                               {
                                                   { "clamp", BorderMode::Clamp },
                                                   { "repeat", BorderMode::Repeat },
                                               });

            m_filter = properties.getEnum<FilterMode>(
                "filter", FilterMode::Bilinear,
                {
                    { "nearest", FilterMode::Nearest },
                    { "bilinear", FilterMode::Bilinear },
                });
        }

        Color evaluate(const Point2& uv) const override {
            Point2 adjustedUV = adjustUV(uv);
            Color sampledColor = (m_filter == FilterMode::Nearest) ? sampleNearest(adjustedUV) : sampleBilinear(adjustedUV);
            return sampledColor * m_exposure;
        }

        Color sampleNearest(const Point2& uv) const {
            // Map UV coordinates to shifted lattice coordinates
            Point2 adjustedUV = Point2(uv.x() * m_image->resolution().x(), uv.y() * m_image->resolution().y()) - Point2(0.5f, 0.5f);
            Point2i latticeCoord = Point2i(std::floor(adjustedUV.x()), std::floor(adjustedUV.y()));

            latticeCoord = applyBorderHandling(latticeCoord);

            return m_image->operator()(Point2(static_cast<float>(latticeCoord.x()) / m_image->resolution().x(),
                                              static_cast<float>(latticeCoord.y()) / m_image->resolution().y()));
        }

        Color sampleBilinear(const Point2& uv) const {
            Point2 adjustedUV = Point2(uv.x() * m_image->resolution().x(), uv.y() * m_image->resolution().y()) - Point2(0.5f, 0.5f);
            
            // Base lattice coordinates
            Point2i baseCoord = Point2i(std::floor(adjustedUV.x()), std::floor(adjustedUV.y()));
            
            // Fractional part for interpolation
            Point2 frac(adjustedUV.x() - baseCoord.x(), adjustedUV.y() - baseCoord.y());

            // Neighboring lattice coordinates for bilinear interpolation
            Point2i coords[4] = {
                applyBorderHandling(baseCoord),
                applyBorderHandling(Point2i(baseCoord.x() + 1, baseCoord.y())),
                applyBorderHandling(Point2i(baseCoord.x(), baseCoord.y() + 1)),
                applyBorderHandling(Point2i(baseCoord.x() + 1, baseCoord.y() + 1))
            };

            // Fetch texel values for each of the four neighboring coordinates
            Color texelValues[4] = {
                m_image->operator()(Point2(static_cast<float>(coords[0].x()) / m_image->resolution().x(),
                                           static_cast<float>(coords[0].y()) / m_image->resolution().y())),
                m_image->operator()(Point2(static_cast<float>(coords[1].x()) / m_image->resolution().x(),
                                           static_cast<float>(coords[1].y()) / m_image->resolution().y())),
                m_image->operator()(Point2(static_cast<float>(coords[2].x()) / m_image->resolution().x(),
                                           static_cast<float>(coords[2].y()) / m_image->resolution().y())),
                m_image->operator()(Point2(static_cast<float>(coords[3].x()) / m_image->resolution().x(),
                                           static_cast<float>(coords[3].y()) / m_image->resolution().y()))
            };

            // Perform bilinear interpolation between the texel values
            Color interpolatedColor = texelValues[0] * (1 - frac.x()) * (1 - frac.y()) +
                                      texelValues[1] * frac.x() * (1 - frac.y()) +
                                      texelValues[2] * (1 - frac.x()) * frac.y() +
                                      texelValues[3] * frac.x() * frac.y();

            return interpolatedColor;
        }

        Point2i applyBorderHandling(const Point2i& latticeCoord) const {
            Point2i adjustedCoord = latticeCoord;
            if (m_border == BorderMode::Clamp) {
                adjustedCoord.x() = std::clamp(latticeCoord.x(), 0, m_image->resolution().x() - 1);
                adjustedCoord.y() = std::clamp(latticeCoord.y(), 0, m_image->resolution().y() - 1);
            } else if (m_border == BorderMode::Repeat) {
                adjustedCoord.x() = latticeCoord.x() % m_image->resolution().x();
                adjustedCoord.y() = latticeCoord.y() % m_image->resolution().y();
            }
            return adjustedCoord;
        }

        Point2 adjustUV(const Point2& uv) const {
            Point2 adjustedUV = uv;
            switch (m_border) {
                case BorderMode::Clamp:
                    adjustedUV.x() = std::clamp(uv.x(), 0.0f, 1.0f);
                    adjustedUV.y() = std::clamp(uv.y(), 0.0f, 1.0f);
                    break;
                case BorderMode::Repeat:
                    adjustedUV.x() = uv.x() - std::floor(uv.x());
                    adjustedUV.y() = uv.y() - std::floor(uv.y());
                    break;
            }
            return adjustedUV;
        }
        
        std::string toString() const override {
            return tfm::format("ImageTexture[\n"
                               "  image = %s,\n"
                               "  exposure = %f,\n"
                               "]",
                               indent(m_image), m_exposure);
        }
    };
} // namespace lightwave

REGISTER_TEXTURE(ImageTexture, "image")