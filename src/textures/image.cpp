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
            // First, adjust the UV coordinates based on the border mode.
            Point2 adjustedUV = uv;
            switch (m_border) {
            case BorderMode::Clamp:
                adjustedUV.x() = std::min(std::max(adjustedUV.x(), 0.0f), 1.0f);
                adjustedUV.y() = std::min(std::max(adjustedUV.y(), 0.0f), 1.0f);
                break;
            case BorderMode::Repeat:
                adjustedUV.x() = adjustedUV.x() - std::floor(adjustedUV.x());
                adjustedUV.y() = adjustedUV.y() - std::floor(adjustedUV.y());
                break;
            }

            // Now sample the color from the image based on the filter mode.
            Color sampledColor;
            switch (m_filter) {
            case FilterMode::Nearest:
                sampledColor = sampleNearest(adjustedUV);
                break;
            case FilterMode::Bilinear:
                sampledColor = sampleBilinear(adjustedUV);
                break;
            }

            // Apply the exposure correction.
            sampledColor *= std::pow(2.0f, m_exposure);

            return sampledColor;
        }

        Color sampleNearest(const Point2& uv) const {
            // Apply border handling and get the color directly using the public interface
            Point2 adjustedUV = applyBorderMode(uv);
            return (*m_image)(adjustedUV);
        }

        Color sampleBilinear(const Point2& uv) const {
            // Apply border handling for the UV coordinates
            Point2 adjustedUV = applyBorderMode(uv);

            // Calculate the base coordinates and fractional parts
            Point2i baseCoord = Point2i(std::floor(adjustedUV.x() * m_image->resolution().x()),
                                        std::floor(adjustedUV.y() * m_image->resolution().y()));
            Point2 frac(adjustedUV.x() * m_image->resolution().x() - baseCoord.x(),
                        adjustedUV.y() * m_image->resolution().y() - baseCoord.y());

            // Retrieve the colors of the four neighboring pixels
            Color c00 = m_image->operator()(Point2(baseCoord.x() / static_cast<float>(m_image->resolution().x()),
                baseCoord.y() / static_cast<float>(m_image->resolution().y())));
            Color c10 = m_image->operator()(Point2((baseCoord.x() + 1) / static_cast<float>(m_image->resolution().x()),
                baseCoord.y() / static_cast<float>(m_image->resolution().y())));
            Color c01 = m_image->operator()(Point2(baseCoord.x() / static_cast<float>(m_image->resolution().x()),
                (baseCoord.y() + 1) / static_cast<float>(m_image->resolution().y())));
            Color c11 = m_image->operator()(Point2((baseCoord.x() + 1) / static_cast<float>(m_image->resolution().x()),
                (baseCoord.y() + 1) / static_cast<float>(m_image->resolution().y())));

            // Perform bilinear interpolation
            Color interpolatedColor = (c00 * (1 - frac.x()) * (1 - frac.y())) +
                (c10 * frac.x() * (1 - frac.y())) +
                (c01 * (1 - frac.x()) * frac.y()) +
                (c11 * frac.x() * frac.y());

            return interpolatedColor;
        }

        // Utility method to apply border mode adjustments to UV coordinates
        Point2 applyBorderMode(const Point2& uv) const {
            Point2 adjustedUV = uv;
            if (m_border == BorderMode::Clamp) {
                adjustedUV.x() = std::clamp(uv.x(), 0.0f, 1.0f);
                adjustedUV.y() = std::clamp(uv.y(), 0.0f, 1.0f);
            }
            else if (m_border == BorderMode::Repeat) {
                adjustedUV.x() = uv.x() - std::floor(uv.x());
                adjustedUV.y() = uv.y() - std::floor(uv.y());
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
