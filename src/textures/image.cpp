#include <lightwave.hpp>

namespace lightwave
{

    class ImageTexture : public Texture
    {
        enum class BorderMode
        {
            Clamp,
            Repeat,
        };

        enum class FilterMode
        {
            Nearest,
            Bilinear,
        };

        ref<Image> m_image;
        float m_exposure;
        BorderMode m_border;
        FilterMode m_filter;

      public:
        ImageTexture(const Properties &properties)
        {
            if (properties.has("filename"))
            {
                m_image = std::make_shared<Image>(properties);
            }
            else
            {
                m_image = properties.getChild<Image>();
            }
            m_exposure = properties.get<float>("exposure", 1);
            m_border = properties.getEnum<BorderMode>("border", BorderMode::Repeat,
                                                      {
                                                          {"clamp", BorderMode::Clamp},
                                                          {"repeat", BorderMode::Repeat},
                                                      });
            m_filter = properties.getEnum<FilterMode>("filter", FilterMode::Bilinear,
                                                      {
                                                          {"nearest", FilterMode::Nearest},
                                                          {"bilinear", FilterMode::Bilinear},
                                                      });
        }

        // Main method to call
        Color evaluate(const Point2 &uv) const override
        {
            Point2 adjustedUV = adjustUV(uv);
            Color sampledColor =
                (m_filter == FilterMode::Nearest) ? sampleNearest(adjustedUV) : sampleBilinear(adjustedUV);
            return sampledColor * m_exposure;
        }

        Color sampleNearest(const Point2 &uv) const
        {
            // Convert UV to lattice coordinates
            Point2i latticeCoord =
                Point2i(std::floor(uv.x() * m_image->resolution().x()), std::floor(uv.y() * m_image->resolution().y()));

            // Apply border handling
            latticeCoord = applyBorderHandling(latticeCoord);

            return m_image->operator()(Point2i(latticeCoord));
        }

        Color sampleBilinear(const Point2 &uv) const
        {
            // Shift UV coordinates to align with pixel centers
            Point2 shiftedUV =
                Point2(uv.x() * m_image->resolution().x() - 0.5f, uv.y() * m_image->resolution().y() - 0.5f);

            // Base lattice coordinates
            Point2i baseCoord = Point2i(std::floor(shiftedUV.x()), std::floor(shiftedUV.y()));

            // Neighboring lattice coordinates for bilinear interpolation
            Point2i coords[4] = {applyBorderHandling(baseCoord),
                                 applyBorderHandling(Point2i(baseCoord.x() + 1, baseCoord.y())),
                                 applyBorderHandling(Point2i(baseCoord.x(), baseCoord.y() + 1)),
                                 applyBorderHandling(Point2i(baseCoord.x() + 1, baseCoord.y() + 1))};

            // Calculate fractional part for interpolation
            Point2 frac = Point2(shiftedUV.x() - baseCoord.x(), shiftedUV.y() - baseCoord.y());

            // Fetch texel values for each coordinate
            Color texelValues[4];
            for (int i = 0; i < 4; ++i)
                texelValues[i] = m_image->operator()(Point2i(coords[i]));

            // Perform bilinear interpolation
            Color interpolatedColor = texelValues[0] * (1 - frac.x()) * (1 - frac.y()) +
                                      texelValues[1] * frac.x() * (1 - frac.y()) +
                                      texelValues[2] * (1 - frac.x()) * frac.y() + texelValues[3] * frac.x() * frac.y();

            return interpolatedColor;
        }

        Point2i applyBorderHandling(const Point2i &latticeCoord) const
        {
            if (m_border == BorderMode::Clamp)
            {
                return Point2i(std::clamp(latticeCoord.x(), 0, m_image->resolution().x() - 1),
                               std::clamp(latticeCoord.y(), 0, m_image->resolution().y() - 1));
            }
            else if (m_border == BorderMode::Repeat)
            {
                return Point2i((latticeCoord.x() % m_image->resolution().x() + m_image->resolution().x()) %
                                   m_image->resolution().x(),
                               (latticeCoord.y() % m_image->resolution().y() + m_image->resolution().y()) %
                                   m_image->resolution().y());
            }
            return latticeCoord; // Default case
        }

        Point2 adjustUV(const Point2 &uv) const
        {
            Point2 adjustedUV = uv;
            switch (m_border)
            {
            case BorderMode::Clamp:
                adjustedUV.x() = std::clamp(uv.x(), 0.0f, 1.0f);
                adjustedUV.y() = std::clamp(uv.y(), 0.0f, 1.0f);
                break;
            case BorderMode::Repeat:
                adjustedUV.x() = uv.x() - std::floor(uv.x());
                adjustedUV.y() = uv.y() - std::floor(uv.y());
                break;
            }
            adjustedUV.y() = 1.0f - adjustedUV.y();
            return adjustedUV;
        }

        std::string toString() const override
        {
            return tfm::format("ImageTexture[\n"
                               "  image = %s,\n"
                               "  exposure = %f,\n"
                               "]",
                               indent(m_image), m_exposure);
        }
    };
} // namespace lightwave

REGISTER_TEXTURE(ImageTexture, "image")