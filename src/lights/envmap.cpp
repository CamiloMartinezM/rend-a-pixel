#include <lightwave.hpp>

namespace lightwave
{
    class EnvironmentMap final : public BackgroundLight
    {
      private:
        /// @brief The texture to use as background
        ref<Texture> m_texture;

        /// @brief An optional transform from local-to-world space
        ref<Transform> m_transform;

        ref<Image> m_image;

        struct Distribution2D
        {
            std::vector<float> func;
            std::vector<float> cdf;
            int count;
        };

        Distribution2D m_rowDistribution;
        std::vector<Distribution2D> m_colDistributions;

        void initializeDistributions()
        {
            int width = m_image->resolution().x();
            int height = m_image->resolution().y();
            std::vector<float> imageLuminance(width * height);

            // Compute luminance for each pixel and row marginal CDF
            m_rowDistribution.func.resize(height);
            m_rowDistribution.cdf.resize(height + 1);
            for (int y = 0; y < height; ++y)
            {
                float rowSum = 0;
                for (int x = 0; x < width; ++x)
                {
                    Color color = m_image->operator()(Point2i(x, y));
                    imageLuminance[y * width + x] = color.luminance();
                    rowSum += color.luminance();
                }
                m_rowDistribution.func[y] = rowSum;
            }
            computeCumulativeDistribution(m_rowDistribution);

            // Compute column CDFs for each row
            m_colDistributions.resize(height);
            for (int y = 0; y < height; ++y)
            {
                m_colDistributions[y].func.resize(width);
                m_colDistributions[y].cdf.resize(width + 1);
                for (int x = 0; x < width; ++x)
                {
                    m_colDistributions[y].func[x] = imageLuminance[y * width + x];
                }
                computeCumulativeDistribution(m_colDistributions[y]);
            }
        }

        void computeCumulativeDistribution(Distribution2D &dist)
        {
            dist.cdf[0] = 0;
            for (size_t i = 1; i < dist.cdf.size(); ++i)
            {
                dist.cdf[i] = dist.cdf[i - 1] + dist.func[i - 1] / dist.func.size();
            }
        }

        float sampleContinuousDistribution(const Distribution2D &dist, float sample, float &pdf) const
        {
            auto it = std::lower_bound(dist.cdf.begin(), dist.cdf.end(), sample);
            int offset = std::max(0, static_cast<int>(it - dist.cdf.begin() - 1));
            float du = sample - dist.cdf[offset];
            if ((dist.cdf[offset + 1] - dist.cdf[offset]) > 0)
            {
                du /= (dist.cdf[offset + 1] - dist.cdf[offset]);
            }
            pdf = dist.func[offset] / (dist.func.size() * (dist.cdf.back() - dist.cdf.front()));
            return (offset + du) / dist.func.size();
        }

      public:
        EnvironmentMap(const Properties &properties)
        {
            if (properties.has("filename"))
            {
                m_image = std::make_shared<Image>(properties);
            }
            else
            {
                m_image = properties.getChild<Image>();
            }
            m_texture = properties.getChild<Texture>();
            m_transform = properties.getOptionalChild<Transform>();
            initializeDistributions();
        }

        BackgroundLightEval evaluate(const Vector &direction) const override
        {
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
            float v = theta / Pi;

            // Return the evaluated texture value at the mapped UV coordinates
            Vector2 warped(u, v);
            Color value = m_texture->evaluate(warped);

            return {
                .value = value,
            };
        }

        DirectLightSample sampleDirect(const Point &origin, Sampler &rng) const override
        {
            float pdfRow, pdfCol;
            float uRow = rng.next();
            float uCol = rng.next();

            float rowSample = sampleContinuousDistribution(m_rowDistribution, uRow, pdfRow);
            float colSample = sampleContinuousDistribution(
                m_colDistributions[static_cast<int>(rowSample * m_image->resolution().y())], uCol, pdfCol);

            float theta = rowSample * Pi;
            float phi = colSample * 2 * Pi;
            float sinTheta = sin(theta);

            Vector direction(cos(phi) * sinTheta, cos(theta), sin(phi) * sinTheta);

            if (sinTheta == 0)
                pdfRow = 0;

            float pdf = Inv2Pi * InvPi * pdfRow * pdfCol / sinTheta;
            Color value = evaluate(direction).value;

            return {
                .wi = direction,
                .weight = value / pdf,
                .distance = Infinity,
            };
        }

        std::string toString() const override
        {
            return tfm::format("EnvironmentMap[\n"
                               "  texture = %s,\n"
                               "  transform = %s\n"
                               "]",
                               indent(m_texture), indent(m_transform));
        }
    };

} // namespace lightwave

REGISTER_LIGHT(EnvironmentMap, "envmap")