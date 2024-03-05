// #include "../util/sampling_utils.h"
// #include "../util/texture_utils.h"
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

        /// @brief Height of the Environment Map, if it exists.
        int m_height = 0;

        /// @brief Width of the Environment Map, if it exists.
        int m_width = 0;

        /// @brief Defines whether an Environment Map was provided. Defaults to false
        bool providedEnvironmentMap = false;

        struct Distribution2D
        {
            std::vector<float> func;
            std::vector<float> cdf;
            int count;
        };

        Distribution2D m_rowDistribution;
        std::vector<Distribution2D> m_colDistributions;

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

        /**
         * @brief Initializes the 2D sampling distribution from the environment map.
         *
         * This function creates a 2D sampling distribution for the environment map on luminance values, taking into
         * account the sine of the latitude to handle distortion at the poles. This distribution is then used for
         * importance sampling the environment map.
         */
        void initializeDistribution2D()
        {
            std::unique_ptr<float[]> imageLuminance(new float[m_width * m_height]);

            // Compute luminance for each pixel and row marginal CDF
            m_rowDistribution.func.resize(m_height);
            m_rowDistribution.cdf.resize(m_width + 1);
            for (int y = 0; y < m_height; ++y)
            {
                float rowSum = 0;
                float sinTheta = sin(Pi * static_cast<float>(y + 0.5f) / (float)m_height);
                for (int x = 0; x < m_width; ++x)
                {
                    Point2 uv((x + 0.5f) / m_width, (y + 0.5f) / m_height);
                    Color c = m_texture->evaluate(uv);
                    imageLuminance[y * m_width + x] = c.luminance();
                    rowSum += c.luminance();
                }
                m_rowDistribution.func[y] = rowSum;
            }
            computeCumulativeDistribution(m_rowDistribution);

            // Compute column CDFs for each row
            m_colDistributions.resize(m_height);
            for (int y = 0; y < m_height; ++y)
            {
                m_colDistributions[y].func.resize(m_width);
                m_colDistributions[y].cdf.resize(m_width + 1);
                for (int x = 0; x < m_width; ++x)
                {
                    m_colDistributions[y].func[x] = imageLuminance[y * m_width + x];
                }
                computeCumulativeDistribution(m_colDistributions[y]);
            }
        }

        /// @brief Convert (u, v) to an (x, y, z) direction vector. It is the opposite conversion done in "evaluate".
        inline std::pair<Vector, float> uvToDirection(const Point2 &uv, float mapPdf) const
        {
            // Convert uv back to spherical coordinates
            float theta = uv.y() * Pi;
            float phi = Pi * (1 - 2 * uv.x());
            float sinTheta = sin(theta);

            // Convert spherical coordinates back to 3D Cartesian coordinates
            Vector direction(sinTheta * cos(phi), cos(theta), sinTheta * sin(phi));

            // Compute Pdf for sampled infinite light direction and convert it to solid-angle measure
            float pdf = (sinTheta == 0) ? 0.0f : InvPi * Inv2Pi * mapPdf / sinTheta;

            // Conversion of the pdf from area measure (in uv space) to solid angle measure
            return {direction.normalized(), pdf};
        }

      public:
        EnvironmentMap(const Properties &properties)
        {
            m_texture = properties.getChild<Texture>();
            m_transform = properties.getOptionalChild<Transform>();
            providedEnvironmentMap = false;
            if (UseImprovedEnvSampling)
            {
                auto imageTexture = dynamic_cast<ImageTexture *>(m_texture.get());
                if (imageTexture)
                {
                    Point2i m_resolution = imageTexture->m_image->resolution();
                    m_width = m_resolution.x();
                    m_height = m_resolution.y();
                    providedEnvironmentMap = true;
                    initializeDistribution2D();
                }
            }
        }

        BackgroundLightEval evaluate(const Vector &direction) const override
        {
            // hints:
            // * if (m_transform) { transform direction vector from world to local coordinates }
            // * find the corresponding pixel coordinate for the given local direction
            // Apply the environment map transform to the input world direction
            Vector transformedDirection = direction;
            if (m_transform)
                transformedDirection = m_transform->inverse(direction);
            transformedDirection = transformedDirection.normalized();

            // Convert 3D Cartesian coordinates to spherical coordinates (θ, φ)
            float theta = acos(transformedDirection.y());                          // θ: polar angle
            float phi = atan2(transformedDirection.z(), transformedDirection.x()); // φ: azimuthal angle

            // Normalize θ and φ to [0, 1] range for texture coordinates
            float u = Inv2Pi * (Pi - phi);
            float v = theta * InvPi;

            // Return the evaluated texture value at the mapped uv coordinates
            Vector2 warped(u, v);
            float pdf = uniformSpherePdf();
            Color colorEval = m_texture->evaluate(warped);
            return {.value = colorEval, .pdf = pdf};
        }

        DirectLightSample sampleDirect(const Point &origin, Sampler &rng) const override
        {
            // Use Uniform Sampling if an HDRI was not provided or if the global constant is set so that the
            // Improved Environment Sampling routine is not used
            if (!providedEnvironmentMap || !UseImprovedEnvSampling)
            {
                Vector direction = squareToUniformSphere(rng.next2D());
                auto E = evaluate(direction);
                return {.wi = direction, .weight = E.value / Inv4Pi, .distance = Infinity, .pdf = E.pdf};
            }

            // Improved Environment Sampling routine
            float pdfRow, pdfCol;
            float uRow = rng.next();
            float uCol = rng.next();

            float u = sampleContinuousDistribution(m_rowDistribution, uRow, pdfRow);
            float v = sampleContinuousDistribution(m_colDistributions[static_cast<int>(u * m_height)], uCol, pdfCol);
            float mapPdf = pdfRow * pdfCol;

            // Return an invalid sample, just in case
            if (mapPdf == 0.0f)
                return DirectLightSample::invalid();

            Point2 uv(u, v);
            auto [direction, pdf] = uvToDirection(uv, mapPdf);

            // Return an invalid sample, just in case
            if (pdf == 0.0f)
                return DirectLightSample::invalid();

            // Lookup the environment map using the uv coordinates
            auto E = m_texture->evaluate(uv);
            return {.wi = direction, .weight = E / pdf, .distance = Infinity, .pdf = pdf};
        }

        LightType getLightType() const override
        {
            return LightType::Envmap;
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