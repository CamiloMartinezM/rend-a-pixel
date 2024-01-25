#include "../util/sampling_utils.h"
#include "../util/texture_utils.h"
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

        /// @brief Raw image with the Environment Map
        ref<Image> m_image;

        /// @brief Defines whether an Environment Map was provided. Defaults to false
        bool providedEnvironmentMap = false;

        /// @brief A unique pointer to the 2D sampling distribution calculated from an environment map image
        std::unique_ptr<Distribution2D> m_distribution;

        /// @brief A unique pointer to a MIPMap holding the environment map in terms of color data
        std::unique_ptr<MIPMap<Color>> Lmap;

        /**
         * @brief Initializes the 2D sampling distribution from the environment map.
         *
         * This function creates a 2D sampling distribution for the environment map stored in the class private member
         * `Lmap`. It converts the environment map into a scalar-valued image based on luminance values, taking into
         * account the sine of the latitude to handle distortion at the poles. This distribution is then used for
         * importance sampling the environment. map.
         */
        void initializeDistribution()
        {
            const int width = m_image->resolution().x();
            const int height = m_image->resolution().y();
            std::unique_ptr<float[]> img(new float[width * height]);

            // Compute scalar-valued image from environment map
            float filter = (float)1.0f / max(width, height);

            for (int v = 0; v < height; ++v)
            {
                float vp = static_cast<float>(v) / height;
                float sinTheta = sin(Pi * static_cast<float>(v + 0.5f) / height);
                for (int u = 0; u < width; ++u)
                {
                    float up = static_cast<float>(u) / width;
                    img[u + v * width] = Lmap->Lookup(Point2(up, vp), filter).luminance();
                    img[u + v * width] *= sinTheta;
                }
            }

            // Compute the 2D sampling distribution of the image
            m_distribution.reset(new Distribution2D(img.get(), width, height));
        }

        /// @brief Convert (u, v) to an (x, y, z) direction vector (optionally in world coordinates).
        inline std::pair<Vector, float> uvToDirection(const Point2 &uv, float mapPdf, bool toWorldTransform) const
        {
            float theta = uv.y() * Pi;
            float phi = uv.x() * 2 * Pi;
            float sinTheta = sin(theta);
            float x = sinTheta * cos(phi);
            float y = cos(theta);
            float z = sinTheta * sin(phi);

            Vector direction(x, y, z);
            float pdf = (sinTheta == 0) ? 0.0f : InvPi * Inv2Pi * mapPdf / sinTheta;

            return {direction, pdf};
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
                    m_image = imageTexture->m_image;

                    Point2i resolution = m_image->resolution();

                    // Convert Image to a suitable format for MIPMap
                    std::unique_ptr<Color[]> texels(new Color[resolution.x() * resolution.y()]);
                    for (int y = 0; y < resolution.y(); ++y)
                        for (int x = 0; x < resolution.x(); ++x)
                            texels[y * resolution.x() + x] = m_image->get(Point2i(x, y));

                    // Initialize MIPMap
                    Lmap.reset(new MIPMap<Color>(resolution, texels.get()));

                    providedEnvironmentMap = true;
                    initializeDistribution();
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

            // Return the evaluated texture value at the mapped UV coordinates
            Vector2 warped(u, v);
            return {.value = m_texture->evaluate(warped)};
        }

        DirectLightSample sampleDirect(const Point &origin, Sampler &rng) const override
        {
            // Use Uniform Sampling if an HDRI was not provided or if the global constant is set so that the
            // Improved Environment Sampling routine is not used
            if (!providedEnvironmentMap || !UseImprovedEnvSampling)
            {
                Vector direction = squareToUniformSphere(rng.next2D());
                auto E = evaluate(direction);
                return {.wi = direction, .weight = E.value / Inv4Pi, .distance = Infinity, .pdf = uniformSpherePdf()};
            }

            // Improved Environment Sampling
            float mapPdf = 0.0f;
            Point2 uv = m_distribution->sampleContinuous(rng.next2D(), &mapPdf);
            if (mapPdf == 0)
                return {.wi = Vector(), .weight = Color::black(), .distance = Infinity, .pdf = 0.0f};

            // Convert the sampled UV to a direction and compute the PDF for spherical mapping
            auto [direction, pdf] = uvToDirection(uv, mapPdf, true);
            auto E = evaluate(direction);
            return {.wi = direction, .weight = E.value / pdf, .distance = Infinity, .pdf = pdf};
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