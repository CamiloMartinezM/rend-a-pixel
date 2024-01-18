/**
 * @file texture.hpp
 * @brief Contains the Texture interface, which models spatially varying properties of materials.
 */

#pragma once

#include <lightwave/color.hpp>
#include <lightwave/core.hpp>
#include <lightwave/math.hpp>

namespace lightwave
{

    /// @brief Models spatially varying material properties (e.g., images or procedural noise).
    class Texture : public Object
    {
      public:
        /**
         * @brief Returns the color at a given texture coordinate.
         * For most applications, the input point will lie in the unit square [0,1)^2, but points outside this
         * domain are also allowed.
         */
        virtual Color evaluate(const Point2 &uv) const = 0;
        /**
         * @brief Returns a scalar value at a given texture coordinate.
         * For most applications, the input point will lie in the unit square [0,1)^2, but points outside this
         * domain are also allowed.
         */
        virtual float scalar(const Point2 &uv) const
        {
            // arbitrary mapping from RGB images to scalar values (typically those will be grayscale anyway and
            // we would ideally have a separate texture interface for scalar values)
            return evaluate(uv).r();
        }
    };

    /// @brief Models a texture that is backed by an image.
    class ImageTexture : public Texture
    {
      private:
        /// @brief Describes how texture coordinates outside the [0, 1] range are handled.
        enum class BorderMode
        {
            Clamp,
            Repeat
        };

        /// @brief Describes the filtering mode used when sampling the texture.
        enum class FilterMode
        {
            Nearest,
            Bilinear,
        };

        /**
         * @brief Samples the texture color using nearest-neighbor interpolation.
         * 
         * This method will select the texel that is closest to the specified UV coordinates
         * and return its color. This is the simplest form of texture filtering, which can
         * result in a pixelated image when scaling up.
         *
         * @param uv The UV coordinates where the color is sampled. These coordinates should
         *           be within the range of the texture's resolution after any necessary adjustments.
         * @return The color of the nearest texel to the specified UV coordinates.
         */
        Color sampleNearest(const Point2 &uv) const;

        /**
         * @brief Samples the texture color using bilinear interpolation.
         * 
         * This method performs texture filtering by averaging the colors of the four
         * texels that are closest to the specified UV coordinates. This results in a
         * smoother image compared to nearest-neighbor interpolation, especially when
         * the texture is magnified.
         *
         * @param uv The UV coordinates where the color is sampled. These coordinates should
         *           be within the range of the texture's resolution after any necessary adjustments.
         * @return The color obtained by bilinearly interpolating the colors of the four surrounding texels.
         */
        Color sampleBilinear(const Point2 &uv) const;

        /**
         * @brief Applies the border handling mode to the given lattice coordinates.
         * 
         * Depending on the border mode set for the texture (clamp or repeat), this method
         * adjusts the provided lattice coordinates to adhere to the texture's border behavior.
         * This ensures that sampling is always performed within the valid range of the texture.
         *
         * @param latticeCoord The lattice coordinates to be adjusted according to the texture's border mode.
         * @return The adjusted lattice coordinates after applying the border handling.
         */
        Point2i applyBorderHandling(const Point2i &latticeCoord) const;

        /**
         * @brief Adjusts the UV coordinates based on the texture's border mode.
         * 
         * This method normalizes the UV coordinates according to the texture's border handling,
         * ensuring that they lie within the texture's valid range. It may clamp the coordinates
         * to the edges or wrap them around depending on the set border mode.
         *
         * @param uv The original UV coordinates to be adjusted.
         * @return The UV coordinates after adjustment.
         */
        Point2 adjustUV(const Point2 &uv) const;

      public:
        BorderMode m_border;
        FilterMode m_filter;

        /// @brief Shared pointer to the image resource.
        ref<Image> m_image;

        /// @brief Exposure level for the texture, affects its brightness.
        float m_exposure;

        /**
         * @brief Constructs an ImageTexture with specified properties.
         * @param properties A structure containing texture initialization properties.
         */
        ImageTexture(const Properties &properties);

        /**
         * @brief Evaluates the color of the texture at a given texture coordinate.
         * @param uv The texture coordinates, within the unit square [0,1)^2.
         * @return The color at the given texture coordinate after applying exposure and sampling.
         */
        Color evaluate(const Point2 &uv) const;

        /**
         * @brief Returns a string representation of the ImageTexture for debugging.
         * @return A string containing the formatted properties of the ImageTexture.
         */
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