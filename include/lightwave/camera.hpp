/**
 * @file camera.hpp
 * @brief Contains the Camera interface and related structures.
 */

#pragma once

#include <lightwave/color.hpp>
#include <lightwave/core.hpp>
#include <lightwave/math.hpp>
#include <lightwave/properties.hpp>
#include <lightwave/transform.hpp>
#include <lightwave/warp.hpp>

namespace lightwave
{
    /// @brief Defines whether Bokeh effects should be enforced on the scene. If it is true, then the required variables
    /// for the Bokeh bokehConfig will be looked for inside <camera type="thinlens"> .
    static const bool UseBokehEffects = true;

    /// @brief Stores the configuration to generate Bokeh shapes on a scene.
    struct Bokeh
    {
        /// @brief How many blades or vertices are used in the lens aperture, with 0 indicating a perfectly round iris.
        /// This parameter simulates the effect of a real camera's aperture blades on the bokeh shape.
        int blades;

        /// @brief Radius of the inner obstructing sphere, simulating the lens' physical constraints.
        /// This parameter affects the size of the bokeh shapes.
        float radius;

        /// @brief Controls how much the projection differs from a perfect circle.
        /// Higher values result in more pronounced deviations from circular bokeh shapes.
        float rateOfChange;

        /// @brief Step angle in radians between vertices of the polygonal aperture shape, calculated as 2 Pi / blades.
        /// It represents the angular distance between consecutive vertices of the aperture polygon, ensuring uniform
        /// spacing.
        float step;

        /// @brief Dynamically allocated array of points defining the vertices of the polygonal aperture shape when
        /// blades > 2. These points are used to construct non-circular bokeh shapes that mimic the effect of aperture
        /// blades.
        std::unique_ptr<Point2[]> corners;

        /// @brief Describes the ratio of the two circles which determine the shape of the bokeh.
        /// This parameter can be used to simulate the effect of an optical vignetting.
        float innerRadius;

        /// @brief Controls the distribution of brightness within the bokeh shape.
        /// 0 even destribution, <0 center is emphasized, >0 edges are emphasized
        float weightDistr;

        /// @brief Determines the extent to which positions not part of the edge are reduced in brightness, effectively
        /// controlling the contrast between the center and the edges of the bokeh shapes.
        float weightStrength;

        /// @brief Precomputed integral value used for normalization in the weighting function, ensuring that the bokeh
        /// effect maintains light conservation.
        float integral;
    };

    /// @brief Initializes the Bokeh shapes based on the configuration stored inside the given parameter.
    void buildBokehShapes(Bokeh &bokehConfig);

    /// @brief The result of sampling a Camera.
    struct CameraSample
    {
        /// @brief The direction vector, pointing away from the camera.
        Ray ray;

        /// @brief The weight of the sample.
        Color weight;
    };

    /// @brief A Camera, representing the relationship between pixel coordinates and rays.
    class Camera : public Object
    {
      protected:
        /// @brief The resolution of the image that is being rendered.
        Vector2i m_resolution;

        /// @brief The transform that leads from local coordinates to world space coordinates.
        ref<Transform> m_transform;

        /// @brief Configuration for simulating the bokeh effect in a scene, based on artificial camera's aperture
        /// characteristics. We use this variable during ray generation to modify rays in a way that simulates the depth
        /// of field and bokeh effects produced by a physical camera lens.
        Bokeh bokehConfig;

        /**
         * @brief Samples a point on the lens and biases it according to the bokeh configuration.
         *
         * @param sample A point sampled from the unit square.
         * @param pFilm The position on the film plane corresponding to the sample, used for vignetting effects.
         * @param weight Reference that will be updated with the weight of the sample based on bokeh characteristics.
         * @return A Point2 representing the biased sample on the lens.
         * 
         * Credits: https://knork.org/realistic-bokeh.html
         *
         * The function operates in three main modes:
         *
         * 1. Polygonal Aperture: If the number of aperture blades is three or more, the function generates a polygonal
         * bokeh shape. It calculates the vertices of the polygon (A, B, C) based on the number of blades and uses these
         * to bias the sample.
         *
         * 2. Optical Vignetting: If the rate of change parameter is significant, the function simulates optical
         * vignetting. This involves adjusting the sample based on its position on the film plane and the specified rate
         * of change.
         *
         * 3. Circular Aperture: If neither of the above conditions is met, the function defaults to using a circular
         * aperture model. This involves simply mapping the square sample to a uniform disk concentric sample.
         *
         * After determining the biased sample point, the function calculates a weight for the sample. This weight is
         * influenced by the edge of the sample (how far it is from the center of the aperture), the distribution of
         * weights across the aperture, and the strength of this weighting, given by weightStrength. The final weight is
         * adjusted to ensure light conservation based on the cached integral.
         */
        inline Point2 biasSampleOnBokeh(const Point2 &sample, const Point2 &pFilm, float &weight) const
        {
            Point2 biasedSample = sample;
            float edge = 0.0f;

            // Polygonal aperture
            if (bokehConfig.blades >= 3)
            {
                // Calculate the position of the sample within the polygonal aperture shape
                float xSample = sample.x();
                int index = static_cast<int>(xSample * bokehConfig.blades);
                index = min(index, bokehConfig.blades - 1); // Ensure index is within bounds
                xSample = (xSample - (index * (1.0f / bokehConfig.blades))) * bokehConfig.blades;

                // Determine the vertices of the current polygon segment
                int previous = index + 1 == bokehConfig.blades ? 0 : index + 1;
                Point2 A(0.0f);                           // Origin
                Point2 B = bokehConfig.corners[index];    // Current vertex
                Point2 C = bokehConfig.corners[previous]; // Next vertex

                // Interpolate within the triangle formed by A, B, and C to find the biased sample
                float r1 = xSample;
                float r2 = sample.y();
                biasedSample = (1 - sqrt(r1)) * A + (sqrt(r1) * (1 - r2)) * B + (sqrt(r1) * r2) * C;
                edge = r1;
            }
            else if (abs(bokehConfig.rateOfChange) > 0.002f) // Optical Vignetting
            {
                // Convert the sample to disk coordinates, accounting for vignetting
                float rateOfChange = bokehConfig.rateOfChange, radius = bokehConfig.innerRadius;
                Point2 s = squareToUniformDiskConcentric(sample);
                float c1 = sqr(s.x()) + sqr(s.y());

                // Convert to [-1, 1], adjusting the sample based on its position on the film plane
                s.setX(s.x() - rateOfChange * pFilm.x());
                s.setY(s.y() - rateOfChange * pFilm.y());

                // Calculate the "circularity" of the adjusted sample
                float c2 = sqr(s.x()) + sqr(s.y());
                c2 *= radius; // Adjust based on the inner radius of the vignetting effect

                // If the sample is outside the valid range, we reject it by setting its weight to 0.0f
                if (c2 > 1)
                    weight = 0.0f;

                biasedSample = s;
                edge = max(c1, c2);
            }
            else // Default Circular Aperture
            {
                // Do the same as the normal thinlens camera model without Bokeh
                biasedSample = squareToUniformDiskConcentric(sample);
                edge = Vector2(biasedSample).length();
            }

            // Calculate the sample weight based on the edge distance and the bokeh weighting configuration
            float power = bokehConfig.weightStrength;
            float distr = bokehConfig.weightDistr;
            float intg = bokehConfig.integral;

            // Adjust edge and distribution parameters if the distribution is centered (distr < 0)
            if (distr < 0)
            {
                edge = 1 - edge;
                distr = -distr;
                power = -power;
            }

            // Calculate the final weight based on the power function, ensuring it is positive
            float w = power * pow(edge, distr) + 1 - intg;
            if (w > 0)
                weight = w;

            return biasedSample;
        }

      public:
        Camera(const Properties &properties)
        {
            m_resolution.x() = properties.get<int>("width");
            m_resolution.y() = properties.get<int>("height");
            m_transform = properties.getChild<Transform>();
        }

        /// @brief Returns the resolution of the image that is being rendered.
        const Vector2i &resolution() const
        {
            return m_resolution;
        }

        /**
         * @brief Helper function to sample the camera model for a given pixel.
         * This function samples a random position within the given pixel, normalizes the pixel coordinates, and
         * then calls the @c CameraSample::sample method for normalized pixel coordinates.
         *
         * @param pixel The pixel coordinates ranging from [0,0] to [resolution().x() - 1, resolution.y() - 1].
         * @param rng A random number generator used to steer the sampling.
         */
        CameraSample sample(const Point2i &pixel, Sampler &rng) const;

        /**
         * @brief Samples a ray according to this camera model in world space coordinates.
         * Sampling begins in local coordinates following the convention that [0,0,1] is the central viewing direction,
         * and then transforms the ray into world coordinates using the supplied @c m_transform object.
         *
         * @param normalized Normalized coordinates ranging from [-1, -1] (bottom left of the image) to [+1, +1] (top
         * right of the image).
         * @param rng A random number generator used to steer the sampling.
         */
        virtual CameraSample sample(const Point2 &normalized, Sampler &rng) const = 0;
    };

} // namespace lightwave
