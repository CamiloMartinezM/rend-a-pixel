#include <lightwave.hpp>

namespace lightwave
{
    /**
     * @brief A thinlens camera model with a given field of view angle, transform, radius and focal distance, based on
     * the perspective camera.
     */
    class Thinlens : public Camera
    {
      private:
        float fov, tanFov, aspectRatio, radius, focalDistance;
        std::string fovAxis;

      public:
        Thinlens(const Properties &properties) : Camera(properties)
        {
            // Perspective variables
            fovAxis = properties.get<std::string>("fovAxis");
            aspectRatio = static_cast<float>(m_resolution.x()) / static_cast<float>(m_resolution.y());
            fov = (properties.get<float>("fov") * Pi) / 180.0f; // fov converted to radians
            tanFov = tan(fov / 2.0f);

            // Depth of Field variables
            focalDistance = properties.get<float>("focalDistance");
            radius = properties.get<float>("radius");

            // Initialize Bokeh effect-related variables from properties, if it is activated by the global variable
            if (UseBokehEffects)
            {
                int edges = properties.get<int>("lenselements", 0);
                bokehConfig.blades = edges;

                // Generate the shapes for the bokeh effect based on the number of blades
                buildBokehShapes(bokehConfig);

                bokehConfig.innerRadius = properties.get<float>("innerRadius", 1.0);
                bokehConfig.rateOfChange = properties.get<float>("rateOfChange", 0.0f);
                bokehConfig.weightDistr = properties.get<float>("weightDistr", 2.0);
                bokehConfig.weightStrength = properties.get<float>("weightStrength", 0.0);

                // Validate the weight strength to ensure light conservation
                if (bokehConfig.weightStrength < 0 || bokehConfig.weightStrength > 1)
                {
                    logger(EError, "For weightStrength values above 1 or below 0 light conservation is violated, using "
                                   "1.0f instead");
                    bokehConfig.weightStrength = 1.0f;
                }

                // Adjust weight strength based on the distribution
                bokehConfig.weightStrength = bokehConfig.weightStrength * (bokehConfig.weightDistr + 1);

                // This variable accounts for the modifications introduced by the weightDistr and weightStrength
                // variables, ensuring that the bokeh effect does not artificially brighten or darken the image
                bokehConfig.integral = bokehConfig.weightStrength / (bokehConfig.weightDistr + 1);
            }
        }

        CameraSample sample(const Point2 &normalized, Sampler &rng) const override
        {
            // Compute the position on the film plane in camera space, based on the normalized coordinates
            float filmX = normalized.x() * tanFov * (fovAxis == "x" ? 1.0f : aspectRatio);
            float filmY = normalized.y() * tanFov * (fovAxis == "x" ? 1.0f / aspectRatio : 1.0f);
            Point2 pFilm(filmX, filmY);

            // Create the ray starting at the camera origin and passing through the film plane position
            Vector rayDirection(filmX, filmY, 1.0f);
            Point rayOrigin(0.0f);
            float weight = 1.0f;

            // If the lens radius is greater than zero, modify the ray for depth of field
            if (radius > 0)
            {
                // Sample a point on the lens, calling the appropiate function depending on it Bokeh effects are on
                Point2 pLens;
                if (UseBokehEffects)
                    pLens = radius * biasSampleOnBokeh(rng.next2D(), pFilm, weight);
                else
                    pLens = radius * squareToUniformDiskConcentric(rng.next2D());

                // Compute the focal plane intersection point
                float ft = focalDistance / rayDirection.z();
                Point focalPoint = rayDirection * ft;

                // Update the ray origin and direction for depth of field effect
                rayOrigin = Point(pLens.x(), pLens.y(), 0.0f);
                rayDirection = (focalPoint - rayOrigin).normalized();
            }

            // Transform the ray to world space
            Ray worldRay = m_transform->apply(Ray(rayOrigin, rayDirection)).normalized();

            return CameraSample{.ray = worldRay, .weight = Color(weight)};
        }

        std::string toString() const override
        {
            return tfm::format("Thinlens[\n"
                               "  width = %d,\n"
                               "  height = %d,\n"
                               "  transform = %s,\n"
                               "  radius = %s,\n"
                               "  focalDistance = %s,\n"
                               "]",
                               m_resolution.x(), m_resolution.y(), indent(m_transform), radius, focalDistance);
        }
    };
} // namespace lightwave

REGISTER_CAMERA(Thinlens, "thinlens")
