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
            fovAxis = properties.get<std::string>("fovAxis");
            radius = properties.get<float>("radius");
            focalDistance = properties.get<float>("focalDistance");

            aspectRatio = static_cast<float>(m_resolution.x()) / static_cast<float>(m_resolution.y());
            fov = (properties.get<float>("fov") * Pi) / 180.0f; // fov converted to radians
            tanFov = tan(fov / 2.0f);
        }

        CameraSample sample(const Point2 &normalized, Sampler &rng) const override
        {
            // Compute the position on the film plane in camera space, based on the normalized coordinates
            float filmX = normalized.x() * tanFov * (fovAxis == "x" ? 1.0f : aspectRatio);
            float filmY = normalized.y() * tanFov * (fovAxis == "x" ? 1.0f / aspectRatio : 1.0f);

            // Create the ray starting at the camera origin and passing through the film plane position
            Vector rayDirection(filmX, filmY, 1.0f);
            Point rayOrigin(0.0f);

            // If the lens radius is greater than zero, modify the ray for depth of field
            if (radius > 0)
            {
                // Sample a point on the lens
                Point2 p = squareToUniformDiskConcentric(rng.next2D());
                Point2 pLens = Point2(radius * p.x(), radius * p.y());

                // Compute the focal plane intersection point
                float ft = focalDistance / rayDirection.z();
                Point focalPoint = rayDirection * ft;

                // Update the ray origin and direction for depth of field effect
                rayOrigin = Point(pLens.x(), pLens.y(), 0.0f);
                rayDirection = (focalPoint - rayOrigin).normalized();
            }

            // Transform the ray to world space
            Ray worldRay = m_transform->apply(Ray(rayOrigin, rayDirection)).normalized();

            return CameraSample{.ray = worldRay, .weight = Color(1.0f)};
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
