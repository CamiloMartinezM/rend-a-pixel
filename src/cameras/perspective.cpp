#include <lightwave.hpp>

namespace lightwave
{
    /**
     * @brief A perspective camera with a given field of view angle and transform.
     *
     * In local coordinates (before applying m_transform), the camera looks in positive z direction [0,0,1].
     * Pixels on the left side of the image ( @code normalized.x < 0 @endcode ) are directed in negative x
     * direction ( @code ray.direction.x < 0 ), and pixels at the bottom of the image ( @code normalized.y < 0 @endcode
     * ) are directed in negative y direction ( @code ray.direction.y < 0 ).
     */
    class Perspective : public Camera
    {
      private:
        float fov, tanFov, aspectRatio;
        std::string fovAxis;

      public:
        Perspective(const Properties &properties) : Camera(properties)
        {
            // hints:
            // * precompute any expensive operations here (most importantly trigonometric functions)
            // * use m_resolution to find the aspect ratio of the image
            fovAxis = properties.get<std::string>("fovAxis");
            aspectRatio = static_cast<float>(m_resolution.x()) / static_cast<float>(m_resolution.y());
            fov = (properties.get<float>("fov") * Pi) / 180.0f; // fov converted to radians
            tanFov = tan(fov / 2.0f);
        }

        CameraSample sample(const Point2 &normalized, Sampler &rng) const override
        {
            // Compute the position on the film plane in camera space, based on the normalized coordinates
            float filmX = normalized.x() * tanFov * (fovAxis == "x" ? 1.0f : aspectRatio);
            float filmY = normalized.y() * tanFov * (fovAxis == "x" ? 1.0f / aspectRatio : 1.0f);

            Vector rayDirection = Vector(filmX, filmY, 1.0f);

            // Transforming the rays to world coordinates via the camera transform specified in the scene description
            // file (stored in m_transform, apply does local to world.
            Ray worldRay = m_transform->apply(Ray(Vector(0.0f), rayDirection)).normalized();

            return CameraSample{.ray = worldRay, .weight = Color(1.0f)};
        }

        std::string toString() const override
        {
            return tfm::format("Perspective[\n"
                               "  width = %d,\n"
                               "  height = %d,\n"
                               "  transform = %s,\n"
                               "]",
                               m_resolution.x(), m_resolution.y(), indent(m_transform));
        }
    };
} // namespace lightwave

REGISTER_CAMERA(Perspective, "perspective")
