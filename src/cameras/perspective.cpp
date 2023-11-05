#include <lightwave.hpp>

namespace lightwave
{

    /**
     * @brief A perspective camera with a given field of view angle and transform.
     *
     * In local coordinates (before applying m_transform), the camera looks in positive z direction [0,0,1].
     * Pixels on the left side of the image ( @code normalized.x < 0 @endcode ) are directed in negative x
     * direction ( @code ray.direction.x < 0 ), and pixels at the bottom of the image ( @code normalized.y < 0 @endcode )
     * are directed in negative y direction ( @code ray.direction.y < 0 ).
     */
    class Perspective : public Camera
    {
    public:
        Perspective(const Properties &properties)
            : Camera(properties)
        {
            // NOT_IMPLEMENTED

            // hints:
            // * precompute any expensive operations here (most importantly trigonometric functions)
            // * use m_resolution to find the aspect ratio of the image
        }

        CameraSample sample(const Point2 &normalized, Sampler &rng) const override
        {
            // NOT_IMPLEMENTED
            return CameraSample{.ray = Ray(Vector(normalized.x(), normalized.x(), 0.f), Vector(0.f, 0.f, 1.f)), .weight = Color(1.0f)};

            // hints:
            // * use m_transform to transform the local camera coordinate system into the world coordinate system
        }

        std::string toString() const override
        {
            return tfm::format(
                "Perspective[\n"
                "  width = %d,\n"
                "  height = %d,\n"
                "  transform = %s,\n"
                "]",
                m_resolution.x(),
                m_resolution.y(),
                indent(m_transform));
        }
    };

}

REGISTER_CAMERA(Perspective, "perspective")
