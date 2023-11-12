#include <lightwave.hpp>
#include <math.h>

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
    private:
        double fov, tan_fov, aspect_ratio;
        
    public:
        Perspective(const Properties &properties)
            : Camera(properties)
        {
            // NOT_IMPLEMENTED

            // hints:
            // * precompute any expensive operations here (most importantly trigonometric functions)
            // * use m_resolution to find the aspect ratio of the image
            aspect_ratio = static_cast<double>(m_resolution.x()) / static_cast<double>(m_resolution.y());
            fov = (properties.get<float>("fov") / 180.0) * M_PI; // fov converted to radians
            tan_fov = tan(fov / 2.0);
        }

        CameraSample sample(const Point2 &normalized, Sampler &rng) const override
        {
            // first transforming the normalized coordinates to the local camera coordinate
            double lc_coord_x = normalized.x() * tan_fov * aspect_ratio;
            double lc_coord_y = normalized.y() * tan_fov;
            Vector lc_coord_point =  Vector(lc_coord_x, lc_coord_y, 1.0f);

            // then transforming the rays to world coordinates via the camera transform specified in
            // the scene description file (stored in m_transform).
            // hints:
            // * use m_transform to transform the local camera coordinate system into the world coordinate system
            Ray world_coord_ray = m_transform->apply(Ray(Point(0.f, 0.f, 0.f), lc_coord_point.normalized()));
            
            return CameraSample{.ray = world_coord_ray, .weight = Color(1.0f)};
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
