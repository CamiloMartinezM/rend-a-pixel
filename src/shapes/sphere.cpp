#include <lightwave.hpp>

namespace lightwave
{
    /// @brief Obtained from:
    /// https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-sphere-intersection.html
    /// For better computing precision
    static bool solveQuadratic(const float &a, const float &b, const float &c, float &x0, float &x1)
    {
        float discr = b * b - 4 * a * c;
        if (discr < 0)
            return false;
        else if (discr == 0)
            x0 = x1 = -0.5 * b / a;
        else
        {
            float q = (b > 0) ? -0.5 * (b + sqrt(discr)) : -0.5 * (b - sqrt(discr));
            x0 = q / a;
            x1 = c / q;
        }
        if (x0 > x1)
            std::swap(x0, x1);

        return true;
    }

    /// @brief A sqhere in the R3-plane with center (0, 0, 0) and radius = 1
    class Sphere : public Shape
    {
        /**
         * @brief Constructs a surface event for a given position, used by @ref
         * intersect to populate the @ref Intersection and by @ref sampleArea to
         * populate the @ref AreaSample .
         * @param surf The surface event to populate with texture coordinates, shading
         * frame and area pdf
         * @param position The hitpoint (i.e., point in [-radius,-radius,-radius] to [+radius,+radius,+radius]), found
         * via intersection or area sampling
         */
        inline void populate(SurfaceEvent &surf, const Point &position) const
        {
            // For assignment_1 the following worked to pass the tests:
            //
            // Map the position onto the sphere's surface
            // surf.uv.x() = (position.x() + radius) / (2 * radius);
            // surf.uv.y() = (position.z() + radius) / (2 * radius);
            // // normal always points in the direction of p - o
            // surf.frame.normal = (position - center_point).normalized();
            // // the tangent always points in positive x direction
            // surf.frame.tangent = radius * surf.frame.normal.cross(Vector(0.f, 0.f, 1.f)).normalized();
            // // the bitagent always points in positive z direction
            // surf.frame.bitangent = surf.frame.normal.cross(surf.frame.tangent).normalized();
            // surf.pdf = 0.f;
            // For assignment_2 the following was made:

            // Normal always points in the direction of p - o
            Vector p_o = (position - centerPoint).normalized();

            surf.position = Point(p_o); // Hit point is normalized

            surf.frame.normal = p_o;

            // Calculate spherical coordinates
            float theta = std::acos(p_o.y());         // inclination
            float phi = std::atan2(p_o.x(), p_o.z()); // azimuth

            // Map the spherical coordinates to UV coordinates
            // U coordinate: phi mapped from [0, 2*PI] to [0, 1]
            // V coordinate: theta mapped from [0, PI] to [0, 1]
            surf.uv.x() = phi * Inv2Pi;
            surf.uv.y() = theta * InvPi;

            // Adjust phi to be in the range [0, 2*PI]
            if (surf.uv.x() < 0)
                surf.uv.x() += 1.0f;

            // Construct the shading frame
            // The tangent is perpendicular to the normal and the up vector (0, 1, 0)
            surf.frame.tangent = Vector(0.f, 1.f, 0.f).cross(surf.frame.normal).normalized();
            // The bitangent is perpendicular to both the normal and the tangent
            surf.frame.bitangent = surf.frame.normal.cross(surf.frame.tangent);
            // Uniform PDF over the sphere's surface
            surf.pdf = Inv4Pi / (radius * radius);
        }

      private:
        float radius;
        Point centerPoint;

        /**
         * @brief Samples a point on the surface of the sphere uniformly using the squareToUniformSphere method.
         * @param sample A 2D point in the unit square.
         * @return A 3D point on the surface of the sphere corresponding to the input sample.
         */
        inline Point SphereUniformSample(const Point2 &sample) const 
        {
            return centerPoint + radius * squareToUniformSphere(sample);
        }

        /// @brief Populates an AreaSample with the given position.
        AreaSample populateAreaSampleWithPosition(const Point &position) const  
        {
            AreaSample sample;
            populate(sample, position);
            return sample;
        }

      public:
        Sphere(const Properties &properties)
        {
            radius = 1;
            centerPoint = Point(0.f, 0.f, 0.f);
        }

        bool intersect(const Ray &ray, Intersection &its, Sampler &rng) const override
        {
            float t0, t1; // Possible intersection points (quadratic roots)

            // Defines the a, b, c values of a quadratic equation of the form: a*x² + b*x + c = 0
            // Where the equation for intersecting a line and a sphere is given by:
            // d²*t² + 2*(o - c)*d*t + ((o - c)² - r²) = 0
            // Where a point p in the line is defined by: p = o + t*d
            // and a point p in the sphere is (p - c)² - r² = 0
            Vector L = ray.origin - centerPoint;
            float a = ray.direction.dot(ray.direction);
            float b = 2 * ray.direction.dot(L);
            float c = L.dot(L) - radius * radius;

            if (!solveQuadratic(a, b, c, t0, t1))
                return false; // if no solution, then no intersection

            if (t0 < Epsilon)
            {
                // if t0 is less than Epsilon, try t1
                if (t1 < Epsilon)
                    return false; // both t0 and t1 are less than Epsilon, hence no intersection
                t0 = t1;          // if t1 > Epsilon, then stay with t1 and discard t0
            }

            // until now we never report an intersection closer than Epsilon (to avoid self-intersections)!
            // we also do not update the intersection if a closer intersection already exists (i.e., its.t is
            // lower than our own t)
            if (t0 > its.t)
                return false;

            // compute the hitpoint
            const Point hit_position = ray(t0);

            // we have determined there was an intersection!
            // we are now free to change the intersection object and return true.
            its.t = t0;
            populate(its,
                     hit_position); // compute the shading frame, texture coordinates and area pdf (same as sampleArea)
            return true;
        }

        Bounds getBoundingBox() const override
        {
            return Bounds(Point{-radius, -radius, -radius}, Point{+radius, +radius, +radius});
        }

        Point getCentroid() const override
        {
            return centerPoint;
        }

        AreaSample sampleArea(Sampler &rng) const override
        {
            // Sample a point on the unit sphere using spherical coordinates, transform the point to the sphere's scale
            // and position, and populate the AreaSample
            return populateAreaSampleWithPosition(SphereUniformSample(rng.next2D()));
        }

        AreaSample sampleArea(Sampler &rng, const Intersection &ref) const override
        {
            // Sample uniformly on sphere if ref.position is inside it
            Point pOrigin = OffsetRayOrigin(ref.position, ref.wo, centerPoint - ref.position);
            if ((pOrigin - centerPoint).lengthSquared() <= radius * radius)
                return sampleArea(rng);

            // Otherwise, sample sphere uniformly inside subtended cone
            Point sampledPoint = subtendedConeUniformSphereSampling(rng.next2D(), centerPoint, radius, ref.position);
            return populateAreaSampleWithPosition(sampledPoint);
        }

        std::string toString() const override
        {
            return "Sphere[]";
        }
    };
} // namespace lightwave

REGISTER_SHAPE(Sphere, "sphere")