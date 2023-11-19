#include <lightwave.hpp>

namespace lightwave {

    // Declaring the solveQuadratic function for visibility within the namespace
    static bool solveQuadratic(const float& a, const float& b, const float& c, float& x0, float& x1);

    /// @brief A sqhere in the R3-plane with center (0, 0, 0) and radius = 1
    class Sphere : public Shape {

      /**
       * @brief Constructs a surface event for a given position, used by @ref
       * intersect to populate the @ref Intersection and by @ref sampleArea to
       * populate the @ref AreaSample .
       * @param surf The surface event to populate with texture coordinates, shading
       * frame and area pdf
       * @param position The hitpoint (i.e., point in [-radius,-radius,-radius] to [+radius,+radius,+radius]), found
       * via intersection or area sampling
       */
        inline void populate(SurfaceEvent& surf, const Point& position) const {
            surf.position = position;

            // Map the position onto the sphere's surface
            // surf.uv.x() = 0.5 + atan2(position.y(), position.x()) / (2 * M_PI);
        
            // surf.uv.y() = 0.5 - asin(position.z()) / M_PI;
            surf.uv.x() = (position.x() + radius) / (2 * radius);
            surf.uv.y() = (position.y() + radius) / (2 * radius);

            // surf.uv.x() = position.x();
            // surf.uv.y() = position.y();

            // normal always points in the direction of p - o
            surf.frame.normal = (position - center_point).normalized();
            // the tangent always points in positive x direction
            surf.frame.tangent = radius * surf.frame.normal.cross(Vector(0.f, 0.f, 1.f));
            // the bitagent always points in positive z direction
            surf.frame.bitangent = surf.frame.normal.cross(surf.frame.tangent);

            surf.pdf = 0.f;
        }

        private:

        float radius;
        Point center_point;

        public:

        Sphere(const Properties& properties) {
            radius = 1;
            center_point = Point(0.f, 0.f, 0.f);
        }

        bool intersect(const Ray& ray, Intersection& its,
                       Sampler& rng) const override
        {
            float t0, t1; // Possible intersection points (quadratic roots)

            // Defines the a, b, c values of a quadratic equation of the form: a*x² + b*x + c = 0
            // Where the equation for intersecting a line and a sphere is given by:
            // d²*t² + 2*(o - c)*d*t + ((o - c)² - r²) = 0
            // Where a point p in the line is defined by: p = o + t*d
            // and a point p in the sphere is (p - c)² - r² = 0
            Vector L = ray.origin - center_point;
            float a = ray.direction.dot(ray.direction);
            float b = 2 * ray.direction.dot(L);
            float c = L.dot(L) - radius * radius;

            if (!solveQuadratic(a, b, c, t0, t1)) return false; // if no solution, then no intersection

            if (t0 < 0) {
                // if t0 is negative, try t1
                if (t1 < 0) return false; // both t0 and t1 are negative, hence no intersection
                t0 = t1; // if t0 > 0, then stay with t1 and discard t0
            }

            // note that we never report an intersection closer than Epsilon (to avoid self-intersections)!
            // we also do not update the intersection if a closer intersection already exists (i.e., its.t is 
            // lower than our own t)
            if (t0 < Epsilon || t0 > its.t)
                return false;

            // compute the hitpoint
            const Point hit_position = ray(t0);

            // dismiss anything outside of the [-radius,-radius,-radius]..[+radius,+radius,+radius] domain.
            if (!getBoundingBox().includes(hit_position))
                return false;

            // we have determined there was an intersection! 
            // we are now free to change the intersection object and return true.
            its.t = t0;
            populate(its, hit_position); // compute the shading frame, texture coordinates and area pdf (same as sampleArea)
            return true;
        }

        Bounds getBoundingBox() const override {
            return Bounds(Point{ -radius, -radius, -radius }, Point{ +radius, +radius, +radius });
        }

        Point getCentroid() const override {
            return center_point;
        }

        AreaSample sampleArea(Sampler& rng) const override {
            NOT_IMPLEMENTED
        }

        std::string toString() const override {
            return "Sphere[]";
        }
    };

    // Obtained from: https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-sphere-intersection.html
    // For better computing precision
    bool solveQuadratic(const float& a, const float& b, const float& c, float& x0, float& x1)
    {
        float discr = b * b - 4 * a * c;
        if (discr < 0) return false;
        else if (discr == 0) x0 = x1 = -0.5 * b / a;
        else {
            float q = (b > 0) ?
                -0.5 * (b + sqrt(discr)) :
                -0.5 * (b - sqrt(discr));
            x0 = q / a;
            x1 = c / q;
        }
        if (x0 > x1) std::swap(x0, x1);

        return true;
    }

} // namespace lightwave

REGISTER_SHAPE(Sphere, "sphere")