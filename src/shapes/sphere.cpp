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

            // Uniform area pdf. updatePdf() overwrites this if it is called after this method
            surf.pdf = Inv4Pi / sqr(radius); 
        }

      private:
        float radius;
        Point centerPoint;

        /**
         * @brief Projects a sampled vector on the surface of the sphere and scales by the sphere radius
         * @param v A previously sampled vector.
         * @return A 3D point on the surface of the sphere.
         */
        inline Point projectBackOnSphere(const Vector &v) const
        {
            Point pObj = centerPoint + radius * v;
            return Point(Vector(pObj) * radius / (pObj - centerPoint).length()); // Scale pObj by the sphere’s radius
        }

        /// @brief Populates an AreaSample with the given position.
        AreaSample populateAreaSampleWithPosition(const Point &position) const
        {
            AreaSample sample;
            populate(sample, position);
            return sample;
        }

        /// @brief Updates the PDF of the sphere based on the sampling method used (uniform, cosine-weighted,
        /// subtended cone).
        inline void updatePdf(SurfaceEvent &surf, const ShapeSamplingMethod usedSamplingMethod,
                              const Vector sampledVector, const Point projectedPoint, const Intersection ref,
                              const bool adjustUniform) const
        {
            if (usedSamplingMethod == ShapeSamplingMethod::Uniform)
            {
                surf.pdf = uniformSpherePdf();
                // if (adjustUniform)
                // {
                //     // Convert area sampling PDF in ss to solid angle measure
                //     Vector n = Vector(projectedPoint).normalized();
                //     surf.pdf /= n.dot(ref.wo) / (ref.position - projectedPoint).lengthSquared();
                //     if (std::isinf(surf.pdf))
                //         surf.pdf = 0;
                // }
            }
            else if (usedSamplingMethod == ShapeSamplingMethod::CosineWeighted)
                surf.pdf = cosineSpherePdf(sampledVector);
            else
                surf.pdf = subtendedConePdf(centerPoint, radius, ref.position);
            surf.pdf /= sqr(radius);
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
            float c = L.dot(L) - sqr(radius);

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
            Vector sampledVector;
            if (SphereSampling == ShapeSamplingMethod::Uniform)
                sampledVector = squareToUniformSphere(rng.next2D());
            else if (SphereSampling == ShapeSamplingMethod::CosineWeighted)
                sampledVector = squareToCosineSphere(rng.next2D());

            // Project the sampled point back on the sphere and populate an AreaSample with it
            Point projectedPoint = projectBackOnSphere(sampledVector);
            AreaSample sampledArea = populateAreaSampleWithPosition(projectedPoint);
            updatePdf(sampledArea, SphereSampling, sampledVector, projectedPoint, Intersection(), false);
            return sampledArea;
        }

        AreaSample sampleArea(Sampler &rng, const Intersection &ref) const override
        {
            // Use default behaviour if the sphere sampling method is to be uniform or cosine-weighted
            if (SphereSampling == ShapeSamplingMethod::Uniform || SphereSampling == ShapeSamplingMethod::CosineWeighted)
                return sampleArea(rng);

            Vector sampledVector;
            ShapeSamplingMethod usedSamplingMethod;
            bool adjustUniform = false;

            // Sample uniformly on sphere if ref.position is inside it
            Point pOrigin = OffsetRayOrigin(ref.position, ref.wo, centerPoint - ref.position);
            if ((pOrigin - centerPoint).lengthSquared() <= radius * radius)
            {
                sampledVector = squareToUniformSphere(rng.next2D());
                usedSamplingMethod = ShapeSamplingMethod::Uniform;
                adjustUniform = true;
            }
            else // Otherwise, sample sphere uniformly inside subtended cone
            {
                sampledVector = subtendedConeSphereSampling4ed(rng.next2D(), centerPoint, radius, ref.position);
                usedSamplingMethod = ShapeSamplingMethod::SubtendedCone;
            }

            // Project the sampled point back on the sphere and populate an AreaSample with it
            Point projectedPoint = projectBackOnSphere(sampledVector);
            AreaSample sampledArea = populateAreaSampleWithPosition(projectedPoint);
            updatePdf(sampledArea, usedSamplingMethod, sampledVector, projectedPoint, ref, adjustUniform);
            return sampledArea;
        }

        std::string toString() const override
        {
            return "Sphere[]";
        }
    };
} // namespace lightwave

REGISTER_SHAPE(Sphere, "sphere")