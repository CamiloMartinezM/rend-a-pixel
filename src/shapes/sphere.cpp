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
        /// @brief Sphere sampling routine
        const ShapeSamplingMethod SphereSampling = ShapeSamplingMethod::SubtendedCone;

        /// @brief Converts a position to a uv coordinate.
        /// @param position The normalized hitpoint as a Vector from the center of the sphere. 
        inline Point2 sphereToUV(const Vector &position) const
        {
            // Calculate spherical coordinates
            float theta = acos(position.y());         // inclination
            float phi = atan2(position.x(), position.z()); // azimuth

            // Map the spherical coordinates to UV coordinates
            // U coordinate: phi mapped from [0, 2*PI] to [0, 1]
            // V coordinate: theta mapped from [0, PI] to [0, 1]
            return Point2(phi * Inv2Pi, theta * InvPi);
        }

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

            // Calculate the uv coordinates based on the hit position
            surf.uv = sphereToUV(p_o);

            // Adjust phi to be in the range [0, 2*PI]
            if (surf.uv.x() < 0.0f)
                surf.uv.x() += 1.0f;

            // Construct the shading frame
            // The tangent is perpendicular to the normal and the up vector (0, 1, 0)
            surf.frame = Frame(surf.frame.normal);
        }

      private:
        Point centerPoint;

        /**
         * @brief Projects a sampled vector on the surface of the sphere and scales by the sphere radius
         * @param v A previously sampled vector.
         * @return A 3D point on the surface of the sphere.
         */
        inline Point projectBackOnSphere(const Vector &v) const
        {
            // Point pObj = centerPoint + v;
            // return Point(Vector(pObj) * radius / (pObj - centerPoint).length());
            // Result is the following after simplification, because radius = 1, centerPoint = (0, 0, 0)
            return v.normalized();
        }

        /// @brief Populates an AreaSample with the given position.
        AreaSample populateAreaSampleWithPosition(const Point &position) const
        {
            AreaSample sample;
            populate(sample, position);
            return sample;
        }

        /// @brief Updates the PDF of the given Surface Event with the Uniform Sphere PDF
        inline void updatePdfUniform(SurfaceEvent &surf) const
        {
            surf.pdf = uniformSpherePdf();
        }

        /// @brief Updates the PDF of the sphere based on the sampling method used (uniform, cosine-weighted,
        /// subtended cone).
        inline void updatePdf(SurfaceEvent &surf, const ShapeSamplingMethod &usedSamplingMethod,
                              const Vector &sampledVector, const Point &refPoint) const
        {
            if (usedSamplingMethod == ShapeSamplingMethod::Uniform)
                updatePdfUniform(surf);
            else if (usedSamplingMethod == ShapeSamplingMethod::CosineWeighted)
                surf.pdf = cosineHemispherePdf(sampledVector);
            else
                surf.pdf = pdfFromSolidAngleMeasure(subtendedConePdf(centerPoint, 1.0f, refPoint),
                                                    (surf.position - refPoint).length(), surf.frame.normal,
                                                    surf.position - refPoint);
        }

      public:
        Sphere(const Properties &properties)
        {
            centerPoint = Point(0.0f);
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
            float c = L.dot(L) - 1; // - sqr(radius), but radius=1

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
            const Point hitPosition = ray(t0);

            // perform the alpha masking test
            if (its.alphaMask)
            {
                // convert hit position to UV coordinates
                Point2 uv = sphereToUV((hitPosition - centerPoint).normalized());             
                float alphaValue = its.alphaMask->scalar(uv);

                // stochastically dismiss the intersection based on alpha value
                if (alphaValue < rng.next())
                    return false; 
            }

            // we have determined there was an intersection!
            // we are now free to change the intersection object and return true.
            its.t = t0;
            populate(its, hitPosition); // compute the shading frame and texture coordinates

            // update the pdf of the intersection based on the sphere sampling
            updatePdf(its, SphereSampling, (hitPosition - centerPoint).normalized(), ray.origin);
            return true;
        }

        Bounds getBoundingBox() const override
        {
            return Bounds(Point(-1.0f), Point(+1.0f));
        }

        Point getCentroid() const override
        {
            return centerPoint;
        }

        AreaSample sampleArea(Sampler &rng) const override
        {
            // Go for Uniform Sphere Sampling as default
            Vector sampledVector = squareToUniformSphere(rng.next2D());

            // Project the sampled point back on the sphere and populate an AreaSample with it
            Point projectedPoint = projectBackOnSphere(sampledVector);
            AreaSample sampledArea = populateAreaSampleWithPosition(projectedPoint);

            // For Uniform Sphere Sampling, sampledArea.pdf = 1/(4 * Pi * r^2)
            updatePdfUniform(sampledArea);
            return sampledArea;
        }

        AreaSample sampleArea(Sampler &rng, const SurfaceEvent &ref) const override
        {
            // Use default behaviour if the sphere sampling method is to be uniform
            if (SphereSampling == ShapeSamplingMethod::Uniform)
                return sampleArea(rng);

            Vector sampledVector, importanceSampledVector;
            ShapeSamplingMethod usedSamplingMethod;
            Point projectedPoint;
            if (SphereSampling == ShapeSamplingMethod::CosineWeighted)
            {
                // Otherwise, sample sphere with cosine-weighted sampling
                usedSamplingMethod = ShapeSamplingMethod::CosineWeighted;

                // Get a vector on the hemisphere
                sampledVector = squareToCosineHemisphere(rng.next2D());

                // Build a normal from the center of the sphere towards the origin (points outside of the surface)
                Vector importantDirection = (ref.position - centerPoint).normalized();

                // Build a frame out of that normal and convert the sampled vector to this coordinate system that
                // points towards the intersection point from the sphere
                Frame importanceSamplingFrame = Frame(importantDirection);
                importanceSampledVector = importanceSamplingFrame.toWorld(sampledVector);
            }
            else
            {
                // Sample uniformly on sphere if ref.position is inside it
                Point pOrigin =
                    offsetRayOrigin(ref.position, ref.frame.normal, (ref.position - centerPoint).normalized());
                if ((pOrigin - centerPoint).lengthSquared() <= sqr(1.0f + Epsilon))
                {
                    importanceSampledVector = squareToUniformSphere(rng.next2D());
                    usedSamplingMethod = ShapeSamplingMethod::Uniform;
                }
                else // Otherwise, sample sphere inside subtended cone
                {
                    usedSamplingMethod = ShapeSamplingMethod::SubtendedCone;
                    importanceSampledVector =
                        subtendedConeSphereSampling4ed(rng.next2D(), centerPoint, 1.0f, ref.position);
                }
            }

            // Project the sampled point back on the sphere and populate an AreaSample with it
            projectedPoint = projectBackOnSphere(importanceSampledVector);
            AreaSample sampledArea = populateAreaSampleWithPosition(projectedPoint);

            // For Uniform, only sampleArea is used; for cosine-weighted, sampleArea and sampledVector are used; and
            // for subtended cone, all parameters are used to calculate the PDF
            updatePdf(sampledArea, usedSamplingMethod, sampledVector, ref.position);
            return sampledArea;
        }

        std::string toString() const override
        {
            return "Sphere[]";
        }
    };
} // namespace lightwave

REGISTER_SHAPE(Sphere, "sphere")