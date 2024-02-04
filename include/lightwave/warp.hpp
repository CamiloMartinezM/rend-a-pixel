/**
 * @file warp.hpp
 * @brief Contains functions that map one domain to another.
 */

#pragma once

#include <lightwave/math.hpp>

namespace lightwave
{

    /**
     * @brief Warps a given point from the unit square ([0,0] to [1,1]) to a unit circle (centered around [0,0] with
     * radius 1), with uniform density given by @code 1 / Pi @endcode .
     * @see Based on http://psgraphics.blogspot.ch/2011/01/improved-code-for-concentric-map.html
     */
    inline Point2 squareToUniformDiskConcentric(const Point2 &sample)
    {
        float r1 = 2 * sample.x() - 1;
        float r2 = 2 * sample.y() - 1;

        float phi, r;
        if (r1 == 0 && r2 == 0)
        {
            r = 0;
            phi = 0;
        }
        else if (r1 * r1 > r2 * r2)
        {
            r = r1;
            phi = Pi4 * (r2 / r1);
        }
        else
        {
            r = r2;
            phi = Pi2 - Pi4 * (r1 / r2);
        }

        float cosPhi = std::cos(phi);
        float sinPhi = std::sin(phi);
        return {r * cosPhi, r * sinPhi};
    }

    /**
     * @brief Warps a given point from the unit square ([0,0] to [1,1]) to a unit sphere (centered around [0,0,0] with
     * radius 1), with uniform density given by @code 1 / (4 * Pi) @endcode .
     */
    inline Vector squareToUniformSphere(const Point2 &sample)
    {
        float z = 1 - 2 * sample.y();
        float r = safe_sqrt(1 - z * z);
        float phi = 2 * Pi * sample.x();
        float cosPhi = cos(phi);
        float sinPhi = sin(phi);
        return {r * cosPhi, r * sinPhi, z};
    }

    /**
     * @brief Warps a given point from the unit square ([0,0] to [1,1]) to a unit hemisphere (centered around [0,0,0]
     * with radius 1, pointing in z direction), with respect to solid angle.
     */
    inline Vector squareToUniformHemisphere(const Point2 &sample)
    {
        Point2 p = squareToUniformDiskConcentric(sample);
        float z = 1.0f - p.x() * p.x() - p.y() * p.y();
        float s = sqrt(z + 1.0f);
        return {s * p.x(), s * p.y(), z};
    }

    /// @brief Returns the density of the @ref squareToUniformHemisphere warping.
    inline float uniformHemispherePdf()
    {
        return Inv2Pi;
    }

    /// @brief Returns the density of the @ref squareToUniformSphere warping.
    inline float uniformSpherePdf()
    {
        return Inv4Pi;
    }

    /**
     * @brief Warps a given point from the unit square ([0,0] to [1,1]) to a unit hemisphere (centered around [0,0,0]
     * with radius 1, pointing in z direction), with density given by @code cos( angle( result, [0,0,1] ) ) @endcode .
     */
    inline Vector squareToCosineHemisphere(const Point2 &sample)
    {
        Point2 p = squareToUniformDiskConcentric(sample);
        float z = safe_sqrt(1.0f - p.x() * p.x() - p.y() * p.y());
        return {p.x(), p.y(), z};
    }

    /**
     * @brief Warps a given point from the unit square ([0,0] to [1,1]) to a unit sphere (centered around [0,0,0] with
     * radius 1, pointing in z direction). It uses the same method in squareToCosineHemisphere(), and a random variable,
     * such as sample.x() to decide whether to flip the z axis, effectively sampling both hemispheres.
     */
    inline Vector squareToCosineSphere(const Point2 &sample)
    {
        // First, perform cosine-weighted sampling on a hemisphere
        Point2 p = squareToUniformDiskConcentric(sample);
        float z = safe_sqrt(1.0f - p.x() * p.x() - p.y() * p.y());

        // Randomly decide whether to use the upper or lower hemisphere
        if (sample.x() < 0.5)
            return {p.x(), p.y(), z}; // Use the upper hemisphere
        else
            return {p.x(), p.y(), -z}; // Use the lower hemisphere (flip the z-coordinate)
    }

    /// @brief Returns the density of the @ref squareToCosineHemisphere warping.
    inline float cosineHemispherePdf(const Vector &vector)
    {
        return InvPi * std::max(vector.z(), float(0));
    }

    /// @brief Returns the density of the @ref squareToCosineSphere warping.
    inline float cosineSpherePdf(const Vector &vector)
    {
        return 2 * InvPi * std::max(vector.z(), float(0));
    }

    /**
     * @brief Returns the Pdf value in solid angle measure, from a value in spatial units.
     * @param pdf Spatial Pdf.
     * @param distance Distance between the shading point and the reference point.
     * @param n Surface normal.
     * @param wi The direction towards the reference point (e.g, a light source). It does not have to be normalized.
     */
    inline float pdfToSolidAngleMeasure(const float &pdf, const float &distance, const Vector &n, const Vector &wi)
    {
        Vector wiNormalized = wi.normalized();
        if (abs(n.dot(-wiNormalized)) <= MachineEpsilon)
            return Infinity;

        float result = pdf * sqr(distance) / abs(n.dot(-wiNormalized));
        return clamp(result, MachineEpsilon, Infinity); 
    }

    /**
     * @brief Returns the Pdf value in spatial units, from a value in solid angle measure.
     * @param pdf Pdf in solid angle measure.
     * @param distance Distance between the shading point and the reference point.
     * @param n Surface normal.
     * @param wi The direction towards the reference point (e.g, a light source). It does not have to be normalized.
     */
    inline float pdfFromSolidAngleMeasure(const float &pdf, const float &distance, const Vector &n, const Vector &wi)
    {
        Vector wiNormalized = wi.normalized();

        if (sqr(distance) <= MachineEpsilon)
            return Infinity;

        float result = pdf * abs(n.dot(-wiNormalized)) / sqr(distance);
        // Clamp the resulting pdf between the machine epsilon and infinity
        return clamp(result, MachineEpsilon, Infinity); 
    }

    /**
     * @brief Converts spherical coordinates to a direction vector in the global coordinate system.
     * From: https://pbr-book.org/3ed-2018/Color_and_Radiometry/Working_with_Radiometric_Integrals#SphericalDirection
     *
     * @param sinTheta The sine of the polar angle theta.
     * @param cosTheta The cosine of the polar angle theta.
     * @param phi The azimuthal angle in radians.
     * @return A direction vector in 3D space corresponding to the spherical coordinates.
     */
    inline Vector sphericalDirection3ed(float sinTheta, float cosTheta, float phi)
    {
        return Vector(sinTheta * cos(phi), sinTheta * sin(phi), cosTheta);
    }

    /**
     * @brief Converts spherical coordinates to a direction vector in a local coordinate system.
     * From: https://pbr-book.org/3ed-2018/Color_and_Radiometry/Working_with_Radiometric_Integrals#SphericalDirection
     *
     * @param sinTheta The sine of the polar angle theta.
     * @param cosTheta The cosine of the polar angle theta.
     * @param phi The azimuthal angle in radians.
     * @param x The x-axis of the local coordinate system.
     * @param y The y-axis of the local coordinate system.
     * @param z The z-axis of the local coordinate system.
     * @return A direction vector in the local coordinate system corresponding to the spherical coordinates.
     */
    inline Vector sphericalDirection3ed(float sinTheta, float cosTheta, float phi, const Vector &x, const Vector &y,
                                        const Vector &z)
    {
        return sinTheta * cos(phi) * x + sinTheta * sin(phi) * y + cosTheta * z;
    }

    /**
     * @brief Converts a theta and phi pair into a unit (x, y, z) vector, applying these equations directly. Notice that
     * the function is given the sine and cosine of theta, rather than theta itself. This is because the sine and cosine
     * of theta are often already available to the caller. This is not normally the case for phi, however, so phi is
     * passed in as is. From:
     * https://pbr-book.org/4ed/Geometry_and_Transformations/Spherical_Geometry#SphericalDirection
     *
     * @param sinTheta The sine of the polar angle theta.
     * @param cosTheta The cosine of the polar angle theta.
     * @param phi The azimuthal angle in radians.
     *
     * @return A (x, y, z) vector in the local coordinate system.
     */
    inline Vector sphericalDirection4ed(float sinTheta, float cosTheta, float phi)
    {
        return Vector(clamp(sinTheta, -1, 1) * cos(phi), clamp(sinTheta, -1, 1) * sin(phi), clamp(cosTheta, -1, 1));
    }

    /**
     * @brief Samples a point on a sphere from a given reference point using uniform sampling within the subtended cone.
     *
     * This function samples points on the sphere that are more likely to be visible from the reference point,
     * based on the solid angle subtended by the sphere. It avoids sampling points on the backside of the sphere
     * that would not be visible. Taken from PBRT's 3rd Edition.
     *
     * From: https://pbr-book.org/3ed-2018/Light_Transport_I_Surface_Reflection/Sampling_Light_Sources
     *
     * @param sample A random 2D point on the unit square used to generate the sample.
     * @param pCenter The center point of the sphere.
     * @param radius The radius of the sphere.
     * @param refPoint The reference point from which the sphere is viewed.
     *
     * @return A 3D point of the sphere in world space (not on the surface of the sphere!).
     */
    inline Vector subtendedConeSphereSampling3ed(const Point2 &sample, const Point &pCenter, const float &radius,
                                                 const Point &refPoint)
    {
        // Compute coordinate system for sphere sampling
        Vector wc = (pCenter - refPoint).normalized();
        Vector wcX, wcY;
        buildOrthonormalBasisPBRT(wc, wcX, wcY);

        // Sample sphere uniformly inside subtended cone
        // Compute theta and phi values for the sample in the cone
        float sinThetaMax2 = radius * radius / (refPoint - pCenter).lengthSquared();
        float cosThetaMax = sqrt(max(0.f, 1 - sinThetaMax2));
        float cosTheta = (1 - sample.x()) + sample.x() * cosThetaMax;
        float sinTheta = sqrt(max(0.f, 1 - cosTheta * cosTheta));
        float phi = sample.y() * 2 * Pi;

        // Compute angle theta from the center of the sphere to the sampled point on the surface
        float dc = (refPoint - pCenter).length();
        float ds = dc * cosTheta - sqrt(max(0.f, radius * radius - dc * dc * sinTheta * sinTheta));
        float cosAlpha = (dc * dc + radius * radius - ds * ds) / (2 * dc * radius);
        float sinAlpha = sqrt(max(0.f, 1 - cosAlpha * cosAlpha));

        // If we wanted to project back on the sphere, compute the surface normal and the sampled point on the sphere
        // Vector nObj = sphericalDirection3ed(sinAlpha, cosAlpha, phi, -wcX, -wcY, -wc);
        // Point p(pCenter + radius * Vector(nObj)); // Sampled point
        // return Point(Vector(p) * radius / (p - pCenter).length()); // Scale the point by the sphere’s radius

        return sphericalDirection3ed(sinAlpha, cosAlpha, phi, -wcX, -wcY, -wc);
    }

    /**
     * @brief Samples a vector on a sphere from a given reference point using uniform sampling within the subtended
     * cone.
     *
     * This function samples points on the sphere that are more likely to be visible from the reference point,
     * based on the solid angle subtended by the sphere. It avoids sampling points on the backside of the sphere
     * that would not be visible. Taken from PBRT's 4th Edition.
     *
     * From: https://pbr-book.org/4ed/Shapes/Spheres#Sampling
     *
     * @param sample A random 2D point on the unit square used to generate the sample.
     * @param pCenter The center point of the sphere.
     * @param radius The radius of the sphere.
     * @param refPoint The reference point from which the sphere is viewed.
     *
     * @return A 3D vector of the sphere in world space (not on the surface of the sphere!).
     */
    inline Vector subtendedConeSphereSampling4ed(const Point2 &sample, const Point &pCenter, const float &radius,
                                                 const Point &refPoint)
    {
        // Compute theta and phi values for sample in cone
        float sinThetaMax = radius / (refPoint - pCenter).length();
        float sin2ThetaMax = sqr(sinThetaMax);
        float cosThetaMax = safe_sqrt(1 - sin2ThetaMax);

        float cosTheta = (cosThetaMax - 1) * sample.x() + 1;
        float sin2Theta = 1 - sqr(cosTheta);

        // Compute cone sample via Taylor series expansion for small angles
        if (sin2ThetaMax < 0.00068523f /* sin^2(1.5 deg) */)
        {
            sin2Theta = sin2ThetaMax * sample.x();
            cosTheta = sqrt(1 - sin2Theta);
        }

        // Compute angle alpha from center of sphere to sampled point on surface
        float cosAlpha = sin2Theta / sinThetaMax + cosTheta * safe_sqrt(1 - sin2Theta / sqr(sinThetaMax));
        float sinAlpha = safe_sqrt(1 - sqr(cosAlpha));
        float phi = sample.y() * 2 * Pi;

        // Compute surface normal and sampled point on sphere
        Vector w = sphericalDirection4ed(sinAlpha, cosAlpha, phi);
        Frame samplingFrame = Frame((pCenter - refPoint).normalized());

        // If we wanted to project back on the sphere
        // Vector n(samplingFrame.toWorld(w));
        // Point p(pCenter + radius * Vector(n)); // Sampled point
        // return Point(Vector(p) * radius / (p - pCenter).length()); // Scale the point by the sphere’s radius

        return samplingFrame.toWorld(-w);
    }

    /// @brief Returns the density of the @ref subtendedConeSphereSampling3ed and subtendedConeSphereSampling4ed
    /// warping.
    inline float subtendedConePdf(const Point &pCenter, const float &radius, const Point &refPoint)
    {
        // Compute general solid angle sphere PDF
        float sin2ThetaMax = sqr(radius) / (refPoint - pCenter).lengthSquared();
        float cosThetaMax = safe_sqrt(1 - sin2ThetaMax);
        float oneMinusCosThetaMax = 1 - cosThetaMax;
        // Compute more accurate oneMinusCosThetaMax for small solid angle
        if (sin2ThetaMax < 0.00068523f /* sin^2(1.5 deg) */)
            oneMinusCosThetaMax = sin2ThetaMax / 2;

        return Inv2Pi / oneMinusCosThetaMax;
    }

    /**
     * @brief Offsets a point slightly along the normal direction to avoid self-intersection based on the geometry
     * normal and the incident direction. From:
     * https://pbr-book.org/3ed-2018/Shapes/Managing_Rounding_Error#OffsetRayOrigin
     *
     * @param p The original point on the surface from which the ray is to be shot.
     * @param n The normal vector at point p on the surface. It should be normalized.
     * @param w The direction vector along which the ray will be shot. It should be normalized.
     *
     * @return The offset point, slightly away from the surface along the normal direction if the
     *         ray is leaving the surface, or slightly into the surface if the ray is entering.
     */
    inline Point offsetRayOrigin(const Point &p, const Vector &n, const Vector &w)
    {
        const Vector pError = Vector(Epsilon);
        float d = Vector(abs(n.x()), abs(n.y()), abs(n.z())).dot(pError);
        Vector offset = d * n;
        if (w.dot(n) < 0)
            offset = -offset;
        Point po = p + offset; // Round offset point po away from p
        return po;
    }

} // namespace lightwave
