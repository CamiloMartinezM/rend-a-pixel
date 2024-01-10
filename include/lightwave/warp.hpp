/**
 * @file warp.hpp
 * @brief Contains functions that map one domain to another.
 */

#pragma once

#include <lightwave/math.hpp>

namespace lightwave {

/**
 * @brief Warps a given point from the unit square ([0,0] to [1,1]) to a unit circle (centered around [0,0] with radius 1),
 * with uniform density given by @code 1 / Pi @endcode .
 * @see Based on http://psgraphics.blogspot.ch/2011/01/improved-code-for-concentric-map.html
 */
inline Point2 squareToUniformDiskConcentric(const Point2 &sample) {
    float r1 = 2 * sample.x() - 1;
    float r2 = 2 * sample.y() - 1;

    float phi, r;
    if (r1 == 0 && r2 == 0) {
        r = 0;
        phi = 0;
    } else if (r1 * r1 > r2 * r2) {
        r = r1;
        phi = Pi4 * (r2 / r1);
    } else {
        r = r2;
        phi = Pi2 - Pi4 * (r1 / r2);
    }

    float cosPhi = std::cos(phi);
    float sinPhi = std::sin(phi);
    return { r * cosPhi, r * sinPhi };
}

/**
 * @brief Warps a given point from the unit square ([0,0] to [1,1]) to a unit sphere (centered around [0,0,0] with radius 1),
 * with uniform density given by @code 1 / (4 * Pi) @endcode .
 */
inline Vector squareToUniformSphere(const Point2 &sample) {
    float z = 1 - 2 * sample.y();
    float r = safe_sqrt(1 - z * z);
    float phi = 2 * Pi * sample.x();
    float cosPhi = std::cos(phi);
    float sinPhi = std::sin(phi);
    return { r * cosPhi, r * sinPhi, z };
}

/**
 * @brief Warps a given point from the unit square ([0,0] to [1,1]) to a unit hemisphere (centered around [0,0,0] with radius 1,
 * pointing in z direction), with respect to solid angle.
 */
inline Vector squareToUniformHemisphere(const Point2 &sample) {
    Point2 p = squareToUniformDiskConcentric(sample);
    float z = 1.0f - p.x() * p.x() - p.y() * p.y();
    float s = sqrt(z + 1.0f);
    return { s * p.x(), s * p.y(), z };
}

/// @brief Returns the density of the @ref squareToUniformHemisphere warping.
inline float uniformHemispherePdf() {
    return Inv2Pi;
}

/**
 * @brief Warps a given point from the unit square ([0,0] to [1,1]) to a unit hemisphere (centered around [0,0,0] with radius 1,
 * pointing in z direction), with density given by @code cos( angle( result, [0,0,1] ) ) @endcode .
 */
inline Vector squareToCosineHemisphere(const Point2 &sample) {
    Point2 p = squareToUniformDiskConcentric(sample);
    float z = safe_sqrt(1.0f - p.x() * p.x() - p.y() * p.y());
    return { p.x(), p.y(), z };
}

/// @brief Returns the density of the @ref squareToCosineHemisphere warping.
inline float cosineHemispherePdf(const Vector &vector) {
    return InvPi * std::max(vector.z(), float(0));
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
inline Vector sphericalDirection(float sinTheta, float cosTheta, float phi)
{
    return Vector(sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta);
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
inline Vector sphericalDirection(float sinTheta, float cosTheta, float phi, const Vector &x, const Vector &y,
                                 const Vector &z)
{
    return sinTheta * cos(phi) * x + sinTheta * sin(phi) * y + cosTheta * z;
}

/**
 * @brief Samples a point on a sphere from a given reference point using uniform sampling within the subtended cone.
 * 
 * This function samples points on the sphere that are more likely to be visible from the reference point,
 * based on the solid angle subtended by the sphere. It avoids sampling points on the backside of the sphere
 * that would not be visible.
 * 
 * From: https://pbr-book.org/3ed-2018/Light_Transport_I_Surface_Reflection/Sampling_Light_Sources
 * 
 * @param sample A random 2D point on the unit square used to generate the sample.
 * @param pCenter The center point of the sphere.
 * @param radius The radius of the sphere.
 * @param refPoint The reference point from which the sphere is viewed.
 * 
 * @return A 3D point on the surface of the sphere in world space.
 */
inline Point subtendedConeUniformSphereSampling(const Point2 &sample, const Point &pCenter, const float &radius,
                                                const Point &refPoint)
{
    // Compute coordinate system for sphere sampling
    Vector wc = (pCenter - refPoint).normalized();
    Vector wcX, wcY;
    buildOrthonormalBasis(wc, wcX, wcY);

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

    // Compute the surface normal and the sampled point on the sphere
    Vector nObj = sphericalDirection(sinAlpha, cosAlpha, phi, -wcX, -wcY, -wc);
    Point pObj = Point(radius * nObj);
    return pObj;
}

/**
 * @brief Offsets a point slightly along the normal direction to avoid self-intersection based on the geometry normal
 * and the incident direction. From: https://pbr-book.org/3ed-2018/Shapes/Managing_Rounding_Error#OffsetRayOrigin
 *
 * @param p The original point on the surface from which the ray is to be shot.
 * @param n The normal vector at point p on the surface. It should be normalized.
 * @param w The direction vector along which the ray will be shot. It should be normalized.
 *
 * @return The offset point, slightly away from the surface along the normal direction if the
 *         ray is leaving the surface, or slightly into the surface if the ray is entering.
 */
inline Point OffsetRayOrigin(const Point &p, const Vector &n, const Vector &w)
{
    const Vector pError = Vector(Epsilon);
    float d = Vector(abs(n.x()), abs(n.y()), abs(n.z())).dot(pError);
    Vector offset = d * n;
    if (w.dot(n) < 0)
        offset = -offset;
    Point po = p + offset; // Round offset point po away from p
    return po;
}
}
