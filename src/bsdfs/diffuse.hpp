/**
 * @brief Functions for dealing with diffuse distribution.
 * @file diffuse.hpp
 */

#pragma once

#include <lightwave/math.hpp>

namespace lightwave::diffuse
{
    /**
     * @brief Evaluates the Pdf of the Bsdf evaluation.
     * @param wo The outgoing direction light is scattered in, pointing away
     * from the surface, in local coordinates.
     * @param wi The incoming direction light comes from, pointing away from the surface, in local coordinates.
     */
    inline float pdf(const Vector &wo, const Vector &wi)
    {
        // Check if wo and wi are in the same hemisphere:
        // * If they are, for a Lambertian reflectance, the PDF is 1 divided by the hemisphere area (2 * Pi).
        // * Otherwise, return 0 as there's no reflection.
        return Frame::sameHemisphere(wo, wi) ? Inv2Pi : 0.0f;
    }
} // namespace lightwave::diffuse
