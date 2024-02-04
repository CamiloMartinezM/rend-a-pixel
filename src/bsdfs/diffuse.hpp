/**
 * @brief Functions for dealing with diffuse distribution.
 * @file diffuse.hpp
 */

#pragma once

#include <lightwave/math.hpp>
#include <lightwave/warp.hpp>

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
        return cosineHemispherePdf(wi);
    }
} // namespace lightwave::diffuse
