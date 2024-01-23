/**
 * @brief Scrambling and Radical Inverse Functions for Quasi-Random Number Generation
 *
 * This header defines a set of functions for computing the radical inverse and its scrambled variants, crucial in
 * generating low-discrepancy sequences for quasi-random number generation.
 *
 * From: https://www.pbr-book.org/4ed/Sampling_and_Reconstruction/Halton_Sampler#RadicalInverse
 *
 * @file scrambling_utils.h
 */

#ifndef SCRAMBLING_UTILS_H
#define SCRAMBLING_UTILS_H

#include "../util/hash_utils.h"
#include "../util/permutation_utils.h"
#include <lightwave.hpp>

namespace lightwave
{
    /**
     * @brief Calculates the Halton Sequence for base 2.
     * Computes the Van der Corput radical inverse h_2(i), using basic bit operations. From:
     * http://holger.dammertz.org/stuff/notes_HammersleyOnHemisphere.html. Hacker's Delight, Henry S. Warren, 2001
     * @param bits The input number for which the radical inverse is to be computed as a 32bit integer.
     * @return The radical inverse of the input number 'bits' using base 2.
     */
    inline float radicalInverse_VdC(uint bits)
    {
        bits = (bits << 16u) | (bits >> 16u);
        bits = ((bits & 0x55555555u) << 1u) | ((bits & 0xAAAAAAAAu) >> 1u);
        bits = ((bits & 0x33333333u) << 2u) | ((bits & 0xCCCCCCCCu) >> 2u);
        bits = ((bits & 0x0F0F0F0Fu) << 4u) | ((bits & 0xF0F0F0F0u) >> 4u);
        bits = ((bits & 0x00FF00FFu) << 8u) | ((bits & 0xFF00FF00u) >> 8u);
        return float(bits) * 2.3283064365386963e-10; // / 0x100000000
    }

    /**
     * @brief Computes the radical inverse of a given integer 'a' using the baseIndex-th prime number as the base.
     *
     * The radical inverse function is a key component in quasi-random number generation, particularly in the
     * context of generating low-discrepancy sequences for use in sample generation. This implementation uses the
     * prime number corresponding to 'baseIndex' as the base for the radical inverse computation. The prime bases
     * are often used in Halton sequences to ensure that the resulting sequence has low-discrepancy properties,
     * which are beneficial for tasks like rendering, where they can help to achieve more uniform sampling. From:
     * https://www.pbr-book.org/4ed/Sampling_and_Reconstruction/Halton_Sampler#RadicalInverse
     *
     * @param baseIndex The index of the prime number to be used as the base for the radical inverse calculation.
     *                  This index is not the prime number itself, but rather the position in an ordered list of
     *                  prime numbers.
     * @param a The integer input for which the radical inverse is computed. Typically the index of the sample.
     * @return The computed radical inverse of 'a' as a floating-point number in the range [0,1).
     */
    inline float radicalInverse(int baseIndex, uint64_t a)
    {
        unsigned int base = Primes[baseIndex];
        if (base == 2) // Use Van der Corput radical inverse
        {
            // Ensure that 'a' fits into a 32-bit unsigned integer.
            // This might truncate 'a' if it's larger than 32 bits.
            uint32_t a32 = static_cast<uint32_t>(a);
            return min(radicalInverse_VdC(a32), OneMinusEpsilon);
        }
        else
        {
            // We have to stop once reversedDigits is >= limit since otherwise the
            // next digit of |a| may cause reversedDigits to overflow.
            uint64_t limit = ~0ull / base - base;
            float invBase = (float)1 / (float)base, invBaseM = 1;
            uint64_t reversedDigits = 0;
            while (a && reversedDigits < limit)
            {
                // Extract least significant digit from _a_ and update _reversedDigits_
                uint64_t next = a / base;
                uint64_t digit = a - next * base;
                reversedDigits = reversedDigits * base + digit;
                invBaseM *= invBase;
                a = next;
            }
            return min(reversedDigits * invBaseM, OneMinusEpsilon);
        }
    }

    /**
     * @brief Computes the inverse of the radical inverse function.
     *
     * This function is designed to reverse the process of the radical inverse function. It takes an integer
     * composed of digits in reverse order (as if they were computed by the RadicalInverse function and then
     * reversed) and reconstructs the original index 'a'. The 'inverse' parameter should be the value obtained by
     * taking the integer part of the radical inverse and reversing its digits. The function accounts for the number
     * of digits (nDigits) that were used in the original radical inverse computation to correctly reconstruct the
     * index, including any leading zeros that were lost during reversal. This is important for ensuring that values
     * like 1234 and 123400, which both reverse to 4321, can be distinguished based on the original number of
     * digits. From: https://www.pbr-book.org/4ed/Sampling_and_Reconstruction/Halton_Sampler#RadicalInverse
     *
     * @param inverse The reversed integer digits obtained from the radical inverse function.
     * @param base The base used for the radical inverse computation.
     * @param nDigits The number of digits that were used in the original radical inverse computation.
     * @return The original index 'a' that corresponds to the reversed digits.
     */
    inline uint64_t InverseRadicalInverse(uint64_t inverse, int base, int nDigits)
    {
        uint64_t index = 0;
        for (int i = 0; i < nDigits; ++i)
        {
            uint64_t digit = inverse % base;
            inverse /= base;
            index = index * base + digit;
        }
        return index;
    }

    /**
     * @brief Computes the scrambled radical inverse of a given number using a specified base and permutation.
     *
     * This method generates low-discrepancy sequences for sampling, by modifying the radical inverse by applying a
     * permutation to each digit: Ψ_b(a) = 0.π_1(d_1)π_2(d_2)...π_m(d_m). This scrambling can lead to better
     * quasi-random sequences by reducing correlation between consecutive numbers in the sequence. From:
     * https://www.pbr-book.org/4ed/Sampling_and_Reconstruction/Halton_Sampler#RadicalInverse
     *
     * @param baseIndex The index of the prime number to be used as the base for the radical inverse calculation.
     *                  This index refers to a sequence of prime numbers (e.g., 2, 3, 5, 7, 11, ...).
     * @param a The input number for which the radical inverse is to be computed. It's often an index of the sample.
     * @param perm The DigitPermutation object used to scramble the digits of the radical inverse.
     *             This object contains the permutation data and methods necessary for the scrambling process.
     *
     * @return The scrambled radical inverse of the input number 'a' using the specified base and digit permutation.
     */
    inline float scrambledRadicalInverse(int baseIndex, uint64_t a, const DigitPermutation &perm)
    {
        unsigned int base = Primes[baseIndex];
        // We have to stop once reversedDigits is >= limit since otherwise the
        // next digit of |a| may cause reversedDigits to overflow.
        uint64_t limit = ~0ull / base - base;
        float invBase = (float)1 / (float)base, invBaseM = 1;
        uint64_t reversedDigits = 0;
        int digitIndex = 0;
        while (1 - (base - 1) * invBaseM < 1 && reversedDigits < limit)
        {
            // Permute least significant digit from _a_ and update _reversedDigits_
            uint64_t next = a / base;
            int digitValue = a - next * base;
            reversedDigits = reversedDigits * base + perm.Permute(digitIndex, digitValue);
            invBaseM *= invBase;
            ++digitIndex;
            a = next;
        }
        return min(invBaseM * reversedDigits, OneMinusEpsilon);
    }

    /**
     * @brief An even more effective scrambling approach defines digit permutations that not only depend on the
     * index of the current digit i, but that also depend on the values of the previous digits d_1 d_2... d_{i-1}.
     */
    inline float OwenScrambledRadicalInverse(int baseIndex, uint64_t a, uint32_t hash)
    {
        unsigned int base = Primes[baseIndex];
        // We have to stop once reversedDigits is >= limit since otherwise the
        // next digit of |a| may cause reversedDigits to overflow.
        uint64_t limit = ~0ull / base - base;
        float invBase = (float)1 / (float)base, invBaseM = 1;
        uint64_t reversedDigits = 0;
        int digitIndex = 0;
        while (1 - invBaseM < 1 && reversedDigits < limit)
        {
            // Compute Owen-scrambled digit for _digitIndex_
            uint64_t next = a / base;
            int digitValue = a - next * base;
            uint32_t digitHash = MixBits(hash ^ reversedDigits);
            digitValue = PermutationElement(digitValue, base, digitHash);
            reversedDigits = reversedDigits * base + digitValue;
            invBaseM *= invBase;
            ++digitIndex;
            a = next;
        }
        return min(invBaseM * reversedDigits, OneMinusEpsilon);
    }
} // namespace lightwave

#endif // SCRAMBLING_UTILS_H