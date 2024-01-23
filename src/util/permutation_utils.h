/**
 * @brief Permutation Utilities for Quasi-Random Number Generation
 *
 * This header defines utilities for generating and managing digit permutations used in scrambling radical inverses for
 * quasi-random number generation. It includes the 'DigitPermutation' class for creating and handling individual
 * permutations and a utility function to generate permutations for each prime number up to a specified maximum.
 *
 * From: https://www.pbr-book.org/4ed/Sampling_and_Reconstruction/Halton_Sampler#RadicalInverse
 *
 * @file permutation_utils.h
 */

#ifndef PERMUTATION_UTILS_H
#define PERMUTATION_UTILS_H

#include "../util/hash_utils.h"
#include "../util/pstd.h"
#include <lightwave.hpp>

// Added to define _Allocator_ used in PBRT
namespace pstd
{
    namespace pmr
    {
        template <typename T> class polymorphic_allocator;
    }
} // namespace pstd

using Allocator = pstd::pmr::polymorphic_allocator<std::byte>;

namespace lightwave
{
    /**
     * @brief Represents a permutation of digits for a given base used in scrambling radical inverses.
     *
     * This class is responsible for generating and storing a permutation of digits for a specific base.
     * The permutations are used to scramble the digits in the computation of radical inverses, which are
     * crucial for generating low-discrepancy sequences for quasi-random number generation. The class uses
     * a provided allocator to handle memory management for the permutations.
     */
    class DigitPermutation
    {
      public:
        /**
         * @brief Default constructor for DigitPermutation.
         */
        DigitPermutation() = default;

        /**
         * @brief Constructs a DigitPermutation for a given base using a seed value.
         *
         * This constructor initializes a DigitPermutation object by calculating the required number of
         * digits for the given base and generating a random permutation for each digit based on the seed.
         *
         * @param base The base to compute the number of digits for permutation.
         * @param seed The seed value used for random number generation in creating the permutations.
         * @param alloc An allocator instance used for memory management.
         */
        DigitPermutation(int base, uint32_t seed, Allocator alloc) : base(base)
        {
            // base, should be uint16_t
            // Compute number of digits needed for _base_
            nDigits = 0;
            float invBase = (float)1 / (float)base, invBaseM = 1;
            while (1 - (base - 1) * invBaseM < 1)
            {
                ++nDigits;
                invBaseM *= invBase;
            }

            permutations = alloc.allocate_object<uint16_t>(nDigits * base);
            // Compute random permutations for all digits
            for (int digitIndex = 0; digitIndex < nDigits; ++digitIndex)
            {
                uint64_t dseed = Hash(base, digitIndex, seed);
                for (int digitValue = 0; digitValue < base; ++digitValue)
                {
                    int index = digitIndex * base + digitValue;
                    permutations[index] = PermutationElement(digitValue, base, dseed);
                }
            }
        }

        /**
         * @brief Permutes a digit based on the stored permutations.
         *
         * Given a digit index and digit value, this method returns the permuted value as per the
         * precomputed permutations. This method is typically called during the scrambling of
         * radical inverse values.
         *
         * @param digitIndex The index of the digit to be permuted.
         * @param digitValue The value of the digit to be permuted.
         * @return The permuted value of the digit.
         */
        int Permute(int digitIndex, int digitValue) const
        {
            // The following conditions should be met: digitIndex < nDigits; digitValue < base
            return permutations[digitIndex * base + digitValue];
        }

        /**
         * @brief Converts the current state of the DigitPermutation object to a string.
         * This is useful for debugging or logging the internal state of the permutation.
         * @return A string representation of the DigitPermutation object.
         */
        std::string toString() const;

      private:
        // DigitPermutation Private Members
        int base;               // The base for the permutations.
        int nDigits;            // The number of digits in the base.
        uint16_t *permutations; // The array of permuted digit values.
    };

    /**
     * @brief Generates a vector of DigitPermutations for each prime number up to a specified maximum.
     *
     * This utility function constructs a new vector and populates it with DigitPermutation objects,
     * one for each prime number up to the maximum prime used in radical inverse calculations. It uses
     * a provided allocator to handle the creation and storage of the permutations.
     *
     * @param seed The seed value for random number generation used in creating the permutations.
     * @param alloc An allocator instance used for memory management of the vector and permutations.
     * @return A pointer to a vector containing the DigitPermutation objects for each base prime number.
     */
    pstd::vector<DigitPermutation> *ComputeRadicalInversePermutations(uint32_t seed, Allocator alloc = {});

} // namespace lightwave

#endif // PERMUTATION_UTILS_H