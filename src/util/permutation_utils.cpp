// pbrt is Copyright(c) 1998-2020 Matt Pharr, Wenzel Jakob, and Greg Humphreys.
// The pbrt source code is licensed under the Apache License, Version 2.0.
// SPDX: Apache-2.0

#include "../util/permutation_utils.h"
#include <sstream>

namespace lightwave
{
    /**
     * @brief Converts the current state of the DigitPermutation object to a string.
     * This is useful for debugging or logging the internal state of the permutation.
     * @return A string representation of the DigitPermutation object.
     */
    std::string DigitPermutation::toString() const
    {
        std::ostringstream stream;

        stream << "[ DigitPermitation base: " << base << " nDigits: " << nDigits << " permutations: ";
        for (int digitIndex = 0; digitIndex < nDigits; ++digitIndex)
        {
            stream << "[" << digitIndex << "] ( ";
            for (int digitValue = 0; digitValue < base; ++digitValue)
            {
                stream << permutations[digitIndex * base + digitValue];
                if (digitValue != base - 1)
                    stream << ", ";
            }
            stream << " ) ";
        }
        stream << "]";

        return stream.str();
    }

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
    pstd::vector<DigitPermutation> *ComputeRadicalInversePermutations(uint32_t seed, Allocator alloc)
    {
        pstd::vector<DigitPermutation> *perms = alloc.new_object<pstd::vector<DigitPermutation>>(alloc);
        perms->resize(NumPrimes);
        for (int i = 0; i < NumPrimes; ++i)
            (*perms)[i] = DigitPermutation(Primes[i], seed, alloc);
        return perms;
    }

} // namespace lightwave
