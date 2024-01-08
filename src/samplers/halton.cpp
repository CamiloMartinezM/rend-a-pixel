#include "../util/scrambling_utils.h"
#include <lightwave.hpp>

namespace lightwave
{
    class Halton : public Sampler
    {
      private:
        pstd::vector<DigitPermutation> *digitPermutations = nullptr;
        static constexpr int MaxHaltonResolution = 128;
        Vector2i fullResolution{(int)Infinity};
        Point2i baseScales, baseExponents;
        int multInverse[2], dimension = 0;
        int64_t haltonIndex = 0;
        uint64_t m_seed;
        Allocator alloc;

        /**
         * @brief Computes the multiplicative inverse of a number modulo n.
         *
         * This function calculates the multiplicative inverse of a given number 'a' modulo 'n' using the Extended
         * Euclidean Algorithm. The multiplicative inverse is a number 'x' such that (a * x) % n == 1.
         *
         * @param a The number for which to find the multiplicative inverse.
         * @param n The modulo.
         * @return The multiplicative inverse of 'a' modulo 'n'.
         */
        inline static uint64_t multiplicativeInverse(int64_t a, int64_t n)
        {
            int64_t x, y;
            extendedGCD(a, n, &x, &y);
            return Mod(x, n);
        }

        /**
         * @brief Computes the greatest common divisor (GCD) of two numbers and the coefficients of Bézout's identity.
         *
         * This function implements the Extended Euclidean Algorithm to find integers x and y such that a * x + b * y =
         * GCD(a, b).
         *
         * @param a The first number.
         * @param b The second number.
         * @param x Pointer to store the coefficient corresponding to 'a'.
         * @param y Pointer to store the coefficient corresponding to 'b'.
         */
        inline static void extendedGCD(uint64_t a, uint64_t b, int64_t *x, int64_t *y)
        {
            if (b == 0)
            {
                *x = 1;
                *y = 0;
                return;
            }
            int64_t d = a / b, xp, yp;
            extendedGCD(b, a % b, &xp, &yp);
            *x = yp;
            *y = xp - (d * yp);
        }

        /**
         * @brief Prepares the sampler for generating samples for a new pixel.
         *
         * This function sets up the Halton sampler for a new pixel. It computes the Halton index based on the pixel
         * position and sample index, considering the base scales and exponents. This method is essential for ensuring
         * that each pixel receives a unique sequence of samples.
         *
         * @param p The pixel coordinates for which to generate samples.
         * @param sampleIndex The index of the sample in the sequence for the given pixel.
         * @param dim The starting dimension for sampling.
         */
        void StartPixelSample(Point2i p, int sampleIndex, int dim)
        {
            haltonIndex = 0;
            int sampleStride = baseScales[0] * baseScales[1];
            if (sampleStride > 1)
            {
                Point2i pm(Mod(p[0], MaxHaltonResolution), Mod(p[1], MaxHaltonResolution));
                for (int i = 0; i < 2; ++i)
                {
                    uint64_t dimOffset = (i == 0) ? InverseRadicalInverse(pm[i], 2, baseExponents[i])
                                                  : InverseRadicalInverse(pm[i], 3, baseExponents[i]);
                    haltonIndex += dimOffset * (sampleStride / baseScales[i]) * multInverse[i];
                }
                haltonIndex %= sampleStride;
            }

            haltonIndex += sampleIndex * sampleStride;
            dimension = max(2, dim);
        }

        /**
         * @brief Sets up the base scales and exponents for the Halton sequence.
         *
         * This function calculates the scales and exponents for the first two dimensions of the Halton sequence based
         * on the resolution of the image and the maximum resolution supported by the Halton sampler. It ensures that
         * the sample values are properly scaled and distributed across the image.
         */
        inline void setUpBaseScalesExponents()
        {
            // Find radical inverse base scales and exponents that cover sampling area
            for (int i = 0; i < 2; ++i)
            {
                int base = (i == 0) ? 2 : 3;
                int scale = 1, exp = 0;
                while (scale < min(fullResolution[i], MaxHaltonResolution))
                {
                    scale *= base;
                    ++exp;
                }
                baseScales[i] = scale;
                baseExponents[i] = exp;
            }

            // Compute multiplicative inverses for _baseScales_
            multInverse[0] = multiplicativeInverse(baseScales[1], baseScales[0]);
            multInverse[1] = multiplicativeInverse(baseScales[0], baseScales[1]);
        }

        /**
         * @brief Generates a sample value for a given dimension using the Halton sequence.
         *
         * This method computes a sample value in the Halton sequence for a specified dimension.
         * The method of computation depends on the randomization strategy set for the Halton sampler.
         * If 'None' is specified, a regular radical inverse is used. If 'PermuteDigits' is specified,
         * the radical inverse is scrambled with a permutation for that dimension. If any other strategy
         * (e.g., 'Owen') is specified, Owen scrambling is applied for further decorrelation.
         *
         * @param dimension The dimension for which to generate the sample. This typically corresponds
         *                  to a dimension in a multi-dimensional integration problem, with each dimension
         *                  using a different prime base in the Halton sequence.
         * @return A sample value in the range [0, 1), computed based on the chosen randomization strategy.
         */
        float SampleDimension(int dimension) const
        {
            if (HaltonRandomizeStrategy == RandomizeStrategy::None)
                return radicalInverse(dimension, haltonIndex);
            else if (HaltonRandomizeStrategy == RandomizeStrategy::PermuteDigits)
                return scrambledRadicalInverse(dimension, haltonIndex, (*digitPermutations)[dimension]);
            else
                return OwenScrambledRadicalInverse(dimension, haltonIndex, MixBits(1 + (dimension << 4)));
        }

      public:
        Halton(const Properties &properties) : Sampler(properties)
        {
            m_seed = properties.get<int>("seed", 1337);
            if (HaltonRandomizeStrategy == RandomizeStrategy::PermuteDigits)
                digitPermutations = ComputeRadicalInversePermutations(m_seed, alloc);
            setUpBaseScalesExponents();
        }

        void seed(int sampleIndex) override
        {
            // No need to re-seed for every sample in Halton as the sequence is deterministic
            dimension = 0;
            haltonIndex = 0;
        }

        void seed(const Point2i &pixel, int sampleIndex) override
        {
            // Starts the pixel sample starting from dimension = 0
            StartPixelSample(pixel, sampleIndex, 0);
        }

        float next() override
        {
            // Compute the next number in the Halton sequence using the current base
            // Corresponds to PBRT's Get1D()
            if (dimension >= NumPrimes)
                dimension = 2;
            return SampleDimension(dimension++);
        }

        Point2 next2D() override
        {
            // Compute the next two numbers in the Halton sequence for 2D sampling
            // Corresponds to PBRT's Get2D()
            if (dimension + 1 >= NumPrimes)
                dimension = 2;
            int dim = dimension;
            dimension += 2;
            return {SampleDimension(dim), SampleDimension(dim + 1)};
        }

        ref<Sampler> clone() const override
        {
            return std::make_shared<Halton>(*this);
        }

        AvailableSampler name() const override
        {
            return AvailableSampler::HaltonSampler;
        }

        void additionalInitialization(const Vector2i &resolution) override
        {
            // Modify the current sampler's resolution
            this->fullResolution = resolution;
            // Re-define base scales and base exponents based on provided image resolution
            setUpBaseScalesExponents();
        }

        std::string toString() const override
        {
            return tfm::format("[ HaltonSampler digitPermutations: %p "
                               "haltonIndex: %d dimension: %d samplesPerPixel: %d "
                               "baseScales: %s baseExponents: %s multInverse: [ %d %d ] ]",
                               indent(digitPermutations), indent(haltonIndex), indent(dimension),
                               indent(m_samplesPerPixel), indent(baseScales), indent(baseExponents),
                               indent(multInverse[0]), indent(multInverse[1]));
        }
    };
} // namespace lightwave

REGISTER_SAMPLER(Halton, "halton")