/**
 * @brief 2D Sampling and Distribution Utilities
 *
 * This header contains the implementation of 1D and 2D distribution structures for Monte Carlo integration.
 * 'Distribution1D' handles one-dimensional probability distribution functions and their sampling methods, while
 * 'Distribution2D' extends these concepts to two dimensions.
 *
 * From:
 * https://pbr-book.org/3ed-2018/Monte_Carlo_Integration/2D_Sampling_with_Multidimensional_Transformations#Distribution2D
 *
 * @file sampling_utils.h
 */

#ifndef SAMPLING_UTILS_H
#define SAMPLING_UTILS_H

#include <lightwave.hpp>

namespace lightwave
{
    /**
     * @brief A 2D vector class template that stores elements in a 1D vector.
     * @tparam T The type of elements stored in the 2D vector.
     */
    template <typename T> class Vector2D
    {
        std::vector<T> data;
        int uRes, vRes;

      public:
        /// @brief Constructor to initialize empty array given the parameters that determine the size.
        Vector2D(int uRes, int vRes) : uRes(uRes), vRes(vRes) { data.resize(uRes * vRes); }
        /// @brief Constructor to initialize from an existing array.
        Vector2D(int uRes, int vRes, const T *d) : uRes(uRes), vRes(vRes)
        {
            data.resize(uRes * vRes);
            if (d)                    
                std::copy(d, d + uRes * vRes, data.begin());
        }

        /// @brief Returns the number of elements along the u-axis.
        int uSize() const { return uRes; }
        /// @brief Returns the number of elements along the v-axis.
        int vSize() const { return vRes; }
        /// @brief Accesses the element at specified coordinates with read-write access.
        T &operator()(int u, int v) { return data[v * uRes + u]; }
        /// @brief Accesses the element at specified coordinates with read-only access.
        const T &operator()(int u, int v) const { return data[v * uRes + u]; }
    };

    /// @brief Represents a 1D probability distribution constructed from an arbitrary function.
    struct Distribution1D
    {
        /**
         * @brief Constructs the 1D distribution from a function defined over a discrete set of points.
         * @param f Array containing the function values at each point.
         * @param n Number of points in the function array.
         */
        Distribution1D(const float *f, int n) : func(f, f + n), cdf(n + 1)
        {
            // Compute integral of step function at
            cdf[0] = 0;
            for (int i = 1; i < n + 1; ++i)
                cdf[i] = cdf[i - 1] + func[i - 1] / n;

            // Transform step function integral into CDF
            funcInt = cdf[n];
            if (funcInt == 0)
            {
                for (int i = 1; i < n + 1; ++i)
                    cdf[i] = float(i) / float(n);
            }
            else
            {
                for (int i = 1; i < n + 1; ++i)
                    cdf[i] /= funcInt;
            }
        }

        /// @brief Returns the number of points in the function array.
        int count() const { return func.size(); }

        /**
         * @brief Samples a value from the continuous distribution.
         * @param u A random variable uniformly distributed over [0, 1).
         * @param pdf Optional pointer to store the computed PDF at the sampled point.
         * @param off Optional pointer to store the offset in the function array corresponding to the sampled value.
         * @return The sampled value, normalized to [0, 1), corresponding to the cumulative distribution function (CDF).
         */
        float sampleContinuous(float u, float *pdf, int *off = nullptr) const
        {
            // Find surrounding CDF segments and offset
            int offset = SpecialMath::FindInterval(cdf.size(), [&](int index) { return cdf[index] <= u; });

            if (off)
                *off = offset;
            // Compute offset along CDF segment
            float du = u - cdf[offset];
            if ((cdf[offset + 1] - cdf[offset]) > 0)
                du /= (cdf[offset + 1] - cdf[offset]);

            // Compute PDF for sampled offset
            if (pdf)
                *pdf = func[offset] / funcInt;

            // Return  corresponding to sample
            return (offset + du) / count();
        }

        /**
         * @brief Samples an index from the discrete distribution.
         * @param u A random variable uniformly distributed over [0, 1).
         * @param pdf Optional pointer to store the computed PDF at the sampled index.
         * @param uRemapped Optional pointer to store the remapped continuous value of u within the sampled segment of
         * the CDF.
         */
        int sampleDiscrete(float u, float *pdf = nullptr, float *uRemapped = nullptr) const
        {
            // Find surrounding CDF segments and offset
            int offset = SpecialMath::FindInterval(cdf.size(), [&](int index) { return cdf[index] <= u; });
            if (pdf)
                *pdf = func[offset] / (funcInt * count());
            if (uRemapped)
                *uRemapped = (u - cdf[offset]) / (cdf[offset + 1] - cdf[offset]);
            return offset;
        }

        /**
         * @brief Computes the PDF for a given index in the discrete distribution.
         * @param index The index of the point in the function array for which to compute the PDF.
         */
        float discretePDF(int index) const
        {
            return func[index] / (funcInt * count());
        }

        std::vector<float> func; // The function values defining the distribution.
        std::vector<float> cdf;  // The cumulative distribution function derived from the function values.
        float funcInt;           // The integral of the function.
    };

    /// @brief Represents a 2D probability distribution constructed from a 2D array of function values.
    class Distribution2D
    {
      public:
        /**
         * @brief Constructs the 2D distribution from a 2D array of function values.
         * @param data Pointer to the array containing the function values in row-major order.
         * @param nu The number of elements along the u-axis (horizontal).
         * @param nv The number of elements along the v-axis (vertical).
         */
        Distribution2D(const float *data, int nu, int nv);

        /**
         * @brief Samples a point from the continuous 2D distribution.
         * @param u A 2D point with components uniformly distributed over [0, 1).
         * @param pdf Pointer to store the computed probability density function value at the sampled point.
         */
        Point2 sampleContinuous(const Point2 &u, float *pdf) const
        {
            float pdfs[2];
            int v;
            float d1 = pMarginal->sampleContinuous(u[1], &pdfs[1], &v);
            float d0 = pConditionalV[v]->sampleContinuous(u[0], &pdfs[0]);
            *pdf = pdfs[0] * pdfs[1];
            return Point2(d0, d1);
        }

        /**
         * @brief Computes the probability density function value for a given 2D point in the distribution.
         * @param p The 2D point for which to compute the PDF, normalized to the [0, 1) range in both dimensions.
         */
        float pdf(const Point2 &p) const
        {
            int iu = SpecialMath::Clamp(int(p[0] * pConditionalV[0]->count()), 0, pConditionalV[0]->count() - 1);
            int iv = SpecialMath::Clamp(int(p[1] * pMarginal->count()), 0, pMarginal->count() - 1);
            return pConditionalV[iv]->func[iu] / pMarginal->funcInt;
        }

      private:
        std::vector<std::unique_ptr<Distribution1D>> pConditionalV; // Conditional distribution for each v.
        std::unique_ptr<Distribution1D> pMarginal;                  // Marginal distribution over v.
    };

} // namespace lightwave

#endif // SCRAMBLING_UTILS_H