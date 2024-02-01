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
    template <typename T> class Vector2D
    {
        std::vector<T> data;
        int uRes, vRes;

      public:
        /// @brief Constructor to initialize empty array given the parameters that determine the size.
        Vector2D(int uRes, int vRes) : uRes(uRes), vRes(vRes)
        {
            data.resize(uRes * vRes);
        }

        /// @brief Constructor to initialize from an existing array
        Vector2D(int uRes, int vRes, const T *d) : uRes(uRes), vRes(vRes)
        {
            data.resize(uRes * vRes); // Allocate space for the data
            if (d)                    // Assign data one element at a time
                for (int v = 0; v < vRes; ++v)
                    for (int u = 0; u < uRes; ++u)
                        (*this)(u, v) = d[v * uRes + u];
        }

        int uSize() const
        {
            return uRes;
        }
        int vSize() const
        {
            return vRes;
        }

        // Accessor that mimics 2D array access
        T &operator()(int u, int v)
        {
            return data[v * uRes + u];
        }

        const T &operator()(int u, int v) const
        {
            return data[v * uRes + u];
        }
    };

    struct Distribution1D
    {
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

        int count() const
        {
            return func.size();
        }

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

        float discretePDF(int index) const
        {
            return func[index] / (funcInt * count());
        }

        std::vector<float> func, cdf;
        float funcInt;
    };
    class Distribution2D
    {
      public:
        Distribution2D(const float *data, int nu, int nv);
        Point2 sampleContinuous(const Point2 &u, float *pdf) const
        {
            float pdfs[2];
            int v;
            float d1 = pMarginal->sampleContinuous(u[1], &pdfs[1], &v);
            float d0 = pConditionalV[v]->sampleContinuous(u[0], &pdfs[0]);
            *pdf = pdfs[0] * pdfs[1];
            return Point2(d0, d1);
        }

        float pdf(const Point2 &p) const
        {
            int iu = SpecialMath::Clamp(int(p[0] * pConditionalV[0]->count()), 0, pConditionalV[0]->count() - 1);
            int iv = SpecialMath::Clamp(int(p[1] * pMarginal->count()), 0, pMarginal->count() - 1);
            return pConditionalV[iv]->func[iu] / pMarginal->funcInt;
        }

      private:
        std::vector<std::unique_ptr<Distribution1D>> pConditionalV;
        std::unique_ptr<Distribution1D> pMarginal;
    };

} // namespace lightwave

#endif // SCRAMBLING_UTILS_H