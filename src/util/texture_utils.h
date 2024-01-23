/**
 * @brief Texture Utilities (particularly MIPMap Implementation)
 *
 * This header defines the MIPMap class and associated utilities. It includes the implementation of mipmapping for
 * efficient texture sampling at various resolutions, supporting operations like EWA (Elliptical Weighted Averaging)
 * filtering, trilinear filtering, and texture lookup optimizations.
 * 
 * From: https://pbr-book.org/3ed-2018/Texture/Image_Texture#MIPMap
 * 
 * @file texture_utils.h
 */

#ifndef TEXTURE_UTILS_H
#define TEXTURE_UTILS_H

#include <lightwave.hpp>
#include <numeric>

namespace lightwave
{
    // MIPMap Helper Declarations
    enum class ImageWrap
    {
        Repeat,
        Black,
        Clamp
    };

    struct ResampleWeight
    {
        int firstTexel;
        float weight[4];
    };

    // MIPMap Declarations
    template <typename T> class MIPMap
    {
      public:
        // MIPMap Public Methods
        MIPMap(const Point2i &resolution, const T *data, bool doTri = false, float maxAniso = 8.f,
               ImageWrap wrapMode = ImageWrap::Repeat);

        int Width() const
        {
            return resolution[0];
        }

        int Height() const
        {
            return resolution[1];
        }

        int Levels() const
        {
            return pyramid.size();
        }

        const T &Texel(int level, int s, int t) const;
        T Lookup(const Point2 &st, float width = 0.f) const;
        T Lookup(const Point2 &st, Vector2 dstdx, Vector2 dstdy) const;

      private:
        // MIPMap Private Methods
        std::unique_ptr<ResampleWeight[]> resampleWeights(int oldRes, int newRes)
        {
            // newRes >= oldRes
            std::unique_ptr<ResampleWeight[]> wt(new ResampleWeight[newRes]);
            float filterwidth = 2.f;
            for (int i = 0; i < newRes; ++i)
            {
                // Compute image resampling weights for _i_th texel
                float center = (i + .5f) * oldRes / newRes;
                wt[i].firstTexel = std::floor((center - filterwidth) + 0.5f);
                for (int j = 0; j < 4; ++j)
                {
                    float pos = wt[i].firstTexel + j + .5f;
                    wt[i].weight[j] = SpecialMath::Lanczos((pos - center) / filterwidth);
                }

                // Normalize filter weights for texel resampling
                float invSumWts = 1 / (wt[i].weight[0] + wt[i].weight[1] + wt[i].weight[2] + wt[i].weight[3]);
                for (int j = 0; j < 4; ++j)
                    wt[i].weight[j] *= invSumWts;
            }
            return wt;
        }

        float clamp(float v)
        {
            return SpecialMath::Clamp(v, 0.f, Infinity);
        }

        Color clamp(const Color &v)
        {
            return v.clamp(0.f, Infinity);
        }

        T triangle(int level, const Point2 &st) const;
        T EWA(int level, Point2 st, Vector2 dst0, Vector2 dst1) const;

        // MIPMap Private Data
        const bool doTrilinear;
        const float maxAnisotropy;
        const ImageWrap wrapMode;
        Point2i resolution;
        std::vector<std::unique_ptr<SpecialMath::Vector2D<T>>> pyramid;
        static constexpr int WeightLUTSize = 128;
        static float weightLut[WeightLUTSize];
    };

    // MIPMap Method Definitions
    template <typename T>
    MIPMap<T>::MIPMap(const Point2i &res, const T *img, bool doTrilinear, float maxAnisotropy, ImageWrap wrapMode)
        : doTrilinear(doTrilinear), maxAnisotropy(maxAnisotropy), wrapMode(wrapMode), resolution(res)
    {
        std::unique_ptr<T[]> resampledImage = nullptr;
        if (!SpecialMath::IsPowerOf2(resolution[0]) || !SpecialMath::IsPowerOf2(resolution[1]))
        {
            // Resample image to power-of-two resolution
            Point2i resPow2(SpecialMath::RoundUpPow2(resolution[0]), SpecialMath::RoundUpPow2(resolution[1]));

            // Resample image in s direction
            std::unique_ptr<ResampleWeight[]> sWeights = resampleWeights(resolution[0], resPow2[0]);
            resampledImage.reset(new T[resPow2[0] * resPow2[1]]);

            // Create a range of integers from 0 to resolution[1] - 1
            std::vector<int> tRange(resolution[1]);
            std::iota(tRange.begin(), tRange.end(), 0);

            // Apply sWeights to zoom in s direction
            lightwave::for_each_parallel(tRange.begin(), tRange.end(), [&](int t) {
                for (int s = 0; s < resPow2[0]; ++s)
                {
                    // Compute texel (s,t) in s-zoomed image
                    resampledImage[t * resPow2[0] + s] = T::black();
                    for (int j = 0; j < 4; ++j)
                    {
                        int origS = sWeights[s].firstTexel + j;
                        if (wrapMode == ImageWrap::Repeat)
                            origS = SpecialMath::Mod(origS, resolution[0]);
                        else if (wrapMode == ImageWrap::Clamp)
                            origS = SpecialMath::Clamp(origS, 0, resolution[0] - 1);
                        if (origS >= 0 && origS < (int)resolution[0])
                            resampledImage[t * resPow2[0] + s] +=
                                sWeights[s].weight[j] * img[t * resolution[0] + origS];
                    }
                }
            });

            // Resample image in t direction
            std::unique_ptr<ResampleWeight[]> tWeights = resampleWeights(resolution[1], resPow2[1]);
            std::vector<T *> resampleBufs;
            int nThreads = lightwave::MaxThreads;
            for (int i = 0; i < nThreads; ++i)
                resampleBufs.push_back(new T[resPow2[1]]);

            // Create a range of integers for 's'
            std::vector<int> sRange(resPow2[0]);
            std::iota(sRange.begin(), sRange.end(), 0);

            // Iterate over the range of 's'
            lightwave::for_each_parallel(sRange.begin(), sRange.end(), [&](int s) {
                T *workData = resampleBufs[ThreadIndex]; // Assuming ThreadIndex can be determined
                for (int t = 0; t < resPow2[1]; ++t)
                {
                    workData[t] = T::black();
                    for (int j = 0; j < 4; ++j)
                    {
                        int offset = tWeights[t].firstTexel + j;
                        if (wrapMode == ImageWrap::Repeat)
                            offset = SpecialMath::Mod(offset, resolution[1]);
                        else if (wrapMode == ImageWrap::Clamp)
                            offset = SpecialMath::Clamp(offset, 0, (int)resolution[1] - 1);
                        if (offset >= 0 && offset < (int)resolution[1])
                            workData[t] += tWeights[t].weight[j] * resampledImage[offset * resPow2[0] + s];
                    }
                }
                for (int t = 0; t < resPow2[1]; ++t)
                    resampledImage[t * resPow2[0] + s] = clamp(workData[t]);
            });

            for (auto ptr : resampleBufs)
                delete[] ptr;

            resolution = resPow2;
        }
        // Initialize levels of MIPMap from image
        int nLevels = 1 + SpecialMath::Log2Int(max(resolution[0], resolution[1]));
        pyramid.resize(nLevels);

        // Initialize most detailed level of MIPMap
        pyramid[0].reset(
            new SpecialMath::Vector2D<T>(resolution[0], resolution[1], resampledImage ? resampledImage.get() : img));
        for (int i = 1; i < nLevels; ++i)
        {
            // Initialize i-th MIPMap level from (i-1) level
            int sRes = max(1, pyramid[i - 1]->uSize() / 2);
            int tRes = max(1, pyramid[i - 1]->vSize() / 2);
            pyramid[i].reset(new SpecialMath::Vector2D<T>(sRes, tRes));

            // Create a range of integers for 't'
            std::vector<int> tRange(tRes);
            std::iota(tRange.begin(), tRange.end(), 0);

            // Filter four texels from finer level of pyramid
            lightwave::for_each_parallel(tRange.begin(), tRange.end(), [&](int t) {
                for (int s = 0; s < sRes; ++s)
                    (*pyramid[i])(s, t) = .25f * (Texel(i - 1, 2 * s, 2 * t) + Texel(i - 1, 2 * s + 1, 2 * t) +
                                                  Texel(i - 1, 2 * s, 2 * t + 1) + Texel(i - 1, 2 * s + 1, 2 * t + 1));
            });
        }

        // Initialize EWA filter weights if needed
        if (weightLut[0] == 0.)
        {
            for (int i = 0; i < WeightLUTSize; ++i)
            {
                float alpha = 2;
                float r2 = float(i) / float(WeightLUTSize - 1);
                weightLut[i] = exp(-alpha * r2) - exp(-alpha);
            }
        }
    }

    template <typename T> const T &MIPMap<T>::Texel(int level, int s, int t) const
    {
        const SpecialMath::Vector2D<T> &l = *pyramid[level];
        // Compute texel (s,t) accounting for boundary conditions
        switch (wrapMode)
        {
        case ImageWrap::Repeat:
            s = SpecialMath::Mod(s, l.uSize());
            t = SpecialMath::Mod(t, l.vSize());
            break;
        case ImageWrap::Clamp:
            s = SpecialMath::Clamp(s, 0, l.uSize() - 1);
            t = SpecialMath::Clamp(t, 0, l.vSize() - 1);
            break;
        case ImageWrap::Black: {
            static const T black = T::black();
            if (s < 0 || s >= (int)l.uSize() || t < 0 || t >= (int)l.vSize())
                return black;
            break;
        }
        }
        return l(s, t);
    }

    template <typename T> T MIPMap<T>::Lookup(const Point2 &st, float width) const
    {
        // Compute MIPMap level for trilinear filtering
        float level = Levels() - 1 + SpecialMath::Log2(std::max(width, (float)1e-8));

        // Perform trilinear interpolation at appropriate MIPMap level
        if (level < 0)
            return triangle(0, st);
        else if (level >= Levels() - 1)
            return Texel(Levels() - 1, 0, 0);
        else
        {
            int iLevel = floor(level);
            float delta = level - iLevel;
            return SpecialMath::Lerp(delta, triangle(iLevel, st), triangle(iLevel + 1, st));
        }
    }

    template <typename T> T MIPMap<T>::triangle(int level, const Point2 &st) const
    {
        level = SpecialMath::Clamp(level, 0, Levels() - 1);
        float s = st[0] * pyramid[level]->uSize() - 0.5f;
        float t = st[1] * pyramid[level]->vSize() - 0.5f;
        int s0 = floor(s), t0 = floor(t);
        float ds = s - s0, dt = t - t0;
        return (1 - ds) * (1 - dt) * Texel(level, s0, t0) + (1 - ds) * dt * Texel(level, s0, t0 + 1) +
               ds * (1 - dt) * Texel(level, s0 + 1, t0) + ds * dt * Texel(level, s0 + 1, t0 + 1);
    }

    template <typename T> T MIPMap<T>::Lookup(const Point2 &st, Vector2 dst0, Vector2 dst1) const
    {
        if (doTrilinear)
        {
            float width = std::max(std::max(std::abs(dst0[0]), std::abs(dst0[1])),
                                   std::max(std::abs(dst1[0]), std::abs(dst1[1])));
            return Lookup(st, width);
        }
        // Compute ellipse minor and major axes
        if (dst0.lengthSquared() < dst1.lengthSquared())
            std::swap(dst0, dst1);
        float majorLength = dst0.length();
        float minorLength = dst1.length();

        // Clamp ellipse eccentricity if too large
        if (minorLength * maxAnisotropy < majorLength && minorLength > 0)
        {
            float scale = majorLength / (minorLength * maxAnisotropy);
            dst1 *= scale;
            minorLength *= scale;
        }
        if (minorLength == 0)
            return triangle(0, st);

        // Choose level of detail for EWA lookup and perform EWA filtering
        float lod = max((float)0, Levels() - (float)1 + SpecialMath::Log2(minorLength));
        int ilod = floor(lod);
        return SpecialMath::Lerp(lod - ilod, EWA(ilod, st, dst0, dst1), EWA(ilod + 1, st, dst0, dst1));
    }

    template <typename T> T MIPMap<T>::EWA(int level, Point2 st, Vector2 dst0, Vector2 dst1) const
    {
        if (level >= Levels())
            return Texel(Levels() - 1, 0, 0);
        // Convert EWA coordinates to appropriate scale for level
        st[0] = st[0] * pyramid[level]->uSize() - 0.5f;
        st[1] = st[1] * pyramid[level]->vSize() - 0.5f;
        dst0[0] *= pyramid[level]->uSize();
        dst0[1] *= pyramid[level]->vSize();
        dst1[0] *= pyramid[level]->uSize();
        dst1[1] *= pyramid[level]->vSize();

        // Compute ellipse coefficients to bound EWA filter region
        float A = dst0[1] * dst0[1] + dst1[1] * dst1[1] + 1;
        float B = -2 * (dst0[0] * dst0[1] + dst1[0] * dst1[1]);
        float C = dst0[0] * dst0[0] + dst1[0] * dst1[0] + 1;
        float invF = 1 / (A * C - B * B * 0.25f);
        A *= invF;
        B *= invF;
        C *= invF;

        // Compute the ellipse's (s,t) bounding box in texture space
        float det = -B * B + 4 * A * C;
        float invDet = 1 / det;
        float uSqrt = sqrt(det * C), vSqrt = sqrt(A * det);
        int s0 = ceil(st[0] - 2 * invDet * uSqrt);
        int s1 = floor(st[0] + 2 * invDet * uSqrt);
        int t0 = ceil(st[1] - 2 * invDet * vSqrt);
        int t1 = floor(st[1] + 2 * invDet * vSqrt);

        // Scan over ellipse bound and compute quadratic equation
        T sum(0.f);
        float sumWts = 0;
        for (int it = t0; it <= t1; ++it)
        {
            float tt = it - st[1];
            for (int is = s0; is <= s1; ++is)
            {
                float ss = is - st[0];
                // Compute squared radius and filter texel if inside ellipse
                float r2 = A * ss * ss + B * ss * tt + C * tt * tt;
                if (r2 < 1)
                {
                    int index = min((int)(r2 * WeightLUTSize), WeightLUTSize - 1);
                    float weight = weightLut[index];
                    sum += Texel(level, is, it) * weight;
                    sumWts += weight;
                }
            }
        }
        return sum / sumWts;
    }

    template <typename T> float MIPMap<T>::weightLut[WeightLUTSize];

} // namespace lightwave

#endif // TEXTURE_UTILS_H
