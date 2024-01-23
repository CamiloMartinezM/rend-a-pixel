/**
 * @brief Template-based Mathematical Functions
 *
 * This file contains a collection of specialized mathematical functions adapted from PBRT. These functions include
 * template-based implementations for operations like modulus, clamping, interval search, etc.
 *
 * @file special_math.hpp
 */

#pragma once

namespace lightwave::SpecialMath
{
    /// @brief Pi
    static constexpr float SM_Pi = 3.14159265358979323846f;
    /// @brief 1 / Pi
    static constexpr float SM_InvPi = 0.31830988618379067154f;
    /// @brief Acceptable deviation from 0 for Sinc function
    static constexpr float SM_Epsilon = 1e-5f;

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

    /**
     * @brief Computes the modulus of two values.
     *
     * This function calculates the modulus of `a` by `b`, effectively returning the remainder of the division of `a` by
     * `b`. Unlike the standard modulus operator `%` in C++, this function is designed to handle floating point numbers
     * as well as integers. It ensures that the result is always non-negative and less than `b`, even when `a` is
     * negative. This is particularly useful in scenarios where a positive modulo result is required irrespective of the
     * sign of `a`.
     *
     * From: https://github.com/mmp/pbrt-v4/blob/master/src/pbrt/util/math.h
     *
     * @tparam T The type of the values. It can be any arithmetic type (integer or floating point).
     *
     * @param a The dividend in the modulus operation.
     * @param b The divisor in the modulus operation. It's assumed to be positive and non-zero.
     *
     * @return Returns the modulus of `a` by `b`. The result is guaranteed to be non-negative
     * and in the range [0, b).
     *
     * Example usage:
     * double result = Mod(-3.5, 2.0);  // result will be 0.5.
     */
    template <typename T> inline T Mod(T a, T b)
    {
        T result = a - (a / b) * b;
        return (T)((result < 0) ? result + b : result);
    }

    /**
     * @brief Clamps a value to a specified range.
     *
     * `Clamp` restricts a given value `val` to lie within the range specified by `low` and `high`. One of the key
     * features of this implementation is its flexibility in handling different types for the clamping range (`low` and
     * `high`) and the value being clamped (`val`). This is particularly useful in scenarios where implicit conversion
     * between these types is legal and desirable. For example, it allows for clamping a float to an integer range, like
     * in Clamp(floatValue, 0, 1), which would otherwise be problematic due to C++’s strict template type resolution
     * rules.
     *
     * From: https://pbr-book.org/3ed-2018/Utilities/Main_Include_File#Clamp
     *
     * @tparam T The type of the value to be clamped. This is also the return type of the function.
     * @tparam U The type of the lower bound of the range. It should be implicitly convertible to T.
     * @tparam V The type of the upper bound of the range. It should be implicitly convertible to T.
     *
     * @param val The value to be clamped.
     * @param low The lower bound of the clamping range. `val` will not be less than this.
     * @param high The upper bound of the clamping range. `val` will not be more than this.
     *
     * @return Returns the clamped value, which is the original value adjusted to be within the range [low, high].
     * If `val` is less than `low`, `low` is returned. If `val` is greater than `high`, `high` is returned. Otherwise,
     * `val` itself is returned.
     *
     * Example usage:
     * float floatValue = 1.5;
     * int clampedValue = Clamp(floatValue, 0, 1);  // clampedValue will be 1.
     */
    template <typename T, typename U, typename V> inline T Clamp(T val, U low, V high)
    {
        if (val < low)
            return low;
        else if (val > high)
            return high;
        else
            return val;
    }

    /**
     * @brief Finds the interval in a procedurally generated array using a predicate.
     *
     * `FindInterval` is a helper function that emulates the behavior of `std::upper_bound`. Unlike `std::upper_bound`,
     * it does not require direct access to an array but instead uses a function object to determine values at various
     * indices. This approach allows for bisection of arrays that are procedurally generated or interpolated from point
     * samples. Additionally, the function incorporates bounds checking to handle corner cases, ensuring a valid
     * interval is selected even if the predicate evaluates to true or false for all entries.
     *
     * From: https://pbr-book.org/3ed-2018/Utilities/Main_Include_File#FindInterval
     *
     * @tparam Predicate A callable type that takes an integer index and returns a boolean. This predicate is used to
     * perform the bisection.
     * @param size The size of the virtual array to bisect. It represents the upper bound of the search space.
     * @param pred The predicate function object used for bisection. It should accept an integer index and return a
     * boolean indicating if the condition is met at that index.
     *
     * @return Returns an integer representing the index of the upper bound of the interval in the virtual array that
     * satisfies the predicate. The return value is clamped to be within the range of [0, size - 2], ensuring validity
     * even in edge cases.
     *
     * Example usage:
     * int size = 10;
     * auto pred = [](int index) { return someCondition; };
     * int intervalIndex = FindInterval(size, pred);
     */
    template <typename Predicate> int FindInterval(int size, const Predicate &pred)
    {
        int first = 0, len = size;
        while (len > 0)
        {
            int half = len >> 1, middle = first + half;
            // Bisect range based on value of pred at middle
            if (pred(middle))
            {
                first = middle + 1;
                len -= half + 1;
            }
            else
                len = half;
        }
        return Clamp(first - 1, 0, size - 2);
    }

    /**
     * @brief Checks if a value is a power of 2.
     * From: https://pbr-book.org/3ed-2018/Utilities/Main_Include_File#x1-Base-2Operations
     * 
     * @tparam T The type of the value to be checked.
     * @param v The value to check.
     * @return True if `v` is a power of 2, false otherwise.
     */
    template <typename T> inline bool IsPowerOf2(T v)
    {
        return v && !(v & (v - 1));
    }

    /**
     * @brief Rounds up an integer to the next highest power of 2.
     *
     * @param v The integer to round up.
     * @return The smallest power of 2 greater than or equal to `v`.
     */
    inline int32_t RoundUpPow2(int32_t v)
    {
        v--;
        v |= v >> 1;
        v |= v >> 2;
        v |= v >> 4;
        v |= v >> 8;
        v |= v >> 16;
        return v + 1;
    }

    /**
     * @brief Calculates the base-2 logarithm of an unsigned integer.
     * From: https://pbr-book.org/3ed-2018/Utilities/Main_Include_File#x1-Base-2Operations
     * 
     * @param v The value to compute the logarithm for.
     * @return The base-2 logarithm of `v`.
     */
    inline int Log2Int(uint32_t v)
    {
        return 31 - __builtin_clz(v);
    }

    /**
     * @brief Calculates the base-2 logarithm of a floating-point number.
     * From: https://pbr-book.org/3ed-2018/Utilities/Main_Include_File#x1-Base-2Operations
     * 
     * @param x The value to compute the logarithm for.
     * @return The base-2 logarithm of `x`.
     */
    inline float Log2(float x)
    {
        const float invLog2 = 1.442695040888963387004650940071;
        return std::log(x) * invLog2;
    }

    /**
     * @brief Performs linear interpolation between two values.
     * From: https://pbr-book.org/3ed-2018/Utilities/Mathematical_Routines#Lerp
     * 
     * @param t The interpolation factor (0.0 to 1.0).
     * @param v1 The start value.
     * @param v2 The end value.
     * @return The interpolated value between `v1` and `v2`.
     */
    inline float Lerp(float t, float v1, float v2)
    {
        return (1 - t) * v1 + t * v2;
    }

    /**
     * @brief Templated version of Lerp for generic types.
     *
     * @tparam T The type of the values.
     * @param t The interpolation factor (0.0 to 1.0).
     * @param v1 The start value of type T.
     * @param v2 The end value of type T.
     * @return The interpolated value between `v1` and `v2`.
     */
    template <typename T> T Lerp(float t, const T &v1, const T &v2)
    {
        return (1 - t) * v1 + t * v2;
    }

    /**
     * @brief Computes the normalized sinc function for a value.
     * From: https://pbr-book.org/3ed-2018/Sampling_and_Reconstruction/Image_Reconstruction#LanczosSincFilter::Sinc
     * 
     * @param x The input value.
     * @return The normalized sinc of `x`.
     */
    inline float Sinc(float x)
    {
        x = std::abs(x);
        if (x < SM_Epsilon)
            return 1;
        return SM_InvPi * std::sin(SM_Pi * x) / x;
    }

    /**
     * @brief Computes the Lanczos function for a value. Currently identical to the Sinc function.
     *
     * @param x The input value.
     * @return The Lanczos of `x`.
     */
    inline float Lanczos(float x)
    {
        return Sinc(x);
    }
} // namespace lightwave::SpecialMath
