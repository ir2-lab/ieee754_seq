#ifndef _IEEE754_SEQ_
#define _IEEE754_SEQ_

#include <cassert>
#include <limits>
#include <cfloat>
#include <array>
#include <cmath>

// get templated function for the # of mantissa digits
template <class T>
struct ieee754_seq_traits
{
    typedef int int_type;
    constexpr static bool is_supported() { return false; }
    constexpr static int mantissa_digits() { return 0; }
    constexpr static int exponent_digits() { return 0; }
};
template <>
struct ieee754_seq_traits<float>
{
    typedef int int_type;
    constexpr static bool is_supported() { return true; }
    constexpr static int mantissa_digits() { return FLT_MANT_DIG; }
    constexpr static int exponent_digits() { return (sizeof(float) << 3) - 1 - FLT_MANT_DIG; }
};
template <>
struct ieee754_seq_traits<double>
{
    typedef long int_type;
    constexpr static bool is_supported() { return true; }
    constexpr static int mantissa_digits() { return DBL_MANT_DIG; }
    constexpr static int exponent_digits() { return (sizeof(double) << 3) - 1 - DBL_MANT_DIG; }
};
template <>
struct ieee754_seq_traits<long double>
{
    typedef __int128 int_type; // only GCC & clang
    constexpr static bool is_supported() { return true; }
    constexpr static int mantissa_digits() { return LDBL_MANT_DIG; }
    constexpr static int exponent_digits()
    {
        return (sizeof(long double) << 3) - 1 - LDBL_MANT_DIG;
    }
};

/**
 * @brief A sequence of IEEE-754 floating point values
 *
 *
 * \code
 *
 * #define N_BIT 4
 * #define EXP_MIN 0
 * #define EXP_MAX 1
 *
 * // float sequence
 * ieee754_seq<float, N_BIT, EXP_MIN, EXP_MAX> sf;
 *
 * // double sequence
 * ieee754_seq<double, N_BIT, EXP_MIN, EXP_MAX> sd;
 *
 * \endcode
 *
 * For example, a typical for-loop over the sequence can be implemented
 * as follows:
 *
 * \code
 *
 * #define N_BIT 4
 * #define EXP_MIN 0
 * #define EXP_MAX 1
 *
 * ieee754_seq<float, N_BIT, EXP_MIN, EXP_MAX> s;
 *
 * for(auto i = s.begin(); i != s.end(); i++) {
 *     float f = *i; // dereference operator returns the floating point number
 *     std::cout << f << std::endl;
 *     A[i] = ... // array indexing OK - implicit conversion to int
 * }
 *
 * \endcode
 *
 *
 * The following can also be used:
 * \code
 * ieee754_seq<float, N_BIT, EXP_MIN, EXP_MAX> R;
 * for(const float& f : R) {
 *     std::cout << f << std::endl;
 * }
 * \endcode
 *
 * The value of std::numeric_limits::is_iec559
 * defined in the <limits> standard library header
 * is checked at compile-time for compatibility
 * of the underlying implementation with the IEEE-754 std.
 *
 * @todo
 * Include a compile-time C++ check for endianess.
 * Currently gcc in Ubuntu 22.04 does not have std::endian
 *
 *
 * @tparam T floating-point type, either float, double or long double
 * @tparam _Nbits # of mantissa bits
 * @tparam _Emin Exponent of the min value = 2^minExp
 * @tparam _Emax Exponent of the max value = 2^maxExp
 *
 */
template <class T, int _Nbits, int _Emin, int _Emax>
struct ieee754_seq
{
    typedef ieee754_seq_traits<T> traits;
    /// \brief Type of the floating-point numbers
    typedef T real_type;
    /// \brief Type of the integral numbers
    typedef typename traits::int_type int_type;
    /// @brief Number of mantissa bits
    constexpr static const int_type mantissa_bits = _Nbits;
    /// @brief Minimum exponent
    constexpr static const int_type min_exponent = _Emin;
    /// @brief Maximum exponent
    constexpr static const int_type max_exponent = _Emax;

    // compile-time checks
    static_assert(std::numeric_limits<real_type>::is_iec559,
                  "floating-point type is not IEC559/IEEE754 conformant");

    static_assert(sizeof(real_type) == sizeof(int_type),
                  "floating-point & integral types must be of the same size");

    static_assert(mantissa_bits >= 0, "mantissa_bits must be a positive integer or zero.");

    static_assert(mantissa_bits <= traits::mantissa_digits(),
                  "mantissa_bits must be less or equal of real_type mantissa bits.");

    static_assert(max_exponent > min_exponent, "_maxExp must be larger than _minExp");

    static_assert(max_exponent <= std::numeric_limits<real_type>::max_exponent,
                  "_maxExp must less or equal to std::numeric_limits<real_type>::max_exponent");

    static_assert(max_exponent > 1 - std::numeric_limits<real_type>::max_exponent,
                  "_maxExp must grater than 1-std::numeric_limits<real_type>::max_exponent");

    static_assert(min_exponent <= std::numeric_limits<real_type>::max_exponent,
                  "_minExp must less or equal to std::numeric_limits<real_type>::max_exponent");

    static_assert(min_exponent > 1 - std::numeric_limits<real_type>::max_exponent,
                  "_minExp must grater than 1-std::numeric_limits<real_type>::max_exponent");

    // these constants are computed at compile time
    // this is ensured by the constexpr

    /// @brief Number of mantissa values
    constexpr static const int_type Nm = (int_type(1) << mantissa_bits);
    /// @brief Minimum real value in the sequence \f$ \mbox{min} = 2^{E_{min}} \f$
    constexpr static const real_type min = (min_exponent >= int_type(0))
            ? int_type(1) << min_exponent
            : real_type(1) / (int_type(1) << -min_exponent);
    /// @brief Maximum real value in the sequence \f$ \mbox{max} = 2^{E_{max}} \f$
    constexpr static const real_type max = (max_exponent >= int_type(0))
            ? int_type(1) << max_exponent
            : real_type(1) / (int_type(1) << -max_exponent);
    /// @brief Number of elements in the sequence
    constexpr static const int_type count =
            (max_exponent - min_exponent) * (int_type(1) << mantissa_bits) + int_type(1);

private:
    typedef ieee754_seq<T, _Nbits, _Emin, _Emax> _seq_t;
    // Exponent bias
    constexpr static const int_type bias =
            (int_type(std::numeric_limits<real_type>::max_exponent - 1) + min_exponent) * Nm;
    // Exponent shift
    constexpr static const int_type shift =
            (int_type(traits::mantissa_digits() - 1) - mantissa_bits);
    // count-1
    constexpr static const int_type dim = (max_exponent - min_exponent) * Nm;

    constexpr static const int_type mantissa_mask = Nm - 1;
    constexpr static const int_type exp_mask = (int_type(1) << traits::exponent_digits()) - 1;

    constexpr static auto log2_m_{ []() constexpr {
        std::array<real_type, Nm> w{};
        real_type f = real_type(1) / Nm;
        for (int_type i = 0; i < Nm; ++i)
            w[i] = std::log2(1 + f * i);
        return w;
    }() };

    /**
     * @brief Convert an index i to log2(xi)
     *
     * @param index
     * @return float
     */
    constexpr static real_type idx2log2(int_type index)
    {

        assert(index >= 0 && index <= dim);

        int m = index & mantissa_mask;
        int e = ((index >> mantissa_bits) & exp_mask) + min_exponent;

        return log2_m_[m] + e;
    }

    /**
     * @brief Convert float value to index
     *
     * The returned value is truncated to within the index range
     *
     * @param val a floating point value
     * @return unsigned int the corresponding index
     */
    constexpr static int_type val2idx(real_type val)
    {
        if (val <= min)
            return int_type(0);
        if (val >= max)
            return dim;
        int_type ll = *reinterpret_cast<int_type *>(&val);
        ll = (ll >> shift) - bias;
        return ll;
    }

    /**
     * @brief Convert an index to a floating point value
     *
     * @param index
     * @return float
     */
    constexpr static real_type idx2val(int_type index)
    {

        assert(index >= 0 && index <= dim);
        int_type ll = (index + bias) << shift;
        return *reinterpret_cast<real_type *>(&ll);
    }

public:
    class iterator
    {
    public:
        /**
         * @brief Constructor with integral initializer.
         * @param i The initial value, defaults to 0
         */
        constexpr iterator(int_type i = int_type(0)) : i_(i) { }

        /**
         * @brief Constructor with floating-point initializer.
         *
         * The ieee754_seq's IntValue is initialized to fromValue() called with argument v.
         *
         * @param v The floating-point number
         */
        constexpr explicit iterator(const real_type &v) : i_(val2idx(v)) { }

        /**
         * @brief Helper function to convert a real number to index
         *
         * The function finds the real number X_I within the corteo range
         * which is closest to v and returns the corresponding index I.
         *
         * See \ref CorteoIdx for more details.
         *
         * If v is smaller than corteo_index::minVal the function returns begin().
         *
         * If v is larger than corteo_index::maxVal the function returns end().
         *
         * @param v is the floating-point number
         * @return the corresponding index
         */
        constexpr static iterator fromValue(real_type v) { return iterator(val2idx(v)); }

        /**
         * @brief Advance the ieee754_seq by one
         * @return A reference to the ieee754_seq
         */
        iterator &operator++()
        {
            i_++;
            return *this;
        }
        iterator operator++(int)
        {
            iterator retval = *this;
            ++(*this);
            return retval;
        }

        /**
         * @brief Reduce the ieee754_seq by one
         * @return A reference to the ieee754_seq
         */
        iterator &operator--()
        {
            i_--;
            return *this;
        }
        iterator operator--(int)
        {
            iterator retval = *this;
            --(*this);
            return retval;
        }

        /**
         * @brief The dereference operator * returns the corresponding real value
         */
        constexpr real_type operator*() { return idx2val(i_); }

        /**
         * @brief Implicit conversion to int_type
         */
        constexpr operator int_type() const { return i_; }

        /**
         * @brief Returns the corresponding real value
         */
        constexpr real_type value() const { return idx2val(i_); }

        /**
         * @brief Returns the corresponding real value
         */
        constexpr real_type log2v() const { return idx2log2(i_); }

        /**
         * @brief Returns an iterator pointing to the 1st point of the range
         */
        constexpr iterator begin() const { return iterator(int_type(0)); }

        /**
         * @brief Returns an iterator pointing to one past the last point of the range
         */
        constexpr iterator end() const { return iterator(count); }

        /**
         * @brief Returns an iterator pointing to the last point of the range
         */
        constexpr iterator rbegin() const { return iterator(dim); }

        /**
         * @brief Returns an iterator pointing to one before the first point
         */
        constexpr iterator rend() const { return iterator(int_type(-1)); }

    private:
        int_type i_;
    };

    /**
     * @brief Returns an iterator pointing to the 1st point of the range
     */
    constexpr int_type size() const { return count; }

    /**
     * @brief Returns an iterator pointing to the 1st point of the range
     */
    constexpr iterator begin() const { return iterator(int_type(0)); }

    /**
     * @brief Returns an iterator pointing to one past the last point of the range
     */
    constexpr iterator end() const { return iterator(count); }

    /**
     * @brief Returns an iterator pointing to the last point of the range
     */
    constexpr iterator rbegin() const { return iterator(dim); }

    /**
     * @brief Returns an iterator pointing to one before the first point
     */
    constexpr iterator rend() const { return iterator(int_type(-1)); }

    /**
     * @brief Return the i-th real value in the sequence
     */
    constexpr static real_type at(const int_type &i) { return idx2val(i); }

    /**
     * @brief Return the index i such that at(i) <= v < at(i+1)
     */
    constexpr static int_type index_of(const real_type &v) { return val2idx(v); }

    /**
     * @brief Return the base-2 logarithm of the i-th real value in the sequence
     */
    constexpr static real_type log2_at(const int_type &i) { return idx2log2(i); }

    /**
     * @brief Return the i-th real value in the sequence
     */
    constexpr real_type operator[](const int_type &i) { return idx2val(i); }

    /**
     * @brief A linear interpolator for a ieee754_seq
     */
    class log_lin_interp
    {
        typedef _seq_t::iterator iterator;

    public:
        /// @brief Default constructor, creates empty object
        log_lin_interp() = default;

        /**
         * @brief Construct a new lin_interp object to interpolate between the data given by y
         *
         * y is a numeric sequence which supports C-style (y[i]) 0-based indexing. The
         * type of the data elements must be convertible to RealType.
         * The size of y must be equal or larger than that of \ref iterator_t.
         *
         * y is copied to an internal buffer.
         *
         * @tparam sequence_t a container type with C-style indexing
         *
         * @param y the data to interpolate
         */
        template <class sequence_t>
        explicit log_lin_interp(const sequence_t &y)
        {
            set(y);
        }

        /**
         * @brief Call this function to set new interpolator data
         *
         * y is a numeric sequence which supports C-style (y[i]) 0-based indexing. The
         * type of the data elements must be convertible to float.
         *
         * The size of y must be equal or larger than that of the corteo index \ref
         * idx_t.
         *
         * y is copied to an internal buffer.
         *
         * @tparam sequence_t a container type with C-style indexing
         *
         * @param y the data to interpolate
         */
        template <class sequence_t>
        void set(const sequence_t &y)
        {
            for (iterator i(0), j(1); i < count - 1; ++i, ++j) {
                y_[i] = y[i];
                d_[i] = (y[j] - y[i]) / (std::log2((*j) / (*i)));
            }
            y_[count - 1] = y[count - 1];
        }

        /// @brief Returns the linearly interpolated value \f$ y = y(x) \f$
        real_type operator()(const real_type &x) const
        {
            if (x <= _seq_t::min)
                return y_.front();
            if (x >= _seq_t::max)
                return y_.back();
            iterator i(x);
            return y_[i] + d_[i] * std::log2(x / (*i));
        }

        const std::array<real_type, count> &data() const { return y_; }

    private:
        std::array<real_type, count> y_, d_;
    };

    /**
     * @brief A log-log interpolator for a ieee754_seq range
     */
    class log_log_interp
    {
        typedef _seq_t::iterator iterator;

    public:
        /// @brief Default constructor, creates empty object
        log_log_interp() = default;

        /**
         * @brief Construct a new log_interp object to perform log-log interpolation with the data
         * given by y
         *
         * y is a numeric sequence which supports C-style (y[i]) 0-based indexing. The type of
         * the data elements must be convertible to float. log(y[i]) must be finite for all i in the
         * range of \ref iterator_type.
         *
         * The size of y must be equal or larger than ieee754_seq::count \ref idx_t.
         * The fisrt count elements of y will be used.
         *
         * y is copied to an internal buffer.
         *
         * @tparam sequence_t a container type with C-style indexing
         * @param y the data to interpolate
         */
        template <class sequence_t>
        explicit log_log_interp(const sequence_t &y)
        {
            set(y);
        }

        /**
         * @brief Call this function to set new interpolator data
         *
         * y is a numeric sequence which supports C-style (y[i]) 0-based indexing. The type of
         * the data elements must be convertible to float.
         *
         * log(y[i]) must be finite for all i in the range
         *
         * The size of y must be equal or larger than ieee754_seq::count \ref idx_t.
         * The fisrt count elements of y will be used.
         *
         * y is copied to an internal buffer.
         *
         * @tparam sequence_t a container type with C-style indexing
         * @param y the data to interpolate
         */
        template <class sequence_t>
        void set(const sequence_t &y)
        {
            for (iterator i(0), j(1); i < count - 1; i++, j++) {
                d_[i] = (std::log2(y[j] / y[i])) / (std::log2((*j) / (*i)));
                y_[i] = y[i];
            }
            y_[count - 1] = y[count - 1];
        }

        /// @brief Returns the log-log interpolated value y(x)
        real_type operator()(const real_type &x) const
        {
            if (x <= min)
                return y_.front();
            if (x >= max)
                return y_.back();
            iterator i(x);
            return y_[i] * std::exp2(d_[i] * std::log2(x / (*i)));
        }

        const std::array<real_type, count> &data() const { return y_; }

    private:
        std::array<real_type, count> y_, d_;
    };
};

#endif
