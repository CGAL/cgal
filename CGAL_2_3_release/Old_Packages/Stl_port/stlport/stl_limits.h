/*
 * Copyright (c) 1997
 * Silicon Graphics Computer Systems, Inc.
 *
 * Copyright (c) 1999 
 * Boris Fomitchev
 *
 * This material is provided "as is", with absolutely no warranty expressed
 * or implied. Any use is at your own risk.
 *
 * Permission to use or copy this software for any purpose is hereby granted 
 * without fee, provided the above notices are retained on all copies.
 * Permission to modify the code and to distribute modified code is granted,
 * provided the above notices are retained, and a notice that the code was
 * modified is included with the above copyright notice.
 *
 */

/* NOTE: This may be not portable code. Parts of numeric_limits<> are
 * inherently machine-dependent, and this file is written to make
 * use of system <limits.h>.  Parts of it (in
 * particular, some of the characteristics of floating-point types)
 * may be incorrect for platform that does not implement IEEE standard
 * for floating types.
 */

#ifndef __STLPORT_LIMITS_H
#define __STLPORT_LIMITS_H

#define __SGI_CPP_LIMITS

# ifndef __STL_CONFIG_H
#  include <stl_config.h>
# endif

#ifndef __STLPORT_CLIMITS
# include <climits>
#endif

#ifndef __STLPORT_CFLOAT
# include <cfloat>
#endif

#ifndef __STL_NO_WCHAR_T
# include <cwchar>
# ifndef WCHAR_MIN
#  define WCHAR_MIN 0
#  define WCHAR_MAX ((wchar_t)~0)
# endif
#endif

__STL_BEGIN_NAMESPACE

enum float_round_style {
  round_indeterminate       = -1,
  round_toward_zero         =  0,
  round_to_nearest          =  1,
  round_toward_infinity     =  2,
  round_toward_neg_infinity =  3
};

enum float_denorm_style {
  denorm_indeterminate = -1,
  denorm_absent        =  0,
  denorm_present       =  1
};

// Base class for all specializations of numeric_limits.

template <class __number>
class _Numeric_limits_base {
public:
  static const bool is_specialized __STL_INLINE_STATIC_INIT ( false );

  static __number min() __STL_NOTHROW { return __number(); }
  static __number max() __STL_NOTHROW { return __number(); }

  static const int digits   __STL_INLINE_STATIC_INIT( 0 );
  static const int digits10 __STL_INLINE_STATIC_INIT ( 0 );

  static const bool is_signed  __STL_INLINE_STATIC_INIT ( false );
  static const bool is_integer __STL_INLINE_STATIC_INIT ( false );
  static const bool is_exact   __STL_INLINE_STATIC_INIT ( false );

  static const int radix __STL_INLINE_STATIC_INIT ( 0 );

  static __number epsilon() __STL_NOTHROW     { return __number(); }
  static __number round_error() __STL_NOTHROW { return __number(); }

  static const int min_exponent  __STL_INLINE_STATIC_INIT ( 0 );
  static const int min_exponent10 __STL_INLINE_STATIC_INIT ( 0 );
  static const int max_exponent   __STL_INLINE_STATIC_INIT ( 0 );
  static const int max_exponent10 __STL_INLINE_STATIC_INIT ( 0 );

  static const bool has_infinity      __STL_INLINE_STATIC_INIT ( false );
  static const bool has_quiet_NaN     __STL_INLINE_STATIC_INIT ( false );
  static const bool has_signaling_NaN __STL_INLINE_STATIC_INIT ( false );
  static const float_denorm_style has_denorm __STL_INLINE_STATIC_INIT ( denorm_absent );
  static const bool has_denorm_loss   __STL_INLINE_STATIC_INIT ( false );

  static __number infinity() __STL_NOTHROW      { return __number(); }
  static __number quiet_NaN() __STL_NOTHROW     { return __number(); }
  static __number signaling_NaN() __STL_NOTHROW { return __number(); }
  static __number denorm_min() __STL_NOTHROW    { return __number(); }

  static const bool is_iec559  __STL_INLINE_STATIC_INIT( false );
  static const bool is_bounded __STL_INLINE_STATIC_INIT( false );
  static const bool is_modulo  __STL_INLINE_STATIC_INIT( false );

  static const bool traps           __STL_INLINE_STATIC_INIT( false);
  static const bool tinyness_before __STL_INLINE_STATIC_INIT( false);
  static const float_round_style round_style __STL_INLINE_STATIC_INIT( round_toward_zero );
};

// Base class for integers.

# ifdef __STL_LIMITED_DEFAULT_TEMPLATES
#  ifdef __STL_LONG_LONG
#   define __STL_LIMITS_MIN_TYPE long long
#   define __STL_LIMITS_MAX_TYPE unsigned long long
#  else
#   define __STL_LIMITS_MIN_TYPE long
#   define __STL_LIMITS_MAX_TYPE unsigned long
#  endif
# else
#   define __STL_LIMITS_MIN_TYPE _Int
#   define __STL_LIMITS_MAX_TYPE _Int
# endif /* __STL_LIMITED_DEFAULT_TEMPLATES */

template <class _Int,
          __STL_LIMITS_MIN_TYPE __imin,
          __STL_LIMITS_MAX_TYPE __imax,
          int __idigits>
class _Integer_limits : public _Numeric_limits_base<_Int> 
{
public:
  static const bool is_specialized __STL_INLINE_STATIC_INIT( true );

  static _Int min() __STL_NOTHROW { return (_Int)__imin; }
  static _Int max() __STL_NOTHROW { return (_Int)__imax; }

  static const int digits __STL_INLINE_STATIC_INIT((__idigits < 0) ? \
			    ((int)((sizeof(_Int) * (CHAR_BIT))) - ((__imin == 0) ? 0 : 1)) \
                            : (__idigits));

  static const int digits10  __STL_INLINE_STATIC_INIT( (digits * 301UL) / 1000 ); 
                                // log 2 = 0.301029995664...

  static const bool is_signed __STL_INLINE_STATIC_INIT(  __imin != 0 );
  static const bool is_integer __STL_INLINE_STATIC_INIT(  true );
  static const bool is_exact __STL_INLINE_STATIC_INIT(  true );
  static const int radix __STL_INLINE_STATIC_INIT(  2 );


  static const bool is_bounded __STL_INLINE_STATIC_INIT(  true );
  static const bool is_modulo __STL_INLINE_STATIC_INIT(  true );
};

// Base class for floating-point numbers.
template <class __number,
         int __Digits, int __Digits10,
         int __MinExp, int __MaxExp,
         int __MinExp10, int __MaxExp10,
         __STL_UINT32_T __InfinityWord,
         __STL_UINT32_T __QNaNWord, __STL_UINT32_T __SNaNWord,
         bool __IsIEC559,
         float_round_style __RoundStyle>
class _Floating_limits : public _Numeric_limits_base<__number>
{
public:
  static const bool is_specialized __STL_INLINE_STATIC_INIT( true );

  static const int digits   __STL_INLINE_STATIC_INIT(  __Digits );
  static const int digits10 __STL_INLINE_STATIC_INIT(  __Digits10 );

  static const bool is_signed __STL_INLINE_STATIC_INIT(  true );

  static const int radix __STL_INLINE_STATIC_INIT(  FLT_RADIX /* 2 */ );

  static const int min_exponent   __STL_INLINE_STATIC_INIT(  __MinExp );
  static const int max_exponent   __STL_INLINE_STATIC_INIT(  __MaxExp );
  static const int min_exponent10 __STL_INLINE_STATIC_INIT(  __MinExp10 );
  static const int max_exponent10 __STL_INLINE_STATIC_INIT(  __MaxExp10 );

  static const bool has_infinity      __STL_INLINE_STATIC_INIT(  true );
  static const bool has_quiet_NaN     __STL_INLINE_STATIC_INIT(  true );
  static const bool has_signaling_NaN __STL_INLINE_STATIC_INIT(  true );
  static const float_denorm_style has_denorm __STL_INLINE_STATIC_INIT(  denorm_indeterminate );
  static const bool has_denorm_loss   __STL_INLINE_STATIC_INIT(  false );

  static __number infinity() __STL_NOTHROW {
    static __STL_UINT32_T _S_inf[sizeof(__number) / sizeof(__STL_UINT32_T)] = 
      { __InfinityWord };
    return *__REINTERPRET_CAST(__number*,&_S_inf);
  }
  static __number quiet_NaN() __STL_NOTHROW {
    static __STL_UINT32_T _S_nan[sizeof(__number) / sizeof(__STL_UINT32_T)] = 
      { __QNaNWord };
    return *__REINTERPRET_CAST(__number*,&_S_nan);
  }
  static __number signaling_NaN() __STL_NOTHROW {
    static __STL_UINT32_T _S_nan[sizeof(__number) / sizeof(__STL_UINT32_T)] = 
      { __SNaNWord };
    return *__REINTERPRET_CAST(__number*,&_S_nan);
  }

  static const bool is_iec559       __STL_INLINE_STATIC_INIT(  __IsIEC559 );
  static const bool is_bounded      __STL_INLINE_STATIC_INIT(  true );
  static const bool traps           __STL_INLINE_STATIC_INIT(  true );
  static const bool tinyness_before __STL_INLINE_STATIC_INIT(  false );

  static const float_round_style round_style __STL_INLINE_STATIC_INIT(  __RoundStyle );
};

// Class numeric_limits

// The unspecialized class.

template<class _Tp> 
class numeric_limits : public _Numeric_limits_base<_Tp> {};

// Specializations for all built-in integral types.

#ifndef __STL_NO_BOOL

__STL_TEMPLATE_NULL
class numeric_limits<bool>
  : public _Integer_limits<bool, false, true, 0>
{};

#endif /* __STL_NO_BOOL */

__STL_TEMPLATE_NULL
class numeric_limits<char>
  : public _Integer_limits<char, CHAR_MIN, CHAR_MAX, -1>
{};

# ifndef __STL_NO_SIGNED_BUILTINS
__STL_TEMPLATE_NULL
class numeric_limits<signed char>
  : public _Integer_limits<signed char, SCHAR_MIN, SCHAR_MAX, -1>
{};
# endif

__STL_TEMPLATE_NULL
class numeric_limits<unsigned char>
  : public _Integer_limits<unsigned char, 0, UCHAR_MAX, -1>
{};

#if !(defined ( __STL_NO_WCHAR_T ) || defined (__STL_WCHAR_T_IS_USHORT))

__STL_TEMPLATE_NULL
class numeric_limits<wchar_t>
  : public _Integer_limits<wchar_t, WCHAR_MIN, WCHAR_MAX, -1>
{};

#endif

__STL_TEMPLATE_NULL
class numeric_limits<short>
  : public _Integer_limits<short, SHRT_MIN, SHRT_MAX, -1>
{};

__STL_TEMPLATE_NULL
class numeric_limits<unsigned short>
  : public _Integer_limits<unsigned short, 0, USHRT_MAX, -1>
{};

__STL_TEMPLATE_NULL
class numeric_limits<int>
  : public _Integer_limits<int, INT_MIN, INT_MAX, -1>
{};

__STL_TEMPLATE_NULL
class numeric_limits<unsigned int>
  : public _Integer_limits<unsigned int, 0, UINT_MAX, -1>
{};

__STL_TEMPLATE_NULL
class numeric_limits<long>
  : public _Integer_limits<long, LONG_MIN, LONG_MAX, -1>
{};

__STL_TEMPLATE_NULL
class numeric_limits<unsigned long>
  : public _Integer_limits<unsigned long, 0, ULONG_MAX, -1>
{};

#ifdef __STL_LONG_LONG

#  ifndef   LONGLONG_MAX
#    define LONGLONG_MAX     0x7fffffffffffffffLL
#  endif
#  ifndef   LONGLONG_MIN
#    define LONGLONG_MIN     (-LONGLONG_MAX-1)
#  endif
#  ifndef   ULONGLONG_MAX
#    define ULONGLONG_MAX    0xffffffffffffffffULL
#  endif

__STL_TEMPLATE_NULL
class numeric_limits<long long>
  : public _Integer_limits<long long, LONGLONG_MIN, LONGLONG_MAX, -1>
{};

__STL_TEMPLATE_NULL
class numeric_limits<unsigned long long>
  : public _Integer_limits<unsigned long long, 0, ULONGLONG_MAX, -1>
{};

#endif /* __STL_LONG_LONG */

// Specializations for all built-in floating-point types.

__STL_TEMPLATE_NULL class numeric_limits<float>
  : public _Floating_limits<float, 
                            FLT_MANT_DIG,   // Binary digits of precision
                            FLT_DIG,        // Decimal digits of precision
                            FLT_MIN_EXP,    // Minimum exponent
                            FLT_MAX_EXP,    // Maximum exponent
                            FLT_MIN_10_EXP, // Minimum base 10 exponent
                            FLT_MAX_10_EXP, // Maximum base 10 exponent
                            0x7f800000ul,    // First word of +infinity
                            0x7f810000ul,    // First word of quiet NaN
                            0x7fc10000ul,    // First word of signaling NaN
                            true,           // conforms to iec559
                            round_to_nearest>
{
public:
  static float min() __STL_NOTHROW { return FLT_MIN; }
  static float denorm_min() __STL_NOTHROW { return FLT_MIN; }
  static float max() __STL_NOTHROW { __STL_USING_VENDOR_CSTD return FLT_MAX; }
  static float epsilon() __STL_NOTHROW { return FLT_EPSILON; }
  static float round_error() __STL_NOTHROW { return 0.5f; } // Units: ulps.
};

__STL_TEMPLATE_NULL class numeric_limits<double>
  : public _Floating_limits<double, 
                            DBL_MANT_DIG,   // Binary digits of precision
                            DBL_DIG,        // Decimal digits of precision
                            DBL_MIN_EXP,    // Minimum exponent
                            DBL_MAX_EXP,    // Maximum exponent
                            DBL_MIN_10_EXP, // Minimum base 10 exponent
                            DBL_MAX_10_EXP, // Maximum base 10 exponent
                            0x7ff00000ul,    // First word of +infinity
                            0x7ff10000ul,    // First word of quiet NaN
                            0x7ff90000ul,    // First word of signaling NaN
                            true,           // conforms to iec559
                            round_to_nearest>
{
public:
  static double min() __STL_NOTHROW { return DBL_MIN; }
  static double denorm_min() __STL_NOTHROW { return DBL_MIN; }
  static double max() __STL_NOTHROW { __STL_USING_VENDOR_CSTD return DBL_MAX; }
  static double epsilon() __STL_NOTHROW { return DBL_EPSILON; }
  static double round_error() __STL_NOTHROW { return 0.5; } // Units: ulps.
};

# ifndef __STL_NO_LONG_DOUBLE

__STL_TEMPLATE_NULL 
class numeric_limits<long double>
  : public _Floating_limits<long double, 
                            LDBL_MANT_DIG,  // Binary digits of precision
                            LDBL_DIG,       // Decimal digits of precision
                            LDBL_MIN_EXP,   // Minimum exponent
                            LDBL_MAX_EXP,   // Maximum exponent
                            LDBL_MIN_10_EXP,// Minimum base 10 exponent
                            LDBL_MAX_10_EXP,// Maximum base 10 exponent
                            0x7ff00000ul,    // First word of +infinity
                            0x7ff10000ul,    // First word of quiet NaN
                            0x7ff90000ul,    // First word of signaling NaN
                            false,          // Doesn't conform to iec559
                            round_to_nearest>
{
public:
  static long double min() __STL_NOTHROW { __STL_USING_VENDOR_CSTD return LDBL_MIN; }
  static long double denorm_min() __STL_NOTHROW { __STL_USING_VENDOR_CSTD return LDBL_MIN; }
  static long double max() __STL_NOTHROW { __STL_USING_VENDOR_CSTD return LDBL_MAX; }
  static long double epsilon() __STL_NOTHROW { return LDBL_EPSILON; }
  static long double round_error() __STL_NOTHROW { return 4; } // Units: ulps.
};

# endif


__STL_END_NAMESPACE

# if !defined (__STL_LINK_TIME_INSTANTIATION)
#  include <stl_limits.c>
# endif

#endif /* __SGI_CPP_LIMITS */

// Local Variables:
// mode:C++
// End:
