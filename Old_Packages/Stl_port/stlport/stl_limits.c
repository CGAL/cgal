/*
 * Copyright (c) 1998,1999
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

# ifndef __STLPORT_LIMITS_C
#  define __STLPORT_LIMITS_C

//==========================================================
//  numeric_limits static members
//==========================================================

__STL_BEGIN_NAMESPACE

# ifdef __STL_NO_STATIC_TEMPLATE_DATA

// gcc seems to allow inline initialization, and
// it is the only compiler without static members, so leave it alone

# define __declare_numeric_base_member(__number, __type, __mem) \
__DECLARE_INSTANCE(const __type,_Numeric_limits_base<__number>:: __mem)

# define __declare_numeric_base_members(__number) \
__declare_numeric_base_member(__number, bool, is_specialized);\
__declare_numeric_base_member(__number, int, digits);\
__declare_numeric_base_member(__number, int, digits10);\
__declare_numeric_base_member(__number, bool, is_signed);\
__declare_numeric_base_member(__number, bool, is_integer);\
__declare_numeric_base_member(__number, bool, is_exact);\
__declare_numeric_base_member(__number, int, radix);\
__declare_numeric_base_member(__number, int, min_exponent);\
__declare_numeric_base_member(__number, int, max_exponent);\
__declare_numeric_base_member(__number, int, min_exponent10);\
__declare_numeric_base_member(__number, int, max_exponent10);\
__declare_numeric_base_member(__number, bool, has_infinity);\
__declare_numeric_base_member(__number, bool, has_quiet_NaN);\
__declare_numeric_base_member(__number, bool, has_signaling_NaN);\
__declare_numeric_base_member(__number, float_denorm_style, has_denorm);\
__declare_numeric_base_member(__number, bool, has_denorm_loss);\
__declare_numeric_base_member(__number, bool, is_iec559);\
__declare_numeric_base_member(__number, bool, is_bounded);\
__declare_numeric_base_member(__number, bool, is_modulo);\
__declare_numeric_base_member(__number, bool, traps);\
__declare_numeric_base_member(__number, bool, tinyness_before);\
__declare_numeric_base_member(__number, float_round_style, round_style);

# else
#define __declare_numeric_base_member(__type, __mem, __init) \
template <class __number> \
  const __type _Numeric_limits_base<__number>:: __mem __STL_OUTLINE_STATIC_INIT(__init)

__declare_numeric_base_member(bool, is_specialized, false);
__declare_numeric_base_member(int, digits, 0);
__declare_numeric_base_member(int, digits10, 0);
__declare_numeric_base_member(bool, is_signed, false);
__declare_numeric_base_member(bool, is_integer, false);
__declare_numeric_base_member(bool, is_exact, false);
__declare_numeric_base_member(int, radix, 0);
__declare_numeric_base_member(int, min_exponent, 0);
__declare_numeric_base_member(int, max_exponent, 0);
__declare_numeric_base_member(int, min_exponent10, 0);
__declare_numeric_base_member(int, max_exponent10, 0);
__declare_numeric_base_member(bool, has_infinity, false);
__declare_numeric_base_member(bool, has_quiet_NaN, false);
__declare_numeric_base_member(bool, has_signaling_NaN, false);
__declare_numeric_base_member(float_denorm_style, has_denorm, denorm_absent);
__declare_numeric_base_member(bool, has_denorm_loss, false);
__declare_numeric_base_member(bool, is_iec559, false);
__declare_numeric_base_member(bool, is_bounded, false);
__declare_numeric_base_member(bool, is_modulo, false);
__declare_numeric_base_member(bool, traps, false);
__declare_numeric_base_member(bool, tinyness_before, false);
__declare_numeric_base_member(float_round_style, round_style, round_toward_zero);
# endif

#undef __declare_numeric_base_member

# ifdef __STL_NO_STATIC_TEMPLATE_DATA

# define __HACK_ILIMITS(_Int, __imin, __imax, __idigits) _Integer_limits<_Int, __imin, __imax, __idigits>
# define __HACK_NOTHING
# define __declare_integer_limits_member(_Int, __imin, __imax, __idigits, __type, __mem) \
  __DECLARE_INSTANCE(const __type, __HACK_ILIMITS(_Int, __imin, __imax, __idigits):: __mem,__HACK_NOTHING)

# define __declare_int_members(_Int, __imin, __imax, __idigits) \
__declare_integer_limits_member(_Int, __imin, __imax, __idigits, bool, is_specialized);\
__declare_integer_limits_member(_Int, __imin, __imax, __idigits, int, digits);\
__declare_integer_limits_member(_Int, __imin, __imax, __idigits, int, digits10);\
__declare_integer_limits_member(_Int, __imin, __imax, __idigits, bool, is_signed);\
__declare_integer_limits_member(_Int, __imin, __imax, __idigits, bool, is_integer);\
__declare_integer_limits_member(_Int, __imin, __imax, __idigits, bool, is_exact);\
__declare_integer_limits_member(_Int, __imin, __imax, __idigits, int, radix);\
__declare_integer_limits_member(_Int, __imin, __imax, __idigits, bool, is_bounded);\
__declare_integer_limits_member(_Int, __imin, __imax, __idigits, bool, is_modulo);


# else
#define __declare_integer_limits_member(__type, __mem, __init) \
template <class _Int, __STL_LIMITS_MIN_TYPE __imin, __STL_LIMITS_MAX_TYPE __imax, int __idigits> \
  const __type _Integer_limits<_Int, __imin, __imax, __idigits>:: __mem __STL_OUTLINE_STATIC_INIT(__init)

__declare_integer_limits_member(bool, is_specialized, true);
__declare_integer_limits_member(int, digits, (__idigits < 0) ? \
			    ((int)((sizeof(_Int) * (CHAR_BIT))) - ((__imin == 0) ? 0 : 1)) \
                            : (__idigits) );
__declare_integer_limits_member(int, digits10, (digits * 301UL) / 1000);
__declare_integer_limits_member(bool, is_signed, __imin != 0);
__declare_integer_limits_member(bool, is_integer, true);
__declare_integer_limits_member(bool, is_exact, true);
__declare_integer_limits_member(int, radix, 2);
__declare_integer_limits_member(bool, is_bounded, true);
__declare_integer_limits_member(bool, is_modulo, true);

#undef __declare_integer_limits_member

# endif

# ifdef __STL_NO_STATIC_TEMPLATE_DATA
# define __HACK_FLIMITS(__Num, __Dig, __Dig10,__MnX, __MxX, __MnX10, __MxX10, __Inf, __QNaN, __SNaN,__IsIEEE, __Sty) \
_Floating_limits<__Num, __Dig, __Dig10, \
                              __MnX, __MxX, __MnX10, __MxX10, \
                              __Inf, __QNaN, __SNaN,__IsIEEE, __Sty>

#define __declare_float_limits_member(__Num, __Dig, __Dig10,__MnX, __MxX, __MnX10, __MxX10, __Inf, __QNaN, __SNaN,__IsIEEE, __Sty, __type, __mem) \
 __DECLARE_INSTANCE(const __type, __HACK_FLIMITS(__Num, __Dig, __Dig10, \
                              __MnX, __MxX, __MnX10, __MxX10, \
                              __Inf, __QNaN, __SNaN,__IsIEEE, __Sty)::__mem, __HACK_NOTHING)

#define __declare_float_members(__Num, __Dig, __Dig10,__MnX, __MxX, __MnX10, __MxX10, __Inf, __QNaN, __SNaN,__IsIEEE, __Sty) \
__declare_float_limits_member(__Num, __Dig, __Dig10,__MnX, __MxX, __MnX10, __MxX10, __Inf, __QNaN, __SNaN,__IsIEEE, __Sty, bool, is_specialized);\
__declare_float_limits_member(__Num, __Dig, __Dig10,__MnX, __MxX, __MnX10, __MxX10, __Inf, __QNaN, __SNaN,__IsIEEE, __Sty, int, digits);\
__declare_float_limits_member(__Num, __Dig, __Dig10,__MnX, __MxX, __MnX10, __MxX10, __Inf, __QNaN, __SNaN,__IsIEEE, __Sty, int, digits10);\
__declare_float_limits_member(__Num, __Dig, __Dig10,__MnX, __MxX, __MnX10, __MxX10, __Inf, __QNaN, __SNaN,__IsIEEE, __Sty, bool, is_signed);\
__declare_float_limits_member(__Num, __Dig, __Dig10,__MnX, __MxX, __MnX10, __MxX10, __Inf, __QNaN, __SNaN,__IsIEEE, __Sty, int, radix);\
__declare_float_limits_member(__Num, __Dig, __Dig10,__MnX, __MxX, __MnX10, __MxX10, __Inf, __QNaN, __SNaN,__IsIEEE, __Sty, int, min_exponent);\
__declare_float_limits_member(__Num, __Dig, __Dig10,__MnX, __MxX, __MnX10, __MxX10, __Inf, __QNaN, __SNaN,__IsIEEE, __Sty, int, max_exponent);\
__declare_float_limits_member(__Num, __Dig, __Dig10,__MnX, __MxX, __MnX10, __MxX10, __Inf, __QNaN, __SNaN,__IsIEEE, __Sty, int, min_exponent10);\
__declare_float_limits_member(__Num, __Dig, __Dig10,__MnX, __MxX, __MnX10, __MxX10, __Inf, __QNaN, __SNaN,__IsIEEE, __Sty, int, max_exponent10);\
__declare_float_limits_member(__Num, __Dig, __Dig10,__MnX, __MxX, __MnX10, __MxX10, __Inf, __QNaN, __SNaN,__IsIEEE, __Sty, bool, has_infinity);\
__declare_float_limits_member(__Num, __Dig, __Dig10,__MnX, __MxX, __MnX10, __MxX10, __Inf, __QNaN, __SNaN,__IsIEEE, __Sty, bool, has_quiet_NaN);\
__declare_float_limits_member(__Num, __Dig, __Dig10,__MnX, __MxX, __MnX10, __MxX10, __Inf, __QNaN, __SNaN,__IsIEEE, __Sty, bool, has_signaling_NaN);\
__declare_float_limits_member(__Num, __Dig, __Dig10,__MnX, __MxX, __MnX10, __MxX10, __Inf, __QNaN, __SNaN,__IsIEEE, __Sty, float_denorm_style, has_denorm);\
__declare_float_limits_member(__Num, __Dig, __Dig10,__MnX, __MxX, __MnX10, __MxX10, __Inf, __QNaN, __SNaN,__IsIEEE, __Sty, bool, has_denorm_loss);\
__declare_float_limits_member(__Num, __Dig, __Dig10,__MnX, __MxX, __MnX10, __MxX10, __Inf, __QNaN, __SNaN,__IsIEEE, __Sty, bool, is_iec559);\
__declare_float_limits_member(__Num, __Dig, __Dig10,__MnX, __MxX, __MnX10, __MxX10, __Inf, __QNaN, __SNaN,__IsIEEE, __Sty, bool, is_bounded);\
__declare_float_limits_member(__Num, __Dig, __Dig10,__MnX, __MxX, __MnX10, __MxX10, __Inf, __QNaN, __SNaN,__IsIEEE, __Sty, bool, traps);\
__declare_float_limits_member(__Num, __Dig, __Dig10,__MnX, __MxX, __MnX10, __MxX10, __Inf, __QNaN, __SNaN,__IsIEEE, __Sty, bool, tinyness_before);\
__declare_float_limits_member(__Num, __Dig, __Dig10,__MnX, __MxX, __MnX10, __MxX10, __Inf, __QNaN, __SNaN,__IsIEEE, __Sty, float_round_style, round_style);

# else

#define __declare_float_limits_member(__type, __mem, __init) \
template <class __number,  \
         int __Digits, int __Digits10,    \
         int __MinExp, int __MaxExp,      \
         int __MinExp10, int __MaxExp10,  \
         __STL_UINT32_T __InfinityWord,    \
         __STL_UINT32_T __QNaNWord, __STL_UINT32_T __SNaNWord, \
         bool __IsIEC559, \
         float_round_style __RoundStyle> \
const __type _Floating_limits< __number, __Digits, __Digits10,    \
         __MinExp, __MaxExp, __MinExp10, __MaxExp10,  \
         __InfinityWord,__QNaNWord, __SNaNWord, __IsIEC559, __RoundStyle>::\
         __mem __STL_OUTLINE_STATIC_INIT(__init)

__declare_float_limits_member(bool, is_specialized, true);  
__declare_float_limits_member(int, digits, __Digits);  
__declare_float_limits_member(int, digits10, __Digits10);  
__declare_float_limits_member(bool, is_signed, true);  
__declare_float_limits_member(int, radix, FLT_RADIX);  
__declare_float_limits_member(int, min_exponent, __MinExp);  
__declare_float_limits_member(int, max_exponent, __MaxExp);  
__declare_float_limits_member(int, min_exponent10, __MinExp10);  
__declare_float_limits_member(int, max_exponent10, __MaxExp10);  
__declare_float_limits_member(bool, has_infinity, true);
__declare_float_limits_member(bool, has_quiet_NaN, true);
__declare_float_limits_member(bool, has_signaling_NaN, true);
__declare_float_limits_member(float_denorm_style, has_denorm, denorm_indeterminate);
__declare_float_limits_member(bool, has_denorm_loss, false);
__declare_float_limits_member(bool, is_iec559, __IsIEC559);
__declare_float_limits_member(bool, is_bounded, true);
__declare_float_limits_member(bool, traps, true);
__declare_float_limits_member(bool, tinyness_before, false);
__declare_float_limits_member(float_round_style, round_style, __RoundStyle);

#  undef __declare_float_limits_member

# endif

# ifdef __STL_NO_STATIC_TEMPLATE_DATA

#ifndef __STL_NO_BOOL
__declare_int_members(bool, false, true, 0)
#endif

__declare_int_members(char, CHAR_MIN, CHAR_MAX, -1)

# ifndef __STL_NO_SIGNED_BUILTINS
__declare_int_members(signed char, SCHAR_MIN, SCHAR_MAX, -1)
# endif

__declare_int_members(unsigned char, 0, UCHAR_MAX, -1)

#if defined (__STL_HAS_WCHAR_T) && !defined ( __STL_WCHAR_T_IS_USHORT)
__declare_int_members(wchar_t, INT_MIN, INT_MAX, -1)
#endif

__declare_int_members(short, SHRT_MIN, SHRT_MAX, -1)
__declare_int_members(unsigned short, 0, USHRT_MAX, -1)
__declare_int_members(int, INT_MIN, INT_MAX, -1)
__declare_int_members(unsigned int, 0, UINT_MAX, -1)
__declare_int_members(long, LONG_MIN, LONG_MAX, -1)
__declare_int_members(unsigned long, 0, ULONG_MAX, -1)

#ifdef __STL_LONG_LONG
__declare_int_members(long long, LONGLONG_MIN, LONG_MAX, -1)
__declare_int_members(unsigned long long, 0, ULONGLONG_MAX, -1)
#endif

__declare_float_members(float, FLT_MANT_DIG,FLT_DIG,
                  	    FLT_MIN_EXP,
                            FLT_MAX_EXP,
                            FLT_MIN_10_EXP,
                            FLT_MAX_10_EXP,
                            0x7f800000u,
                            0x7f810000u,
                            0x7fc10000u,
                            true,
                            round_to_nearest)

__declare_float_members(double, 
                            DBL_MANT_DIG,
                            DBL_DIG,
                            DBL_MIN_EXP,
                            DBL_MAX_EXP,
                            DBL_MIN_10_EXP,
                            DBL_MAX_10_EXP,
                            0x7ff00000u,
                            0x7ff10000u,
                            0x7ff90000u,
                            true,
                            round_to_nearest)

# if !defined (__STL_NO_LONG_DOUBLE)
__declare_float_members(long double, 
                            LDBL_MANT_DIG,
                            LDBL_DIG,
                            LDBL_MIN_EXP,
                            LDBL_MAX_EXP,
                            LDBL_MIN_10_EXP,
                            LDBL_MAX_10_EXP,
                            0x7ff00000u,
                            0x7ff10000u,
                            0x7ff90000u,
                            false,
                            round_to_nearest)
# endif

# endif


# undef __HACK_ILIMITS
# undef __HACK_NOTHING
# undef __declare_int_members
# undef __declare_float_members
# undef __STL_LIMITS_MIN_TYPE
# undef __STL_LIMITS_MAX_TYPE
__STL_END_NAMESPACE

#endif /* __STLPORT_LIMITS_C_INCLUDED */
