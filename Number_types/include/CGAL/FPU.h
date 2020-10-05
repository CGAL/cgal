// Copyright (c) 1998-2019
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sylvain Pion, Marc Glisse

#ifndef CGAL_FPU_H
#define CGAL_FPU_H

#include <CGAL/assertions.h>
#include <CGAL/use.h>
#include <cstring> // std::memcpy

#ifndef __INTEL_COMPILER
#include <cmath> // for HUGE_VAL
#endif

// This file specifies some platform dependant functions, regarding the FPU
// directed rounding modes.  There is only support for double precision.
//
// It also contains the definition of the Protect_FPU_rounding<> class,
// which helps to protect blocks of code needing a particular rounding mode.

#if defined __alpha__  && defined __linux__
extern "C" {
#  include <fenv.h>
}
#elif defined __SUNPRO_CC && defined __sun
#  include <ieeefp.h>
#elif defined __osf || defined __osf__
#  ifdef __GNUG__
     // GCC seems to remove (fixincludes) read_rnd/write_rnd...
#    include "/usr/include/float.h"
#  else
#    include <cfloat>
#  endif
#elif defined _MSC_VER || defined __sparc__ || \
     (defined __i386__ && !defined __PGI && !defined __SUNPRO_CC \
      && !defined __SSE2__)
   // Nothing to include.
#else
   // By default we use the ISO C99 version.
#  include <fenv.h>
#endif


// Some useful constants

#if defined CGAL_CFG_NO_LIMITS
#  if defined CGAL_CFG_DENORMALS_COMPILE_BUG
     // For compilers crashing when dealing with denormalized values.
     // So we have to generate it at run time instead.
#    define CGAL_IA_MIN_DOUBLE (CGAL::internal::get_static_minimin())
#  else
#    define CGAL_IA_MIN_DOUBLE (5e-324)
#  endif
#  define CGAL_IA_MAX_DOUBLE (1.7976931348623157081e+308)
#else
#  include <limits>
#  define CGAL_IA_MIN_DOUBLE std::numeric_limits<double>::denorm_min()
#  define CGAL_IA_MAX_DOUBLE (std::numeric_limits<double>::max)()
#endif


// Pure and safe SSE2 mode (g++ -mfpmath=sse && (-msse2 || -march=pentium4))
// can be detected by :
// TODO : see what Intel and VC++ have to say about this.
#if defined _M_X64 || ( defined __FLT_EVAL_METHOD__ && defined __SSE2_MATH__ && \
                        (__FLT_EVAL_METHOD__ == 0 || __FLT_EVAL_METHOD__ == 1) )
#  define CGAL_SAFE_SSE2
#  include <xmmintrin.h>
#endif

// The CGAL_FPU_HAS_EXCESS_PRECISION macro is defined if some computations with
// double can use more than the 53bits of precision of IEEE754, and/or if the
// exponent has a wider range.  This can produce double rounding effects and
// other bad things that we need to protect against.
// The typical offender is the traditional FPU of x86 (SSE2-only mode is not affected).
// Are there others, besides itanium and m68k?
#if !defined CGAL_IA_NO_X86_OVER_UNDER_FLOW_PROTECT && \
  (((defined __i386__ || defined __x86_64__) && !defined CGAL_SAFE_SSE2) \
   || defined __ia64__ \
   || defined _M_IX86 || defined _M_IA64 \
   || (defined FLT_EVAL_METHOD && FLT_EVAL_METHOD != 0 && FLT_EVAL_METHOD != 1))
#  define CGAL_FPU_HAS_EXCESS_PRECISION
#endif

// Presence of SSE2 (for explicit use)
#if  defined(__SSE2__) \
  || (defined(_M_IX86_FP) && _M_IX86_FP >= 2) \
  || defined(_M_X64)
#  include <emmintrin.h>
#  if defined __SSE3__
#    include <pmmintrin.h>
#  endif
#  if defined __SSE4_1__
#    include <smmintrin.h>
#  endif
#  if defined __AVX__
#    include <immintrin.h>
#  endif
#  define CGAL_HAS_SSE2 1
#endif

// Only define CGAL_USE_SSE2 for 64 bits where malloc has a suitable
// alignment, 32 bits is too dangerous.
#if defined(CGAL_HAS_SSE2) && (defined(__x86_64__) || defined(_M_X64))
#  define CGAL_USE_SSE2 1
#endif

#ifdef CGAL_CFG_DENORMALS_COMPILE_BUG
double& get_static_minimin(); // Defined in Interval_arithmetic_impl.h
#endif

namespace  CGAL {

#ifdef CGAL_HEADER_ONLY
// Defined in test_FPU_rounding_mode_impl.h
struct Check_FPU_rounding_mode_is_restored;
inline const Check_FPU_rounding_mode_is_restored&
get_static_check_fpu_rounding_mode_is_restored();
#endif

// Inline function to stop compiler optimizations that shouldn't happen with
// pragma fenv on.
// - constant propagation
// - migration of fesetround across floating point operations
// - (-a)-b -> -(a+b)
// - (-a)*b -> -(a*b)
// etc
inline double IA_opacify(double x)
{
#ifdef __llvm__
  // LLVM's support for inline asm is completely messed up:
  // http://llvm.org/bugs/show_bug.cgi?id=17958
  // http://llvm.org/bugs/show_bug.cgi?id=17959
  // etc.
  // This seems to produce code that is ok (not optimal but better than
  // volatile). In case of trouble, use volatile instead.
# ifdef CGAL_HAS_SSE2
  asm volatile ("" : "+x"(x) );
# elif (defined __VFP_FP__ && !defined __SOFTFP__) || defined __aarch64__
  // ARM
  asm volatile ("" : "+w"(x) );
# else
  asm volatile ("" : "+m"(x) );
# endif
  return x;
#elif defined __xlC__
  // PowerPC - XL C++ (the z/OS version supposedly does not define this macro)
  // If we give it an alternative "+fm", it gets confused and generates worse code.
  asm volatile ("" : "+f"(x) );
  return x;
#elif defined __GNUG__
  // Intel used not to emulate this perfectly, we'll see.
  // If we create a version of IA_opacify for vectors, note that gcc < 4.8
  // fails with "+g" and we need to use "+mx" instead.
  // "+X" ICEs ( http://gcc.gnu.org/bugzilla/show_bug.cgi?id=59155 ) and
  // may not be safe?
  // The constraint 'g' doesn't include floating point registers ???
  // Intel has a bug where -mno-sse still defines __SSE__ and __SSE2__
  // (-mno-sse2 works though), no work-around for now.
# if defined __SSE2_MATH__ || (defined __INTEL_COMPILER && defined __SSE2__)
#  if __GNUC__ * 100 + __GNUC_MINOR__ >= 409
  // ICEs in reload/LRA with older versions.
  asm volatile ("" : "+gx"(x) );
#  else
  asm volatile ("" : "+mx"(x) );
#  endif
# elif (defined __i386__ || defined __x86_64__)
  // "+f" doesn't compile on x86(_64)
  // ( http://gcc.gnu.org/bugzilla/show_bug.cgi?id=59157 )
  // Don't mix "t" with "g": http://gcc.gnu.org/bugzilla/show_bug.cgi?id=59180
  // We can't put "t" with "x" either, prefer "x" for -mfpmath=sse,387.
  // ( http://gcc.gnu.org/bugzilla/show_bug.cgi?id=59181 )
  asm volatile ("" : "+mt"(x) );
# elif (defined __VFP_FP__ && !defined __SOFTFP__) || defined __aarch64__
  // ARM
  asm volatile ("" : "+gw"(x) );
# elif defined __powerpc__ || defined __POWERPC__
  // PowerPC
  asm volatile ("" : "+gd"(x) );
# elif defined __sparc
  // Sparc
  asm volatile ("" : "+ge"(x) );
# elif defined __ia64
  // Itanium
  asm volatile ("" : "+gf"(x) );
# else
  asm volatile ("" : "+g"(x) );
# endif
  return x;
#else
  volatile double e = x;
  return e;
#endif
}

// Inline function to drop excess precision before we forget the rounding mode,
// and stop compiler optimizations at the same time.
inline double IA_force_to_double(double x)
{
#ifndef CGAL_FPU_HAS_EXCESS_PRECISION
  return IA_opacify (x);
#else
#if defined __GNUG__
#  ifdef CGAL_HAS_SSE2
  // For an explanation of volatile:
  // http://gcc.gnu.org/bugzilla/show_bug.cgi?id=56027
  asm volatile ("" : "+mx"(x) );
#  else
  // Similar to writing to a volatile and reading back, except that calling
  // it k times in a row only goes through memory once.
  asm volatile ("" : "+m"(x) );
#  endif
  return x;
#else
  volatile double e = x;
  return e;
#endif
#endif
}

#ifdef CGAL_USE_SSE2
// Vector version of IA_opacify
inline __m128d IA_opacify128(__m128d x)
{
# ifdef __GNUG__
#  ifdef __llvm__
  asm volatile ("":"+x"(x));
#  else
  asm volatile ("":"+xm"(x));
#  endif
  return x;
# else
  volatile __m128d e = x;
#  ifdef _MSC_VER
  // With VS, __m128d is a union, where volatile doesn't disappear automatically
  // However, this version generates wrong code with clang, check before enabling it for more compilers.
  std::memcpy(&x, (void*)&e, 16);
  return x;
#  else
  return e;
#  endif
# endif
}

// Weaker version. It still blocks transformations like -(a-b) to b-a, but does not prevent migrating across fesetround. When an operation has 2 arguments and one uses IA_opacify128, the other one can stick to IA_opacify128_weak
inline __m128d IA_opacify128_weak(__m128d x)
{
# ifdef __GNUG__
#  ifdef __llvm__
  asm ("":"+x"(x));
#  else
  asm ("":"+xm"(x));
#  endif
  return x;
# else
  return IA_opacify128(x);
# endif
}

// _mm_shuffle_pd would work everywhere, but it is too opaque for some optimizations (yes, this is a thin line)
inline __m128d swap_m128d(__m128d x){
# ifdef __llvm__
  return __builtin_shufflevector(x, x, 1, 0);
# elif defined __GNUC__ && !defined __INTEL_COMPILER \
  && __GNUC__ * 100 + __GNUC_MINOR__ >= 407
  return __builtin_shuffle(x, (__m128i){ 1, 0 });
# else
  return _mm_shuffle_pd(x, x, 1);
# endif
}
#endif

// Interval arithmetic needs to protect against double-rounding effects
// caused by excess FPU precision, even if it forces the 53bit mantissa
// precision, because there is no way to fix the problem for the exponent
// which has the same problem.  This affects underflow and overflow cases.
// In case one does not care about such "extreme" situations, one can
// set CGAL_IA_NO_X86_OVER_UNDER_FLOW_PROTECT to pretend there is no excess
// precision.
#if defined CGAL_FPU_HAS_EXCESS_PRECISION
#  define CGAL_IA_FORCE_TO_DOUBLE(x) CGAL::IA_force_to_double(x)
#elif 1
// LLVM doesn't have -frounding-math so needs extra protection.
// GCC also migrates fesetround calls over FP instructions, so protect
// everyone.
#  define CGAL_IA_FORCE_TO_DOUBLE(x) CGAL::IA_opacify(x)
#else
// Unused, reserved to compilers without excess precision and pragma fenv on.
// ??? Should we trust Visual Studio not to optimize too much and let it use
// this when CGAL_IA_NO_X86_OVER_UNDER_FLOW_PROTECT?
#  define CGAL_IA_FORCE_TO_DOUBLE(x) (x)
#endif

// We sometimes need to stop constant propagation,
// because operations are done with a wrong rounding mode at compile time.
#ifndef CGAL_IA_DONT_STOP_CONSTANT_PROPAGATION
#  define CGAL_IA_STOP_CPROP(x)    CGAL::IA_opacify(x)
#else
#  define CGAL_IA_STOP_CPROP(x)    (x)
#endif

// std::sqrt(double) on VC++ and CygWin is buggy when not optimizing.
#if defined ( _MSC_VER ) && ! defined ( _WIN64 )
inline double IA_bug_sqrt(double d)
{
  _asm
  {
    fld d
    fsqrt
    fstp d
  }
  return d;
}

#  define CGAL_BUG_SQRT(d) IA_bug_sqrt(d)


#elif defined __SSE2_MATH__
// For SSE2, we need to call __builtin_sqrt() instead of libc's sqrt().
#  define CGAL_BUG_SQRT(d) __builtin_sqrt(d)
#elif defined __CYGWIN__
inline double IA_bug_sqrt(double d)
{
  double r;
  asm volatile ("fsqrt" : "=t"(r) : "0"(d));
  return r;
}
#  define CGAL_BUG_SQRT(d) IA_bug_sqrt(d)
#else
#  define CGAL_BUG_SQRT(d) std::sqrt(d)
#endif

// Here are the operator macros that make use of the above.
// With GCC, we can do slightly better : test with __builtin_constant_p()
// that both arguments are constant before stopping one of them.
// Use inline functions instead ?
#define CGAL_IA_ADD(a,b) CGAL_IA_FORCE_TO_DOUBLE((a)+CGAL_IA_STOP_CPROP(b))
#define CGAL_IA_SUB(a,b) CGAL_IA_FORCE_TO_DOUBLE(CGAL_IA_STOP_CPROP(a)-(b))
#define CGAL_IA_MUL(a,b) CGAL_IA_FORCE_TO_DOUBLE(CGAL_IA_STOP_CPROP(a)*CGAL_IA_STOP_CPROP(b))
#define CGAL_IA_DIV(a,b) CGAL_IA_FORCE_TO_DOUBLE(CGAL_IA_STOP_CPROP(a)/CGAL_IA_STOP_CPROP(b))
inline double CGAL_IA_SQUARE(double a){
  double b = CGAL_IA_STOP_CPROP(a); // only once
  return CGAL_IA_FORCE_TO_DOUBLE(b*b);
}
#define CGAL_IA_SQRT(a) \
        CGAL_IA_FORCE_TO_DOUBLE(CGAL_BUG_SQRT(CGAL_IA_STOP_CPROP(a)))


#if defined CGAL_SAFE_SSE2

#define CGAL_IA_SETFPCW(CW) _MM_SET_ROUNDING_MODE(CW)
#define CGAL_IA_GETFPCW(CW) CW = _MM_GET_ROUNDING_MODE()
typedef unsigned int FPU_CW_t;
#define CGAL_FE_TONEAREST    _MM_ROUND_NEAREST
#define CGAL_FE_TOWARDZERO   _MM_ROUND_TOWARD_ZERO
#define CGAL_FE_UPWARD       _MM_ROUND_UP
#define CGAL_FE_DOWNWARD     _MM_ROUND_DOWN

#elif defined __i386__ && !defined __PGI && !defined __SUNPRO_CC \
      && !defined CGAL_HAS_SSE2
// If we use both 387 and sse2, be safe and drop to fe[gs]etround.
// Can we test CGAL_USE_SSE2 instead?

// The GNU libc version (cf powerpc) is nicer, but doesn't work on libc 5 :(
// This one also works with CygWin.
// Note that the ISO C99 version may not be enough because of the extended
// mantissa issue on x86 (may be required by some kinds of computation, but
// as far as CGAL::Interval_nt<> is concerned, the double-rounding issues
// are taking care of there).
#define CGAL_IA_SETFPCW(CW) asm volatile ("fldcw %0" : :"m" (CW))
#define CGAL_IA_GETFPCW(CW) asm volatile ("fnstcw %0" : "=m" (CW))
typedef unsigned short FPU_CW_t;
#define CGAL_FE_TONEAREST    (0x000 | 0x127f)
#define CGAL_FE_TOWARDZERO   (0xc00 | 0x127f)
#define CGAL_FE_UPWARD       (0x800 | 0x127f)
#define CGAL_FE_DOWNWARD     (0x400 | 0x127f)

#elif defined __SUNPRO_CC && defined __sun
#define CGAL_IA_SETFPCW(CW) fpsetround(fp_rnd(CW))
#define CGAL_IA_GETFPCW(CW) CW = fpgetround()
typedef unsigned int FPU_CW_t;
#define CGAL_FE_TONEAREST    FP_RN
#define CGAL_FE_TOWARDZERO   FP_RZ
#define CGAL_FE_UPWARD       FP_RP
#define CGAL_FE_DOWNWARD     FP_RM

#elif defined __sparc__
#define CGAL_IA_SETFPCW(CW) asm volatile ("ld %0,%%fsr" : :"m" (CW))
#define CGAL_IA_GETFPCW(CW) asm volatile ("st %%fsr,%0" : "=m" (CW))
typedef unsigned int FPU_CW_t;
#define CGAL_FE_TONEAREST    (0x0        | 0x20000000 | 0x1f)
#define CGAL_FE_TOWARDZERO   (0x40000000 | 0x20000000 | 0x1f)
#define CGAL_FE_UPWARD       (0x80000000 | 0x20000000 | 0x1f)
#define CGAL_FE_DOWNWARD     (0xc0000000 | 0x20000000 | 0x1f)

#elif defined __mips__
#define CGAL_IA_SETFPCW(CW) asm volatile ("ctc1 %0,$31" : :"r" (CW))
#define CGAL_IA_GETFPCW(CW) asm volatile ("cfc1 %0,$31" : "=r" (CW)); CW &= 3
typedef unsigned int FPU_CW_t;
#define CGAL_FE_TONEAREST    (0x0)
#define CGAL_FE_TOWARDZERO   (0x1)
#define CGAL_FE_UPWARD       (0x2)
#define CGAL_FE_DOWNWARD     (0x3)

#elif defined __osf || defined __osf__  // Not yet supported.
#define CGAL_IA_SETFPCW(CW) write_rnd(CW)
#define CGAL_IA_GETFPCW(CW) CW = read_rnd()
typedef unsigned int FPU_CW_t;
#define CGAL_FE_TONEAREST    FP_RND_RN
#define CGAL_FE_TOWARDZERO   FP_RND_RZ
#define CGAL_FE_UPWARD       FP_RND_RP
#define CGAL_FE_DOWNWARD     FP_RND_RM

#elif defined ( _MSC_VER )
#if ( _MSC_VER < 1400)
#define CGAL_IA_SETFPCW(CW) _controlfp (CW, _MCW_RC )
#define CGAL_IA_GETFPCW(CW) CW = _controlfp (0, 0 ) &  _MCW_RC
typedef unsigned short FPU_CW_t;
#else
#define CGAL_IA_SETFPCW(CW) unsigned int dummy; _controlfp_s (&dummy, CW, _MCW_RC )
#define CGAL_IA_GETFPCW(CW)_controlfp_s (&CW, 0, 0 ); CW  &=  _MCW_RC
typedef unsigned int FPU_CW_t;
#endif

#define CGAL_FE_TONEAREST    _RC_NEAR
#define CGAL_FE_TOWARDZERO   _RC_CHOP
#define CGAL_FE_UPWARD       _RC_UP
#define CGAL_FE_DOWNWARD     _RC_DOWN
# elif defined __VFP_FP__ && !defined __SOFTFP__
#define CGAL_IA_SETFPCW(CW) asm volatile ("VMSR FPSCR, %0" : :"r" (CW))
#define CGAL_IA_GETFPCW(CW) asm volatile ("VMRS %0, FPSCR" : "=r" (CW)); CW &= CGAL_FE_TOWARDZERO
typedef unsigned int FPU_CW_t;
#define CGAL_FE_TONEAREST    (0x0)
#define CGAL_FE_TOWARDZERO   (0xC00000)
#define CGAL_FE_UPWARD       (0x400000)
#define CGAL_FE_DOWNWARD     (0x800000)
# elif defined  __aarch64__
#define CGAL_IA_SETFPCW(CW) asm volatile ("MSR FPCR, %0" : :"r" (CW))
#define CGAL_IA_GETFPCW(CW) asm volatile ("MRS %0, FPCR" : "=r" (CW))
typedef unsigned int FPU_CW_t;
#define CGAL_FE_TONEAREST    (0x0)
#define CGAL_FE_TOWARDZERO   (0xC00000)
#define CGAL_FE_UPWARD       (0x400000)
#define CGAL_FE_DOWNWARD     (0x800000)

#else
// This is a version following the ISO C99 standard, which aims at portability.
// The drawbacks are speed on one hand, and also, on x86, it doesn't fix the
// extended mantissa issue (this is not a problem for IA, but it is one for
// some future modular computations).
#define CGAL_IA_SETFPCW(CW)  fesetround(CW)
#define CGAL_IA_GETFPCW(CW)  CW = fegetround()
typedef int FPU_CW_t;
#define CGAL_FE_TONEAREST    FE_TONEAREST
#define CGAL_FE_TOWARDZERO   FE_TOWARDZERO
#define CGAL_FE_UPWARD       FE_UPWARD
#define CGAL_FE_DOWNWARD     FE_DOWNWARD
#endif

// User interface:

inline
FPU_CW_t
FPU_get_cw (void)
{
    FPU_CW_t cw;
    CGAL_IA_GETFPCW(cw);
    return cw;
}

// User interface (cont):

inline
void
FPU_set_cw (FPU_CW_t cw)
{
  CGAL_IA_SETFPCW(cw);
}

inline
FPU_CW_t
FPU_get_and_set_cw (FPU_CW_t cw)
{
    FPU_CW_t old = FPU_get_cw();
    FPU_set_cw(cw);
    return old;
}


// A class whose constructor sets the FPU mode to +inf, saves a backup of it,
// and whose destructor resets it back to the saved state.

template <bool Protected = true> struct Protect_FPU_rounding;

template <>
struct Protect_FPU_rounding<true>
{
  Protect_FPU_rounding(FPU_CW_t r = CGAL_FE_UPWARD)
    : backup( FPU_get_and_set_cw(r) ) {}

  ~Protect_FPU_rounding()
  {
     FPU_set_cw(backup);
  }

private:
  FPU_CW_t backup;
};

template <>
struct Protect_FPU_rounding<false>
{
  Protect_FPU_rounding() {}
  Protect_FPU_rounding(FPU_CW_t /*= CGAL_FE_UPWARD*/) {}
};


// A wrapper on top of the Protect_FPU_rounding to add "expensive" checks
// of the rounding mode.  It is used internally, to benefit from the
// protector declarations to add checks in non-protected mode.

template <bool Protected = true>
struct Checked_protect_FPU_rounding
  : Protect_FPU_rounding<Protected>
{
  Checked_protect_FPU_rounding()
  {
    CGAL_expensive_assertion(FPU_get_cw() == CGAL_FE_UPWARD);
  }

  Checked_protect_FPU_rounding(FPU_CW_t r)
    : Protect_FPU_rounding<Protected>(r)
  {
    CGAL_expensive_assertion(FPU_get_cw() == CGAL_FE_UPWARD);
  }
};


// The class Set_ieee_double_precision forces the double precision (53bit mantissa),
// to protect from double rounding effects on x86 FPU.
// ( Note that it also sets the rounding mode to nearest. )
// Its destructor restores the FPU state as it was previously.
// Note that this affects "long double" as well, and other potential side effects.
// And note that it does not (cannot) "fix" the same problem for the exponent.

struct Set_ieee_double_precision
#ifdef CGAL_FPU_HAS_EXCESS_PRECISION
  : public Protect_FPU_rounding<>
{
  Set_ieee_double_precision()
    : Protect_FPU_rounding<>(CGAL_FE_TONEAREST) {}
};
#else
{
  Set_ieee_double_precision() {} // only to kill compiler warnings.
};
#endif


// The following function serves the same goal as Set_ieee_double_precision but
// does the change globally (no destructor resets the previous behavior).
inline void force_ieee_double_precision()
{
#ifdef CGAL_FPU_HAS_EXCESS_PRECISION
    FPU_set_cw(CGAL_FE_TONEAREST);
#endif
}

} //namespace CGAL

#ifdef CGAL_HEADER_ONLY
#include <CGAL/test_FPU_rounding_mode_impl.h>
#endif // CGAL_HEADER_ONLY

#endif // CGAL_FPU_H
