// ============================================================================
//
// Copyright (c) 1998,1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Interval_arithmetic/_FPU.h
// revision      : $Revision$
// revision_date : $Date$
// package       : Interval Arithmetic
// author(s)     : Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_FPU_H
#define CGAL_FPU_H

// This file specifies some platform dependant functions, regarding the FPU
// directed rounding modes.  There is only support for double precision.


// Some useful constants

#define CGAL_IA_MIN_DOUBLE (5e-324) // subnormal
#define CGAL_IA_MAX_DOUBLE (1.7976931348623157081e+308)


#if defined(__osf__)
#undef __osf
#define __osf
#endif

#if (	!defined(__i386)  && \
	!defined(__sparc) && \
	!defined(__alpha) && \
	!defined(__mips) \
	)
#error "Architecture not supported."
#endif


// The test-suite fails for Irix 5.3, because its assembler is BAD:
// /usr/tmp/cca004Mh.s: Assembler messages:
// /usr/tmp/cca004Mh.s:6396: Error: ERROR: Illegal operands `ctc1'
// /usr/tmp/cca004Mh.s:6439: Error: ERROR: Illegal operands `ctc1'
// /usr/tmp/cca004Mh.s:6482: Error: ERROR: Illegal operands `ctc1'
// *** Error code 1 (bu21)
//
// Does it need passing some options ? (like -Wa,... on Sun ?)

#ifndef CGAL_IA_DONT_USE_ASSEMBLY
#if ( defined(__GNUC__) && \
    ( defined(__i386)   || \
      defined(__sparc)  || \
      defined(__mips)   || \
      defined(__alpha)  ) )
#define CGAL_IA_USE_ASSEMBLY
#endif
#endif // CGAL_IA_DONT_USE_ASSEMBLY

#ifdef __linux
#include <fpu_control.h>
#endif

#ifndef CGAL_IA_USE_ASSEMBLY
#ifdef __sun
#include <ieeefp.h>
#endif
#ifdef __osf
#include <float.h>
#endif
#ifdef __sgi
  // This #define is forced for backward compatibility with Irix 5.3.
  // I think it slows down Irix 6.2...
#define __SGI_ieeefph__
#ifdef __SGI_ieeefph__
#include <ieeefp.h>
#else
#include <sys/fpu.h>
#endif
#endif
#endif // CGAL_IA_USE_ASSEMBLY

CGAL_BEGIN_NAMESPACE

#ifdef __i386
#ifdef CGAL_IA_USE_ASSEMBLY
#define CGAL_IA_SETFPCW(CW) asm volatile ("fldcw %0" : :"m" (CW))
#define CGAL_IA_GETFPCW(CW) asm volatile ("fstcw %0" : "=m" (CW))
#endif
typedef unsigned short FPU_CW_t;
enum {  //               rounding | def. mask
    FPU_cw_near = _FPU_RC_NEAREST | 0x127f,
    FPU_cw_zero = _FPU_RC_ZERO    | 0x127f,
    FPU_cw_up   = _FPU_RC_UP      | 0x127f,
    FPU_cw_down = _FPU_RC_DOWN    | 0x127f
};
#endif // __i386

#ifdef __sparc
#ifdef CGAL_IA_USE_ASSEMBLY
#define CGAL_IA_SETFPCW(CW) asm volatile ("ld %0,%%fsr" : :"m" (CW))
#define CGAL_IA_GETFPCW(CW) asm volatile ("st %%fsr,%0" : "=m" (CW))
typedef unsigned int FPU_CW_t;
enum {  //        rounding   | precision  | def.mask
    FPU_cw_near = 0x0        | 0x20000000 | 0x1f,
    FPU_cw_zero = 0x40000000 | 0x20000000 | 0x1f,
    FPU_cw_up   = 0x80000000 | 0x20000000 | 0x1f,
    FPU_cw_down = 0xc0000000 | 0x20000000 | 0x1f
};
#else
typedef fp_rnd FPU_CW_t;
enum {
    FPU_cw_near = FP_RN,
    FPU_cw_zero = FP_RZ,
    FPU_cw_up   = FP_RP,
    FPU_cw_down = FP_RM
};
#endif
#endif // __sparc

#ifdef __mips
#ifdef CGAL_IA_USE_ASSEMBLY
#define CGAL_IA_SETFPCW(CW) asm volatile ("ctc1 %0,$31" : :"r" (CW))
#define CGAL_IA_GETFPCW(CW) asm volatile ("cfc1 %0,$31" : "=r" (CW))
typedef unsigned int FPU_CW_t;
enum {
    FPU_cw_near = 0x0,
    FPU_cw_zero = 0x1,
    FPU_cw_up   = 0x2,
    FPU_cw_down = 0x3
};
#else
#ifdef __SGI_ieeefph__ // 2 cases for C... I love IRIX !!!
typedef fp_rnd FPU_CW_t;
enum {
    FPU_cw_near = FP_RN,
    FPU_cw_zero = FP_RZ,
    FPU_cw_up   = FP_RP,
    FPU_cw_down = FP_RM
};
#else
typedef int FPU_CW_t;
enum {
    FPU_cw_near = ROUND_TO_NEAREST,
    FPU_cw_zero = ROUND_TO_ZERO,
    FPU_cw_up   = ROUND_TO_PLUS_INFINITY,
    FPU_cw_down = ROUND_TO_MINUS_INFINITY
};
#endif
#endif
#endif // __mips

#ifdef __alpha // This one is not really supported.
#ifdef CGAL_IA_USE_ASSEMBLY
#define CGAL_IA_SETFPCW(CW) asm volatile ("mt_fpcr %0; excb" : :"f" (CW))
#define CGAL_IA_GETFPCW(CW) asm volatile ("excb; mf_fpcr %0" : "=f" (CW))
typedef unsigned long FPU_CW_t;
enum { //         rounding
    // I guess it won't work, because enum == int.
    FPU_cw_zero = 0x0000000000000000UL,
    FPU_cw_near = 0x0800000000000000UL,
    FPU_cw_up   = 0x0c00000000000000UL,
    FPU_cw_down = 0x0400000000000000UL
};
#else
typedef unsigned int FPU_CW_t;
enum {
    FPU_cw_zero = FP_RND_RZ,
    FPU_cw_near = FP_RND_RN,
    FPU_cw_up   = FP_RND_RP,
    FPU_cw_down = FP_RND_RM
};
#endif
#endif // __alpha


// Main functions:
static inline FPU_CW_t FPU_get_cw (void);
static inline void FPU_set_cw (FPU_CW_t);

static inline FPU_CW_t FPU_get_cw (void)
{
    FPU_CW_t cw;
#ifdef CGAL_IA_USE_ASSEMBLY
    CGAL_IA_GETFPCW(cw);
#else

#ifdef __linux
#error  // It seems there's no C function in libc5 !!!
#endif

#ifdef __osf
    cw = read_rnd();
#endif

#ifdef __sun
    cw = fpgetround();
#endif

#ifdef __sgi
#ifdef __SGI_ieeefph__
    cw = fpgetround();
#else
    union fpc_csr fpu_ctl;
    fpu_ctl.fc_word = get_fpc_csr();
    cw = fpu_ctl.fc_struct.rounding_mode;
#endif
#endif

#endif
    return cw;
}

static inline void FPU_set_cw (FPU_CW_t cw)
{
#ifdef CGAL_IA_USE_ASSEMBLY
    CGAL_IA_SETFPCW(cw);
#else

#ifdef __linux
    __setfpucw(cw);
#endif

#ifdef __osf
    write_rnd(cw);
#endif

#ifdef __sun
    fpsetround(cw);
#endif

#ifdef __sgi
#ifdef __SGI_ieeefph__
    fpsetround(cw);
#else
    union fpc_csr fpu_ctl;
    fpu_ctl.fc_word = get_fpc_csr();
    fpu_ctl.fc_struct.rounding_mode = cw;
    set_fpc_csr (fpu_ctl.fc_word);
#endif
#endif

#endif
}

// Obscolete: wrappers for the old interface.

#if 1
static inline void FPU_set_rounding_to_zero (void)
{ FPU_set_cw(FPU_cw_zero); }

static inline void FPU_set_rounding_to_nearest (void)
{ FPU_set_cw(FPU_cw_near); }

static inline void FPU_set_rounding_to_infinity (void)
{ FPU_set_cw(FPU_cw_up); }

static inline void FPU_set_rounding_to_minus_infinity (void)
{ FPU_set_cw(FPU_cw_down); }
#endif

CGAL_END_NAMESPACE

#endif // CGAL_FPU_H
