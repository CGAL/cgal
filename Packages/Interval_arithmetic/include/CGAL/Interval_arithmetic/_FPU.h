// ============================================================================
//
// Copyright (c) 1998 The CGAL Consortium
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

/*
 * Specify some *very* platform dependant functions, regarding the FPU
 * directed rounding modes. There is only support for double precision.
 *
 * I try to make it work with gcc/g++/cc/CC.
 * And I also try to provide the equivalent assembly code for GNU C.
 */

#if defined(__osf__)
#undef __osf
#define __osf
#endif

#if defined(__i386__) || defined(i386)	/* Needed by some old egcs versions */
#undef __i386
#define __i386
#endif

#if (	!defined(__i386)  && \
	!defined(__sparc) && \
	!defined(__alpha) && \
	!defined(__mips) \
	)
#error "Architecture not recognized."
#endif


// The test-suite fails for Mips.  Maybe an egcs-1.0.3 bug (1.1 is ok !):
// /usr/tmp/cca004Mh.s: Assembler messages:
// /usr/tmp/cca004Mh.s:6396: Error: ERROR: Illegal operands `ctc1'
// /usr/tmp/cca004Mh.s:6439: Error: ERROR: Illegal operands `ctc1'
// /usr/tmp/cca004Mh.s:6482: Error: ERROR: Illegal operands `ctc1'
// *** Error code 1 (bu21)

#if ( defined(__GNUC__) && \
    ( defined(__i386) || \
      defined(__mips) || \
      defined(__sparc) || \
      defined(__alpha) ) )
#define CGAL_IA_USE_ASSEMBLY
#endif

#ifndef CGAL_IA_USE_ASSEMBLY
#ifdef __linux
#include <fpu_control.h>
#endif
#ifdef __osf
#include <float.h>
#endif
#ifdef __sgi
  /* This #define is forced for backward compatibility with Irix 5.3. */
  /* I think it slows down Irix 6.2... */
#define __SGI_ieeefph__
#ifdef __SGI_ieeefph__
#include <ieeefp.h>
#else
#include <sys/fpu.h>
#endif
#endif
#ifdef __sun
#include <ieeefp.h>
#endif
#else	// CGAL_IA_USE_ASSEMBLY
#ifdef __i386
#define CGAL_IA_SETFPCW(CW) asm volatile ("fldcw %0" : : "g" (CW))
#define CGAL_IA_GETFPCW(CW) asm volatile ("fstcw %0" : "=m" (CW) : )
// x86:                                       rounding | precision | def. mask
static const unsigned short CGAL_FPU_cw_zero = 0xC00   |  0x200    | 0x107f; 
static const unsigned short CGAL_FPU_cw_near = 0x0     |  0x200    | 0x107f;
static const unsigned short CGAL_FPU_cw_up   = 0x800   |  0x200    | 0x107f; 
static const unsigned short CGAL_FPU_cw_down = 0x400   |  0x200    | 0x107f;
#endif
#ifdef __sparc
#define CGAL_IA_SETFPCW(CW) asm volatile ("ld %0,%%fsr" : : "g" (CW))
// Sparc:                                    rounding   | precision  |def. mask
static const unsigned int CGAL_FPU_cw_zero = 0x40000000 | 0x20000000 | 0x1f;
static const unsigned int CGAL_FPU_cw_near = 0x0        | 0x20000000 | 0x1f;
static const unsigned int CGAL_FPU_cw_up   = 0x80000000 | 0x20000000 | 0x1f;
static const unsigned int CGAL_FPU_cw_down = 0xc0000000 | 0x20000000 | 0x1f;
#endif
#ifdef __mips
#define CGAL_IA_SETFPCW(CW) asm volatile ("lw $8,%0; ctc1 $8 $31": :"m" (CW))
#define CGAL_IA_GETFPCW(CW) asm volatile ("cfc1 $8 $31; sw $8,%0": "=m" (CW) :)
static const unsigned int CGAL_FPU_cw_zero = 0x1;
static const unsigned int CGAL_FPU_cw_near = 0x0;
static const unsigned int CGAL_FPU_cw_up   = 0x2;
static const unsigned int CGAL_FPU_cw_down = 0x3;
#endif
#ifdef __alpha
#define CGAL_IA_SETFPCW(CW) asm volatile ("mt_fpcr %0; excb" : : "f" (CW))
// Alpha:                      rounding
static const unsigned long CGAL_FPU_cw_zero = 0x0000000000000000UL;
static const unsigned long CGAL_FPU_cw_near = 0x0800000000000000UL;
static const unsigned long CGAL_FPU_cw_up   = 0x0c00000000000000UL;
static const unsigned long CGAL_FPU_cw_down = 0x0400000000000000UL;
#endif
#endif

static inline void CGAL_FPU_set_rounding_to_zero (void);
static inline void CGAL_FPU_set_rounding_to_nearest (void);
static inline void CGAL_FPU_set_rounding_to_infinity (void);
static inline void CGAL_FPU_set_rounding_to_minus_infinity (void);


/* Code of those inline functions: */

static inline void CGAL_FPU_set_rounding_to_zero (void)
{
#ifdef CGAL_IA_USE_ASSEMBLY
	CGAL_IA_SETFPCW(CGAL_FPU_cw_zero);
#else
#ifdef __linux
#ifdef __i386 
        __setfpucw(0x1e72);
#else
        __setfpucw(0x1f72);	/* FIX ME */
#endif
#endif

#ifdef __osf
	write_rnd(FP_RND_RZ);
#endif

#ifdef __sun
        fpsetround(FP_RZ);
#endif

#ifdef __sgi
#ifdef __SGI_ieeefph__
	fp_rnd mode = FP_RZ;
	fpsetround(mode);
#else
	union fpc_csr fpu_ctl;
	fpu_ctl.fc_word = get_fpc_csr();
	fpu_ctl.fc_struct.rounding_mode = ROUND_TO_ZERO;
	set_fpc_csr (fpu_ctl.fc_word);
#endif
#endif
#endif
}

static inline void CGAL_FPU_set_rounding_to_nearest (void)
{
#ifdef CGAL_IA_USE_ASSEMBLY
	CGAL_IA_SETFPCW(CGAL_FPU_cw_near);
#else
#ifdef __osf
	write_rnd(FP_RND_RN);
#endif

#ifdef __sun
	fpsetround(FP_RN);
#endif

#ifdef __linux
#ifdef __i386
	__setfpucw(0x1272);
#else
	__setfpucw(0x1372);
#endif
#endif

#ifdef __sgi
#ifdef __SGI_ieeefph__
	fp_rnd mode = FP_RN;
	fpsetround(mode);
#else
	union fpc_csr fpu_ctl;
	fpu_ctl.fc_word = get_fpc_csr();
	fpu_ctl.fc_struct.rounding_mode = ROUND_TO_NEAREST;
	set_fpc_csr (fpu_ctl.fc_word);
#endif
#endif
#endif
}

static inline void CGAL_FPU_set_rounding_to_infinity (void)
{
#ifdef CGAL_IA_USE_ASSEMBLY
	CGAL_IA_SETFPCW(CGAL_FPU_cw_up);
#else
#ifdef __osf
	write_rnd(FP_RND_RP);
#endif

#ifdef __sun
        fpsetround(FP_RP);
#endif

#ifdef __linux
#ifdef __i386
	__setfpucw(0x1a72);
#else
	__setfpucw(0x1b72);
#endif
#endif

#ifdef __sgi
#ifdef __SGI_ieeefph__
	fp_rnd mode = FP_RP;
	fpsetround(mode);
#else
	union fpc_csr fpu_ctl;
	fpu_ctl.fc_word = get_fpc_csr();
	fpu_ctl.fc_struct.rounding_mode = ROUND_TO_PLUS_INFINITY;
	set_fpc_csr (fpu_ctl.fc_word);
#endif
#endif
#endif
}

static inline void CGAL_FPU_set_rounding_to_minus_infinity (void)
{
#ifdef CGAL_IA_USE_ASSEMBLY
	CGAL_IA_SETFPCW(CGAL_FPU_cw_down);
#else
#ifdef __osf
	write_rnd(FP_RND_RM);
#endif

#ifdef __sun
        fpsetround(FP_RM);
#endif

#ifdef __linux
#ifdef __i386
        __setfpucw(0x1672);
#else
        __setfpucw(0x1772);
#endif
#endif

#ifdef __sgi
#ifdef __SGI_ieeefph__
	fp_rnd mode = FP_RM;
	fpsetround(mode);
#else
	union fpc_csr fpu_ctl;
	fpu_ctl.fc_word = get_fpc_csr();
	fpu_ctl.fc_struct.rounding_mode = ROUND_TO_MINUS_INFINITY;
	set_fpc_csr (fpu_ctl.fc_word);
#endif
#endif
#endif
}

// Rouding mode empiric testing.

enum CGAL_FPU_rounding_mode
{
    CGAL_FPU_MINUS_INFINITY,
    CGAL_FPU_ZERO,
    CGAL_FPU_NEAREST,
    CGAL_FPU_PLUS_INFINITY
};


// The result of 1+epsilon, 1-epsilon, -1+epsilon is enough
// to detect exactly the rounding mode.
// (epsilon = MIN_DOUBLE, ulp = 2^-52 or 2^-53).
// ----------------------------------------------------
// rounding mode:	 +inf	 -inf	 0	 nearest
// ----------------------------------------------------
//  1+epsilon		 1+ulp	 1	 1	 1
//  1-epsilon		 1	 1-ulp	 1-ulp	 1
// -1+epsilon		-1+ulp	-1	-1+ulp	-1
// ----------------------------------------------------

static inline CGAL_FPU_rounding_mode CGAL_FPU_get_rounding_mode ()
{
    const double m = 5e-324; // CGAL_Interval_nt_advanced::min_double;
    // If not marked "volatile", the result is false when optimizing
    // because the constants are pre-computed at compile time !!!
    volatile double x, y, z, xe, ye, ze;
    x =  1.0; xe = x + m;
    y =  1.0; ye = y - m;
    z = -1.0; ze = z + m;
    if ((x == xe) && (y == ye) && (z == ze)) return CGAL_FPU_NEAREST;
    if (y == ye) return CGAL_FPU_PLUS_INFINITY;
    if (z == ze) return CGAL_FPU_MINUS_INFINITY;
    return CGAL_FPU_ZERO;
}

#endif // CGAL_FPU_H
