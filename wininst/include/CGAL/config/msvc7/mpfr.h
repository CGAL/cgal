#ifndef __MPFR_LIBGMP_DLL
#define __MPFR_LIBGMP_DLL 1
#endif
/* mpfr.h -- Include file for mpfr.

Copyright 1999, 2000, 2001, 2002, 2003 Free Software Foundation, Inc.

This file is part of the MPFR Library.

The MPFR Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 2.1 of the License, or (at your
option) any later version.

The MPFR Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the MPFR Library; see the file COPYING.LIB.  If not, write to
the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
MA 02111-1307, USA. */

#ifndef __MPFR_H
#define __MPFR_H

/* check if stdio.h is included */
#if defined (FILE)                                              \
  || defined (H_STDIO)                                          \
  || defined (_H_STDIO)               /* AIX */                 \
  || defined (_STDIO_H)               /* glibc, Sun, SCO */     \
  || defined (_STDIO_H_)              /* BSD, OSF */            \
  || defined (__STDIO_H)              /* Borland */             \
  || defined (__STDIO_H__)            /* IRIX */                \
  || defined (_STDIO_INCLUDED)        /* HPUX */                \
  || defined (__dj_include_stdio_h_)  /* DJGPP */               \
  || defined (_FILE_DEFINED)          /* Microsoft */           \
  || defined (__STDIO__)              /* Apple MPW MrC */
#define _MPFR_H_HAVE_FILE 1
#endif

/* Definition of rounding modes */

#define GMP_RNDN 0
#define GMP_RNDZ 1
#define GMP_RNDU 2
#define GMP_RNDD 3

/* Definition of exponent limits */

#define MPFR_EMAX_DEFAULT ((mp_exp_t) (((unsigned long) 1 << 30) - 1))
#define MPFR_EMIN_DEFAULT (-(MPFR_EMAX_DEFAULT))

#define MPFR_EMIN_MIN MPFR_EMIN_DEFAULT
#define MPFR_EMIN_MAX MPFR_EMAX_DEFAULT
#define MPFR_EMAX_MIN MPFR_EMIN_DEFAULT
#define MPFR_EMAX_MAX MPFR_EMAX_DEFAULT

/* Flags */

#define MPFR_FLAGS_UNDERFLOW 1
#define MPFR_FLAGS_OVERFLOW 2
#define MPFR_FLAGS_NAN 4
#define MPFR_FLAGS_INEXACT 8
#define MPFR_FLAGS_ALL 15

/* Definitions of types and their semantics */

typedef unsigned long int mp_prec_t; /* easy to change if necessary */
#define MPFR_PREC_MIN 2
#define MPFR_PREC_MAX (ULONG_MAX >> 1)
/* Limit mainly due to the multiplication code. */

typedef int mp_rnd_t;

typedef struct {
  mp_prec_t _mpfr_prec; /* WARNING : for the mpfr type, the precision */
                              /* should be understood as the number of BITS,*/
			      /* not the number of mp_limb_t's. This means  */
			      /* that the corresponding number of allocated
				 limbs is >= ceil(_mp_prec/BITS_PER_MP_LIMB) */
  mp_size_t _mpfr_size;         /* MPFR_ABSSIZE(.) is the number of allocated
				 limbs the field _mp_d points to.
				 The sign is that of _mpfr_size.
				 The number 0 is such that _mp_d[k-1]=0
				 where k = ceil(_mp_prec/BITS_PER_MP_LIMB) */
  mp_exp_t _mpfr_exp;
  mp_limb_t *_mpfr_d;
}
__mpfr_struct;

/*
   The number represented is

    sign(_mpfr_size)*(_mpfr_d[k-1]/B+_mpfr_d[k-2]/B^2+...+_mpfr_d[0]/B^k)*2^_mpfr_exp

   where k=ceil(_mp_prec/BITS_PER_MP_LIMB) and B=2^BITS_PER_MP_LIMB.

   For the msb (most significant bit) normalized representation, we must have
   _mpfr_d[k-1]>=B/2, unless the number is zero (in that case its sign is still
   given by sign(_mpfr_size)).

   We must also have the last k*BITS_PER_MP_LIMB-_mp_prec bits set to zero.
*/

typedef __mpfr_struct mpfr_t[1];
typedef __mpfr_struct *mpfr_ptr;
typedef __gmp_const __mpfr_struct *mpfr_srcptr;

#define MPFR_SIGN(x) (((x)->_mpfr_size >> 31) ? -1 : 1)

/* Prototypes */

#ifndef _PROTO
#if defined (__STDC__) || defined (__cplusplus)
#define _PROTO(x) x
#else
#define _PROTO(x) ()
#endif
#endif

/* _PROTO will be renamed __GMP_PROTO in gmp 4.1 */
#ifndef __GMP_PROTO
#define __GMP_PROTO(x) _PROTO(x)
#endif

#if defined (__GNUC__)
#define __MPFR_DECLSPEC_EXPORT  __declspec(__dllexport__)
#define __MPFR_DECLSPEC_IMPORT  __declspec(__dllimport__)
#endif
#if defined (_MSC_VER) || defined (__BORLANDC__)
#define __MPFR_DECLSPEC_EXPORT  __declspec(dllexport)
#define __MPFR_DECLSPEC_IMPORT  __declspec(dllimport)
#endif
#ifdef __WATCOMC__
#define __MPFR_DECLSPEC_EXPORT  __export
#define __MPFR_DECLSPEC_IMPORT  __import
#endif
#ifdef __IBMC__
#define __MPFR_DECLSPEC_EXPORT  _Export
#define __MPFR_DECLSPEC_IMPORT  _Import
#endif

#if __MPFR_LIBGMP_DLL
#if __MPFR_WITHIN_MPFR
/* compiling to go into a DLL libgmp */
#define __MPFR_DECLSPEC  __MPFR_DECLSPEC_EXPORT
#else
/* compiling to go into an application which will link to a DLL libgmp */
#define __MPFR_DECLSPEC  __MPFR_DECLSPEC_IMPORT
#endif
#else
/* all other cases */
#define __MPFR_DECLSPEC
#endif

#if __MPFR_LIBGMP_DLL
#if __GMP_WITHIN_MPFRXX
/* compiling to go into a DLL libgmpxx */
#define __MPFR_DECLSPEC_XX  __MPFR_DECLSPEC_EXPORT
#else
/* compiling to go into a application which will link to a DLL libgmpxx */
#define __MPFR_DECLSPEC_XX  __MPFR_DECLSPEC_IMPORT
#endif
#else
/* all other cases */
#define __MPFR_DECLSPEC_XX
#endif

#if defined (__cplusplus)
extern "C" {
#endif

__MPFR_DECLSPEC extern unsigned int __gmpfr_flags;
__MPFR_DECLSPEC extern mp_exp_t __gmpfr_emin;
__MPFR_DECLSPEC extern mp_exp_t __gmpfr_emax;
__MPFR_DECLSPEC mp_exp_t mpfr_get_emin _PROTO ((void));
__MPFR_DECLSPEC int mpfr_set_emin _PROTO ((mp_exp_t));
__MPFR_DECLSPEC mp_exp_t mpfr_get_emax _PROTO ((void));
__MPFR_DECLSPEC int mpfr_set_emax _PROTO ((mp_exp_t));
__MPFR_DECLSPEC void mpfr_clear_flags _PROTO ((void));
__MPFR_DECLSPEC void mpfr_clear_underflow _PROTO ((void));
__MPFR_DECLSPEC void mpfr_clear_overflow _PROTO ((void));
__MPFR_DECLSPEC void mpfr_clear_nanflag _PROTO ((void));
__MPFR_DECLSPEC void mpfr_clear_inexflag _PROTO ((void));
__MPFR_DECLSPEC int mpfr_check_range _PROTO ((mpfr_ptr, int, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_underflow_p _PROTO ((void));
__MPFR_DECLSPEC int mpfr_overflow_p _PROTO ((void));
__MPFR_DECLSPEC int mpfr_nanflag_p _PROTO ((void));
__MPFR_DECLSPEC int mpfr_inexflag_p _PROTO ((void));
__MPFR_DECLSPEC void mpfr_init2 _PROTO ((mpfr_ptr, mp_prec_t));

__MPFR_DECLSPEC void mpfr_init _PROTO ((mpfr_ptr));
__MPFR_DECLSPEC int mpfr_prec_round _PROTO ((mpfr_ptr, mp_prec_t, mp_rnd_t));
#define mpfr_round_prec(x,r,p) mpfr_prec_round(x,p,r) /* compatibility 2.0.1 */
__MPFR_DECLSPEC int mpfr_can_round _PROTO ((mpfr_ptr, mp_exp_t, mp_rnd_t, mp_rnd_t, mp_prec_t));
__MPFR_DECLSPEC mp_exp_t mpfr_get_exp _PROTO ((mpfr_srcptr));
__MPFR_DECLSPEC int mpfr_set_exp _PROTO ((mpfr_ptr, mp_exp_t));
__MPFR_DECLSPEC int mpfr_set_d _PROTO ((mpfr_ptr, double, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_set_ld _PROTO ((mpfr_ptr, long double, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_set_z _PROTO ((mpfr_ptr, mpz_srcptr, mp_rnd_t));
__MPFR_DECLSPEC void mpfr_set_nan _PROTO ((mpfr_ptr));
__MPFR_DECLSPEC void mpfr_set_inf _PROTO ((mpfr_ptr, int));
__MPFR_DECLSPEC mp_exp_t mpfr_get_z_exp _PROTO ((mpz_ptr, mpfr_srcptr));
__MPFR_DECLSPEC int mpfr_set_q _PROTO ((mpfr_ptr, mpq_srcptr, mp_rnd_t));
__MPFR_DECLSPEC double mpfr_get_d _PROTO ((mpfr_srcptr, mp_rnd_t));
__MPFR_DECLSPEC long double mpfr_get_ld _PROTO ((mpfr_srcptr, mp_rnd_t));
__MPFR_DECLSPEC double mpfr_get_d1 _PROTO ((mpfr_srcptr));
__MPFR_DECLSPEC double mpfr_get_d_2exp _PROTO ((long *, mpfr_srcptr, mp_rnd_t));
__MPFR_DECLSPEC long mpfr_get_si _PROTO ((mpfr_srcptr, mp_rnd_t));
__MPFR_DECLSPEC unsigned long mpfr_get_ui _PROTO ((mpfr_srcptr, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_set_f _PROTO ((mpfr_ptr, mpf_srcptr, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_set_si _PROTO ((mpfr_ptr, long, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_set_ui _PROTO ((mpfr_ptr, unsigned long, mp_rnd_t));
__MPFR_DECLSPEC void mpfr_random _PROTO ((mpfr_ptr));
__MPFR_DECLSPEC void mpfr_random2 _PROTO ((mpfr_ptr, mp_size_t, mp_exp_t));
__MPFR_DECLSPEC int mpfr_urandomb _PROTO ((mpfr_ptr, gmp_randstate_t));
__MPFR_DECLSPEC void mpfr_clear _PROTO ((mpfr_ptr));
__MPFR_DECLSPEC void mpfr_nextabove _PROTO ((mpfr_ptr));
__MPFR_DECLSPEC void mpfr_nextbelow _PROTO ((mpfr_ptr));
__MPFR_DECLSPEC void mpfr_nexttoward _PROTO ((mpfr_ptr, mpfr_srcptr));
__MPFR_DECLSPEC int mpfr_set_str _PROTO ((mpfr_ptr, __gmp_const char *, int, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_init_set_str _PROTO ((mpfr_ptr, __gmp_const char *, int, mp_rnd_t));
__MPFR_DECLSPEC char* mpfr_get_str _PROTO ((char *, mp_exp_t *, int, size_t, mpfr_srcptr, mp_rnd_t));

#ifdef _MPFR_H_HAVE_FILE
__MPFR_DECLSPEC size_t mpfr_inp_str _PROTO ((mpfr_ptr, FILE *, int, mp_rnd_t));
__MPFR_DECLSPEC size_t mpfr_out_str _PROTO ((FILE *, int, size_t, mpfr_srcptr, mp_rnd_t));
#endif

__MPFR_DECLSPEC int mpfr_mul _PROTO ((mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_pow_ui _PROTO ((mpfr_ptr, mpfr_srcptr, unsigned long int, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_ui_pow_ui _PROTO ((mpfr_ptr, unsigned long int, unsigned long int, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_div _PROTO ((mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_agm _PROTO ((mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_sqrt _PROTO ((mpfr_ptr, mpfr_srcptr, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_sqrt_ui _PROTO ((mpfr_ptr, unsigned long, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_add _PROTO ((mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_add_ui _PROTO ((mpfr_ptr, mpfr_srcptr, unsigned long, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_sub_ui _PROTO ((mpfr_ptr, mpfr_srcptr, unsigned long, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_add_one_ulp _PROTO ((mpfr_ptr, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_sub _PROTO ((mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_ui_sub _PROTO ((mpfr_ptr, unsigned long, mpfr_srcptr, mp_rnd_t));
__MPFR_DECLSPEC void mpfr_reldiff _PROTO ((mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_const_pi _PROTO ((mpfr_ptr, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_const_log2 _PROTO ((mpfr_ptr, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_const_euler _PROTO ((mpfr_ptr, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_log _PROTO ((mpfr_ptr, mpfr_srcptr, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_exp _PROTO ((mpfr_ptr, mpfr_srcptr, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_exp2 _PROTO ((mpfr_ptr, mpfr_srcptr, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_sin _PROTO ((mpfr_ptr, mpfr_srcptr, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_sin_cos _PROTO ((mpfr_ptr, mpfr_ptr, mpfr_srcptr, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_cos _PROTO ((mpfr_ptr, mpfr_srcptr, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_tan _PROTO ((mpfr_ptr, mpfr_srcptr, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_mul_ui _PROTO((mpfr_ptr, mpfr_srcptr, unsigned long int, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_cmp_ui_2exp _PROTO ((mpfr_srcptr, unsigned long int, mp_exp_t));
__MPFR_DECLSPEC int mpfr_cmp_si_2exp _PROTO ((mpfr_srcptr, long int, mp_exp_t));
__MPFR_DECLSPEC int mpfr_mul_2exp _PROTO((mpfr_ptr, mpfr_srcptr, unsigned long int, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_div_2exp _PROTO((mpfr_ptr, mpfr_srcptr, unsigned long int, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_mul_2ui _PROTO((mpfr_ptr, mpfr_srcptr, unsigned long int, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_div_2ui _PROTO((mpfr_ptr, mpfr_srcptr, unsigned long int, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_mul_2si _PROTO((mpfr_ptr, mpfr_srcptr, long int, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_div_2si _PROTO((mpfr_ptr, mpfr_srcptr, long int, mp_rnd_t));
__MPFR_DECLSPEC void mpfr_set_prec _PROTO((mpfr_ptr, mp_prec_t));
__MPFR_DECLSPEC void mpfr_set_prec_raw _PROTO((mpfr_ptr, mp_prec_t));
__MPFR_DECLSPEC void mpfr_set_default_prec _PROTO((mp_prec_t));
__MPFR_DECLSPEC mp_prec_t mpfr_get_default_prec _PROTO((void));
__MPFR_DECLSPEC extern mp_prec_t __gmpfr_default_fp_bit_precision;
__MPFR_DECLSPEC extern mp_rnd_t __gmpfr_default_rounding_mode;
__MPFR_DECLSPEC __gmp_const char * mpfr_print_rnd_mode _PROTO((mp_rnd_t));
__MPFR_DECLSPEC int mpfr_neg _PROTO((mpfr_ptr, mpfr_srcptr, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_sub_one_ulp _PROTO((mpfr_ptr, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_div_ui _PROTO((mpfr_ptr, mpfr_srcptr, unsigned long int, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_ui_div _PROTO((mpfr_ptr, unsigned long int, mpfr_srcptr, mp_rnd_t));
__MPFR_DECLSPEC mp_prec_t mpfr_get_prec _PROTO((mpfr_srcptr));
__MPFR_DECLSPEC void mpfr_set_default_rounding_mode _PROTO((mp_rnd_t));
__MPFR_DECLSPEC int mpfr_eq _PROTO((mpfr_srcptr, mpfr_srcptr, unsigned long));
__MPFR_DECLSPEC int mpfr_rint _PROTO((mpfr_ptr, mpfr_srcptr, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_round _PROTO((mpfr_ptr, mpfr_srcptr));
__MPFR_DECLSPEC int mpfr_trunc _PROTO((mpfr_ptr, mpfr_srcptr));
__MPFR_DECLSPEC int mpfr_ceil _PROTO((mpfr_ptr, mpfr_srcptr));
__MPFR_DECLSPEC int mpfr_floor _PROTO((mpfr_ptr, mpfr_srcptr));
__MPFR_DECLSPEC int mpfr_frac _PROTO((mpfr_ptr, mpfr_srcptr, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_fits_ulong_p _PROTO((mpfr_srcptr, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_fits_slong_p _PROTO((mpfr_srcptr, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_fits_uint_p _PROTO((mpfr_srcptr, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_fits_sint_p _PROTO((mpfr_srcptr, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_fits_ushort_p _PROTO((mpfr_srcptr, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_fits_sshort_p _PROTO((mpfr_srcptr, mp_rnd_t));
__MPFR_DECLSPEC void mpfr_extract _PROTO((mpz_ptr, mpfr_srcptr, unsigned int));
__MPFR_DECLSPEC void mpfr_swap _PROTO((mpfr_ptr, mpfr_ptr));
__MPFR_DECLSPEC void mpfr_dump _PROTO((mpfr_srcptr, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_set4 _PROTO ((mpfr_ptr, mpfr_srcptr, mp_rnd_t, int));
__MPFR_DECLSPEC int mpfr_cmp3 _PROTO ((mpfr_srcptr, mpfr_srcptr, int));
__MPFR_DECLSPEC int mpfr_cmp_d _PROTO ((mpfr_srcptr, double));
__MPFR_DECLSPEC int mpfr_cmpabs _PROTO ((mpfr_srcptr, mpfr_srcptr));

#define mpfr_cmp_abs mpfr_cmpabs /* keep for compatibility with mpfr-2.0.1 */

__MPFR_DECLSPEC int mpfr_nan_p _PROTO((mpfr_srcptr));
__MPFR_DECLSPEC int mpfr_inf_p _PROTO((mpfr_srcptr));
__MPFR_DECLSPEC int mpfr_number_p _PROTO((mpfr_srcptr));
__MPFR_DECLSPEC int mpfr_acos _PROTO ((mpfr_ptr, mpfr_srcptr, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_asin _PROTO ((mpfr_ptr, mpfr_srcptr, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_atan _PROTO ((mpfr_ptr, mpfr_srcptr, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_erf _PROTO ((mpfr_ptr, mpfr_srcptr, mp_rnd_t));

__MPFR_DECLSPEC int mpfr_sinh _PROTO((mpfr_ptr, mpfr_srcptr, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_tanh _PROTO((mpfr_ptr, mpfr_srcptr, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_factorial _PROTO ((mpfr_ptr, unsigned long int, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_ui_pow _PROTO ((mpfr_ptr, unsigned long int, mpfr_srcptr, mp_rnd_t));

__MPFR_DECLSPEC int mpfr_atanh _PROTO((mpfr_ptr, mpfr_srcptr, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_acosh _PROTO((mpfr_ptr, mpfr_srcptr, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_asinh _PROTO((mpfr_ptr, mpfr_srcptr, mp_rnd_t));

__MPFR_DECLSPEC int mpfr_cosh _PROTO((mpfr_ptr, mpfr_srcptr, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_sinh _PROTO((mpfr_ptr, mpfr_srcptr, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_tanh _PROTO((mpfr_ptr, mpfr_srcptr, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_asin _PROTO ((mpfr_ptr, mpfr_srcptr, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_atan _PROTO ((mpfr_ptr, mpfr_srcptr, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_fac_ui _PROTO ((mpfr_ptr, unsigned long int, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_fma _PROTO ((mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mpfr_srcptr, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_hypot _PROTO ((mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_pow _PROTO ((mpfr_ptr, mpfr_srcptr,mpfr_srcptr, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_pow_si _PROTO ((mpfr_ptr, mpfr_srcptr, long int, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_integer_p _PROTO ((mpfr_srcptr));
__MPFR_DECLSPEC int mpfr_log2 _PROTO ((mpfr_ptr, mpfr_srcptr, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_log10 _PROTO ((mpfr_ptr, mpfr_srcptr, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_log1p _PROTO ((mpfr_ptr, mpfr_srcptr, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_expm1 _PROTO ((mpfr_ptr, mpfr_srcptr, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_cbrt _PROTO ((mpfr_ptr, mpfr_srcptr, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_gamma _PROTO ((mpfr_ptr, mpfr_srcptr, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_zeta _PROTO ((mpfr_ptr, mpfr_srcptr, mp_rnd_t));

__MPFR_DECLSPEC int mpfr_min _PROTO ((mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_max _PROTO ((mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_dim _PROTO ((mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_copysign _PROTO ((mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mp_rnd_t));

__MPFR_DECLSPEC int mpfr_mul_z _PROTO ((mpfr_ptr, mpfr_srcptr, mpz_srcptr, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_div_z _PROTO ((mpfr_ptr, mpfr_srcptr, mpz_srcptr, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_add_z _PROTO ((mpfr_ptr, mpfr_srcptr, mpz_srcptr, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_sub_z _PROTO ((mpfr_ptr, mpfr_srcptr, mpz_srcptr, mp_rnd_t));

__MPFR_DECLSPEC int mpfr_mul_q _PROTO ((mpfr_ptr, mpfr_srcptr, mpq_srcptr, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_div_q _PROTO ((mpfr_ptr, mpfr_srcptr, mpq_srcptr, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_add_q _PROTO ((mpfr_ptr, mpfr_srcptr, mpq_srcptr, mp_rnd_t));
__MPFR_DECLSPEC int mpfr_sub_q _PROTO ((mpfr_ptr, mpfr_srcptr, mpq_srcptr, mp_rnd_t));

__MPFR_DECLSPEC int mpfr_greater_p _PROTO ((mpfr_srcptr, mpfr_srcptr));
__MPFR_DECLSPEC int mpfr_greaterequal_p _PROTO ((mpfr_srcptr, mpfr_srcptr));
__MPFR_DECLSPEC int mpfr_less_p _PROTO ((mpfr_srcptr, mpfr_srcptr));
__MPFR_DECLSPEC int mpfr_lessequal_p _PROTO ((mpfr_srcptr, mpfr_srcptr));
__MPFR_DECLSPEC int mpfr_lessgreater_p _PROTO ((mpfr_srcptr, mpfr_srcptr));
__MPFR_DECLSPEC int mpfr_equal_p _PROTO ((mpfr_srcptr, mpfr_srcptr));
__MPFR_DECLSPEC int mpfr_unordered_p _PROTO ((mpfr_srcptr, mpfr_srcptr));

#if defined (__cplusplus)
}
#endif

/* prevent from using mpfr_get_e{min,max} as lvalues */
#define mpfr_get_emin() (__gmpfr_emin + 0)
#define mpfr_get_emax() (__gmpfr_emax + 0)

#define mpfr_clear_flags() \
  ((void) (__gmpfr_flags = 0))
#define mpfr_clear_underflow() \
  ((void) (__gmpfr_flags &= MPFR_FLAGS_ALL ^ MPFR_FLAGS_UNDERFLOW))
#define mpfr_clear_overflow() \
  ((void) (__gmpfr_flags &= MPFR_FLAGS_ALL ^ MPFR_FLAGS_OVERFLOW))
#define mpfr_clear_nanflag() \
  ((void) (__gmpfr_flags &= MPFR_FLAGS_ALL ^ MPFR_FLAGS_NAN))
#define mpfr_clear_inexflag() \
  ((void) (__gmpfr_flags &= MPFR_FLAGS_ALL ^ MPFR_FLAGS_INEXACT))
#define mpfr_underflow_p() \
  ((int) (__gmpfr_flags & MPFR_FLAGS_UNDERFLOW))
#define mpfr_overflow_p() \
  ((int) (__gmpfr_flags & MPFR_FLAGS_OVERFLOW))
#define mpfr_nanflag_p() \
  ((int) (__gmpfr_flags & MPFR_FLAGS_NAN))
#define mpfr_inexflag_p() \
  ((int) (__gmpfr_flags & MPFR_FLAGS_INEXACT))

#define mpfr_round(a,b) mpfr_rint((a), (b), -1)
#define mpfr_trunc(a,b) mpfr_rint((a), (b), GMP_RNDZ)
#define mpfr_ceil(a,b)  mpfr_rint((a), (b), GMP_RNDU)
#define mpfr_floor(a,b) mpfr_rint((a), (b), GMP_RNDD)

#define mpfr_cmp_ui(b,i) mpfr_cmp_ui_2exp((b),(i),0)
#define mpfr_cmp_si(b,i) mpfr_cmp_si_2exp((b),(i),0)
#define mpfr_set(a,b,r) mpfr_set4(a,b,r,MPFR_SIGN(b))
#define mpfr_abs(a,b,r) mpfr_set4(a,b,r,1)
#define mpfr_cmp(b, c) mpfr_cmp3(b, c, 1)
#define mpfr_sgn(x) mpfr_cmp_ui(x,0)
#define mpfr_mul_2exp(y,x,n,r) mpfr_mul_2ui((y),(x),(n),(r))
#define mpfr_div_2exp(y,x,n,r) mpfr_div_2ui((y),(x),(n),(r))

#define mpfr_init_set_si(x, i, rnd) \
 ( mpfr_init(x), mpfr_set_si((x), (i), (rnd)) )
#define mpfr_init_set_ui(x, i, rnd) \
 ( mpfr_init(x), mpfr_set_ui((x), (i), (rnd)) )
#define mpfr_init_set_d(x, d, rnd) \
 ( mpfr_init(x), mpfr_set_d((x), (d), (rnd)) )
#define mpfr_init_set_z(x, i, rnd) \
 ( mpfr_init(x), mpfr_set_z((x), (i), (rnd)) )
#define mpfr_init_set_q(x, i, rnd) \
 ( mpfr_init(x), mpfr_set_q((x), (i), (rnd)) )
#define mpfr_init_set(x, y, rnd) \
 ( mpfr_init(x), mpfr_set((x), (y), (rnd)) )
#define mpfr_init_set_f(x, y, rnd) \
 ( mpfr_init(x), mpfr_set_f((x), (y), (rnd)) )

#endif
