/* mpf2mpfr.h -- Compatibility include file with mpf.

Copyright 1999, 2000, 2001, 2002, 2004 Free Software Foundation, Inc.

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

/* types */
#define mpf_t mpfr_t
#define mpf_srcptr mpfr_srcptr
#define mpf_ptr mpfr_ptr

/* functions which don't take as argument the rounding mode */
#undef mpf_ceil
#define mpf_ceil mpfr_ceil
#undef mpf_clear
#define mpf_clear mpfr_clear
#undef mpf_cmp
#define mpf_cmp mpfr_cmp
#undef mpf_cmp_si
#define mpf_cmp_si mpfr_cmp_si
#undef mpf_cmp_ui
#define mpf_cmp_ui mpfr_cmp_ui
#undef mpf_cmp_d
#define mpf_cmp_d mpfr_cmp_d
#undef mpf_eq
#define mpf_eq mpfr_eq
#undef mpf_floor
#define mpf_floor mpfr_floor
#undef mpf_get_prec
#define mpf_get_prec mpfr_get_prec
#undef mpf_init
#define mpf_init mpfr_init
#undef mpf_init2
#define mpf_init2 mpfr_init2
#undef mpf_integer_p
#define mpf_integer_p mpfr_integer_p
#undef mpf_random2
#define mpf_random2 mpfr_random2
#undef mpf_set_default_prec
#define mpf_set_default_prec mpfr_set_default_prec
#undef mpf_get_default_prec
#define mpf_get_default_prec mpfr_get_default_prec
#undef mpf_set_prec
#define mpf_set_prec(x,p) mpfr_set_prec(x, p)
#undef mpf_set_prec_raw
#define mpf_set_prec_raw mpfr_set_prec_raw
#undef mpf_trunc
#define mpf_trunc mpfr_trunc
#undef mpf_sgn
#define mpf_sgn mpfr_sgn
#undef mpf_swap
#define mpf_swap mpfr_swap

/* functions which take as argument the rounding mode */
#undef mpf_abs
#define mpf_abs(x,y) mpfr_abs(x,y,__gmpfr_default_rounding_mode)
#undef mpf_add
#define mpf_add(x,y,z) mpfr_add(x,y,z,__gmpfr_default_rounding_mode)
#undef mpf_add_ui
#define mpf_add_ui(x,y,z) \
             mpfr_add_ui(x,y,z,__gmpfr_default_rounding_mode)
#undef mpf_div
#define mpf_div(x,y,z) mpfr_div(x,y,z,__gmpfr_default_rounding_mode)
#undef mpf_div_ui
#define mpf_div_ui(x,y,z) \
                          mpfr_div_ui(x,y,z,__gmpfr_default_rounding_mode)
#undef mpf_div_2exp
#define mpf_div_2exp(x,y,z) \
                         mpfr_div_2exp(x,y,z,__gmpfr_default_rounding_mode)
#undef mpf_dump
#define mpf_dump(x) \
                mpfr_dump(x,__gmpfr_default_rounding_mode)
#undef mpf_fits_slong_p
#define mpf_fits_slong_p(x) mpfr_fits_ulong_p(x,__gmpfr_default_rounding_mode)
#undef mpf_fits_ulong_p
#define mpf_fits_ulong_p(x) mpfr_fits_ulong_p(x,__gmpfr_default_rounding_mode)
#undef mpf_fits_sint_p
#define mpf_fits_sint_p(x) mpfr_fits_uint_p(x,__gmpfr_default_rounding_mode)
#undef mpf_fits_uint_p
#define mpf_fits_uint_p(x) mpfr_fits_uint_p(x,__gmpfr_default_rounding_mode)
#undef mpf_fits_sshort_p
#define mpf_fits_sshort_p(x) mpfr_fits_ushort_p(x,__gmpfr_default_rounding_mode)
#undef mpf_fits_ushort_p
#define mpf_fits_ushort_p(x) mpfr_fits_ushort_p(x,__gmpfr_default_rounding_mode)
#undef mpf_get_str
#define mpf_get_str(x,y,z,t,u) \
               mpfr_get_str(x,y,z,t,u,__gmpfr_default_rounding_mode)
#undef mpf_get_d
#define mpf_get_d(x) mpfr_get_d(x,__gmpfr_default_rounding_mode)
#undef mpf_get_d_2exp
#define mpf_get_d_2exp(e,x) mpfr_get_d_2exp(e,x,__gmpfr_default_rounding_mode)
#undef mpf_get_ui
#define mpf_get_ui(x) mpfr_get_ui(x,__gmpfr_default_rounding_mode)
#undef mpf_get_si
#define mpf_get_si(x) mpfr_get_ui(x,__gmpfr_default_rounding_mode)
#undef mpf_inp_str
#define mpf_inp_str(x,y,z) mpfr_inp_str(x,y,z,__gmpfr_default_rounding_mode)
#undef mpf_set_str
#define mpf_set_str(x,y,z) mpfr_set_str(x,y,z,__gmpfr_default_rounding_mode)
#undef mpf_init_set
#define mpf_init_set(x,y) mpfr_init_set(x,y,__gmpfr_default_rounding_mode)
#undef mpf_init_set_d
#define mpf_init_set_d(x,y) mpfr_init_set_d(x,y,__gmpfr_default_rounding_mode)
#undef mpf_init_set_si
#define mpf_init_set_si(x,y) mpfr_init_set_si(x,y,__gmpfr_default_rounding_mode)
#undef mpf_init_set_str
#define mpf_init_set_str(x,y,z) mpfr_init_set_str(x,y,z,__gmpfr_default_rounding_mode)
#undef mpf_init_set_ui
#define mpf_init_set_ui(x,y) mpfr_init_set_ui(x,y,__gmpfr_default_rounding_mode)
#undef mpf_mul
#define mpf_mul(x,y,z) mpfr_mul(x,y,z,__gmpfr_default_rounding_mode)
#undef mpf_mul_2exp
#define mpf_mul_2exp(x,y,z) mpfr_mul_2exp(x,y,z,__gmpfr_default_rounding_mode)
#undef mpf_mul_ui
#define mpf_mul_ui(x,y,z) mpfr_mul_ui(x,y,z,__gmpfr_default_rounding_mode)
#undef mpf_neg
#define mpf_neg(x,y) mpfr_neg(x,y,__gmpfr_default_rounding_mode)
#undef mpf_out_str
#define mpf_out_str(x,y,z,t) mpfr_out_str(x,y,z,t,__gmpfr_default_rounding_mode)
#undef mpf_pow_ui
#define mpf_pow_ui(x,y,z) mpfr_pow_ui(x,y,z,__gmpfr_default_rounding_mode)
#undef mpf_reldiff
#define mpf_reldiff(x,y,z) mpfr_reldiff(x,y,z,__gmpfr_default_rounding_mode)
#undef mpf_set
#define mpf_set(x,y) mpfr_set(x,y,__gmpfr_default_rounding_mode)
#undef mpf_set_d
#define mpf_set_d(x,y) mpfr_set_d(x,y,__gmpfr_default_rounding_mode)
#undef mpf_set_q
#define mpf_set_q(x,y) mpfr_set_q(x,y,__gmpfr_default_rounding_mode)
#undef mpf_set_si
#define mpf_set_si(x,y) mpfr_set_si(x,y,__gmpfr_default_rounding_mode)
#undef mpf_set_ui
#define mpf_set_ui(x,y) mpfr_set_ui(x,y,__gmpfr_default_rounding_mode)
#undef mpf_set_z
#define mpf_set_z(x,y) mpfr_set_z(x,y,__gmpfr_default_rounding_mode)
#undef mpf_sqrt
#define mpf_sqrt(x,y) mpfr_sqrt(x,y,__gmpfr_default_rounding_mode)
#undef mpf_sqrt_ui
#define mpf_sqrt_ui(x,y) mpfr_sqrt_ui(x,y,__gmpfr_default_rounding_mode)
#undef mpf_sub
#define mpf_sub(x,y,z) mpfr_sub(x,y,z,__gmpfr_default_rounding_mode)
#undef mpf_sub_ui
#define mpf_sub_ui(x,y,z) mpfr_sub_ui(x,y,z,__gmpfr_default_rounding_mode)
#undef mpf_ui_div
#define mpf_ui_div(x,y,z) mpfr_ui_div(x,y,z,__gmpfr_default_rounding_mode)
#undef mpf_ui_sub
#define mpf_ui_sub(x,y,z) mpfr_ui_sub(x,y,z,__gmpfr_default_rounding_mode)
#undef mpf_urandomb
#define mpf_urandomb(x,y,n) mpfr_urandomb(x,y)
