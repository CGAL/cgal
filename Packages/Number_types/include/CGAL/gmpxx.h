// Copyright (c) 2002,2003  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Sylvain Pion
 
#ifndef CGAL_GMPXX_H
#define CGAL_GMPXX_H

#include <CGAL/basic.h>
#include <CGAL/Number_type_traits.h>
#include <CGAL/Interval_arithmetic.h>
#include <utility>

#include <gmpxx.h>
#include <mpfr.h>

// This file gathers the necessary adaptors so that the following
// C++ number types that come with GMP can be used by CGAL :
// - mpz_class
// - mpq_class

// - mpf_class support is commented out until to_interval() is implemented.
//   It is probably not very useful with CGAL anyway.

// Note that GMP++ use the expression template mechanism, which makes things
// a little bit complicated in order to make square(x+y) work for example.
// Reading gmpxx.h shows that ::__gmp_expr<T, T> is the mp[zqf]_class proper,
// while ::__gmp_expr<T, U> is the others "expressions".

CGAL_BEGIN_NAMESPACE

#ifdef CGAL_NEW_NT_TRAITS

template<>
struct Number_type_traits<mpz_class>
  : public CGALi::Default_euclidean_ring_number_type_traits<mpz_class>
{
  typedef Tag_true   Has_gcd;
  typedef Tag_false  Has_division;
  typedef Tag_false  Has_sqrt;

  typedef Tag_true   Has_exact_ring_operations;
  typedef Tag_false  Has_exact_division;
  typedef Tag_false  Has_exact_sqrt;

  static inline std::pair<double,double>
  to_interval (const mpz_class & z) {
    // GMP returns the closest double (seen in the code).
    Protect_FPU_rounding<true> P(CGAL_FE_TONEAREST);
    double app = CGAL::to_double(z);
    // If it's lower than 2^53, then it's exact.
    if (CGAL_CLIB_STD::fabs(app) < double(1<<26)*double(1<<27))
      return to_interval(app);
  
    // If the double approximation is infinite, then adding a small +/-
    // epsilon is not correct enough.
    if (! CGAL::is_finite(app)) {
      if (app > 0)
	return std::pair<double, double>(CGAL_IA_MAX_DOUBLE,CGALi::infinity);
      return std::pair<double, double>(-CGALi::infinity, -CGAL_IA_MAX_DOUBLE);
    }

    FPU_set_cw(CGAL_FE_UPWARD);
    Interval_nt<false> approx(app);
    approx += Interval_nt<false>::smallest();
    return approx.pair();
  }
};

template<>
struct Number_type_traits<mpq_class>
  : public CGALi::Default_field_number_type_traits<mpz_class>
{
  typedef Tag_false  Has_gcd;
  typedef Tag_true   Has_division;
  typedef Tag_false  Has_sqrt;

  typedef Tag_true   Has_exact_ring_operations;
  typedef Tag_true   Has_exact_division;
  typedef Tag_false  Has_exact_sqrt;

  typedef Tag_true   Has_rational_traits;

  static inline std::pair<double, double>
  to_interval (const mpq_class & q) {
    Interval_nt<> quot = Interval_nt<>(CGAL::to_interval(q.get_num())) /
      Interval_nt<>(CGAL::to_interval(q.get_den()));
    return  quot.pair();
  }
};

template <>
struct Rational_traits<mpq_class> {
  typedef mpz_class RT;
  RT numerator   (const mpq_class & r) const { return r.get_num(); }
  RT denominator (const mpq_class & r) const { return r.get_den(); }

  mpq_class make_rational(const RT & n, const RT & d) const
  { return mpq_class(n, d); } 
};

#else

template <>
struct Number_type_traits<mpz_class> {
  typedef Tag_false Has_gcd;
  typedef Tag_false Has_division;
  typedef Tag_true  Has_sqrt;


};

template <>
struct Number_type_traits<mpq_class> {
  typedef Tag_false Has_gcd;
  typedef Tag_true  Has_division;
  typedef Tag_false Has_sqrt;
};

template <>
struct Rational_traits<mpq_class> {
  typedef mpz_class RT;
  RT numerator   (const mpq_class & r) const { return r.get_num(); }
  RT denominator (const mpq_class & r) const { return r.get_den(); }

  mpq_class make_rational(const RT & n, const RT & d) const
  { return mpq_class(n, d); } 
};

#ifndef CGAL_USE_ADL_FOR_NT
template < typename T, typename U >
inline
::__gmp_expr<T, T>
sqrt(const ::__gmp_expr<T, U> &e)
{
    return ::sqrt(e);
}
#endif

template < typename T, typename U >
inline
double
to_double(const ::__gmp_expr<T, U> & e)
{ return ::__gmp_expr<T, T>(e).get_d(); }

template < typename T, typename U >
inline
bool
is_finite(const ::__gmp_expr<T, U> &)
{ return true; }

template < typename T, typename U >
inline
bool
is_valid(const ::__gmp_expr<T, U> &)
{ return true; }

template < typename T, typename U >
inline
io_Operator
io_tag(const ::__gmp_expr<T, U> &)
{ return io_Operator(); }

template < typename T, typename U >
std::pair<double,double>
to_interval (const ::__gmp_expr<T, U> & z)
{
  // Calls the functions below after dealing with the expression template.
  return to_interval(::__gmp_expr<T, T>(z));
}

inline
std::pair<double, double>
to_interval (const mpz_class & z)
{
  mpfr_t x;
  mpfr_init2 (x, 53); /* Assume IEEE-754 */
  mpfr_set_z (x, z.get_mpz_t(), GMP_RNDD);
  double i = mpfr_get_d (x, GMP_RNDD); /* EXACT but can overflow */
  mpfr_set_z (x, z.get_mpz_t(), GMP_RNDU);
  double s = mpfr_get_d (x, GMP_RNDU); /* EXACT but can overflow */
  mpfr_clear (x);
  return std::pair<double, double>(i, s);
}

inline
std::pair<double, double>
to_interval (const mpq_class & q)
{
  mpfr_t x;
  mpfr_init2 (x, 53); /* Assume IEEE-754 */
  mpfr_set_q (x, q.get_mpq_t(), GMP_RNDD);
  double i = mpfr_get_d (x, GMP_RNDD); /* EXACT but can overflow */
  mpfr_set_q (x, q.get_mpq_t(), GMP_RNDU);
  double s = mpfr_get_d (x, GMP_RNDU); /* EXACT but can overflow */
  mpfr_clear (x);
  return std::pair<double, double>(i, s);
}

#ifndef CGAL_USE_ADL_FOR_NT
// These are necessary due to expression-templates.
template < typename T, typename U >
inline
::__gmp_expr<T, T>
abs(const ::__gmp_expr<T, U>& x) { return ::abs(x); }
#endif

template < typename T, typename U >
inline
::__gmp_expr<T, T>
square(const ::__gmp_expr<T, U>& x) { return x*x; }

template < typename T, typename U >
inline
Sign
sign(const ::__gmp_expr<T, U> & e)
{ return (Sign) ::sgn(e); }

template < typename T, typename U1, typename U2 >
inline
Comparison_result
compare(const ::__gmp_expr<T, U1> & e1,
        const ::__gmp_expr<T, U2> & e2)
{
  // cmp returns any int value, not just -1/0/1...
  return (Comparison_result) CGAL_NTS sign(::cmp(e1, e2));
}

template < typename T, typename U >
inline
bool
is_zero(const ::__gmp_expr<T, U> & e)
{ return ::sgn(e) == 0; }

template < typename T, typename U >
inline
bool
is_one(const ::__gmp_expr<T, U> & e)
{ return e == 1; }

template < typename T, typename U >
inline
bool
is_positive(const ::__gmp_expr<T, U> & e)
{ return ::sgn(e) > 0; }

template < typename T, typename U >
inline
bool
is_negative(const ::__gmp_expr<T, U> & e)
{ return ::sgn(e) < 0; }

#endif

CGAL_END_NAMESPACE

#endif // CGAL_GMPXX_H
