// Copyright (c) 2001-2004  Utrecht University (The Netherlands),
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

#ifndef CGAL_MP_FLOAT_H
#define CGAL_MP_FLOAT_H

#include <CGAL/basic.h>
#include <CGAL/Interval_arithmetic.h>
#ifndef CGAL_NEW_NT_TRAITS
#  include <CGAL/Quotient.h>
#endif
#include <iostream>
#include <vector>
#include <algorithm>
#include <CGAL/Number_type_traits.h>

// MP_Float : multiprecision scaled integers.

// Some invariants on the internal representation :
// - zero is represented by an empty vector, and whatever exp.
// - no leading or trailing zero in the vector => unique

// The main algorithms are :
// - Addition/Subtraction
// - Multiplication
// - Comparison
// - to_double() / to_interval()
// - Construction from a double.
// - IOs

// TODO :
// - The exponent really overflows sometimes -> make it multiprecision.
// - Write a generic wrapper that adds an exponent to be used by MP integers.
// - Karatsuba (or other) ?  Would be fun to implement at least.
// - Division, sqrt... : different options :
//   - nothing
//   - convert to double, take approximation, compute over double, reconstruct
//   - exact division as separate function (for Bareiss...)
//   - returns the quotient of the division, forgetting the rest.

CGAL_BEGIN_NAMESPACE

class MP_Float;

MP_Float operator+(const MP_Float &a, const MP_Float &b);
MP_Float operator-(const MP_Float &a, const MP_Float &b);
MP_Float operator*(const MP_Float &a, const MP_Float &b);
MP_Float operator/(const MP_Float &a, const MP_Float &b);

#ifndef CGAL_NEW_NT_TRAITS
Comparison_result
compare (const MP_Float & a, const MP_Float & b);
#endif

class MP_Float
{
public:
#ifndef CGAL_NEW_NT_TRAITS
  typedef Tag_false  Has_gcd;
  typedef Tag_true   Has_division;
  typedef Tag_false  Has_sqrt;
#endif

  typedef short limb;
  typedef int   limb2;
  typedef double exponent_type;

  typedef std::vector<limb>  V;
  typedef V::const_iterator  const_iterator;
  typedef V::iterator        iterator;

private:

#ifdef CGAL_NEW_NT_TRAITS
public:
  static const unsigned log_limb = 8 * sizeof(limb);
private:
  static const limb2 base = 1 << log_limb;
  static const V::size_type limbs_per_double = 2 + 53/log_limb;

  static double trunc_max() {
    double value = double(base)*(base/2-1)/double(base-1);
    return value;
  }

  static double trunc_min() {
    double value = double(-base)*(base/2)/double(base-1);
    return value;
  }

public:
  Comparison_result compare(const MP_Float& b) const
  {
    if (is_zero())
      return (Comparison_result) - b.sign();
    if (b.is_zero())
      return (Comparison_result) sign();

    for (exponent_type i = std::max(max_exp(), b.max_exp()) - 1;
	 i >= std::min(min_exp(), b.min_exp()); i--)
      {
	if (of_exp(i) > b.of_exp(i))
	  return LARGER;
	if (of_exp(i) < b.of_exp(i))
	  return SMALLER;
      }
    return EQUAL;
  }

  MP_Float square() const
  {
    if (is_zero())
      return MP_Float();
    
    MP_Float r;
    r.exp = 2*exp;
    r.v.assign(2*v.size(), 0);
    for(unsigned i=0; i<v.size(); i++)
      {
	unsigned j;
	limb2 carry = 0;
	limb carry2 = 0;
	for(j=0; j<i; j++)
	  {
	    // There is a risk of overflow here :(
	    // It can only happen when a.v[i] == a.v[j] == -2^15 (log_limb...)
	    limb2 tmp0 = std::multiplies<limb2>()(v[i], v[j]);
	    limb2 tmp1 = carry + (limb2) r.v[i+j] + tmp0;
	    limb2 tmp = tmp0 + tmp1;

	    limb tmpcarry;
	    split(tmp, tmpcarry, r.v[i+j]);
	    carry = tmpcarry + (limb2) carry2;

	    // Is there a more efficient way to handle this carry ?
	    if (tmp > 0 && tmp0 < 0 && tmp1 < 0)
	      {
		// If my calculations are correct, this case should
		// never happen.
		CGAL_assertion(false);
	      }
	    else if (tmp < 0 && tmp0 > 0 && tmp1 > 0)
	      carry2 = 1;
	    else
	      carry2 = 0;
	  }
	// last round for j=i :
	limb2 tmp0 = carry + (limb2) r.v[i+i]
	  + std::multiplies<limb2>()(v[i], v[i]);
	split(tmp0, r.v[i+i+1], r.v[i+i]);
	r.v[i+i+1] += carry2;
      }
    r.canonicalize();
    return r;
  }

  MP_Float sqrt() const {
    return MP_Float( CGAL::sqrt(to_double()) );
  }

private:
  // Returns (first * 2^second), an approximation of b.
  inline
  std::pair<double, int> to_double_exp() const
  {
    if (is_zero())
      return std::make_pair(0.0, 0);

    exponent_type exp = max_exp();
    int steps = std::min(limbs_per_double, v.size());
    double d_exp_1 =
      CGAL_CLIB_STD::ldexp(1.0, - static_cast<int>(log_limb));
    double d_exp   = 1.0;
    double d = 0;

    for (exponent_type i = exp - 1; i > exp - 1 - steps; i--) {
      d_exp *= d_exp_1;
      d += d_exp * of_exp(i);
    }

    // The cast is necessary for SunPro.
    return std::make_pair(d, static_cast<int>(exp * log_limb));
  }

  // Returns (first * 2^second), an interval surrounding b.
  inline
  std::pair<std::pair<double, double>, int>
  to_interval_exp() const
  {
    if (is_zero())
      return std::make_pair(std::pair<double, double>(0, 0), 0);

    exponent_type exp = max_exp();
    int steps = std::min(limbs_per_double, v.size());
    double d_exp_1 = CGAL_CLIB_STD::ldexp(1.0, - (int) log_limb);
    double d_exp   = 1.0;

    Protect_FPU_rounding<true> P;
    Interval_nt_advanced d = 0;

    exponent_type i;
    for (i = exp - 1; i > exp - 1 - steps; i--) {
      d_exp *= d_exp_1;
      if (d_exp == 0) // Take care of underflow.
	d_exp = CGAL_IA_MIN_DOUBLE;
      d += d_exp * of_exp(i);
    }

    if (i >= min_exp() && d.is_point()) {
      if (of_exp(i) > 0)
	d += Interval_nt_advanced(0, d_exp);
      else if (of_exp(i) < 0)
	d += Interval_nt_advanced(-d_exp, 0);
      else
	d += Interval_nt_advanced(-d_exp, d_exp);
    }
  
#ifdef CGAL_EXPENSIVE_ASSERTION // force it always in early debugging
    if (d.is_point())
      CGAL_assertion(MP_Float(d.inf()) == (*this));
    else
      CGAL_assertion(MP_Float(d.inf()) <= (*this) &&
		     MP_Float(d.sup()) >= (*this));
#endif

    return std::make_pair(d.pair(), static_cast<int>(exp * log_limb));
  }

public:
  double to_double() const
  {
    std::pair<double, int> ap = to_double_exp();
    return ap.first * CGAL_CLIB_STD::ldexp(1.0, ap.second);
  }

  // FIXME : This function deserves proper testing...
  std::pair<double,double> to_interval() const
  {
    std::pair<std::pair<double, double>, int> ap = to_interval_exp();
    double scale = CGAL_CLIB_STD::ldexp(1.0, ap.second);
    return (Interval_nt<>(ap.first) * scale).pair();
  }

public:
#endif
  void remove_leading_zeros()
  {
    while (!v.empty() && v.back() == 0)
      v.pop_back();
  }

  void remove_trailing_zeros()
  {
    if (v.empty() || v.front() != 0)
      return;

    iterator i = v.begin();
    for (++i; *i == 0; ++i)
      ;
    exp += i-v.begin();
    v.erase(v.begin(), i);
  }

  // This union is used to convert an unsigned short to a short with
  // the same binary representation, without invoking implementation-defined
  // behavior (standard 4.7.3).
  // It is needed by PGCC, which behaves differently from the others.
  union to_signed {
      unsigned short us;
      short s;
  };

public:

  // Splits a limb2 into 2 limbs (high and low).
  static
  void split(limb2 l, limb & high, limb & low)
  {
    to_signed l2 = {l};
    low = l2.s;
    high = (l - low) >> (8*sizeof(limb));
  }

  // Given a limb2, returns the higher limb.
  static
  limb higher_limb(limb2 l)
  {
      limb high, low;
      split(l, high, low);
      return high;
  }

  void canonicalize()
  {
    remove_leading_zeros();
    remove_trailing_zeros();
  }

  MP_Float()
      : exp(0)
  {
    CGAL_assertion(sizeof(limb2) == 2*sizeof(limb));
    CGAL_assertion(v.empty());
    // Creates zero.
  }

#if 0
  // Causes ambiguities
  MP_Float(limb i)
  : v(1,i), exp(0)
  {
    remove_leading_zeros();
  }
#endif

  MP_Float(limb2 i)
  : v(2), exp(0)
  {
    split(i, v[1], v[0]);
    canonicalize();
  }

  MP_Float(double d);

  MP_Float operator-() const
  {
    return MP_Float() - *this;
  }

  MP_Float& operator+=(const MP_Float &a) { return *this = *this + a; }
  MP_Float& operator-=(const MP_Float &a) { return *this = *this - a; }
  MP_Float& operator*=(const MP_Float &a) { return *this = *this * a; }
  MP_Float& operator/=(const MP_Float &a) { return *this = *this / a; }

  exponent_type max_exp() const
  {
    return v.size() + exp;
  }

  exponent_type min_exp() const
  {
    return exp;
  }

  limb of_exp(exponent_type i) const
  {
    if (i < exp || i >= max_exp())
      return 0;
    return v[static_cast<int>(i-exp)];
  }

  bool is_zero() const
  {
    return v.empty();
  }

  Sign sign() const
  {
    if (v.empty())
      return ZERO;
    if (v.back()>0)
      return POSITIVE;
    CGAL_assertion(v.back()<0);
    return NEGATIVE;
  }

  void swap(MP_Float &m)
  {
    std::swap(v, m.v);
    std::swap(exp, m.exp);
  }

  V v;
  exponent_type exp;
};

inline
void swap(MP_Float &m, MP_Float &n)
{ m.swap(n); }

inline
bool operator<(const MP_Float &a, const MP_Float &b)
{ return CGAL_NTS compare(a, b) == SMALLER; }

inline
bool operator>(const MP_Float &a, const MP_Float &b)
{ return b < a; }

inline
bool operator>=(const MP_Float &a, const MP_Float &b)
{ return ! (a < b); }

inline
bool operator<=(const MP_Float &a, const MP_Float &b)
{ return ! (a > b); }

inline
bool operator==(const MP_Float &a, const MP_Float &b)
{ return (a.v == b.v) && (a.v.empty() || (a.exp == b.exp)); }

inline
bool operator!=(const MP_Float &a, const MP_Float &b)
{ return ! (a == b); }

std::ostream &
operator<< (std::ostream & os, const MP_Float &b);

// This one is for debug.
std::ostream &
print (std::ostream & os, const MP_Float &b);

std::istream &
operator>> (std::istream & is, MP_Float &b);


#ifdef CGAL_NEW_NT_TRAITS

template<>
struct Number_type_traits<MP_Float>
  : public CGALi::Default_ring_number_type_traits<MP_Float>
{
  typedef Tag_true   Has_division;
  typedef Tag_true   Has_sqrt;

  typedef Tag_true   Has_exact_ring_operations;
  typedef Tag_false  Has_exact_division;
  typedef Tag_false  Has_exact_sqrt;

  static inline Sign sign(const MP_Float& a) { return a.sign(); }

  static inline MP_Float square(const MP_Float& a) {
    return a.square();
    //    return CGALi::square(a);
  }

  static inline MP_Float sqrt(const MP_Float& a) {
    return a.sqrt();
    //    return CGALi::sqrt(a);
  }

  // to_double() returns, not the closest double, but a one
  // bit error is allowed.
  // We guarantee : to_double(MP_Float(double d)) == d.
  static inline double to_double(const MP_Float& a) {
    return a.to_double();
    //    return CGALi::to_double(a);
  }

  static inline std::pair<double,double>
  to_interval(const MP_Float &b) {
    return b.to_interval();
    //    return CGALi::to_interval(b);
  }
  
  static inline Comparison_result compare(const MP_Float& a,
					  const MP_Float& b) {
    return a.compare(b);
    //    return CGALi::compare(a, b);
  }

  static inline bool is_valid(const MP_Float&)  { return true; }
  static inline bool is_finite(const MP_Float&) { return true; }
};

#else // CGAL_NEW_NT_TRAITS

inline
Sign
sign (const MP_Float &a)
{
  return a.sign();
}

MP_Float
square(const MP_Float&);

MP_Float
sqrt(const MP_Float &d);

// to_double() returns, not the closest double, but a one bit error is allowed.
// We guarantee : to_double(MP_Float(double d)) == d.
double
to_double(const MP_Float &b);

// Overloaded in order to protect against overflow.
double
to_double(const Quotient<MP_Float> &b);

std::pair<double,double>
to_interval(const MP_Float &b);

std::pair<double, double>
to_interval(const Quotient<MP_Float> &b);

inline
void
simplify_quotient(MP_Float & numerator, MP_Float & denominator)
{
  // Currently only simplifies the two exponents.
  numerator.exp -= denominator.exp;
  denominator.exp = 0;
}

inline
bool
is_finite(const MP_Float &)
{
  return true;
}

inline
bool
is_valid(const MP_Float &)
{
  return true;
}

inline
io_Operator
io_tag(const MP_Float &)
{
  return io_Operator();
}

#endif

CGAL_END_NAMESPACE

#endif // CGAL_MP_FLOAT_H
