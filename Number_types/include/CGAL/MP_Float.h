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
// $URL$
// $Id$
// 
//
// Author(s)     : Sylvain Pion

#ifndef CGAL_MP_FLOAT_H
#define CGAL_MP_FLOAT_H

#include <CGAL/basic.h>
#include <CGAL/Interval_nt.h>
#include <CGAL/Root_of_2.h> // needed because of the [] operation on Root_of_2 
#include <CGAL/MP_Float_fwd.h>
#include <iostream>
#include <vector>
#include <algorithm>

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
#ifdef CGAL_MP_FLOAT_ALLOW_INEXACT
MP_Float operator/(const MP_Float &a, const MP_Float &b);
#endif

Comparison_result
compare (const MP_Float & a, const MP_Float & b);

class MP_Float
{
public:
  typedef Tag_false  Has_gcd;
#ifdef CGAL_MP_FLOAT_ALLOW_INEXACT
  typedef Tag_true   Has_division;
  typedef Tag_true   Has_sqrt;
#else
  typedef Tag_false  Has_division;
  typedef Tag_false  Has_sqrt;
#endif

  typedef Tag_true   Has_exact_ring_operations;
  typedef Tag_false  Has_exact_division;
  typedef Tag_false  Has_exact_sqrt;

  typedef short      limb;
  typedef int        limb2;
  typedef double     exponent_type;

  typedef std::vector<limb>  V;
  typedef V::const_iterator  const_iterator;
  typedef V::iterator        iterator;

private:

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

  // The constructors from float/double/long_double are factorized in the
  // following template :
  template < typename T >
  void construct_from_builtin_fp_type(T d);

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

  MP_Float(float d);

  MP_Float(double d);

  MP_Float(long double d);

  MP_Float operator-() const
  {
    return MP_Float() - *this;
  }

  MP_Float& operator+=(const MP_Float &a) { return *this = *this + a; }
  MP_Float& operator-=(const MP_Float &a) { return *this = *this - a; }
  MP_Float& operator*=(const MP_Float &a) { return *this = *this * a; }
#ifdef CGAL_MP_FLOAT_ALLOW_INEXACT
  MP_Float& operator/=(const MP_Float &a) { return *this = *this / a; }
#endif

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

  // Converts to a rational type (e.g. Gmpq).
  template < typename T >
  T to_rational() const
  {
    const unsigned log_limb = 8 * sizeof(MP_Float::limb);

    if (is_zero())
      return 0;

    MP_Float::const_iterator i;
    exponent_type exp = min_exp() * log_limb;
    T res = 0;

    for (i = v.begin(); i != v.end(); i++)
    {
      res += CGAL_CLIB_STD::ldexp(static_cast<double>(*i),
                                  static_cast<int>(exp));
      exp += log_limb;
    }

    return res;
  }

  std::size_t size() const
  {
    return v.size();
  }

  // Returns a scaling factor (in limbs) which would be good to extract to get
  // a value with an exponent close to 0.
  exponent_type find_scale() const
  {
    return exp + v.size();
  }

  // Rescale the value by some factor (in limbs).  (substract the exponent)
  void rescale(exponent_type scale)
  {
    if (v.size() != 0)
      exp -= scale;
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

inline
Sign
sign (const MP_Float &a)
{
  return a.sign();
}

MP_Float
square(const MP_Float&);

MP_Float
approximate_sqrt(const MP_Float &d);

MP_Float
approximate_division(const MP_Float &n, const MP_Float &d);

#ifdef CGAL_MP_FLOAT_ALLOW_INEXACT
inline
MP_Float
operator/(const MP_Float &a, const MP_Float &b)
{
  return approximate_division(a, b);
}

inline
MP_Float
sqrt(const MP_Float &d)
{
  return approximate_sqrt(d);
}
#endif

// to_double() returns, not the closest double, but a one bit error is allowed.
// We guarantee : to_double(MP_Float(double d)) == d.
double
to_double(const MP_Float &b);

std::pair<double,double>
to_interval(const MP_Float &b);

template < typename > class Quotient;

// Overloaded in order to protect against overflow.
double
to_double(const Quotient<MP_Float> &b);

std::pair<double, double>
to_interval(const Quotient<MP_Float> &b);

std::pair<double, int>
to_double_exp(const MP_Float &b);

// Returns (first * 2^second), an interval surrounding b.
std::pair<std::pair<double, double>, int>
to_interval_exp(const MP_Float &b);

inline
void
simplify_quotient(MP_Float & numerator, MP_Float & denominator)
{
  // Currently only simplifies the two exponents.
#if 0
  // This better version causes problems as the I/O is changed for
  // Quotient<MP_Float>, which then does not appear as rational 123/345,
  // 1.23/3.45, this causes problems in the T2 test-suite (to be investigated).
  numerator.exp -= denominator.exp 
                    + (MP_Float::exponent_type) denominator.v.size();
  denominator.exp = - (MP_Float::exponent_type) denominator.v.size();
#else
  numerator.exp -= denominator.exp;
  denominator.exp = 0;
#endif
}

inline void simplify_root_of_2(MP_Float &/*a*/, MP_Float &/*b*/, MP_Float&/*c*/) {
#if 0
  if(CGAL::is_zero(a)) {
  	simplify_quotient(b,c); return;
  } else if(CGAL::is_zero(b)) {
  	simplify_quotient(a,c); return;
  } else if(CGAL::is_zero(c)) {
  	simplify_quotient(a,b); return;
  }  	
  MP_Float::exponent_type va = a.exp + 
    (MP_Float::exponent_type) a.v.size();	
  MP_Float::exponent_type vb = b.exp + 
    (MP_Float::exponent_type) b.v.size();
  MP_Float::exponent_type vc = c.exp + 
    (MP_Float::exponent_type) c.v.size();
  MP_Float::exponent_type min = (std::min)((std::min)(va,vb),vc);	
  MP_Float::exponent_type max = (std::max)((std::max)(va,vb),vc);
  MP_Float::exponent_type med = (min+max)/2.0;
  a.exp -= med;
  b.exp -= med;
  c.exp -= med;  
#endif    	
}

namespace CGALi {
  inline void simplify_3_exp(int &a, int &b, int &c) {
    int min = (std::min)((std::min)(a,b),c);	
    int max = (std::max)((std::max)(a,b),c);
    int med = (min+max)/2;
    a -= med;
    b -= med;
    c -= med;
  }
}

inline
double
to_double(const Root_of_2<MP_Float> &x)
{
  typedef MP_Float RT;
  typedef Quotient<RT> FT;
  typedef CGAL::Rational_traits< FT > Rational;
  Rational r;
  const RT r1 = r.numerator(x.alpha());
  const RT d1 = r.denominator(x.alpha());

  if(x.is_rational()) {
    std::pair<double, int> n = to_double_exp(r1);
    std::pair<double, int> d = to_double_exp(d1);
    double scale = CGAL_CLIB_STD::ldexp(1.0, n.second - d.second);
    return (n.first / d.first) * scale;
  }

  const RT r2 = r.numerator(x.beta());
  const RT d2 = r.denominator(x.beta());
  const RT r3 = r.numerator(x.gamma());
  const RT d3 = r.denominator(x.gamma());

  std::pair<double, int> n1 = to_double_exp(r1);
  std::pair<double, int> v1 = to_double_exp(d1);
  double scale1 = CGAL_CLIB_STD::ldexp(1.0, n1.second - v1.second);

  std::pair<double, int> n2 = to_double_exp(r2);
  std::pair<double, int> v2 = to_double_exp(d2);
  double scale2 = CGAL_CLIB_STD::ldexp(1.0, n2.second - v2.second);

  std::pair<double, int> n3 = to_double_exp(r3);
  std::pair<double, int> v3 = to_double_exp(d3);
  double scale3 = CGAL_CLIB_STD::ldexp(1.0, n3.second - v3.second);

  return ((n1.first / v1.first) * scale1) + 
         ((n2.first / v2.first) * scale2) *
         std::sqrt((n3.first / v3.first) * scale3);
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

std::ostream &
operator<< (std::ostream & os, const MP_Float &b);

// This one is for debug.
std::ostream &
print (std::ostream & os, const MP_Float &b);

std::istream &
operator>> (std::istream & is, MP_Float &b);

// TODO : needs function "bool divides(n, d)" for validity checking, and maybe other things.
// TODO : needs function div(), with remainder.

namespace CGALi {

inline // Move it to libCGAL once it's stable.
MP_Float
exact_division_internal(MP_Float n, MP_Float d, bool & divides)
{
  CGAL_assertion(d != 0);

  typedef MP_Float::exponent_type  exponent_type;
  // Rescale them to have to_double() values with reasonnable exponents.
  exponent_type scale_n = n.find_scale();
  exponent_type scale_d = d.find_scale();
  n.rescale(scale_n);
  d.rescale(scale_d);

  // A simple criteria for detecting when the division is not exact
  // is that if it is exact, then the result must have smaller or
  // equal bit length than "n".
  // Which we approximate with size() and a confortable margin.
  exponent_type max_size_if_exact = n.size() + d.size() + 200;

  // School division algorithm.
  MP_Float res = to_double(n) / to_double(d);
  MP_Float remainder = n - res * d;
  while ( remainder != 0 )
  {
    // We also have to rescale here, since remainder can diminish towards 0.
    exponent_type scale_rem = remainder.find_scale();
    remainder.rescale(scale_rem);
    res.rescale(scale_rem);
    scale_n += scale_rem;

    // A double approximation of the quotient
    // (imagine school division with base ~2^53).
    double approx = to_double(remainder) / to_double(d);
    CGAL_assertion(approx != 0);
    approx = (approx + (4*approx)) - (4*approx); // chop-off the last bit.
    res += approx;
    remainder -= approx * d;
    if (res.size() > max_size_if_exact)
    {
      divides = false;
      return MP_Float();
    }
  }

  divides = true;
  // Scale back the result.
  res.rescale(scale_d - scale_n);
  return res;
}

} // namespace CGALi

inline // Move it to libCGAL once it's stable.
MP_Float
exact_division(const MP_Float & n, const MP_Float & d)
{
  bool exact_div;
  MP_Float res = CGALi::exact_division_internal(n, d, exact_div);
  CGAL_assertion_msg(exact_div, "exact_division() called with operands which do not divide");
  return res;
}

inline // Move it to libCGAL once it's stable.
bool
divides(const MP_Float & n, const MP_Float & d)
{
  bool exact_div;
  CGALi::exact_division_internal(n, d, exact_div);
  return exact_div;
}

CGAL_END_NAMESPACE

#endif // CGAL_MP_FLOAT_H
