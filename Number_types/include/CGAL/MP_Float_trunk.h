// Copyright (c) 2001-2006  INRIA Sophia-Antipolis (France).
// All rights reserved.
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
// $URL: svn+ssh://hemmer@scm.gforge.inria.fr/svn/cgal/trunk/Number_types/include/CGAL/MP_Float.h $
// $Id: MP_Float.h 34597 2006-09-29 08:45:07Z pmachado $
// 
//
// Author(s)     : Sylvain Pion

#ifndef CGAL_MP_FLOAT_H
#define CGAL_MP_FLOAT_H

#include <CGAL/basic.h>
#include <CGAL/Interval_nt.h>
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
// - Integral division div(), exact_division(), gcd(), operator%().
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

CGAL_BEGIN_NAMESPACE

class MP_Float;

MP_Float operator+(const MP_Float &a, const MP_Float &b);
MP_Float operator-(const MP_Float &a, const MP_Float &b);
MP_Float operator*(const MP_Float &a, const MP_Float &b);
#ifdef CGAL_MP_FLOAT_ALLOW_INEXACT
MP_Float operator/(const MP_Float &a, const MP_Float &b);
#endif
MP_Float operator%(const MP_Float &a, const MP_Float &b);

Comparison_result
compare (const MP_Float & a, const MP_Float & b);

class MP_Float
{
public:
  typedef Tag_true   Has_gcd;
#ifdef CGAL_MP_FLOAT_ALLOW_INEXACT
  typedef Tag_true   Has_division;
  typedef Tag_true   Has_sqrt;
#else
  typedef Tag_false  Has_division;
  typedef Tag_false  Has_sqrt;
#endif

  typedef Tag_true   Has_exact_ring_operations;
  typedef Tag_true   Has_exact_division;
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

#ifdef CGAL_ROOT_OF_2_ENABLE_HISTOGRAM_OF_NUMBER_OF_DIGIT_ON_THE_COMPLEX_CONSTRUCTOR
  int tam() const { return v.size(); }
#endif

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
  MP_Float& operator%=(const MP_Float &a) { return *this = *this % a; }

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
#elif 1
  numerator.exp -= denominator.exp;
  denominator.exp = 0;
#else
  if (numerator != 0 && denominator != 0) {
    numerator.exp -= denominator.exp;
    denominator.exp = 0;
    const MP_Float g = gcd(numerator, denominator);
    numerator = exact_division(numerator, g);
    denominator = exact_division(denominator, g);
  }
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

// TODO : needs function div(), with remainder.

namespace CGALi {

// This compares the absolute values of the odd-mantissa.
// (take the mantissas, get rid of all powers of 2, compare
// the absolute values)
inline
Sign
compare_bitlength(const MP_Float &a, const MP_Float &b)
{
  if (a.is_zero())
    return b.is_zero() ? EQUAL : SMALLER;
  if (b.is_zero())
    return LARGER;

  MP_Float aa = abs(a);
  MP_Float bb = abs(b);

  if (aa.size() > (bb.size() + 2)) return LARGER;
  if (bb.size() > (aa.size() + 2)) return SMALLER;

  // multiply by 2 till last bit is 1.
  while (((aa.v[0]) & 1) == 0) // last bit is zero
    aa = aa + aa;

  while (((bb.v[0]) & 1) == 0) // last bit is zero
    bb = bb + bb;

  // sizes might have changed
  if (aa.size() > bb.size()) return LARGER;
  if (aa.size() < bb.size()) return SMALLER;

  for (int i = aa.size()-1; i >= 0; --i)
  {
    if (aa.v[i] > bb.v[i]) return LARGER;
    if (aa.v[i] < bb.v[i]) return SMALLER;
  }
  return EQUAL;
}

inline // Move it to libCGAL once it's stable.
std::pair<MP_Float, MP_Float> // <quotient, remainder>
division(const MP_Float & n, const MP_Float & d)
{
  typedef MP_Float::exponent_type  exponent_type;

  MP_Float remainder = n, divisor = d;

  CGAL_precondition(divisor != 0);

  // Rescale d to have a to_double() value with reasonnable exponent.
  exponent_type scale_d = divisor.find_scale();
  divisor.rescale(scale_d);
  const double dd = to_double(divisor);

  MP_Float res = 0;
  exponent_type scale_remainder = 0;

  bool first_time_smaller_than_divisor = true;

  // School division algorithm.

  while ( remainder != 0 )
  {
    // We have to rescale, since remainder can diminish towards 0.
    exponent_type tmp_scale = remainder.find_scale();
    remainder.rescale(tmp_scale);
    res.rescale(tmp_scale);
    scale_remainder += tmp_scale;

    // Compute a double approximation of the quotient
    // (imagine school division with base ~2^53).
    double approx = to_double(remainder) / dd;
    CGAL_assertion(approx != 0);
    res += approx;
    remainder -= approx * divisor;

    if (remainder == 0)
      break;

    // Then we need to fix it up by checking if neighboring double values
    // are closer to the exact result.
    // There should not be too many iterations, because approx is only a few ulps
    // away from the optimal.
    // If we don't do the fixup, then spurious bits can be introduced, which
    // will require an unbounded amount of additional iterations to be eliminated.

    // The direction towards which we need to try to move from "approx".
    double direction = (sign(remainder) == sign(dd))
                     ?  std::numeric_limits<double>::infinity()
                     : -std::numeric_limits<double>::infinity();

    while (true)
    {
      const double approx2 = nextafter(approx, direction);
      const double delta = approx2 - approx;
      MP_Float new_remainder = remainder - delta * divisor;
      if (abs(new_remainder) < abs(remainder)) {
        remainder = new_remainder;
        res += delta;
        approx = approx2;
      }
      else {
        break;
      }
    }

    if (remainder == 0)
      break;

    // Test condition for non-exact division (with remainder).
    if (compare_bitlength(remainder, divisor) == SMALLER)
    {
      if (! first_time_smaller_than_divisor)
      {
        // Scale back.
        res.rescale(scale_d - scale_remainder);
        remainder.rescale(- scale_remainder);
        CGAL_postcondition(res * d  + remainder == n);
        return std::make_pair(res, remainder);
      }
      first_time_smaller_than_divisor = false;
    }
  }

  // Scale back the result.
  res.rescale(scale_d - scale_remainder);
  CGAL_postcondition(res * d == n);
  return std::make_pair(res, MP_Float(0));
}

} // namespace CGALi

inline // Move it to libCGAL once it's stable.
MP_Float
exact_division(const MP_Float & n, const MP_Float & d)
{
  std::pair<MP_Float, MP_Float> res = CGALi::division(n, d);
  CGAL_assertion_msg(res.second == 0,
                  "exact_division() called with operands which do not divide");
  return res.first;
}

inline // Move it to libCGAL once it's stable.
bool
divides(const MP_Float & n, const MP_Float & d)
{
  return CGALi::division(n, d).second == 0;
}

inline
MP_Float
operator%(const MP_Float& n1, const MP_Float& n2)
{
  return CGALi::division(n1, n2).second;
}

inline
MP_Float
div(const MP_Float& n1, const MP_Float& n2)
{
  return CGALi::division(n1, n2).first;
}

inline
MP_Float
gcd( const MP_Float& a, const MP_Float& b)
{
  CGAL_precondition( a != 0 );
  CGAL_precondition( b != 0 );
  MP_Float x = a, y = b;
  while (true) {
    x = x % y;
    if (x == 0) {
      CGAL_postcondition(divides(a, y) && divides(b, y));
      return y;
    }
    swap(x, y);
  }
}

CGAL_END_NAMESPACE

#endif // CGAL_MP_FLOAT_H
