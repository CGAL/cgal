// ============================================================================
//
// Copyright (c) 2001 The CGAL Consortium
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
// file          : include/CGAL/MP_Float.h
// revision      : $Revision$
// revision_date : $Date$
// package       : Number_types
// author(s)     : Sylvain Pion
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_MP_FLOAT_H
#define CGAL_MP_FLOAT_H

#include <CGAL/basic.h>
#include <CGAL/Interval_arithmetic.h>
#include <iostream>
#include <vector>
#include <functional>
#include <cmath>

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
// - Documentation.
// - Karatsuba (or other) ?  Would be fun to implement at least.
// - Division, sqrt... : different options :
//   - nothing
//   - convert to double, take approximation, compute over double, reconstruct
//   - exact division as separate function (for Bareiss...)
//   - returns the quotient of the division, forgetting the rest.

CGAL_BEGIN_NAMESPACE

class MP_Float;

Comparison_result
compare_noinline (const MP_Float & a, const MP_Float & b);

namespace NTS {
inline
Comparison_result
compare (const MP_Float & a, const MP_Float & b)
{
  return compare_noinline (a, b);
}
}

class MP_Float
{
public:

  typedef short limb;
  typedef int   limb2;

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

public:

  void canonicalize()
  {
    remove_leading_zeros();
    remove_trailing_zeros();
  }

  // Accessory function : Given a limb2, this returns the higher limb.
  // The lower limb is simply obtained by casting to a limb.
  static
  limb higher_limb(limb2 l)
  {
    return (l - (limb) l) >> (8*sizeof(limb));
  }

  MP_Float()
  {
    CGAL_assertion(sizeof(limb2) == 2*sizeof(limb));
    CGAL_assertion(v.empty());
    // Creates zero.
  }

  MP_Float(limb i)
  : v(1,i), exp(0)
  {
    remove_leading_zeros();
  }

  MP_Float(limb2 i)
  : v(2), exp(0)
  {
    v[0] = i;
    v[1] = higher_limb(i);
    canonicalize();
  }

  MP_Float(double d);

  MP_Float operator-() const
  {
    return MP_Float() - *this;
  }

  MP_Float operator+(const MP_Float &) const;
  MP_Float operator-(const MP_Float &) const;
  MP_Float operator*(const MP_Float &) const;
  MP_Float operator/(const MP_Float &) const;

  MP_Float& operator+=(const MP_Float &a)
  {
    return *this = *this + a;
  }

  MP_Float& operator-=(const MP_Float &a)
  {
    return *this = *this - a;
  }

  MP_Float& operator*=(const MP_Float &a)
  {
    return *this = *this * a;
  }

  MP_Float& operator/=(const MP_Float &a)
  {
    return *this = *this / a;
  }

  bool operator<(const MP_Float &b) const
  {
      return CGAL_NTS compare(*this, b) == SMALLER;
  }

  bool operator>(const MP_Float &f) const
  {
      return f<*this;
  }

  bool operator>=(const MP_Float &f) const
  {
      return !(*this<f);
  }

  bool operator<=(const MP_Float &f) const
  {
      return !(*this>f);
  }

  bool operator==(const MP_Float &b) const
  {
    return (v == b.v) && ((exp == b.exp) || v.empty());
  }

  bool operator!=(const MP_Float &b) const
  {
    return ! (*this == b);
  }

  int max_exp() const
  {
    return v.size() + exp;
  }

  int min_exp() const
  {
    return exp;
  }

  limb of_exp(int i) const
  {
    if (i < exp || i >= max_exp())
      return 0;
    return v[i-exp];
  }

  bool is_zero() const
  {
    return v.empty();
  }

  bool is_one() const
  {
    return v.size() == 1 && v.back() == 1;
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

  V v;
  int exp;
};


namespace NTS {

inline
Sign
sign (const MP_Float &a)
{
  return a.sign();
}

} // namespace NTS

MP_Float
sqrt(const MP_Float &d);

// to_double() returns, not the closest double, but a one bit error is allowed.
// We guarantee : to_double(MPI(double d)) == d.
double
to_double(const MP_Float &b);

Interval_base
to_interval(const MP_Float &b);

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

inline
Number_tag
number_type_tag(const MP_Float &)
{
  return Number_tag();
}

std::ostream &
operator<< (std::ostream & os, const MP_Float &b);

// This one is for debug.
std::ostream &
print (std::ostream & os, const MP_Float &b);

std::istream &
operator>> (std::istream & is, MP_Float &b);

CGAL_END_NAMESPACE

#endif // CGAL_MP_FLOAT_H
