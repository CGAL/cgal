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
// package       : Interval Arithmetic
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
#include <utility>
#include <functional>
#include <cmath>

// MP_Float : multiprecision scaled integer.
// Not optimized for large operands at the moment.

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
// - Should the MP_Float unconditionaly replace MP_Integer ?
// - the non-inline stuff must go in src/MP_Float.C
//   (when the code will be stabilized/ported).
// - implement missing CGAL requirements.
// - Documentation.
// - Karatsuba (or other) ?  Would be fun to implement at least.
// - Division, sqrt... : different options :
//   - nothing
//   - convert to double, take approximation, compute over double, reconstruct
//   - exact division as separate function (for Bareiss...)
//   - returns the quotient of the division, forgetting the rest.

CGAL_BEGIN_NAMESPACE

class MP_Float;
double to_double(const MP_Float &);

class MP_Float
{
public:

  typedef short limb;
  typedef int   limb2;

  static const unsigned log_limb = 8*sizeof(limb);
  static const limb2    base     = 1<<log_limb;
  static const unsigned limbs_per_double = 2 + 53/log_limb;
  static const double trunc_max = double(base)*(base/2-1)/double(base-1);
  static const double trunc_min = double(-base)*(base/2)/double(base-1);

  typedef std::vector<limb>  V;
  typedef V::const_iterator  const_iterator;
  typedef V::iterator        iterator;

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

  MP_Float(double d)
  {
    if (d == 0)
      return;

    CGAL_expensive_assertion(is_finite(d) && is_valid(d));
    CGAL_expensive_assertion_code(double bak = d;)

    // This is subtle, because ints are not symetric against 0.

    // First, find the exponent.
    exp = 1 - limbs_per_double;
    while (d < trunc_min || d > trunc_max) {
      exp++;
      d *= 1.0/base;
    }

    while (d >= trunc_min/base && d <= trunc_max/base) {
      exp--;
      d *= base;
    }

    // Then, compute the limbs.
    v.resize(limbs_per_double);
    for (int i = limbs_per_double - 1; i > 0; i--) {
      v[i] = (limb) std::rint(d);
      if (d-v[i] >= double(base/2-1)/(base-1))
        v[i]++;
      d -= v[i];
      d *= base;
    }

    // The last limb fits directly.
    v[0] = (limb) d;

    remove_trailing_zeros();

    CGAL_expensive_assertion(d == (limb) std::rint(d));
    CGAL_assertion(v.back() != 0);
    CGAL_expensive_assertion(CGAL::to_double(*this) == bak);
  }

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

  bool operator<(const MP_Float &) const;

  bool operator==(const MP_Float &b) const
  {
    return (v == b.v) && ((exp == b.exp) || v.empty());
  }

  bool operator!=(const MP_Float &b) const
  {
    return ! (*this == b);
  }

  void canonicalize()
  {
    remove_leading_zeros();
    remove_trailing_zeros();
  }

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
    for (i++; *i == 0; i++)
      ;
    exp += i-v.begin();
    v.erase(v.begin(), i);
  }

public:

  limb of_exp(int i) const
  {
    if (i < exp || i >= max_exp())
      return 0;
    return v[i-exp];
  }

  int max_exp() const
  {
    return v.size() + exp;
  }

  int min_exp() const
  {
    return exp;
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

  // Accessory function : Given a limb2, this returns the higher limb.
  // The lower limb is simply obtained by casting to a limb.
  static
  limb higher_limb(limb2 l)
  {
    return (l - (limb) l) >> log_limb;
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

inline
Comparison_result
compare (const MP_Float & a, const MP_Float & b)
{
  if (a.is_zero())
    return (Comparison_result) - b.sign();
  if (b.is_zero())
    return (Comparison_result) a.sign();

  for (int i = std::max(a.max_exp(), b.max_exp()) - 1;
          i >= std::min(a.min_exp(), b.min_exp()); i--)
  {
    if (a.of_exp(i) > b.of_exp(i))
      return LARGER;
    if (a.of_exp(i) < b.of_exp(i))
      return SMALLER;
  }
  return EQUAL;
}

} // namespace NTS

inline
bool
MP_Float::operator<(const MP_Float & b) const
{
  return CGAL_NTS compare(*this, b) == SMALLER;
}

// Common code for operator+ and operator-.
template <class BinOp>
inline
void
Add_Sub(MP_Float &r, const MP_Float &a, const MP_Float &b)
{
  CGAL_assertion(!a.is_zero() && !b.is_zero());

  int min_exp = std::min(a.min_exp(), b.min_exp());
  int max_exp = std::max(a.max_exp(), b.max_exp());
  r.exp = min_exp;
  r.v.resize(max_exp - min_exp + 1); // One more for the carry.
  r.v[0] = 0;
  for(int i = 0; i < max_exp - min_exp; i++)
  {
    MP_Float::limb2 tmp = r.v[i] + BinOp()(a.of_exp(i+min_exp),
                                             b.of_exp(i+min_exp));
    r.v[i] = tmp;
    r.v[i+1] = MP_Float::higher_limb(tmp);
  }
  r.canonicalize();
}

// a few of the following should not be inline.
inline
MP_Float
MP_Float::operator+(const MP_Float &b) const
{
  if (is_zero())
    return b;
  if (b.is_zero())
    return *this;

  MP_Float r;
  Add_Sub<std::plus<limb2> >(r, *this, b);
  return r;
}

inline
MP_Float
MP_Float::operator-(const MP_Float &b) const
{
  if (is_zero())
    return -b;
  if (b.is_zero())
    return *this;

  MP_Float r;
  Add_Sub<std::minus<limb2> >(r, *this, b);
  return r;
}

inline
MP_Float
MP_Float::operator*(const MP_Float &b) const
{
  MP_Float r;
  if (is_zero() || b.is_zero())
    return r;

  r.exp = exp + b.exp;
  r.v.assign(v.size() + b.v.size(), 0);
  for(unsigned i=0; i<v.size(); i++)
  {
    unsigned j;
    limb2 carry = 0;
    for(j=0; j<b.v.size(); j++)
    {
      limb2 tmp = carry + (limb2) r.v[i+j]
                        + std::multiplies<limb2>()(v[i], b.v[j]);
      r.v[i+j] = tmp;
      carry = higher_limb(tmp);
    }
    r.v[i+j] = carry;
  }
  r.canonicalize();
  return r;
}

inline
MP_Float
MP_Float::operator/(const MP_Float &d) const
{
  return MP_Float(CGAL::to_double(*this)/CGAL::to_double(d));
}

inline
MP_Float
sqrt(const MP_Float &d)
{
  return MP_Float(CGAL_NTS sqrt(CGAL::to_double(d)));
}

// to_double() returns, not the closest double, but a one bit error is allowed.
// We guarantee : to_double(MPI(double d)) == d.
inline
double
to_double(const MP_Float &b)
{
  if (b.is_zero())
    return 0;

  int exp = b.max_exp();
  int steps = std::min(MP_Float::limbs_per_double, b.v.size());
  double d_exp_1 = std::ldexp(1.0,     - MP_Float::log_limb);
  double d_exp   = std::ldexp(1.0, exp * MP_Float::log_limb);
  double d = 0;

  for (int i = exp - 1; i > exp - 1 - steps; i--) {
    d_exp *= d_exp_1;
    d += d_exp * b.of_exp(i);
  }

  return d;
}

// FIXME : This function deserves proper testing...
inline
Interval_base
to_interval(const MP_Float &b)
{
  if (b.is_zero())
    return 0;

  int exp = b.max_exp();
  int steps = std::min(MP_Float::limbs_per_double, b.v.size());
  double d_exp_1 = std::ldexp(1.0,     - MP_Float::log_limb);
  double d_exp   = std::ldexp(1.0, exp * MP_Float::log_limb);

  // We take care of overflow.  The following should be enough.
  if (!CGAL_NTS is_finite(d_exp))
    return Interval_base::Largest;

  Protect_FPU_rounding<true> P;
  Interval_nt_advanced d = 0;

  int i;
  for (i = exp - 1; i > exp - 1 - steps; i--) {
    d_exp *= d_exp_1;
    if (d_exp == 0) // Take care of underflow.
      d_exp = CGAL_IA_MIN_DOUBLE;
    d += d_exp * b.of_exp(i);
  }

  if (i >= b.min_exp() && d.is_point()) {
    if (b.of_exp(i) > 0)
      d += Interval_nt_advanced(0, d_exp);
    else if (b.of_exp(i) < 0)
      d += Interval_nt_advanced(-d_exp, 0);
    else
      d += Interval_nt_advanced(-d_exp, d_exp);
  }
  
// #ifdef CGAL_EXPENSIVE_ASSERTION
#ifdef CGAL_ASSERTION // force it always in early debugging
  if (d.is_point())
    CGAL_assertion(MP_Float(d.inf()) == b);
  else
    CGAL_assertion(MP_Float(d.inf()) <= b && MP_Float(d.sup()) >= b);
#endif

  return d;
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

inline
Number_tag
number_type_tag(const MP_Float &)
{
  return Number_tag();
}

inline
std::ostream &
operator<< (std::ostream & os, const MP_Float &b)
{
  // Binary format would be nice and not hard to have too (useful ?).
  if (b.is_zero())
    return os << 0;

  MP_Float::const_iterator i;
  int exp = b.min_exp() * MP_Float::log_limb;
  double approx = 0; // only for giving an idea.

  for (i = b.v.begin(); i != b.v.end(); i++)
  {
    os << ((*i > 0) ? " +" : " ") << *i;

    if (exp != 0)
      os << " * 2^" << exp;

    approx += std::ldexp(*i, exp);

    exp += MP_Float::log_limb;
  }

  os << "  [ double approx == " << approx << " ]";

  return os;
}

inline
std::istream &
operator>> (std::istream & is, MP_Float &b)
{
  double i;
  is >> i;
  b = MP_Float(i);
  return is;
}

CGAL_END_NAMESPACE

#endif // CGAL_MP_FLOAT_H
