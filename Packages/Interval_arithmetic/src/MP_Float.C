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
// file          : src/MP_Float.C
// revision      : $Revision$
// revision_date : $Date$
// package       : Interval Arithmetic
// author(s)     : Sylvain Pion
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#include <CGAL/basic.h>
#include <CGAL/MP_Float.h>

CGAL_BEGIN_NAMESPACE

// I use macros because otherwise it's not STD compliant...
#define CGAL_MP_FLOAT_TRUNC_MAX (double(base)*(base/2-1)/double(base-1))
#define CGAL_MP_FLOAT_TRUNC_MIN (double(-base)*(base/2)/double(base-1))

MP_Float::MP_Float(double d)
{
    if (d == 0)
      return;

    CGAL_expensive_assertion(is_finite(d) && is_valid(d));
    CGAL_expensive_assertion_code(double bak = d;)

    // This is subtle, because ints are not symetric against 0.

    // First, find the exponent.
    exp = 1 - limbs_per_double;
    while (d < CGAL_MP_FLOAT_TRUNC_MIN || d > CGAL_MP_FLOAT_TRUNC_MAX) {
      exp++;
      d *= 1.0/base;
    }

    while (d >= CGAL_MP_FLOAT_TRUNC_MIN/base
  	&& d <= CGAL_MP_FLOAT_TRUNC_MAX/base) {
      exp--;
      d *= base;
    }

    // Then, compute the limbs.
    v.resize(limbs_per_double);
    for (int i = limbs_per_double - 1; i > 0; i--) {
      v[i] = (limb) ::rint(d);
      if (d-v[i] >= double(base/2-1)/(base-1))
        v[i]++;
      d -= v[i];
      d *= base;
    }

    // The last limb fits directly.
    v[0] = (limb) d;

    remove_trailing_zeros();

    CGAL_expensive_assertion(d == (limb) ::rint(d));
    CGAL_assertion(v.back() != 0);
    CGAL_expensive_assertion(CGAL::to_double(*this) == bak);
}

namespace NTS {

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

MP_Float
MP_Float::operator/(const MP_Float &d) const
{
  return MP_Float(CGAL::to_double(*this)/CGAL::to_double(d));
}

MP_Float
sqrt(const MP_Float &d)
{
  return MP_Float(CGAL_NTS sqrt(CGAL::to_double(d)));
}

// to_double() returns, not the closest double, but a one bit error is allowed.
// We guarantee : to_double(MPI(double d)) == d.
double
to_double(const MP_Float &b)
{
  if (b.is_zero())
    return 0;

  int exp = b.max_exp();
  int steps = std::min(MP_Float::limbs_per_double, b.v.size());
  double d_exp_1 = std::ldexp(1.0, (int) - MP_Float::log_limb);
  double d_exp   = std::ldexp(1.0, exp * MP_Float::log_limb);
  double d = 0;

  for (int i = exp - 1; i > exp - 1 - steps; i--) {
    d_exp *= d_exp_1;
    d += d_exp * b.of_exp(i);
  }

  return d;
}

// FIXME : This function deserves proper testing...
Interval_base
to_interval(const MP_Float &b)
{
  if (b.is_zero())
    return 0;

  int exp = b.max_exp();
  int steps = std::min(MP_Float::limbs_per_double, b.v.size());
  double d_exp_1 = std::ldexp(1.0, (int) - MP_Float::log_limb);
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
  
#ifdef CGAL_EXPENSIVE_ASSERTION // force it always in early debugging
  if (d.is_point())
    CGAL_assertion(MP_Float(d.inf()) == b);
  else
    CGAL_assertion(MP_Float(d.inf()) <= b && MP_Float(d.sup()) >= b);
#endif

  return d;
}

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

    approx += std::ldexp(double(*i), exp);

    exp += MP_Float::log_limb;
  }

  os << "  [ double approx == " << approx << " ]";

  return os;
}

std::istream &
operator>> (std::istream & is, MP_Float &b)
{
  double i;
  is >> i;
  b = MP_Float(i);
  return is;
}

CGAL_END_NAMESPACE
