// ============================================================================
//
// Copyright (c) 2000,2001 The CGAL Consortium
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
// file          : include/CGAL/MP_Integer.h
// revision      : $Revision$
// revision_date : $Date$
// package       : Interval Arithmetic
// author(s)     : Sylvain Pion
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_MP_INTEGER_H
#define CGAL_MP_INTEGER_H

// MP_Integer : very simple multiprecision integer.
// It is not optimized for speed, but for simplicity.

#include <CGAL/basic.h>
#include <CGAL/Interval_arithmetic.h>
#include <iostream>
#include <vector>
#include <functional>

// TODO :
// - Have a typedef Exact_Integer defined to Lazy<MP_Integer/leda_integer/Gmpz>
// - Idem for rational.
// - add "explicit" for the ctors ?
// - the non-inline stuff must go in src/MP_Integer.C.
// - Use iterators instead of operator[], it should be faster and cleaner.
// - implement missing CGAL requirements.
// - Documentation.
// - Implement Karatsuba in order to have less than quadratic multiplication
//   for large operands ?  Maybe not worth it, but maybe worth a bench.

CGAL_BEGIN_NAMESPACE

class MP_Integer
{
public:

  typedef short limb;
  typedef int   limb2;

  static const unsigned log_limb = 8*sizeof(limb);
  static const limb2    base     = 1<<log_limb;

  typedef std::vector<limb>  V;
  typedef V::const_iterator  const_iterator;
  typedef V::iterator        iterator;

  MP_Integer()
  {
    CGAL_assertion(sizeof(limb2) == 2*sizeof(limb));
  }

  MP_Integer(limb i)
  : v(1,i)
  {
    remove_leading_zeros();
  }

  MP_Integer(limb2 i)
  : v(2)
  {
    v[0] = i;
    v[1] = higher_limb(i);
    remove_leading_zeros();
  }

  MP_Integer operator-() const
  {
    return MP_Integer() - (*this);
  }
  MP_Integer operator+(const MP_Integer &) const;
  MP_Integer operator-(const MP_Integer &) const;
  MP_Integer operator*(const MP_Integer &) const;
  MP_Integer operator/(const MP_Integer &) const;

  MP_Integer& operator+=(const MP_Integer &a)
  {
    *this = *this + a;
    return *this;
  }

  MP_Integer& operator-=(const MP_Integer &a)
  {
    *this = *this - a;
    return *this;
  }

  MP_Integer& operator*=(const MP_Integer &a)
  {
    *this = *this * a;
    return *this;
  }

  MP_Integer& operator/=(const MP_Integer &a)
  {
    *this = *this / a;
    return *this;
  }

  bool operator<(const MP_Integer &) const;

  
  bool operator>(const MP_Integer &f) const
  {
      return f<*this;
  }

  bool operator>=(const MP_Integer &f) const
  {
      return !(*this<f);
  }

  bool operator<=(const MP_Integer &f) const
  {
      return !(*this>f);
  }

  bool operator==(const MP_Integer &b) const
  {
    return v == b.v;
  }

  bool operator!=(const MP_Integer &b) const
  {
    return v != b.v;
  }

  void remove_leading_zeros()
  {
    while (!v.empty() && v.back()==0)
      v.pop_back();
  }

  // Accessory function : Given a limb2, this returns the higher limb.
  // The lower limb is simply obtained by casting to a limb.
  static
  limb higher_limb(limb2 l)
  {
    return (l - (limb) l) >> log_limb;
  }

  V v;
};


namespace NTS {

Comparison_result
compare (const MP_Integer & a, const MP_Integer & b)
{
  if (a.v.size() > b.v.size())
    return (a.v.back() > 0) ? LARGER : SMALLER;
  if (a.v.size() < b.v.size())
    return (b.v.back() < 0) ? LARGER : SMALLER;

  for (int i=a.v.size()-1; i>=0; i--)
  {
      if (a.v[i] > b.v[i])
	  return LARGER;
      if (a.v[i] < b.v[i])
	  return SMALLER;
  }
  return EQUAL;
}

} // namespace NTS

bool
MP_Integer::operator<(const MP_Integer & b) const
{
  return CGAL_NTS compare(*this, b) == SMALLER;
}

// Common code for operator+ and operator-.
template <class BinOp>
void
Add_Sub(MP_Integer &r, const MP_Integer &a, const MP_Integer &b)
{
  unsigned nsize = std::max(a.v.size(), b.v.size());
  r.v.resize(nsize);
  MP_Integer::limb2 carry = 0;
  for(unsigned i=0; i<nsize; i++)
  {
    MP_Integer::limb2 tmp = carry + BinOp()(i<a.v.size() ? a.v[i] : 0,
                                            i<b.v.size() ? b.v[i] : 0);
    r.v[i] = tmp;
    carry = MP_Integer::higher_limb(tmp);
  }
  if (carry != 0)
    r.v.push_back(carry);
  else
    r.remove_leading_zeros();
}

MP_Integer
MP_Integer::operator+(const MP_Integer &b) const
{
  MP_Integer r;
  Add_Sub<std::plus<limb2> >(r, *this, b);
  return r;
}

MP_Integer
MP_Integer::operator-(const MP_Integer &b) const
{
  MP_Integer r;
  Add_Sub<std::minus<limb2> >(r, *this, b);
  return r;
}

MP_Integer
MP_Integer::operator*(const MP_Integer &b) const
{
  MP_Integer r;
  r.v.assign(v.size() + b.v.size(),0);
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
  r.remove_leading_zeros();
  return r;
}

#if 1 // Should I define the following or not ?
MP_Integer
MP_Integer::operator/(const MP_Integer &) const
{
  CGAL_assertion_msg(false, "Sorry, MP_Integer division not implemented");
  abort();
  return MP_Integer();
}

MP_Integer
sqrt(const MP_Integer &)
{
  CGAL_assertion_msg(false, "Sorry, MP_Integer square root not implemented");
  abort();
  return MP_Integer();
}
#endif

double
to_double(const MP_Integer &b)
{
  if (b == 0)
    return 0;
  // Linear complexity.  Could be made constant time.
  double d=0;
  for (MP_Integer::const_iterator i = b.v.end()-1; i >= b.v.begin(); i--)
    d = d*MP_Integer::base + *i;
  return d;
}

Interval_base
to_interval(const MP_Integer &b)
{
  if (b == 0)
    return 0;
  // Linear complexity.  Could be made constant time.
  Protect_FPU_rounding<true> P;
  Interval_nt_advanced d=0;
  for (MP_Integer::const_iterator i = b.v.end()-1; i >= b.v.begin(); i--)
    d = d*MP_Integer::base + *i;
  return d;
}

bool is_finite(const MP_Integer &)
{
  return true;
}

bool is_valid(const MP_Integer &)
{
  return true;
}

inline
io_Operator
io_tag(const MP_Integer&)
{ return io_Operator(); }

inline
Number_tag
number_type_tag(const MP_Integer& )
{ return Number_tag(); }

std::ostream &
operator<< (std::ostream & os, const MP_Integer &b)
{
  // Binary format would be nice and not hard to have too.
  if (b == 0)
    return os << 0;

  MP_Integer::const_iterator i = b.v.begin();
  unsigned exp = 0;

  os << *i;

  for (i++; i != b.v.end(); i++)
  {
    exp += MP_Integer::log_limb;

    os << ((*i > 0) ? " +" : " ") << *i << " * 2^" << exp;
  }

  return os;
}

std::istream &
operator>> (std::istream & is, MP_Integer &b)
{
  int i;
  is >> i;
  b = MP_Integer(i);
  return is;
}

CGAL_END_NAMESPACE

#endif // CGAL_MP_INTEGER_H
