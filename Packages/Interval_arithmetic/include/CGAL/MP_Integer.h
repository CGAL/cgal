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
#include <iostream>
#include <vector>
#include <utility>

// TODO :
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
  friend std::ostream & operator<< (std::ostream &, const MP_Integer &);

  // limb2 needs to be twice the sizeof limb.
  typedef short limb;
  typedef int   limb2;

public:
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

  bool operator<(const MP_Integer &) const;
  bool operator==(const MP_Integer &b) const
  {
    return b.v == v;
  }

private:

  void remove_leading_zeros()
  {
    while (!v.empty() && v.back()==0)
      v.pop_back();
  }

  unsigned int size() const
  {
    return v.size();
  }

  limb operator[](unsigned i) const
  {
    return i<size() ? v[i] : 0;
  }

  // Accessory function : Given a limb2, this returns the higher limb.
  // The lower limb is simply obtained by casting to a limb.
  static
  limb higher_limb(limb2 l)
  {
    return (l - (limb) l) >> (8*sizeof(limb));
  }

  std::vector<limb> v;
};

// The following must go in src/MP_Integer.C

bool
MP_Integer::operator<(const MP_Integer & b) const
{
  for (int i=std::max(size(), b.size()); i>=0; i--)
  {
      if ((*this)[i] < b[i])
	  return true;
      if ((*this)[i] > b[i])
	  return false;
  }
  return false;
}

// Maybe we can merge the 2 following functions using STL's Add/Sub function
// objects and an accessory function ?  Maybe even the function above ?
MP_Integer
MP_Integer::operator+(const MP_Integer &b) const
{
  unsigned nsize = std::max(size(), b.size());
  MP_Integer r;
  r.v.resize(nsize);
  limb2 carry = 0;
  for(unsigned i=0; i<nsize; i++)
  {
    limb2 tmp = carry + (limb2) (*this)[i] + (limb2) b[i];
    r.v[i] = tmp;
    carry = higher_limb(tmp);
  }
  if (carry != 0)
    r.v.push_back(carry);
  return r;
}

MP_Integer
MP_Integer::operator-(const MP_Integer &b) const
{
  unsigned nsize = std::max(size(), b.size());
  MP_Integer r;
  r.v.resize(nsize);
  limb2 carry = 0;
  for(unsigned i=0; i<nsize; i++)
  {
    limb2 tmp = carry + (limb2) (*this)[i] - (limb2) b[i];
    r.v[i] = tmp;
    carry = higher_limb(tmp);
  }
  if (carry != 0)
    r.v.push_back(carry);
  return r;
}

MP_Integer
MP_Integer::operator*(const MP_Integer &b) const
{
  MP_Integer r;
  r.v.assign(size() + b.size(),0);
  for(unsigned i=0; i<size(); i++)
  {
    unsigned j;
    limb2 carry = 0;
    for(j=0; j<b.size(); j++)
    {
      limb2 tmp = carry + (limb2) r.v[i+j] + (limb2) v[i] * (limb2) b.v[j];
      r.v[i+j] = tmp;
      carry = higher_limb(tmp);
    }
    r.v[i+j] = carry;
  }
  r.remove_leading_zeros();
  return r;
}

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
  CGAL_assertion(false, "Sorry, MP_Integer square root not implemented");
  abort();
  return MP_Integer();
}

std::ostream &
operator<< (std::ostream & os, const MP_Integer &b)
{
  if (b.size() == 0)
    return os << 0;

  for (int i=b.size()-1; i>=0; i--)
  {
    os << "(" << b[i] << ")"; // parenthesis are here for negatives...
    if (i!=0)
      os << "* 2^" << (8*sizeof(MP_Integer::limb))*i << " + ";
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
