// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// 
// release       :
// release_date  :
// 
// file          : include/CGAL/Double_IEEE.h
// revision      : $Revision$
// revision_date : $Date$
// package       : Interval Arithmetic
// author(s)     : Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ======================================================================
 
#ifndef CGAL_DOUBLE_IEEE_H
#define CGAL_DOUBLE_IEEE_H

#include <CGAL/basic.h>
#include <CGAL/double.h>
#include <CGAL/number_utils.h>

CGAL_BEGIN_NAMESPACE

class Double_IEEE
{
  double _d;

  typedef Double_IEEE Self;

  friend std::istream& operator>>(std::istream&, Self& );
  friend std::ostream& operator<<(std::ostream&, const Self);

public:
  Double_IEEE() {}

  Double_IEEE(const double d) : _d(d) {}

  double d() const
  {
    return _d;
  }
#if 1
    // The copy constructors/assignment: useless.
    // The default ones are ok, but these are faster...
    // [check that]
  Double_IEEE(const Self d)
    : _d(d._d) {}
    
  Self & operator=(const Self & d)
  {
      _d = d._d;
      return *this;
  }
#endif

  operator double() const
  {
    return _d;
  }

  Self operator-() const
  {
    return Self(-_d);
  }

  Self operator+(const Self d) const
  {
    return Self(_d + d._d);
  }

  Self operator-(const Self d) const
  {
    return Self(_d - d._d);
  }

  Self operator*(const Self d) const
  {
    return Self(_d * d._d);
  }

  Self operator/(const Self d) const
  {
    return Self(_d / d._d);
  }

  Self&  operator+=(const Self d)
  {
    _d += d._d;
    return *this;
  }

  Self&  operator-=(const Self d)
  {
    _d -= d._d;
    return *this;
  }

  Self&  operator*=(const Self d)
  {
    _d *= d._d;
    return *this;
  }

  Self&  operator/=(const Self d)
  {
    _d /= d._d;
    return *this;
  }

  bool operator==(const Self d) const
  {
    return _d = d._d;
  }

  bool operator!=(const Self d) const
  {
    return !(*this == d);
  }

  bool operator<(const Self d) const
  {
    return _d < d._d;
  }

  bool operator>(const Self d) const
  {
    return _d > d._d;
  }

  bool operator<=(const Self d) const
  {
    return _d <= d._d;
  }

  bool operator>=(const Self d) const
  {
    return _d >= d._d;
  }
};

inline
bool
is_valid(const Double_IEEE &d)
{
    return is_valid(d.d());
}

inline
bool
is_finite(const Double_IEEE &d)
{
    return is_finite(d.d());
}

inline
double
to_double(const Double_IEEE &d)
{
    return d.d();
}

inline
Number_tag
number_type_tag(const Double_IEEE &)
{
  return Number_tag();
}

inline
io_Operator
io_tag(const Double_IEEE &)
{
  return io_Operator();
}

extern std::ostream  &operator<<(std::ostream& os, const Self);

CGAL_END_NAMESPACE

#endif  // CGAL_DOUBLE_IEEE_H

