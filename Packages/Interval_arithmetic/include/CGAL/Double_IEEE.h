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
friend std::istream& operator>>(std::istream&, Double_IEEE& );
friend std::ostream& operator<<(std::ostream&, const Double_IEEE);

public:
  Double_IEEE()
//    : _d(0.0)
  {}

  Double_IEEE(double d)
    : _d(d)
  {}

  double d() const
  {
    return _d;
  }

  operator double() const
  {
    return _d;
  }

  Double_IEEE operator+(const Double_IEEE de) const
  {
    return Double_IEEE(_d + de._d);
  }

  Double_IEEE operator-(const Double_IEEE de) const
  {
    return Double_IEEE(_d - de._d);
  }

  Double_IEEE operator-() const
  {
    return Double_IEEE(-_d);
  }

  Double_IEEE operator*(const Double_IEEE de) const
  {
    return Double_IEEE(_d * de._d);
  }

  Double_IEEE operator/(const Double_IEEE de) const
  {
    return Double_IEEE(_d / de._d);
  }

  Double_IEEE&  operator+=(const Double_IEEE de)
  {
    _d += de._d;
    return *this;
  }

  Double_IEEE&  operator-=(const Double_IEEE de)
  {
    _d -= de._d;
    return *this;
  }

  Double_IEEE&  operator*=(const Double_IEEE de)
  {
    _d *= de._d;
    return *this;
  }

  Double_IEEE&  operator/=(const Double_IEEE de)
  {
    _d /= de._d;
    return *this;
  }

  bool operator==(const Double_IEEE de) const
  {
    return _d = de._d;
  }

  bool operator!=(const Double_IEEE de) const
  {
    return !(*this == de);
  }

  bool operator<(const Double_IEEE de) const
  {
    return _d < de._d;
  }

  bool operator>(const Double_IEEE de) const
  {
    return _d > de._d;
  }

  bool operator<=(const Double_IEEE de) const
  {
    return _d <= de._d;
  }

  bool operator>=(const Double_IEEE de) const
  {
    return _d >= de._d;
  }

private:
  double _d;
};

inline
bool
is_valid(const Double_IEEE &de)
{
    return is_valid(de.d());
}

inline
bool
is_finite(const Double_IEEE &de)
{
    return is_finite(de.d());
}

inline
double
to_double(const Double_IEEE &de)
{
    return de.d();
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

extern std::ostream  &operator<<(std::ostream& os, const Double_IEEE &);

CGAL_END_NAMESPACE

#endif  // CGAL_DOUBLE_IEEE_H

