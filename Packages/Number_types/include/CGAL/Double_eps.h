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
// file          : Double_eps.h
// package       : Number_types
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_DOUBLE_EPS_H
#define CGAL_DOUBLE_EPS_H

#ifndef CGAL_DOUBLE_H
#include <CGAL/double.h>
#endif // CGAL_DOUBLE_H
#ifndef CGAL_NUMBER_UTILS_H
#include <CGAL/number_utils.h>
#endif // CGAL_NUMBER_UTILS_H

CGAL_BEGIN_NAMESPACE


class Double_eps
{
friend double set_eps(double eps);
friend std::istream& operator>>(std::istream&, Double_eps& );
friend std::ostream& operator<<(std::ostream&, const Double_eps& );

public:
  Double_eps()
    : _d(0.0)
  {}

  Double_eps(double d)
    : _d(d)
  {}

  double d() const
  {
    return _d;
  }

  double eps() const
  {
    return _eps;
  }

  operator double() const
  {
    return _d;
  }

  Double_eps operator+(const Double_eps &de) const
  {
    return Double_eps(_d + de._d);
  }

  Double_eps operator-(const Double_eps &de) const
  {
    return Double_eps(_d - de._d);
  }

  Double_eps operator-() const
  {
    return Double_eps(-_d);
  }

  Double_eps operator*(const Double_eps &de) const
  {
    return Double_eps(_d * de._d);
  }

  Double_eps operator/(const Double_eps &de) const
  {
    return Double_eps(_d / de._d);
  }

  Double_eps&  operator+=(const Double_eps &de)
  {
    _d += de._d;
    return *this;
  }

  Double_eps&  operator-=(const Double_eps &de)
  {
    _d -= de._d;
    return *this;
  }

  Double_eps&  operator*=(const Double_eps &de)
  {
    _d *= de._d;
    return *this;
  }

  Double_eps&  operator/=(const Double_eps &de)
  {
    _d /= de._d;
    return *this;
  }

  bool operator==(const Double_eps &de) const
  {
    return CGAL_NTS abs(_d - de._d) <= _eps;
  }

  bool operator!=(const Double_eps &de) const
  {
    return !(*this == de);
  }

  bool operator<(const Double_eps &de) const
  {
    if (*this == de) {
      return false;
    }
    return _d < de._d;
  }

  bool operator>(const Double_eps &de) const
  {
    if (*this == de) {
      return false;
    }
    return _d > de._d;
  }

  bool operator<=(const Double_eps &de) const
  {
    if (*this == de) {
      return true;
    }
    return _d <= de._d;
  }

  bool operator>=(const Double_eps &de) const
  {
    if (*this == de) {
      return true;
    }
    return _d >= de._d;
  }

private:
  double _d;
  static double _eps;
};

inline bool is_valid(const Double_eps &de)
  {
    return is_valid(de.d());
  }

inline bool is_finite(const Double_eps &de)
  {
    return is_finite(de.d());
  }

inline double to_double(const Double_eps &de)
  {
    return de.d();
  }

inline Double_eps numerator(const Double_eps &d)
{
  return d;
}

inline Double_eps denominator(const Double_eps &)
{
  return Double_eps(1.0);
}

inline Number_tag number_type_tag(const Double_eps &)
{
  return Number_tag();
}

inline io_Operator io_tag(const Double_eps &)
{
  return io_Operator();
}

extern double set_eps(double eps);

extern std::ostream  &operator<<(std::ostream& os, const Double_eps &);

CGAL_END_NAMESPACE

#endif  // CGAL_DOUBLE_EPS_H

