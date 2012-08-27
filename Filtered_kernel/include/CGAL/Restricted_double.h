// Copyright (c) 1999,2000  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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

#ifndef CGAL_RESTRICTED_DOUBLE_H
#define CGAL_RESTRICTED_DOUBLE_H

// This file contains the description of the class Restricted_double.
// The goal of this class is to be run by some overloaded predicates,
// This is just a wrapper for a double, but with limited functionality,
// to detect non-intented use.
// It's the type over which is run the epsilon variant, and which throws
// the exceptions.

// TODO: I need to add some missing operators and functions, min/max...

#include <CGAL/basic.h>
// #include <CGAL/double.h>

namespace CGAL {

#ifdef CGAL_IA_CHECK_RESTRICT
struct Restricted_double
{
  typedef Restricted_double Self;
  struct unsafe_comparison {};          // Exception class.

  Restricted_double () {}
  Restricted_double (const double &d) : _d(d) {}
  Restricted_double (const int &i) : _d(double(i)) {}
// Add operator= for efficiency.
  Self operator+ (const Self &f) const { return _d+f._d; }
  Self operator- (const Self &f) const { return _d-f._d; }
  Self operator* (const Self &f) const { return _d*f._d; }
  Self operator/ (const Self &f) const { return _d/f._d; }
  Self operator- ()              const { return -_d; }

  Self& operator+=(const Self &f) { return *this = _d + f._d; }
  Self& operator-=(const Self &f) { return *this = _d - f._d; }
  Self& operator*=(const Self &f) { return *this = _d * f._d; }
  Self& operator/=(const Self &f) { return *this = _d / f._d; }

#if 0
  bool operator< (const Self &f) const { return _d < f._d; }
#endif
  double dbl () const { return _d; }

private:
  double _d;
};

inline
Restricted_double
sqrt(const Restricted_double &f)
{ return std::sqrt(f.dbl()); }

inline
Restricted_double
abs(const Restricted_double &f)
{ return CGAL::abs(f.dbl()); }

inline
double
to_double(const Restricted_double &f)
{ return f.dbl(); }

#else
typedef double Restricted_double;
#endif

} //namespace CGAL

#endif // CGAL_RESTRICTED_DOUBLE_H
