// ============================================================================
//
// Copyright (c) 1999 The CGAL Consortium
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
// file          : include/CGAL/Restricted_double.h
// revision      : $Revision$
// revision_date : $Date$
// package       : Interval Arithmetic
// author(s)     : Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_RESTRICTED_DOUBLE_H
#define CGAL_RESTRICTED_DOUBLE_H

// This file contains the description of the class Restricted_double.
// The goal of this class is to be run by some overloaded predicates,
// This is just a wrapper for a double, but with limited functionality,
// to detect non-intented use.
// It's the type over which is run the _SAF epsilon variant, and which throws
// the exceptions.

// TODO: I need to add some missing operators and functions, min/max...

#include <CGAL/basic.h>
// #include <CGAL/double.h>

CGAL_BEGIN_NAMESPACE

struct Restricted_double
{
  typedef Restricted_double Self;
  struct unsafe_comparison {};          // Exception class.

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
{
    return std::sqrt(f.dbl());
}

inline
Restricted_double
abs(const Restricted_double &f)
{
    return std::fabs(f.dbl());
}

// Now the epsilon predicates, which might throw the exception.

inline
Sign
lexicographical_sign_SAF(const Restricted_double &,
       	    const Restricted_double &,
       	    const double &)
{
    // Not finished...
    throw Restricted_double::unsafe_comparison();
}

inline
Comparison_result
compare_SAF(const Restricted_double &a,
       	    const Restricted_double &b,
       	    const double &epsilon)
{
    if (a.dbl() > b.dbl()+epsilon) return LARGER;
    if (a.dbl() < b.dbl()-epsilon) return SMALLER;
    if (a.dbl()==b.dbl() && epsilon==0) return EQUAL;
    throw Restricted_double::unsafe_comparison();
}

inline
Sign
sign_SAF(const Restricted_double &a, const double &epsilon)
{
    // return compare_SAF(a,0,epsilon);
    if (a.dbl()> epsilon) return POSITIVE;
    if (a.dbl()<-epsilon) return NEGATIVE;
    if (a.dbl()==0 && epsilon==0) return ZERO;
    throw Restricted_double::unsafe_comparison();
}

CGAL_END_NAMESPACE

#endif // CGAL_RESTRICTED_DOUBLE_H
