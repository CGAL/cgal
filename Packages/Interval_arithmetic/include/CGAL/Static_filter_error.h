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
// file          : include/CGAL/Static_filter_error.h
// revision      : $Revision$
// revision_date : $Date$
// package       : Interval Arithmetic
// author(s)     : Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_STATIC_FILTER_ERROR_H
#define CGAL_STATIC_FILTER_ERROR_H

// This file contains the description of the class Static_filter_error.
// The goal of this class is to be run by some overloaded predicates,
// to compute error bound done in these functions.
// 
// The original idea is from Olivier Devillers.

// TODO: I need to add some missing operators and functions, min/max...

#include <CGAL/basic.h>
#include <CGAL/Interval_arithmetic/_FPU.h>

CGAL_BEGIN_NAMESPACE

struct Static_filter_error
{
  typedef Static_filter_error Sfe;

  Static_filter_error (const double &b, const double &e, const int &d)
      : _b(b), _e(e), _d(d) {}

  Sfe operator+ (const Sfe &f) const
  {
      CGAL_warning_msg(_d == f._d,
	      "you are adding variables of different homogeneous degree");
      double b = _b + f._b;
      FPU_CW_t backup = FPU_get_and_set_cw(FPU_cw_up);
      double e = (ulp(b)/2 + _e) + f._e;
      FPU_set_cw(backup);
      return Sfe(b, e, _d);
  }

  Sfe operator* (const Sfe &f) const
  {
      double b = _b * f._b;
      FPU_CW_t backup = FPU_get_and_set_cw(FPU_cw_up);
      double e = (ulp(b)/2 + _e * f._e) +  _e * f._b + _b * f._e;
      FPU_set_cw(backup);
      return Sfe(b, e, _d+f._d);
  }

  Sfe operator- (const Sfe &f) const { return *this + f; }
  Sfe operator- ()             const { return *this; }
  Sfe operator/ (const Sfe &f) const { abort(); } // Division not supported.

  double error()  const { return _e; }
  double bound()  const { return _b; }
  int    degree() const { return _d; }

  double ulp const (const double &d)
  {
      // You are supposed to call this function with rounding towards
      // +infinity, and on a positive number.
      CGAL_assertion(d>=0);
      double u = (d + CGAL_IA_MIN_DOUBLE) - d;
      CGAL_assertion(u!=0);
      return u;
  }

private:
  // _b is a bound on the absolute value of the _double_ value of the
  //    variable.
  // _e is a bound on the absolute error (difference between _b and the
  //    _real_ value of the variable.
  // _d is the degree of the variable, it allows some additionnal checks.
  double _b, _e;
  int _d;
};

inline
Static_filter_error
sqrt(const Static_filter_error &f)
{
  CGAL_warning_msg(f.degree() & 1 == 0,
	  "you really want a non integer degree ???");
  double b = sqrt(f.bound());
  FPU_CW_t backup = FPU_get_and_set_cw(FPU_cw_up);
  double e = sqrt(f.error()) + f.ulp(b)/2;
  FPU_set_cw(backup);
  return Static_filter_error(b, e, f.degree()/2);
}


// Now the overloaded built-in predicates, that will _set_ the epsilons.

inline
Comparison_result
compare(const Static_filter_error &a,
       	const Static_filter_error &b,
       	double &epsilon)
{
  Static_filter_error c = a-b;
  epsilon = c.error();
  return EQUAL; // ??
}

inline
Sign
sign(const Static_filter_error &a, double &epsilon)
{
  epsilon = a.error();
  return ZERO; // ??
}

CGAL_END_NAMESPACE

#endif // CGAL_STATIC_FILTER_ERROR_H
