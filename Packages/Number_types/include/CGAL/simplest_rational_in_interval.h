// ======================================================================
//
// Copyright (c) 2002 The CGAL Consortium
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
// file          : include/CGAL/simplest_rational_in_interval.h
// package       : Number_types
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri, Susan Hert, Sylvain Pion
//
// coordinator   : MPI, Saarbruecken
// ======================================================================
 
#ifndef CGAL_SIMPLEST_RATIONAL_IN_INTERVAL_H
#define CGAL_SIMPLEST_RATIONAL_IN_INTERVAL_H


#include <CGAL/basic.h>
#include <CGAL/to_rational.h>
#include <cassert>
#include <climits>
#include <cmath>


CGAL_BEGIN_NAMESPACE

/* simplest_rational_in_interval(x,y) returns the rational number with
     the smallest denominator in the interval [x,y].  See Knuth,
     "Seminumerical algorithms", page 654, answer to exercise
     4.53-39. */

template <class Rational>
Rational 
simplest_rational_in_interval(double x, double y) {

  if(x == y){
    return to_rational<Rational>(x);
  }

  if(x > y){
    std::swap(x,y);
  }

  Rational r;  // Return value. 
  typename Rational_traits<Rational>::RT r_numerator, r_denominator;
  Rational_traits<Rational> t;
  // Deal with negative arguments.  We only have to deal with the case
  // where both x and y are negative -- when exactly one is negative
  // the best rational in the interval [x,y] is 0.
  if (x < 0 && y < 0) {
    // Both arguments are negative: solve positive case and negate
    return  - simplest_rational_in_interval<Rational>(fabs(x),fabs(y));
  } else if (x <= 0 || y <= 0) {
    // One argument is 0, or arguments are on opposite sides of 0:
    // simplest rational in interval is 0 exactly. 
    r_numerator = 0;
    r_denominator = 1;
  } else { // x > 0 && y > 0
    double xc = CGAL_CLIB_STD::floor(1/x); // First coefficient of cf for x.
    double xr = CGAL_CLIB_STD::fmod(1/x,1); // Remaining fractional part of x.
    double yc = CGAL_CLIB_STD::floor(1/y); // First coefficient of cf for y.
    double yr = CGAL_CLIB_STD::fmod(1/y,1); // Remaining fractional part of y.

    if (xc < yc) {
      // Return 1/(xc+1).
      r_numerator = 1;
      r_denominator = xc + 1;
    } else if (yc < xc) {
      // Return 1/(yc+1).
      r_numerator = 1;
      r_denominator = yc + 1;
    } else  {  // xc == yc
      // Recurse to find s, the rational with the lowest denominator
      //  between xr and yr.
      Rational  s(simplest_rational_in_interval<Rational>(xr,yr));

      // Return 1/(xc + s).

      r_numerator = t.denominator(s);
      typename Rational_traits<Rational>::RT  xc_rt(xc);
      r_denominator = t.numerator(s) + xc_rt * t.denominator(s);
    }
  }

  return t.make_rational(r_numerator, r_denominator);
}

CGAL_END_NAMESPACE

#endif //CGAL_SIMPLEST_RATIONAL_IN_INTERVAL_H
