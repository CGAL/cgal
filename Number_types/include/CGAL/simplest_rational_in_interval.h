// Copyright (c) 2002  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
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
// Author(s)     : Andreas Fabri, Susan Hert, Sylvain Pion

#ifndef CGAL_SIMPLEST_RATIONAL_IN_INTERVAL_H
#define CGAL_SIMPLEST_RATIONAL_IN_INTERVAL_H

#include <CGAL/number_type_basic.h>
#include <CGAL/to_rational.h>
#include <CGAL/use.h>
#include <climits>
#include <cmath>

namespace CGAL {

/* simplest_rational_in_interval(x,y) returns the rational number with
     the smallest denominator in the interval [x,y].  See Knuth,
     "Seminumerical algorithms", page 654, answer to exercise
     4.53-39. */

template <class Rational>
Rational
simplest_rational_in_interval(double x, double y) {

    typedef Fraction_traits<Rational> FT;
    typedef typename FT::Is_fraction Is_fraction;
    typedef typename FT::Numerator_type Numerator_type;
    typedef typename FT::Denominator_type Denominator_type;
    typedef typename FT::Decompose Decompose;
    typedef typename FT::Compose Compose;

    // Must be a fraction
    CGAL_USE_TYPE(Is_fraction);
    CGAL_static_assertion((::boost::is_same<Is_fraction, Tag_true>::value));
    // Numerator_type,Denominator_type must be the same
    CGAL_USE_TYPE(Denominator_type);
    CGAL_static_assertion((::boost::is_same<Numerator_type, Denominator_type>::value));


  if(x == y){
    return to_rational<Rational>(x);
  }

  if(x > y){
    std::swap(x,y);
  }

  Rational r;  // Return value.
  Numerator_type r_numerator, r_denominator;
  // Deal with negative arguments.  We only have to deal with the case
  // where both x and y are negative -- when exactly one is negative
  // the best rational in the interval [x,y] is 0.
  if (x < 0 && y < 0) {
    // Both arguments are negative: solve positive case and negate
    return  - simplest_rational_in_interval<Rational>(CGAL::abs(x),CGAL::abs(y));
  } else if (x <= 0 || y <= 0) {
    // One argument is 0, or arguments are on opposite sides of 0:
    // simplest rational in interval is 0 exactly.
    r_numerator = 0;
    r_denominator = 1;
  } else { // x > 0 && y > 0
    double xc = std::floor(1/x); // First coefficient of cf for x.
    double xr = std::fmod(1/x,1); // Remaining fractional part of x.
    double yc = std::floor(1/y); // First coefficient of cf for y.
    double yr = std::fmod(1/y,1); // Remaining fractional part of y.

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

      Numerator_type xc_rt(xc);
      Numerator_type s_num,s_den;
      Decompose()(s,s_num,s_den);
      r_numerator = s_den;
      r_denominator = s_num + xc_rt * s_den;
    }
  }

  return Compose()(r_numerator, r_denominator);
}

} //namespace CGAL

#endif //CGAL_SIMPLEST_RATIONAL_IN_INTERVAL_H
