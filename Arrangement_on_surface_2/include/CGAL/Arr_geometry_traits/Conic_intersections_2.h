// Copyright (c) 2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s)     : Ron Wein <wein@post.tau.ac.il>

#ifndef CGAL_CONIC_INTERSECTIONS_2_H
#define CGAL_CONIC_INTERSECTIONS_2_H

/*! \file
 * Implementation of functions related to the intersection of conics curves.
 */

namespace CGAL {

/*!
 * Compute the roots of the resultants of the two bivariate polynomials:
 *   C1:  r1*x^2 + s1*y^2 + t1*xy + u1*x + v1*y + w1 = 0
 *   C2:  r2*x^2 + s2*y^2 + t2*xy + u2*x + v2*y + w2 = 0
 * \param deg1 The degree of the first curve.
 * \param deg2 The degree of the second curve.
 * \param xs Output: The real-valued roots of the polynomial, sorted in an
 *                   ascending order.
 * \pre xs must be a vector of size 4.
 * \return The number of distinct roots found.
 */
template <class Nt_traits>
int
 _compute_resultant_roots (Nt_traits& nt_traits,
                           const typename Nt_traits::Integer& r1,
                           const typename Nt_traits::Integer& s1,
                           const typename Nt_traits::Integer& t1,
                           const typename Nt_traits::Integer& u1,
                           const typename Nt_traits::Integer& v1,
                           const typename Nt_traits::Integer& w1,
                           const int& deg1,
                           const typename Nt_traits::Integer& r2,
                           const typename Nt_traits::Integer& s2,
                           const typename Nt_traits::Integer& t2,
                           const typename Nt_traits::Integer& u2,
                           const typename Nt_traits::Integer& v2,
                           const typename Nt_traits::Integer& w2,
                           const int& deg2,
                           typename Nt_traits::Algebraic *xs)
{ 
  if (deg1 == 2 && deg2 == 1)
  {
    // If necessary, swap roles between the two curves, so that the first
    // curve always has the minimal degree.
    return (_compute_resultant_roots (nt_traits,
				      r2, s2, t2, u2, v2, w2, 
				      deg2,
				      r1, s1, t1, u1, v1, w1, 
				      deg1,
				      xs));
  }
  
  // Act according to the degree of the first conic curve.
  const typename Nt_traits::Integer  _two = 2;
  typename Nt_traits::Integer        c[5];
  unsigned int                       degree = 4;
  typename Nt_traits::Algebraic     *xs_end;

  if (deg1 == 1)
  {
    // The first curve has no quadratic coefficients, and represents a line.
    if (CGAL::sign (v1) == ZERO)
    {
      // The first line is u1*x + w1 = 0, therefore:
      xs[0] = nt_traits.convert(-w1) / nt_traits.convert(u1);
      return (1);
    }
    
    // We can write the first curve as: y = -(u1*x + w1) / v1.
    if (deg2 == 1)
    {
      // The second curve is also a line. We therefore get the linear
      // equation c[1]*x + c[0] = 0:
      c[1] = v1*u2 - u1*v2;
      c[0] = v1*w2 - w1*v2;

      if (CGAL::sign (c[1]) == ZERO)
	// The two lines are parallel:
	return (0);

      xs[0] =  nt_traits.convert(-c[0]) /  nt_traits.convert(c[1]);
      return (1);
    }

    // We substitute this expression into the equation of the second
    // conic, and get the quadratic equation c[2]*x^2 + c[1]*x + c[0] = 0:
    c[2] = u1*u1*s2 - u1*v1*t2 + v1*v1*r2;
    c[1] = _two*u1*w1*s2 - u1*v1*v2 - v1*w1*t2 + v1*v1*u2;
    c[0] = w1*w1*s2 - v1*w1*v2 + v1*v1*w2;
    
    xs_end = nt_traits.solve_quadratic_equation (c[2], c[1], c[0],
						 xs);
    return static_cast<int>(xs_end - xs);
  }

  // At this stage, both curves have degree 2. We obtain a qaurtic polynomial
  // whose roots are the x-coordinates of the intersection points.
  if (CGAL::sign (s1) == ZERO && CGAL::sign (s2) == ZERO)
  {
    // If both s1 and s2 are zero, we can write the two curves as:
    //   C1: (t1*x + v1)*y + (r1*x^2 + u1*x + w1) = 0
    //   C2: (t2*x + v2)*y + (r2*x^2 + u2*x + w2) = 0
    // By writing the resultant of these two polynomials we get a cubic
    // polynomial, whose coefficients are given by:
    c[3] = r2*t1 - r1*t2;
    c[2] = t1*u2 - t2*u1 + r2*v1 - r1*v2;
    c[1] = t1*w2 - t2*w1 + u2*v1 - u1*v2;
    c[0] = v1*w2 - v2*w1;
    
    degree = 3;
  }
  else
  {
    // We can write the two curves as:
    //   C1: (s1)*y^2 + (t1*x + v1)*y + (r1*x^2 + u1*x + w1) = 0
    //   C2: (s2)*y^2 + (t2*x + v2)*y + (r2*x^2 + u2*x + w2) = 0
    // By writing the resultant of these two polynomials we get a quartic
    // polynomial, whose coefficients are given by:
    c[4] = -_two*s1*s2*r1*r2 + s1*t2*t2*r1 - s1*t2*t1*r2 +
      s1*s1*r2*r2 - s2*t1*r1*t2 + s2*t1*t1*r2 + s2*s2*r1*r1;

    c[3] = -t2*r1*v1*s2 - u2*t1*t2*s1 - v2*r1*t1*s2 -
      r2*t1*v2*s1 - _two*s1*s2*r1*u2 - t2*u1*t1*s2 + u2*t1*t1*s2 -
      r2*v1*t2*s1 + u1*t2*t2*s1 + _two*v2*r1*t2*s1 + _two*u2*r2*s1*s1 + 
      _two*r2*v1*t1*s2 + _two*u1*r1*s2*s2 - _two*s1*s2*u1*r2;
    
    c[2] = -r2*v1*v2*s1 + u2*u2*s1*s1 + _two*w2*r2*s1*s1 +
      _two*u2*v1*t1*s2 - u2*v1*t2*s1 + w2*t1*t1*s2 - _two*s1*s2*u1*u2 - 
      w2*t1*t2*s1 + v2*v2*r1*s1 + u1*u1*s2*s2 - v2*r1*v1*s2 +
      _two*w1*r1*s2*s2 - u2*t1*v2*s1 - t2*u1*v1*s2 - _two*s1*s2*r1*w2 -
      _two*s1*s2*w1*r2 + r2*v1*v1*s2 + w1*t2*t2*s1 - v2*u1*t1*s2 -
      t2*w1*t1*s2 + _two*v2*u1*t2*s1;

    c[1] = _two*w2*u2*s1*s1 + _two*w2*v1*t1*s2 - w2*v1*t2*s1 +
      _two*v2*w1*t2*s1 + _two*w1*u1*s2*s2 - v2*u1*v1*s2 - _two*s1*s2*u1*w2 -
      v2*w1*t1*s2 + u2*v1*v1*s2 - t2*w1*v1*s2 - w2*t1*v2*s1 + 
      v2*v2*u1*s1 - u2*v1*v2*s1 - _two*s1*s2*w1*u2;
    
    c[0] = s2*v1*v1*w2 - s1*v2*v1*w2 - s2*v1*w1*v2 + s2*s2*w1*w1 -
      _two*s1*s2*w1*w2 + s1*w1*v2*v2 + s1*s1*w2*w2;

    degree = 4;
  }
  
  // Compute the roots of the resultant polynomial.
  typename Nt_traits::Polynomial  poly = 
                                    nt_traits.construct_polynomial (c, degree);

  xs_end = nt_traits.compute_polynomial_roots (poly,
					       xs);
  return static_cast<int>(xs_end - xs);
}

/*!
 * Compute the roots of the resultants of the two bivariate polynomials:
 *   C1:  r*x^2 + s*y^2 + t*xy + u*x + v*y + w = 0
 *   C2:  A*x + B*y + C = 0
 * \param deg1 The degree of the first curve.
 * \param xs Output: The real-valued roots of the polynomial, sorted in an
 *                   ascending order.
 * \pre xs must be a vector of size 4.
 * \return The number of distinct roots found.
 */
template <class Nt_traits>
int
_compute_resultant_roots (Nt_traits& nt_traits,
                          const typename Nt_traits::Algebraic& r,
			      const typename Nt_traits::Algebraic& s,
			      const typename Nt_traits::Algebraic& t,
			      const typename Nt_traits::Algebraic& u,
			      const typename Nt_traits::Algebraic& v,
			      const typename Nt_traits::Algebraic& w,
			      const int& deg1,
			      const typename Nt_traits::Algebraic& A,
			      const typename Nt_traits::Algebraic& B,
			      const typename Nt_traits::Algebraic& C,
			      typename Nt_traits::Algebraic *xs)
{
  if (deg1 == 1)
  {
    // We should actually compute the intersection of two line:
    // (u*x + v*y + w = 0) and (A*x + B*y + C = 0):
    const typename Nt_traits::Algebraic   denom = A*v - B*u;

    if (CGAL::sign (denom) == CGAL::ZERO)
      // The two lines are parallel and do not intersect.
      return (0);

    xs[0] = (B*w - C*v) / denom;
    return (1);
  }

  if (CGAL::sign (B) == CGAL::ZERO)
  {
    // The first line is A*x + C = 0, therefore:
    xs[0] = -C / A;
    return (1);
  }

  // We can write the first curve as: y = -(A*x + C) / B.
  // We substitute this expression into the equation of the conic, and get
  // the quadratic equation c[2]*x^2 + c[1]*x + c[0] = 0:
  const typename Nt_traits::Algebraic  _two = 2;
  typename Nt_traits::Algebraic        c[3];
  typename Nt_traits::Algebraic       *xs_end;

  c[2] = A*A*s - A*B*t + B*B*r;
  c[1] = _two*A*C*s - A*B*v - B*C*t + B*B*u;
  c[0] = C*C*s - B*C*v + B*B*w;

  xs_end = nt_traits.solve_quadratic_equation (c[2], c[1], c[0],
                                               xs);
  return static_cast<int>(xs_end - xs);
}

} //namespace CGAL

#endif
