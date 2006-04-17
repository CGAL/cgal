// Copyright (c) 2006  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source: $
// $Revision$ $Date$
// $Name:  $
//
// Author(s)     : Ron Wein   <wein@post.tau.ac.il>

#ifndef CGAL_RSA_MINKOWSKI_SUM_H
#define CGAL_RSA_MINKOWSKI_SUM_H

#include <CGAL/Minkowski_sum_2/RSA_Minkowski_sum_2.h>

CGAL_BEGIN_NAMESPACE

/*!
 * Compute the Minkowski sum of a linear polygon and a polygon that contains
 * circular arcs, which represents the rotational swept-area of a linear
 * polygon.
 * Note that as the input polygons may not be convex, the Minkowski sum may
 * not be a simple polygon. The result is therefore represented as
 * polygon with holes whose arcs are conic curves (which are either line
 * segments with irrational coefficients, or circular arcs with rational
 * coefficients).
 * \param pgn1 The linear polygon.
 * \param pgn2 The polygon with circular arcs.
 * \pre Both input polygons are simple.
 * \return A polygon with holes that reprsents the sum.
 */
template <class ConicTraits, class Container>
typename Gps_traits_2<ConicTraits>::Polygon_with_holes_2
rsa_minkowski_sum_2
    (const ConicTraits& ,
     const Polygon_2<typename ConicTraits::Rat_kernel, Container>& pgn1,
     const typename Gps_circle_segment_traits_2<typename 
                              ConicTraits::Rat_kernel>::Polygon_2& pgn2)
{
  typedef RSA_Minkowski_sum_2<ConicTraits, Container>        RSA_sum_2;
  typedef typename RSA_sum_2::Sum_polygon_2                  Sum_polygon_2;

  RSA_sum_2                    rsa_sum;
  Sum_polygon_2                sum_bound;
  std::list<Sum_polygon_2>     sum_holes;

  rsa_sum (pgn1, pgn2, 
           sum_bound, std::back_inserter(sum_holes));

  return (typename Gps_traits_2<ConicTraits>::Polygon_with_holes_2
          (sum_bound, sum_holes.begin(), sum_holes.end()));
}

CGAL_END_NAMESPACE

#endif
