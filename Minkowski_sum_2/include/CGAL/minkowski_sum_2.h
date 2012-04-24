// Copyright (c) 2006  Tel-Aviv University (Israel).
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
// Author(s)     : Ron Wein   <wein@post.tau.ac.il>

#ifndef CGAL_MINKOWSKI_SUM_2_H
#define CGAL_MINKOWSKI_SUM_2_H

#include <CGAL/basic.h>
#include <CGAL/Polygon_with_holes_2.h>

#include <CGAL/Minkowski_sum_2/Minkowski_sum_conv_2.h>
#include <CGAL/Minkowski_sum_2/Minkowski_sum_decomp_2.h>
#include <list>

namespace CGAL {

/*!
 * Compute the Minkowski sum of two simple polygons using the convolution
 * method.
 * Note that as the input polygons may not be convex, their Minkowski sum may
 * not be a simple polygon. The result is therefore represented as a polygon
 * with holes.
 * \param pgn1 The first polygon.
 * \param pgn2 The second polygon.
 * \return The resulting polygon with holes, representing the sum.
 */

template <class Kernel, class Container>
Polygon_with_holes_2<Kernel,Container>
minkowski_sum_2 (const Polygon_2<Kernel,Container>& pgn1,
                 const Polygon_2<Kernel,Container>& pgn2)
{
  Minkowski_sum_by_convolution_2<Kernel, Container>  mink_sum;
  Polygon_2<Kernel,Container>                        sum_bound;
  std::list<Polygon_2<Kernel,Container> >            sum_holes;

  if (pgn1.size() > pgn2.size())
    mink_sum (pgn1, pgn2, sum_bound, std::back_inserter(sum_holes));
  else
    mink_sum (pgn2, pgn1, sum_bound, std::back_inserter(sum_holes));

  return (Polygon_with_holes_2<Kernel,Container> (sum_bound,
                                                  sum_holes.begin(),
                                                  sum_holes.end()));
}

/*!
 * Compute the Minkowski sum of two simple polygons by decomposing each
 * polygon to convex sub-polygons and computing the union of the pairwise
 * Minkowski sums of the sub-polygons.
 * Note that as the input polygons may not be convex, their Minkowski sum may
 * not be a simple polygon. The result is therefore represented as a polygon
 * with holes.
 * \param pgn1 The first polygon.
 * \param pgn2 The second polygon.
 * \param decomp A functor for decomposing polygons.
 * \param sum Output: The resulting polygon with holes, representing the sum.
 */
template <class Kernel, class Container, class DecompositionStrategy>
Polygon_with_holes_2<Kernel,Container>
minkowski_sum_2 (const Polygon_2<Kernel,Container>& pgn1,
                 const Polygon_2<Kernel,Container>& pgn2,
                 const DecompositionStrategy&)
{
  Minkowski_sum_by_decomposition_2<DecompositionStrategy,Container>        mink_sum;

  typedef Polygon_with_holes_2<Kernel,Container>                 Polygon_with_holes_2;
  
  Polygon_with_holes_2 sum;

  sum = mink_sum (pgn1, pgn2);
  
  return (sum);
}

} //namespace CGAL

#endif
