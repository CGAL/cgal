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

#include <CGAL/Minkowski_sum_2/Minkowski_sum_by_reduced_convolution_2.h>
#include <CGAL/Minkowski_sum_2/Minkowski_sum_conv_2.h>
#include <CGAL/Minkowski_sum_2/Minkowski_sum_decomp_2.h>
#include <list>

namespace CGAL {

/*!
\ingroup PkgMinkowskiSum2

Computes the Minkowski sum \f$ P \oplus Q\f$ of the two given polygons.
The function computes the reduced convolution of the two polygons and
extracts those loops of the convolution which are part of the Minkowsi
sum. This method works very efficiently, regardless of whether `P` and
`Q` are convex or non-convex.
Note that as the input polygons may not be convex, their Minkowski
sum may not be a simple polygon. The result is therefore represented
as a polygon with holes.

\pre Both `P` and `Q` are simple, counterclockwise-oriented polygons.
*/

template <class Kernel, class Container>
Polygon_with_holes_2<Kernel,Container>
minkowski_sum_reduced_convolution_2 (const Polygon_2<Kernel,Container>& pgn1,
                                        const Polygon_2<Kernel,Container>& pgn2)
{
  Minkowski_sum_by_reduced_convolution_2<Kernel, Container> mink_sum;
  Polygon_2<Kernel,Container>                              sum_bound;
  std::list<Polygon_2<Kernel,Container> >                  sum_holes;

  if (pgn1.size() > pgn2.size())
    mink_sum (pgn1, pgn2, sum_bound, std::back_inserter(sum_holes));
  else
    mink_sum (pgn2, pgn1, sum_bound, std::back_inserter(sum_holes));

  return (Polygon_with_holes_2<Kernel,Container> (sum_bound,
                                                  sum_holes.begin(),
                                                  sum_holes.end()));
}

/*!
\ingroup PkgMinkowskiSum2

Computes the Minkowski sum \f$ P \oplus Q\f$ of the two given polygons.
The function computes the convolution cycles of the two polygons and
extract the regions having positive winding number with respect to these
cycles. This method work very efficiently, regardless of whether `P`
and `Q` are convex or non-convex.
Note that as the input polygons may not be convex, their Minkowski
sum may not be a simple polygon. The result is therefore represented
as a polygon with holes.

\pre Both `P` and `Q` are simple polygons.
*/

template <class Kernel, class Container>
Polygon_with_holes_2<Kernel,Container>
minkowski_sum_full_convolution_2 (const Polygon_2<Kernel,Container>& P,
                                     const Polygon_2<Kernel,Container>& Q)
{
  Minkowski_sum_by_convolution_2<Kernel, Container>  mink_sum;
  Polygon_2<Kernel,Container>                        sum_bound;
  std::list<Polygon_2<Kernel,Container> >            sum_holes;

  if (P.size() > Q.size())
    mink_sum (P, Q, sum_bound, std::back_inserter(sum_holes));
  else
    mink_sum (Q, P, sum_bound, std::back_inserter(sum_holes));

  return (Polygon_with_holes_2<Kernel,Container> (sum_bound,
                                                  sum_holes.begin(),
                                                  sum_holes.end()));
}

/*!
\ingroup PkgMinkowskiSum2

Computes the Minkowski sum \f$ P \oplus Q\f$ of the two given polygons
using the reduced convolution method.
Note that as the input polygons may not be convex, their Minkowski
sum may not be a simple polygon. The result is therefore represented
as a polygon with holes.

\pre Both `P` and `Q` are simple, counterclockwise-oriented polygons.

\sa `CGAL::minkowski_sum_reduced_convolution_2()`
\sa `CGAL::minkowski_sum_full_convolution_2()`
 */

template <class Kernel, class Container>
Polygon_with_holes_2<Kernel,Container>
minkowski_sum_2 (const Polygon_2<Kernel,Container>& P,
                 const Polygon_2<Kernel,Container>& Q)
{
  return minkowski_sum_reduced_convolution_2(P, Q);
}

/*!
\ingroup PkgMinkowskiSum2

Computes the Minkowski sum \f$ P \oplus Q\f$ of the two given polygons.
If the input polygons `P` and `Q` are not convex, the function
decomposes them into convex sub-polygons \f$ P_1, \ldots, P_k\f$ and
\f$ Q_1, \ldots, Q_{\ell}\f$ and computes the union of pairwise sub-sums
(namely \f$ \bigcup_{i,j}{(P_i \oplus Q_j)}\f$).
The decomposition is performed using the given decomposition strategy
`decomp`, which must be an instance of a class that models the
concept `PolygonConvexDecomposition`.
Note that as the input polygons may not be convex, their Minkowski
sum may not be a simple polygon. The result is therefore represented
as a polygon with holes.
\pre Both `P` and `Q` are simple polygons.
*/

template <class Kernel, class Container, class DecompositionStrategy>
Polygon_with_holes_2<Kernel,Container>
minkowski_sum_2 (const Polygon_2<Kernel,Container>& P,
                 const Polygon_2<Kernel,Container>& Q,
                 const DecompositionStrategy&)
{
  Minkowski_sum_by_decomposition_2<DecompositionStrategy,Container>        mink_sum;

  typedef Polygon_with_holes_2<Kernel,Container>                 Polygon_with_holes_2;
  
  Polygon_with_holes_2 sum;

  sum = mink_sum (P, Q);
  
  return (sum);
}

} //namespace CGAL

#endif
