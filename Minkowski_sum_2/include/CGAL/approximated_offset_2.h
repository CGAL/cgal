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

#ifndef CGAL_APPROXIMATED_OFFSET_H
#define CGAL_APPROXIMATED_OFFSET_H

#include <CGAL/Minkowski_sum_2/Approx_offset_base_2.h>
#include <CGAL/Minkowski_sum_2/Offset_conv_2.h>
#include <CGAL/Minkowski_sum_2/Offset_decomp_2.h>

namespace CGAL {

/*!
\ingroup PkgMinkowskiSum2

Provides a guaranteed approximation of the offset of the given polygon
`P` by a given radius `r` - namely, the function computes the
Minkowski sum \f$ P \oplus B_r\f$, where \f$ B_r\f$ is a disc of radius
`r` centered at the origin.
The function actually outputs a set \f$ S\f$ that contains the Minkowski sum,
such that the approximation error is bounded by `eps`.
Note that as the input polygon may not be convex, its offset may not be a
simple polygon. The result is therefore represented as a polygon with
holes, whose edges are either line segments or circular arcs.
\pre `P` is a simple polygon.
*/
template <class Kernel, class Container>
typename Gps_circle_segment_traits_2<Kernel>::Polygon_with_holes_2
approximated_offset_2 (const Polygon_2<Kernel, Container>& P,
                       const typename Kernel::FT& r,
                       const double& eps)
{
  typedef Approx_offset_base_2<Kernel, Container>            Base;
  typedef Offset_by_convolution_2<Base>                      Approx_offset_2;
  typedef typename Approx_offset_2::Offset_polygon_2         Offset_polygon_2;

  Base                                               base (eps);
  Approx_offset_2                                    approx_offset (base);
  Offset_polygon_2                                   offset_bound;
  std::list<Offset_polygon_2>                        offset_holes;

  approx_offset (P, r,
                 offset_bound, std::back_inserter(offset_holes));

  return (typename Gps_circle_segment_traits_2<Kernel>::Polygon_with_holes_2
          (offset_bound, offset_holes.begin(), offset_holes.end()));
}

/*!
\ingroup PkgMinkowskiSum2

Provides a guaranteed approximation of offset the given polygon with holes
`pwh` by a given radius `r`, such that the approximation error is bounded
by `eps`. It does so by offsetting outer boundary of `pwh` and insetting
its holes.
The result is represented as a generalized polygon with holes, such that the edges
of the polygon correspond to line segment and circular arcs.
\pre `pwh` is <I>not</I> unbounded (it has a valid outer boundary).
*/
template <class Kernel, class Container>
typename Gps_circle_segment_traits_2<Kernel>::Polygon_with_holes_2
approximated_offset_2 (const Polygon_with_holes_2<Kernel, Container>& pwh,
                       const typename Kernel::FT& r,
                       const double& eps)
{
  typedef Approx_offset_base_2<Kernel, Container>            Base;
  typedef Offset_by_convolution_2<Base>                      Approx_offset_2;
  typedef typename Approx_offset_2::Offset_polygon_2         Offset_polygon_2;

  Base                                               base (eps);
  Approx_offset_2                                    approx_offset (base);
  Offset_polygon_2                                   offset_bound;
  std::list<Offset_polygon_2>                        offset_holes;

  approx_offset (pwh, r,
                 offset_bound, std::back_inserter(offset_holes));

  return (typename Gps_circle_segment_traits_2<Kernel>::Polygon_with_holes_2
          (offset_bound, offset_holes.begin(), offset_holes.end()));
}

/*!
\ingroup PkgMinkowskiSum2

Provides a guaranteed approximation of the offset of the given polygon
`P` by a radius `r`, as described above.
If the input polygon `P` is not convex, the function
decomposes it into convex sub-polygons \f$ P_1, \ldots, P_k\f$ and computes
the union of the sub-offsets (namely \f$ \bigcup_{i}{(P_i \oplus B_r)}\f$).
The decomposition is performed using the given decomposition strategy
`decomp`, which must be an instance of a class that models the
concept `PolygonConvexDecomposition`.
\pre `P` is a simple polygon.
*/
template <class Kernel, class Container, class DecompositionStrategy>
typename Gps_circle_segment_traits_2<Kernel>::Polygon_with_holes_2
approximated_offset_2 (const Polygon_2<Kernel, Container>& P,
                       const typename Kernel::FT& r,
                       const double& eps,
                       const DecompositionStrategy& decomp)
{
  typedef Approx_offset_base_2<Kernel, Container>            Base;
  typedef Offset_by_decomposition_2<Base, DecompositionStrategy>
                                                             Approx_offset_2;
  typedef typename Approx_offset_2::Offset_polygon_2         Offset_polygon_2;

  Base                                               base (eps);
  Approx_offset_2                                    approx_offset (base);
  Offset_polygon_2                                   offset_bound;
  std::list<Offset_polygon_2>                        offset_holes;

  approx_offset (P, r,
                 offset_bound, std::back_inserter(offset_holes));

  return (typename Gps_circle_segment_traits_2<Kernel>::Polygon_with_holes_2
          (offset_bound, offset_holes.begin(), offset_holes.end()));
}

/*!
\ingroup PkgMinkowskiSum2

Provides a guaranteed approximation of the inset, or inner offset, of
the given polygon `P` by a given radius `r`. Namely, the
function computes the set of points inside the polygon whose distance
from \f$ P\f$'s boundary is at least \f$ r\f$:
\f$ \{ p \in P \;|\; {\rm dist}(p, \partial P) \geq r \}\f$,
with the approximation error bounded by `eps`.
Note that as the input polygon may not be convex, its inset may comprise
several disconnected components. The result is therefore represented as a
sequence of generalized polygons, whose edges are either line segments or
circular arcs.
The output sequence is returned via the output iterator `oi`, whose
value-type must be `Gps_circle_segment_traits_2::Polygon_2`.
\pre `P` is a simple polygon.
*/
template <class Kernel, class Container, class OutputIterator>
OutputIterator
approximated_inset_2 (const Polygon_2<Kernel, Container>& P,
                      const typename Kernel::FT& r,
                      const double& eps,
                      OutputIterator oi)
{
  typedef Approx_offset_base_2<Kernel, Container>            Base;
  typedef Offset_by_convolution_2<Base>                      Approx_offset_2;
  typedef typename Approx_offset_2::Offset_polygon_2         Offset_polygon_2;

  Base                                               base (eps);
  Approx_offset_2                                    approx_offset (base);
  Offset_polygon_2                                   offset_bound;
  std::list<Offset_polygon_2>                        offset_holes;

  oi = approx_offset.inset (P, r,
                            oi);

  return (oi);
}

} //namespace CGAL

#endif
