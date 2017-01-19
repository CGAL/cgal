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
// Author(s)     : Ron Wein   <wein_r@yahoo.com>

#ifndef CGAL_APPROXIMATED_OFFSET_H
#define CGAL_APPROXIMATED_OFFSET_H

#include <CGAL/license/Minkowski_sum_2.h>


#include <CGAL/Minkowski_sum_2/Approx_offset_base_2.h>
#include <CGAL/Minkowski_sum_2/Offset_conv_2.h>
#include <CGAL/Minkowski_sum_2/Offset_decomp_2.h>

namespace CGAL {

/*!
 * Approximate the offset of a given simple polygon by a given radius,
 * using the convolution method.
 * Note that as the input polygon may not be convex, its offset may not be
 * simply connected. The result is therefore represented as a polygon with
 * holes.
 * \param pgn The polygon.
 * \param r The offset radius.
 * \param eps The approximation-error bound.
 * \return The approximated offset polygon.
 */
template <class Kernel, class Container>
typename Gps_circle_segment_traits_2<Kernel>::Polygon_with_holes_2
approximated_offset_2 (const Polygon_2<Kernel, Container>& pgn,
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

  approx_offset (pgn, r,
                 offset_bound, std::back_inserter(offset_holes));

  return (typename Gps_circle_segment_traits_2<Kernel>::Polygon_with_holes_2
          (offset_bound, offset_holes.begin(), offset_holes.end()));
}

/*!
 * Approximate the offset of a given polygon with holes by a given radius,
 * using the convolution method.
 * The result is represented as a polygon with holes whose edges are line
 * segments and circular arcs.
 * \param pwh The polygon with holes.
 * \param r The offset radius.
 * \param eps The approximation-error bound.
 * \pre The polygon is bounded (has a valid outer boundary).
 * \return The approximated offset polygon.
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
 * Approximate the offset of a given simple polygon by a given radius,
 * by decomposing it to convex sub-polygons and computing the union of their
 * offsets.
 * Note that as the input polygon may not be convex, its offset may not be
 * simply connected. The result is therefore represented as a polygon with
 * holes.
 * \param pgn The polygon.
 * \param r The offset radius.
 * \param eps The approximation-error bound.
 * \param decomp A functor for decomposing polygons.
 * \return The approximated offset polygon.
 */
template <class Kernel, class Container, class DecompositionStrategy>
typename Gps_circle_segment_traits_2<Kernel>::Polygon_with_holes_2
approximated_offset_2 (const Polygon_2<Kernel, Container>& pgn,
                       const typename Kernel::FT& r,
                       const double& eps,
                       const DecompositionStrategy&)
{
  typedef Approx_offset_base_2<Kernel, Container>            Base;
  typedef Offset_by_decomposition_2<Base, DecompositionStrategy>
                                                             Approx_offset_2;
  typedef typename Approx_offset_2::Offset_polygon_2         Offset_polygon_2;

  Base                                               base (eps);
  Approx_offset_2                                    approx_offset (base);
  Offset_polygon_2                                   offset_bound;
  std::list<Offset_polygon_2>                        offset_holes;

  approx_offset (pgn, r,
                 offset_bound, std::back_inserter(offset_holes));

  return (typename Gps_circle_segment_traits_2<Kernel>::Polygon_with_holes_2
          (offset_bound, offset_holes.begin(), offset_holes.end()));
}

/*!
 * Approximate the inset of a given simple polygon by a given radius, using
 * the convolution method.
 * Note that as the input polygon may not be convex, its inset may not be
 * simply connected. The result is therefore represented as a set of polygons.
 * \param pgn The polygon.
 * \param r The inset radius.
 * \param eps The approximation-error bound.
 * \param oi An output iterator for the inset polygons.
 *           Its value-type must be
 *           Gps_circle_segment_traits_2<Kernel>::Polygon_2.
 * \return A past-the-end iterator for the inset polygons.
 */
template <class Kernel, class Container, class OutputIterator>
OutputIterator
approximated_inset_2 (const Polygon_2<Kernel, Container>& pgn,
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

  oi = approx_offset.inset (pgn, r,
                            oi);

  return (oi);
}

} //namespace CGAL

#endif
