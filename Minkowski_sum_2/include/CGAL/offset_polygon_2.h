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

#ifndef CGAL_OFFSET_POLYGON_H
#define CGAL_OFFSET_POLYGON_H

#include <CGAL/license/Minkowski_sum_2.h>


#include <CGAL/Minkowski_sum_2/Exact_offset_base_2.h>
#include <CGAL/Minkowski_sum_2/Offset_conv_2.h>
#include <CGAL/Minkowski_sum_2/Offset_decomp_2.h>

namespace CGAL {

/*!
 * Compute the offset of a given simple polygon by a given radius,
 * using the convolution method.
 * Note that as the input polygon may not be convex, its offset may not be
 * simply connected. The result is therefore represented as a polygon with
 * holes.
 * \param pgn The polygon.
 * \param r The offset radius.
 * \return The offset polygon.
 */
template <class ConicTraits, class Container>
typename Gps_traits_2<ConicTraits>::Polygon_with_holes_2
offset_polygon_2 (const Polygon_2<typename ConicTraits::Rat_kernel,
                                  Container>& pgn,
                  const typename ConicTraits::Rat_kernel::FT& r,
                  const ConicTraits& )
{
  typedef Exact_offset_base_2<ConicTraits, Container>        Base;
  typedef Offset_by_convolution_2<Base>                      Exact_offset_2;
  typedef typename Exact_offset_2::Offset_polygon_2          Offset_polygon_2;

  Base                                               base;
  Exact_offset_2                                     exact_offset (base);
  Offset_polygon_2                                   offset_bound;
  std::list<Offset_polygon_2>                        offset_holes;

  exact_offset (pgn, r,
                offset_bound, std::back_inserter(offset_holes));

  return (typename Gps_traits_2<ConicTraits>::Polygon_with_holes_2
          (offset_bound, offset_holes.begin(), offset_holes.end()));
}

/*!
 * Compute the offset of a given polygon with holes by a given radius,
 * using the convolution method.
 * The result is represented as a polygon with holes whose edges are line
 * segments and circular arcs.
 * \param pwh The polygon with holes.
 * \param r The offset radius.
 * \pre The polygon is bounded (has a valid outer boundary).
 * \return The offset polygon.
 */
template <class ConicTraits, class Container>
typename Gps_traits_2<ConicTraits>::Polygon_with_holes_2
offset_polygon_2 (const Polygon_with_holes_2<typename ConicTraits::Rat_kernel,
                                             Container>& pwh,
                  const typename ConicTraits::Rat_kernel::FT& r,
                  const ConicTraits& )
{
  typedef Exact_offset_base_2<ConicTraits, Container>        Base;
  typedef Offset_by_convolution_2<Base>                      Exact_offset_2;
  typedef typename Exact_offset_2::Offset_polygon_2          Offset_polygon_2;

  Base                                               base;
  Exact_offset_2                                     exact_offset (base);
  Offset_polygon_2                                   offset_bound;
  std::list<Offset_polygon_2>                        offset_holes;

  exact_offset (pwh, r,
                offset_bound, std::back_inserter(offset_holes));

  return (typename Gps_traits_2<ConicTraits>::Polygon_with_holes_2
          (offset_bound, offset_holes.begin(), offset_holes.end()));
}

/*!
 * Compute the offset of a given simple polygon by a given radius,
 * by decomposing it to convex sub-polygons and computing the union of their
 * offsets.
 * Note that as the input polygon may not be convex, its offset may not be
 * simply connected. The result is therefore represented as a polygon with
 * holes.
 * \param pgn The polygon.
 * \param r The offset radius.
 * \param decomp A functor for decomposing polygons.
 * \return The offset polygon.
 */
template <class ConicTraits, class Container, class DecompositionStrategy>
typename Gps_traits_2<ConicTraits>::Polygon_with_holes_2
offset_polygon_2 (const Polygon_2<typename ConicTraits::Rat_kernel,
                                  Container>& pgn,
                  const typename ConicTraits::Rat_kernel::FT& r,
                  const DecompositionStrategy&,
                  const ConicTraits& )
{
  typedef Exact_offset_base_2<ConicTraits, Container>        Base;
  typedef Offset_by_decomposition_2<Base, DecompositionStrategy>
                                                             Exact_offset_2;
  typedef typename Exact_offset_2::Offset_polygon_2          Offset_polygon_2;

  Base                                               base;
  Exact_offset_2                                     exact_offset (base);
  Offset_polygon_2                                   offset_bound;
  std::list<Offset_polygon_2>                        offset_holes;

  exact_offset (pgn, r,
                offset_bound, std::back_inserter(offset_holes));

  return (typename Gps_traits_2<ConicTraits>::Polygon_with_holes_2
          (offset_bound, offset_holes.begin(), offset_holes.end()));
}

/*!
 * Compute the inset of a given simple polygon by a given radius, using the
 * convolution method.
 * Note that as the input polygon may not be convex, its inset may not be
 * simply connected. The result is therefore represented as a set of polygons.
 * \param pgn The polygon.
 * \param r The inset radius.
 * \param oi An output iterator for the inset polygons.
 *           Its value-type must be Gps_traits_2<ConicTraits>::Polygon_2.
 * \return A past-the-end iterator for the inset polygons.
 */
template <class ConicTraits, class Container, class OutputIterator>
OutputIterator
inset_polygon_2 (const Polygon_2<typename ConicTraits::Rat_kernel,
                                 Container>& pgn,
                 const typename ConicTraits::Rat_kernel::FT& r,
                 const ConicTraits& ,
                 OutputIterator oi)
{
  typedef Exact_offset_base_2<ConicTraits, Container>        Base;
  typedef Offset_by_convolution_2<Base>                      Exact_offset_2;

  Base                                               base;
  Exact_offset_2                                     exact_offset (base);

  oi = exact_offset.inset (pgn, r,
                           oi);

  return (oi);
}

} //namespace CGAL

#endif
