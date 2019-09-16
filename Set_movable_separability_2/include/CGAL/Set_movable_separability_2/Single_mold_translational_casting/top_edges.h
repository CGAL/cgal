// Copyright (c) 2016 Tel-Aviv University (Israel).
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
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s): Shahar    <shasha94@gmail.com>
//            Efi Fogel <efif@gmail.com>

#ifndef CGAL_SMS_2_TOP_EDGES_SINGLE_MOLD_TRANSLATIONAL_CASTING_H
#define CGAL_SMS_2_TOP_EDGES_SINGLE_MOLD_TRANSLATIONAL_CASTING_H

#include <CGAL/license/Set_movable_separability_2.h>


#include <iostream>
#include <list>

#include <CGAL/Polygon_2.h>
#include <CGAL/enum.h>
#include <CGAL/Set_movable_separability_2/internal/Circle_arrangment.h>
#include <CGAL/Set_movable_separability_2/internal/Utils.h>

namespace CGAL {
namespace Set_movable_separability_2 {
namespace Single_mold_translational_casting {

/* Legend:
 * point = Represented as Direction_2. It is the intersection between the
 *   fitting Direction_2 and the unit circle
 *
 * Arc = Represented as A pair of point. clockwise arc between the first
 *   point and the second point. (each of its sides might be open or closed)
 *
 * SegmentOuterCircle  = Arc that represent all the directions that points
 *   out from the polygon if it start from the
 *   fitting segment. This arc is always open half circle.
 */

/*! \ingroup PkgSetMovableSeparability2Funcs
 *
 * Same as above with the additional `orientation` argument.
 * If the orientation of the polygon is known upon invocation, specify it.
 * Otherwise, it has to be computed.  Note that finding the orientation of a
 * polygon requires time linear in the number of edges.
 *
 * \param polygon the input polygon.
 * \param oi the output iterator. Its value type is a pair, where
 *         (i) the first element in the pair is an iterator to a top edge, and
 *         (ii) the second element is a closed range of pullout directions
 *              represented as a pair of the extreme directions in the range
 *              of type `CastingTraits::Direction_2`.
 * \param orientation the orientation of `polygon`.
 * \param traits the traits to use.
 * \return the past-the-end iterator of the output container.
 * \pre `polygon` must be non-degenerate (has at least 3 vertices), simple, and
 * does not have three consecutive collinear vertices.
 */
template <typename CastingTraits_2, typename OutputIterator>
OutputIterator top_edges(const CGAL::Polygon_2<CastingTraits_2>& polygon,
                         OutputIterator oi,
                         CGAL::Orientation orientation,
                         CastingTraits_2& traits)
{
  /* Legend
   * point = Represented as  Direction_2. It is the intersection between the
   *   fitting Direction_2 and the unit circle
   *
   * arc = Represented as A pair of point. clockwise arc between the first
   *   point and the second point. (each of its sides might be open or closed)
   */
  typedef CastingTraits_2                               Traits;

  CGAL_precondition(polygon.is_simple());
  CGAL_precondition(!internal::is_any_edge_colinear(polygon, traits));

  auto e_it = polygon.edges_begin();
  auto segment_outer_circle =
    internal::get_segment_outer_circle<Traits>(*e_it++, orientation);
  typedef internal::Circle_arrangment<Traits> Circle_arrangment;
  Circle_arrangment circle_arrangment(traits,
                                      segment_outer_circle,polygon.edges_begin());

  for (; e_it != polygon.edges_end(); ++e_it) {
    segment_outer_circle =
      internal::get_segment_outer_circle<Traits>(*e_it, orientation);
    circle_arrangment.add_segment_outer_circle(segment_outer_circle, e_it);
    if (circle_arrangment.all_is_covered_twice()) return oi;
  }
  circle_arrangment.get_all_1_edges(oi);
  return oi;
}

/*! \fn OutputIterator top_edges(const CGAL::Polygon_2<CastingTraits_2>& polygon,
 *                               OutputIterator oi, CastingTraits_2& traits)
 * \param polygon the input polygon that we want to check if is castable or not.
 * \param oi the output iterator to put the top edges in
 * \param traits the traits to use.
 * \return all the possible top edges of the polygon and there pullout direction
 *  a pair of Directions is build this way [firstClockwise,secondClockwise]
 *   (with no rotation)
 */
template <typename CastingTraits_2, typename OutputIterator>
OutputIterator top_edges(const CGAL::Polygon_2<CastingTraits_2>& polygon,
                         OutputIterator oi, CastingTraits_2& traits)
{
  CGAL::Orientation orientation = polygon.orientation();
  return top_edges(polygon, oi, orientation, traits);
}

/*! \fn OutputIterator top_edges(const CGAL::Polygon_2<CastingTraits_2>& polygon,
 *                               OutputIterator oi)
 * \param polygon the input polygon that we want to check if is castable or not.
 * \param oi the output iterator to put the top edges in
 * \param orientation the orientation of `polygon`.
 * \return all the possible top edges of the polygon and there pullout direction
 *  a pair of Directions is build this way [firstClockwise,secondClockwise]
 *   (with no rotation)
 */
template <typename CastingTraits_2, typename OutputIterator>
OutputIterator top_edges(const CGAL::Polygon_2<CastingTraits_2>& polygon,
                         OutputIterator oi, CGAL::Orientation orientation)
{
  CastingTraits_2 traits;
  return top_edges(polygon, oi, orientation, traits);
}

/*! \fn OutputIterator top_edges(const CGAL::Polygon_2<CastingTraits_2>& polygon,
 *                               OutputIterator oi)
 * \param polygon the input polygon that we want to check if is castable or not.
 * \param oi the output iterator to put the top edges in
 * \return all the possible top edges of the polygon and there pullout direction
 *  a pair of Directions is build this way [firstClockwise,secondClockwise]
 *   (with no rotation)
 */
template <typename CastingTraits_2, typename OutputIterator>
OutputIterator top_edges(const CGAL::Polygon_2<CastingTraits_2>& polygon,
                         OutputIterator oi)
{
  CGAL::Orientation orientation = polygon.orientation();
  CastingTraits_2 traits;
  return top_edges(polygon, oi, orientation, traits);
}

} // namespace Single_mold_translational_casting
} // namespace Set_movable_separability_2
} // namespace CGAL

#endif
