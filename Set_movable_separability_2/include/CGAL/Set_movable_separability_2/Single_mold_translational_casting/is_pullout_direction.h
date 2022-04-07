// Copyright (c) 2016 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Shahar    <shasha94@gmail.com>
//            Efi Fogel <efif@gmail.com>

#ifndef CGAL_SMS_2_IS_PULLOUT_DIRECTION_SINGLE_MOLD_TRANSLATIONAL_CASTING_H
#define CGAL_SMS_2_IS_PULLOUT_DIRECTION_SINGLE_MOLD_TRANSLATIONAL_CASTING_H

#include <CGAL/license/Set_movable_separability_2.h>


#include <limits>

#include <CGAL/Polygon_2.h>
#include <CGAL/enum.h>
#include <CGAL/Set_movable_separability_2/internal/Utils.h>
#include <CGAL/Set_movable_separability_2/internal/Circle_arrangment.h>

namespace CGAL {
namespace Set_movable_separability_2 {
namespace Single_mold_translational_casting {

/*! Given a simple polygon, an edge of the polygon and a pullout direction (not
 * rotated) this function determines whether a cavity (of a mold in the plane)
 * that has the shape of the polygon can be used so that the polygon could be
 * casted in the mold with the input edge and being the top edge and then pulled
 * out in the input direction (without rotation) of the mold without colliding
 * into the mold (but possibly sliding along the mold surface).
 *
 * The type that substitutes the template parameter `%CastingTraits_2` must be
 * a model of the concept `CastingTraits_2`.
 *
 * \param[in] pgn the input polygon.
 * \param[in] it an iterator to an edge in pgn.
 * \param[in] orientation the orientation of `pgn`.
 * \param[in] d the pullout direction
 * \return if the polygon can be pullout through edge i with direction d
 *
 * \pre `png` must be non-degenerate (has at least 3 vertices), simple, and
 * does not have three consecutive collinear vertices.
 */
template <typename CastingTraits_2>
bool is_pullout_direction
(const CGAL::Polygon_2<CastingTraits_2>& pgn,
 const typename CGAL::Polygon_2<CastingTraits_2>::Edge_const_iterator& i,
 const typename CastingTraits_2::Direction_2& d,
 CGAL::Orientation orientation, CastingTraits_2& traits)
{
  //NOT CHECKED AT ALL
  CGAL_precondition(pgn.is_simple());
  CGAL_precondition(!internal::is_any_edge_colinear(pgn, traits));

  auto e_it = pgn.edges_begin();
  auto cc_in_between = traits.counterclockwise_in_between_2_object();

  for (; e_it != pgn.edges_end(); ++e_it) {
    auto segment_outer_circle =
      internal::get_segment_outer_circle<CastingTraits_2>(*e_it, orientation);
    bool isordered = !cc_in_between(d, segment_outer_circle.second,
                                    segment_outer_circle.first);
    if (isordered == (e_it == i)) return false;
  }

  return true;
}

/*! Same as above without the traits argument.
 */
template <typename CastingTraits_2>
bool is_pullout_direction
(const CGAL::Polygon_2<CastingTraits_2>& pgn,
 const typename CGAL::Polygon_2<CastingTraits_2>::Edge_const_iterator& i,
 const typename CastingTraits_2::Direction_2& d, CGAL::Orientation orientation)
{
  CastingTraits_2 traits;
  return is_pullout_direction(pgn, i, d, orientation, traits);
}

/*! Same as above without the orientation and traits arguments.
 */
template <typename CastingTraits_2>
bool is_pullout_direction
(const CGAL::Polygon_2<CastingTraits_2>& pgn,
 const typename CGAL::Polygon_2<CastingTraits_2>::Edge_const_iterator& i,
 const typename CastingTraits_2::Direction_2& d)
{
  CGAL::Orientation orientation = pgn.orientation();
  CastingTraits_2 traits;
  return is_pullout_direction(pgn, i, d, orientation, traits);
}

/*! Same as above without the orientation argument.
 */
template <typename CastingTraits_2>
bool is_pullout_direction
(const CGAL::Polygon_2<CastingTraits_2>& pgn,
 const typename CGAL::Polygon_2<CastingTraits_2>::Edge_const_iterator& i,
 const typename CastingTraits_2::Direction_2& d, CastingTraits_2& traits)
{
  CGAL::Orientation orientation = pgn.orientation();
  return is_pullout_direction(pgn, i, d, orientation, traits);
}

/*! Given a simple polygon, and a pullout direction (not rotated)
 * this function determines whether a cavity (of a mold in the plane)
 * that has the shape of the polygon can be used so that the polygon could be
 * casted in the mold with the input edge and being the top edge and then pulled
 * out in the input direction (without rotation) of the mold without colliding
 * into the mold (but possibly sliding along the mold surface).
 *
 * The type that substitutes the template parameter `%CastingTraits_2` must be
 * a model of the concept `CastingTraits_2`.
 *
 * \param[in] pgn the input polygon.
 * \param[in] d the pullout direction
 * \return if the polygon can be pullout through some edge with direction d
 *         the top edge, otherwise, pgn.edges_end()
 *
 * \pre `png` must be non-degenerate (has at least 3 vertices),simple, and
 * does not have three consecutive collinear vertices.
 */
template <typename CastingTraits_2>
typename CGAL::Polygon_2<CastingTraits_2>::Edge_const_iterator
is_pullout_direction(const CGAL::Polygon_2<CastingTraits_2>& pgn,
                     typename CastingTraits_2::Direction_2& d,
                     CGAL::Orientation orientation, CastingTraits_2& traits)
{
  //NOT CHECKED AT ALL
  typedef CGAL::Polygon_2<CastingTraits_2>              Polygon_2;
  typedef typename Polygon_2::Edge_const_iterator       Edge_iter;

  CGAL_precondition(pgn.is_simple());
  CGAL_precondition(!internal::is_any_edge_colinear(pgn, traits));

  Edge_iter e_it = pgn.edges_begin();
  auto segment_outer_circle =
    internal::get_segment_outer_circle<CastingTraits_2>(*e_it++, orientation);
  auto cc_in_between = traits.counterclockwise_in_between_2_object();
  Edge_iter top_edge= pgn.edges_end();
  for (; e_it != pgn.edges_end(); ++e_it) {
    segment_outer_circle =
      internal::get_segment_outer_circle<CastingTraits_2>(*e_it, orientation);
    bool isordered = !cc_in_between(d,
                                    segment_outer_circle.second,
                                    segment_outer_circle.first);
    if (!isordered) {
      // unlikely, this if must be true atleast once for any polygon - add ref
      // to paper
      if (top_edge== pgn.edges_end()) top_edge=e_it;
      else return pgn.edges_end();
    }
  }
  CGAL_postcondition(top_edge!=pgn.edges_end());
  return top_edge;
}

/*! Same as above without the traits argument.
 */
template <typename CastingTraits_2>
typename CGAL::Polygon_2<CastingTraits_2>::Edge_const_iterator
is_pullout_direction(const CGAL::Polygon_2<CastingTraits_2>& pgn,
                     typename CastingTraits_2::Direction_2& d,
                     CGAL::Orientation orientation)
{
  CastingTraits_2 traits;
  return is_pullout_direction(pgn, d, orientation, traits);
}

/*! Same as above without the orientation argument.
 */
template <typename CastingTraits_2>
typename CGAL::Polygon_2<CastingTraits_2>::Edge_const_iterator
is_pullout_direction(const CGAL::Polygon_2<CastingTraits_2>& pgn,
                     typename CastingTraits_2::Direction_2& d,
                     CastingTraits_2& traits)
{
  CGAL::Orientation orientation = pgn.orientation();
  return is_pullout_direction(pgn, d, orientation, traits);
}

/*! Same as above without the orientation and traits arguments.
 */
template <typename CastingTraits_2>
typename CGAL::Polygon_2<CastingTraits_2>::Edge_const_iterator
is_pullout_direction(const CGAL::Polygon_2<CastingTraits_2>& pgn,
                     typename CastingTraits_2::Direction_2& d)
{
  CGAL::Orientation orientation = pgn.orientation();
  CastingTraits_2 traits;
  return is_pullout_direction(pgn, d, orientation, traits);
}

} // namespace Single_mold_translational_casting
} // namesapce Set_movable_separability_2
} // namesapce CGAL

#endif
