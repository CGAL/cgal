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
// Author(s): Shahar    <shasha94@gmail.com>
//            Efi Fogel <efif@gmail.com>

#ifndef CGAL_PULLOUT_DIRECTIONS_SINGLE_MOLD_TRANSLATIONAL_CASTING_2_H
#define CGAL_PULLOUT_DIRECTIONS_SINGLE_MOLD_TRANSLATIONAL_CASTING_2_H

#include <CGAL/Polygon_2.h>

namespace CGAL {

namespace Set_movable_separability_2 {

/*! Given a simple polygon and an edge of the polygon, this function determines
 * whether a cavity (of a mold in the plane) that has the shape of the polygon
 * can be used so that the polygon could be casted in the mold with the input
 * edge being the top edge and then pulled out of the mold without colliding
 * into the mold (but possibly sliding along the mold surface). If the polygon
 * is <em>castable</em> this way, the function computes the closed range of pull
 * directions.
 *
 * The type that substitutes the template parameter `%CastingTraits_2` must be
 * a model of the concept `CastingTraits_2`.
 *
 * \param[in] pgn the input polygon.
 * \param[in] i the index of an edge in pgn.
 * \return a pair of elements, where the first is a Boolean that indicates
 *         whether the input edge is a valid top edge, and the second
 *         is a closed range of pull-out directions represented as a pair
 *         of the extreme directions in the range. If the input edge is not
 *         a valid top edge, the range is nondeterministic.
 *
 * \pre `png` must be non-degenerate (has at least 3 vertices), simple, and
 * does not have three consecutive collinear vertices.
 */
template <typename CastingTraits_2>
std::pair<bool, std::pair<typename CastingTraits_2::Direction_2,
                          typename CastingTraits_2::Direction_2> >
pullout_directions_single_mold_translational_casting_2
(const CGAL::Polygon_2<CastingTraits_2>& pgn, size_t i, Kernel& kernel)
{
  typedef CastingTraits_2               Casting_traits_2;
  typename Casting_traits_2::Direction_2 d1, d2;
  return std::make_pair(false, std::make_pair(d1, d2));
}

template <typename CastingTraits_2>
std::pair<bool, std::pair<typename CastingTraits_2::Direction_2,
                          typename CastingTraits_2::Direction_2> >
/*! Same as above with the additional traits argument.
 * \param[in] pgn the input polygon.
 * \param[in] i the index of an edge in pgn.
 * \param[in] traits the traits to use.
 */
pullout_directions_single_mold_translational_casting_2
(const CGAL::Polygon_2<CastingTraits_2>& pgn, size_t i)
{
  CastingTraits_2 traits;
  return pullout_directions_single_mold_translational_casting_2(pgn, i, traits);
}

} // end of namespace Set_movable_separability_2
} // end of namespace CGAL

#endif
