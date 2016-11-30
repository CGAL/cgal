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

#ifndef CGAL_IS_PULLOUT_DIRECTION_SINGLE_MOLD_TRANSLATIONAL_CASTING_2_H
#define CGAL_IS_PULLOUT_DIRECTION_SINGLE_MOLD_TRANSLATIONAL_CASTING_2_H

namespace CGAL {

namespace Set_movable_separability_2 {

/*!
 */
template <typename CastingTraits_2>
bool is_pullout_direction_single_mold_translational_casting_2
(const CGAL::Polygon_2<CastingTraits_2>& pgn, size_t i,
 typename CastingTraits_2::Direction_2& d, CastingTraits_2& traits)
{
  return false;
}

/*!
 */
template <typename CastingTraits_2>
bool is_pullout_direction_single_mold_translational_casting_2
(const CGAL::Polygon_2<CastingTraits_2>& pgn, size_t i,
 typename CastingTraits_2::Direction_2& d)
{
  CastingTraits_2 traits;
  return is_pullout_direction_single_mold_translational_casting_2(pgn, i, d,
                                                                  traits);
}

/*!
 */
template <typename CastingTraits_2>
bool is_pullout_direction_single_mold_translational_casting_2
(const CGAL::Polygon_2<CastingTraits_2>& pgn,
 typename CastingTraits_2::Direction_2& d, CastingTraits_2& traits)
{
  return false;
}

/*!
 */
template <typename CastingTraits_2>
bool is_pullout_direction_single_mold_translational_casting_2
(const CGAL::Polygon_2<CastingTraits_2>& pgn, size_t i,
 typename CastingTraits_2::Direction_2& d)
{
  CastingTraits_2 traits;
  return is_pullout_direction_single_mold_translational_casting_2(pgn, d, traits);
}

} /* end namesapce Set_movable_separability_2 */
} /* end namesapce CGAL */

#endif
