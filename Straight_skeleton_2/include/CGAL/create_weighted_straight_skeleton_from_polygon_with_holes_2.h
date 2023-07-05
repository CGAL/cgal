// Copyright (c) 2006-2008 Fernando Luis Cacciola Carballal. All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//

// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé
//
#ifndef CGAL_CREATE_WEIGHTED_STRAIGHT_SKELETON_FROM_POLYGON_WITH_HOLES_2_H
#define CGAL_CREATE_WEIGHTED_STRAIGHT_SKELETON_FROM_POLYGON_WITH_HOLES_2_H

#include <CGAL/license/Straight_skeleton_2.h>

#include <CGAL/create_weighted_straight_skeleton_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Straight_skeleton_2/Polygon_iterators.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <type_traits>

namespace CGAL {

template <typename Polygon,
          typename Weights,
          typename K>
std::shared_ptr< Straight_skeleton_2<K> >
inline
create_interior_weighted_straight_skeleton_2(const Polygon& poly_with_holes,
                                             const Weights& weights,
                                             const K& k,
                                             std::enable_if_t<
                                               CGAL_SS_i::has_Hole_const_iterator<Polygon>::value>* = nullptr)
{
  return create_interior_weighted_straight_skeleton_2(poly_with_holes.outer_boundary().vertices_begin(),
                                                      poly_with_holes.outer_boundary().vertices_end(),
                                                      poly_with_holes.holes_begin(),
                                                      poly_with_holes.holes_end(),
                                                      std::begin(weights[0]), std::end(weights[0]),
                                                      std::next(std::begin(weights)), std::end(weights),
                                                      k);
}

// create_exterior_weightd_straight_skeleton_2() for polygon with holes is simply in create_straight_skeleton_2.h
// as the holes do not matter: call create_exterior_straight_skeleton_2 for each boundary (outer & holes).

} // namespace CGAL

#endif // CGAL_CREATE_WEIGHTED_STRAIGHT_SKELETON_FROM_POLYGON_WITH_HOLES_2_H
