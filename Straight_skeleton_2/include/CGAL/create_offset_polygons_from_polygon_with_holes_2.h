// Copyright (c) 2006-2008 Fernando Luis Cacciola Carballal. All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//

// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//
#ifndef CGAL_CREATE_OFFSET_POLYGONS_FROM_POLYGON_WITH_HOLES_2_H
#define CGAL_CREATE_OFFSET_POLYGONS_FROM_POLYGON_WITH_HOLES_2_H

#include <CGAL/license/Straight_skeleton_2.h>

#include <CGAL/Straight_skeleton_2/Straight_skeleton_aux.h>
#include <CGAL/arrange_offset_polygons_2.h>
#include <CGAL/create_offset_polygons_2.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <boost/shared_ptr.hpp>

#include <type_traits>
#include <vector>

namespace CGAL {

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/// INTERIOR

/*! create_interior_skeleton_and_offset_polygons_2 (no sorting of the result) */

// overload where PolygonWithHoles actually is a type of Polygon that supports holes
template <class OutPolygon_ = CGAL::Default,
          class FT, class PolygonWithHoles,
          class OfK = Exact_predicates_inexact_constructions_kernel,
          class SsK = Exact_predicates_inexact_constructions_kernel>
std::vector<boost::shared_ptr<CGAL_SS_i::Polygon_return_type<OutPolygon_, PolygonWithHoles, OfK> > >
inline
create_interior_skeleton_and_offset_polygons_2(const FT& aOffset,
                                               const PolygonWithHoles& aPoly,
                                               const OfK& ofk = OfK(),
                                               const SsK& ssk = SsK(),
                                               std::enable_if_t<
                                                 CGAL_SS_i::has_Hole_const_iterator<PolygonWithHoles>::value>* = nullptr)
{
  using OutPolygon = CGAL_SS_i::Polygon_return_type<OutPolygon_, PolygonWithHoles, OfK>;

  return create_interior_skeleton_and_offset_polygons_2<OutPolygon>(aOffset, aPoly.outer_boundary(),
                                                                    aPoly.holes_begin(), aPoly.holes_end(),
                                                                    ofk, ssk);
}

/*! create_interior_skeleton_and_offset_polygons_with_holes_2 (orders the resulting polygons) */

// 'Polygon' might be a polygon with holes or not, but it returns a polygon with holes
template <class OutPolygonWithHoles_ = CGAL::Default,
          class FT, class Polygon,
          class OfK = Exact_predicates_inexact_constructions_kernel,
          class SsK = Exact_predicates_inexact_constructions_kernel>
std::vector<boost::shared_ptr<CGAL_SS_i::Polygon_with_holes_return_type<OutPolygonWithHoles_, Polygon, OfK> > >
inline
create_interior_skeleton_and_offset_polygons_with_holes_2(const FT& aOffset,
                                                          const Polygon& aPoly,
                                                          const OfK& ofk = OfK(),
                                                          const SsK& ssk = SsK())
{
  using OutPolygon = typename CGAL_SS_i::Default_return_polygon_type<Polygon, OfK>::type;
  using OutPolygonWithHoles = CGAL_SS_i::Polygon_with_holes_return_type<OutPolygonWithHoles_, Polygon, OfK>;

  return arrange_offset_polygons_2<OutPolygonWithHoles>(
           create_interior_skeleton_and_offset_polygons_2<OutPolygon>(aOffset, aPoly, ofk, ssk));
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/// EXTERIOR

/*! create_exterior_skeleton_and_offset_polygons_with_holes_2 (orders the resulting polygons) */

// Polygon might be a Polygon with holes or not, but it returns a Polygon with holes
template <class OutPolygonWithHoles_ = CGAL::Default,
          class FT, class Polygon,
          class OfK = Exact_predicates_inexact_constructions_kernel,
          class SsK = Exact_predicates_inexact_constructions_kernel>
std::vector<boost::shared_ptr<CGAL_SS_i::Polygon_with_holes_return_type<OutPolygonWithHoles_, Polygon, OfK> > >
inline
create_exterior_skeleton_and_offset_polygons_with_holes_2(const FT& aOffset,
                                                          const Polygon& aPoly,
                                                          const OfK& ofk = OfK(),
                                                          const SsK& ssk = SsK())
{
  using OutPolygon = typename CGAL_SS_i::Default_return_polygon_type<Polygon, OfK>::type;
  using OutPolygonWithHoles = CGAL_SS_i::Polygon_with_holes_return_type<OutPolygonWithHoles_, Polygon, OfK>;

  std::vector<boost::shared_ptr<OutPolygon> > raw_output =
    create_exterior_skeleton_and_offset_polygons_2<OutPolygon>(aOffset, aPoly, ofk, ssk);

  CGAL_postcondition(raw_output.size() >= 2);

  // filter offset of the outer frame
  std::swap(raw_output[0], raw_output.back());
  raw_output.pop_back();

  for (boost::shared_ptr<OutPolygon> ptr : raw_output) {
    if (ptr->size() > 1) {
      // keep the first in place is just to get the same behavior as for Polygon_2
      auto first = std::next(ptr->begin());
      std::reverse(first, ptr->end());
    }
  }

  return arrange_offset_polygons_2<OutPolygonWithHoles>(raw_output);
}

/*! create_interior_skeleton_and_offset_polygons_2 with a polygon with holes */

// overload where PolygonWithHoles actually is a type of Polygon that supports holes
template <class OutPolygon_ = CGAL::Default,
          class FT, class PolygonWithHoles,
          class OfK = Exact_predicates_inexact_constructions_kernel,
          class SsK = Exact_predicates_inexact_constructions_kernel>
std::vector<boost::shared_ptr<CGAL_SS_i::Polygon_return_type<OutPolygon_, PolygonWithHoles, OfK> > >
inline
create_exterior_skeleton_and_offset_polygons_2(const FT& aOffset,
                                               const PolygonWithHoles& aPoly,
                                               const OfK& ofk = OfK(),
                                               const SsK& ssk = SsK(),
                                               std::enable_if_t<
                                                 CGAL_SS_i::has_Hole_const_iterator<PolygonWithHoles>::value>* = nullptr)
{
  using OutPolygon = CGAL_SS_i::Polygon_return_type<OutPolygon_, PolygonWithHoles, OfK>;

  std::vector<boost::shared_ptr<OutPolygon> > polygons =
    create_exterior_skeleton_and_offset_polygons_2<OutPolygon>(aOffset, aPoly.outer_boundary(), ofk, ssk);

  for (typename PolygonWithHoles::Hole_const_iterator hit=aPoly.holes_begin(); hit!=aPoly.holes_end(); ++hit)
  {
    typename PolygonWithHoles::Polygon_2 hole = *hit;
    hole.reverse_orientation();
    std::vector<boost::shared_ptr<OutPolygon> > hole_polygons =
        create_interior_skeleton_and_offset_polygons_2<OutPolygon>(aOffset,
                                                                   hole,
                                                                   ofk, ssk);
    polygons.insert(polygons.end(), hole_polygons.begin(), hole_polygons.end());
  }

  return polygons;
}

} // namespace CGAL

#endif // CGAL_CREATE_OFFSET_POLYGONS_FROM_POLYGON_WITH_HOLES_2_H
