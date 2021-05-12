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

#include <CGAL/arrange_offset_polygons_2.h>
#include <CGAL/create_offset_polygons_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>

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
template<class FT, class PolygonWithHoles, class OfK, class SsK,
         class OutPolygon = typename CGAL_SS_i::Default_return_polygon_type<PolygonWithHoles, OfK>::type> // Hole-less polygon type
std::vector<boost::shared_ptr<OutPolygon> >
inline
create_interior_skeleton_and_offset_polygons_2(const FT& aOffset,
                                               const PolygonWithHoles& aPoly,
                                               const OfK& ofk,
                                               const SsK& ssk,
                                               typename std::enable_if<
                                                 CGAL_SS_i::has_Hole_const_iterator<PolygonWithHoles>::value>::type* = nullptr)
{
  return create_interior_skeleton_and_offset_polygons_2(aOffset, aPoly.outer_boundary(),
                                                        aPoly.holes_begin(), aPoly.holes_end(),
                                                        ofk, ssk);
}

/*! create_interior_skeleton_and_offset_polygons_with_holes_2 (orders the resulting polygons) */

// Polygon might be a Polygon with holes or not, but it returns a Polygon with holes
template<class FT, class Polygon, class OfK, class SsK,
         class OutPolygonWithHoles = typename CGAL_SS_i::Default_return_polygon_with_holes_type<Polygon, OfK>::type>
std::vector<boost::shared_ptr<OutPolygonWithHoles> >
inline
create_interior_skeleton_and_offset_polygons_with_holes_2(const FT& aOffset,
                                                          const Polygon& aPoly,
                                                          const OfK& ofk,
                                                          const SsK& ssk)
{
  return arrange_offset_polygons_2<OutPolygonWithHoles>(
           create_interior_skeleton_and_offset_polygons_2(aOffset, aPoly, ofk, ssk));
}

template<class FT, class Polygon, class OfK,
         class OutPolygonWithHoles = typename CGAL_SS_i::Default_return_polygon_with_holes_type<Polygon, OfK>::type>
std::vector<boost::shared_ptr<OutPolygonWithHoles> >
inline
create_interior_skeleton_and_offset_polygons_with_holes_2(const FT& aOffset,
                                                          const Polygon& aPoly,
                                                          const OfK& ofk)
{
  return create_interior_skeleton_and_offset_polygons_with_holes_2(aOffset, aPoly, ofk,
                                                                   Exact_predicates_inexact_constructions_kernel());
}

template<class FT, class Polygon,
         class OutPolygonWithHoles = typename CGAL_SS_i::Default_return_polygon_with_holes_type<
                                       Polygon, Exact_predicates_inexact_constructions_kernel>::type>
std::vector<boost::shared_ptr<OutPolygonWithHoles> >
inline
create_interior_skeleton_and_offset_polygons_with_holes_2(const FT& aOffset,
                                                          const Polygon& aPoly)
{
  return create_interior_skeleton_and_offset_polygons_with_holes_2(aOffset, aPoly,
                                                                   Exact_predicates_inexact_constructions_kernel());
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/// EXTERIOR

/*! create_interior_skeleton_and_offset_polygons_2 with a polygon with holes */

// overload where PolygonWithHoles actually is a type of Polygon that supports holes
template<class FT, class PolygonWithHoles, class OfK, class SsK,
         class OutPolygon = typename CGAL_SS_i::Default_return_polygon_type<PolygonWithHoles, OfK>::type>
std::vector<boost::shared_ptr<OutPolygon> >
inline
create_exterior_skeleton_and_offset_polygons_2(const FT& aOffset,
                                               const PolygonWithHoles& aPoly,
                                               const OfK& ofk,
                                               const SsK& ssk,
                                               typename std::enable_if<
                                                 CGAL_SS_i::has_Hole_const_iterator<PolygonWithHoles>::value>::type* = nullptr)
{
  return create_exterior_skeleton_and_offset_polygons_2(aOffset, aPoly.outer_boundary(), ofk, ssk);
}

/*! create_exterior_skeleton_and_offset_polygons_with_holes_2 (orders the resulting polygons) */

// Polygon might be a Polygon with holes or not, but it returns a Polygon with holes
template<class FT, class Polygon, class OfK, class SsK,
         class OutPolygonWithHoles = typename CGAL_SS_i::Default_return_polygon_with_holes_type<Polygon, OfK>::type>
std::vector<boost::shared_ptr<OutPolygonWithHoles> >
inline
create_exterior_skeleton_and_offset_polygons_with_holes_2(const FT& aOffset,
                                                          const Polygon& aPoly,
                                                          const OfK& ofk,
                                                          const SsK& ssk)
{
  return arrange_offset_polygons_2<OutPolygonWithHoles>(
           create_exterior_skeleton_and_offset_polygons_2(aOffset, aPoly, ofk, ssk));
}

template<class FT, class Polygon, class OfK,
         class OutPolygonWithHoles = typename CGAL_SS_i::Default_return_polygon_with_holes_type<Polygon, OfK>::type>
std::vector<boost::shared_ptr<OutPolygonWithHoles> >
inline
create_exterior_skeleton_and_offset_polygons_with_holes_2(const FT& aOffset,
                                                          const Polygon& aPoly,
                                                          const OfK& ofk)
{
  return create_exterior_skeleton_and_offset_polygons_with_holes_2(aOffset, aPoly, ofk,
                                                                   Exact_predicates_inexact_constructions_kernel());
}

template<class FT, class Polygon,
         class OutPolygonWithHoles = typename CGAL_SS_i::Default_return_polygon_with_holes_type<
                                    Polygon, Exact_predicates_inexact_constructions_kernel>::type>
std::vector<boost::shared_ptr<OutPolygonWithHoles> >
inline
create_exterior_skeleton_and_offset_polygons_with_holes_2(const FT& aOffset,
                                                          const Polygon& aPoly)
{
  return create_exterior_skeleton_and_offset_polygons_with_holes_2(aOffset, aPoly,
                                                                   Exact_predicates_inexact_constructions_kernel());
}

} // end namespace CGAL

#endif
