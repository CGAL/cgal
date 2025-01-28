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
#ifndef CGAL_CREATE_WEIGHTED_OFFSET_POLYGONS_FROM_POLYGON_WITH_HOLES_2_H
#define CGAL_CREATE_WEIGHTED_OFFSET_POLYGONS_FROM_POLYGON_WITH_HOLES_2_H

#include <CGAL/license/Straight_skeleton_2.h>

#include <CGAL/arrange_offset_polygons_2.h>
#include <CGAL/create_offset_polygons_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <type_traits>
#include <vector>
#include <memory>

namespace CGAL {

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/// INTERIOR

/*! create_interior_skeleton_and_offset_polygons_2 (no sorting of the result) */

// overload where PolygonWithHoles actually is a type of Polygon that supports holes
template<class FT, class PolygonWithHoles, class Weights, class OfK, class SsK,
         class OutPolygon = typename CGAL_SS_i::Default_return_polygon_type<PolygonWithHoles, OfK>::type> // Hole-less polygon type
std::vector<std::shared_ptr<OutPolygon> >
inline
create_interior_weighted_skeleton_and_offset_polygons_2(const FT& aOffset,
                                                        const PolygonWithHoles& aPoly,
                                                        const Weights& aWeights,
                                                        const OfK& ofk,
                                                        const SsK& ssk,
                                                        std::enable_if_t<
                                                          CGAL_SS_i::has_Hole_const_iterator<PolygonWithHoles>::value>* = nullptr)
{
  return create_interior_weighted_skeleton_and_offset_polygons_2(aOffset, aPoly.outer_boundary(),
                                                                 aPoly.holes_begin(), aPoly.holes_end(),
                                                                 aWeights[0],
                                                                 std::next(std::begin(aWeights)),
                                                                 std::end(aWeights),
                                                                 ofk, ssk);
}

/*! create_interior_weighted_skeleton_and_offset_polygons_with_holes_2 (orders the resulting polygons) */

// Polygon might be a Polygon with holes or not, but it returns a Polygon with holes
template<class FT, class Polygon, class Weights, class OfK, class SsK,
         class OutPolygonWithHoles = typename CGAL_SS_i::Default_return_polygon_with_holes_type<Polygon, OfK>::type>
std::vector<std::shared_ptr<OutPolygonWithHoles> >
inline
create_interior_weighted_skeleton_and_offset_polygons_with_holes_2(const FT& aOffset,
                                                                   const Polygon& aPoly,
                                                                   const Weights& aWeights,
                                                                   const OfK& ofk,
                                                                   const SsK& ssk)
{
  return arrange_offset_polygons_2<OutPolygonWithHoles>(
           create_interior_weighted_skeleton_and_offset_polygons_2(aOffset, aPoly, aWeights, ofk, ssk));
}

template<class FT, class Polygon, class Weights, class OfK,
         class OutPolygonWithHoles = typename CGAL_SS_i::Default_return_polygon_with_holes_type<Polygon, OfK>::type>
std::vector<std::shared_ptr<OutPolygonWithHoles> >
inline
create_interior_weighted_skeleton_and_offset_polygons_with_holes_2(const FT& aOffset,
                                                                   const Polygon& aPoly,
                                                                   const Weights& aWeights,
                                                                   const OfK& ofk)
{
  return create_interior_weighted_skeleton_and_offset_polygons_with_holes_2(aOffset, aPoly, aWeights, ofk,
                                                                            Exact_predicates_inexact_constructions_kernel());
}

template<class FT, class Polygon, class Weights,
         class OutPolygonWithHoles = typename CGAL_SS_i::Default_return_polygon_with_holes_type<
                                       Polygon, Exact_predicates_inexact_constructions_kernel>::type>
std::vector<std::shared_ptr<OutPolygonWithHoles> >
inline
create_interior_weighted_skeleton_and_offset_polygons_with_holes_2(const FT& aOffset,
                                                                   const Polygon& aPoly,
                                                                   const Weights& aWeights)
{
  return create_interior_weighted_skeleton_and_offset_polygons_with_holes_2(aOffset, aPoly, aWeights,
                                                                            Exact_predicates_inexact_constructions_kernel());
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/// EXTERIOR

/*! create_exterior_skeleton_and_offset_polygons_with_holes_2 (orders the resulting polygons) */

// Polygon might be a Polygon with holes or not, but it returns a Polygon with holes
template<class FT, class Polygon, class Weights, class OfK, class SsK,
         class OutPolygonWithHoles = typename CGAL_SS_i::Default_return_polygon_with_holes_type<Polygon, OfK>::type>
std::vector<std::shared_ptr<OutPolygonWithHoles> >
create_exterior_weighted_skeleton_and_offset_polygons_with_holes_2(const FT& aOffset,
                                                                   const Polygon& aPoly,
                                                                   const Weights& aWeights,
                                                                   const OfK& ofk,
                                                                   const SsK& ssk)
{
  typedef typename CGAL_SS_i::Default_return_polygon_type<Polygon, OfK>::type Polygon_;
  std::vector<std::shared_ptr<Polygon_> > raw_output =
    create_exterior_weighted_skeleton_and_offset_polygons_2(aOffset, aPoly, aWeights, ofk, ssk);

  // filter offset of the outer frame
  std::swap(raw_output[0], raw_output.back());
  raw_output.pop_back();

  for(std::shared_ptr<Polygon_> ptr : raw_output)
    ptr->reverse_orientation();

  return arrange_offset_polygons_2<OutPolygonWithHoles>(raw_output);
}

/*! create_interior_skeleton_and_offset_polygons_2 with a polygon with holes */

// overload where PolygonWithHoles actually is a type of Polygon that supports holes
template<class FT, class PolygonWithHoles, class Weights, class OfK, class SsK,
         class OutPolygon = typename CGAL_SS_i::Default_return_polygon_type<PolygonWithHoles, OfK>::type>
std::vector<std::shared_ptr<OutPolygon> >
inline
create_exterior_weighted_skeleton_and_offset_polygons_2(const FT& aOffset,
                                                        const PolygonWithHoles& aPoly,
                                                        const Weights& aWeights,
                                                        const OfK& ofk,
                                                        const SsK& ssk,
                                                        std::enable_if_t<
                                                          CGAL_SS_i::has_Hole_const_iterator<PolygonWithHoles>::value>* = nullptr)
{
  CGAL_precondition(aWeights.size() == aPoly.number_of_holes() + 1);

  std::vector<std::shared_ptr<OutPolygon> > polygons =
    create_exterior_weighted_skeleton_and_offset_polygons_2(aOffset, aPoly.outer_boundary(), {aWeights[0]}, ofk, ssk);

  std::size_t weight_pos = 1;
  for(typename PolygonWithHoles::Hole_const_iterator hit=aPoly.holes_begin(); hit!=aPoly.holes_end(); ++hit, ++weight_pos)
  {
    typename PolygonWithHoles::Polygon_2 hole = *hit;
    hole.reverse_orientation();
    std::vector<std::shared_ptr<OutPolygon> > hole_polygons =
        create_interior_skeleton_and_offset_polygons_2(aOffset,
                                                       hole,
                                                       {aWeights[weight_pos]},
                                                       ofk, ssk);
    polygons.insert(polygons.end(), hole_polygons.begin(), hole_polygons.end());
  }

  return polygons;
}

template<class FT, class Polygon, class Weights, class OfK,
         class OutPolygonWithHoles = typename CGAL_SS_i::Default_return_polygon_with_holes_type<Polygon, OfK>::type>
std::vector<std::shared_ptr<OutPolygonWithHoles> >
inline
create_exterior_weighted_skeleton_and_offset_polygons_with_holes_2(const FT& aOffset,
                                                                   const Polygon& aPoly,
                                                                   const Weights& aWeights,
                                                                   const OfK& ofk)
{
  return create_exterior_weighted_skeleton_and_offset_polygons_with_holes_2(aOffset, aPoly, aWeights, ofk,
                                                                            Exact_predicates_inexact_constructions_kernel());
}

template<class FT, class Polygon, class Weights,
         class OutPolygonWithHoles = typename CGAL_SS_i::Default_return_polygon_with_holes_type<
                                    Polygon, Exact_predicates_inexact_constructions_kernel>::type>
std::vector<std::shared_ptr<OutPolygonWithHoles> >
inline
create_exterior_weighted_skeleton_and_offset_polygons_with_holes_2(const FT& aOffset,
                                                                   const Polygon& aPoly,
                                                                   const Weights& aWeights)
{
  return create_exterior_weighted_skeleton_and_offset_polygons_with_holes_2(aOffset, aPoly, aWeights,
                                                                            Exact_predicates_inexact_constructions_kernel());
}

} // namespace CGAL

#endif // CGAL_CREATE_WEIGHTED_OFFSET_POLYGONS_FROM_POLYGON_WITH_HOLES_2_H
