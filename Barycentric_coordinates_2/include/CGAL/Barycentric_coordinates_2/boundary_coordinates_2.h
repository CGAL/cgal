// Copyright (c) 2014 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Dmitry Anisimov, David Bommes, Kai Hormann, Pierre Alliez
//

#ifndef CGAL_BARYCENTRIC_BOUNDARY_COORDINATES_2_H
#define CGAL_BARYCENTRIC_BOUNDARY_COORDINATES_2_H

#include <CGAL/license/Barycentric_coordinates_2.h>

// Internal includes.
#include <CGAL/Barycentric_coordinates_2/internal/utils_2.h>

namespace CGAL {
namespace Barycentric_coordinates {

  /*!
    \ingroup PkgBarycentricCoordinates2RefFunctions

    \brief computes 2D boundary coordinates.

    This function computes boundary barycentric coordinates at a given `query` point
    with respect to the vertices of a simple `polygon`, that is one
    coordinate per vertex. The coordinates are stored in a destination range
    beginning at `c_begin`.

    If `query` is at the vertex, the corresponding coordinate is set to one, while
    all other coordinates are zero. If `query` is on the edge, the two corresponding
    coordinates are segment coordinates, while all other coordinates are set to zero.
    If `query` is not on the boundary, all the coordinates are set to zero.

    Internally, `segment_coordinates_2()` are used.

    \tparam VertexRange
    a model of `ConstRange` whose iterator type is `RandomAccessIterator`

    \tparam OutIterator
    a model of `OutputIterator` that accepts values of type `GeomTraits::FT`

    \tparam GeomTraits
    a model of `BarycentricTraits_2`

    \tparam PointMap
    a model of `ReadablePropertyMap` whose key type is `VertexRange::value_type` and
    value type is `GeomTraits::Point_2`

    \param polygon
    an instance of `VertexRange` with 2D points, which form a simple polygon

    \param query
    a query point

    \param c_begin
    the beginning of the destination range with the computed coordinates

    \param traits
    a traits class with geometric objects, predicates, and constructions;
    this parameter can be omitted if the traits class can be deduced from the point type

    \param point_map
    an instance of `PointMap` that maps a vertex from `polygon` to `Point_2`

    \return an output iterator to the element in the destination range,
    one past the last coordinate stored + the flag indicating whether the
    query point belongs to the polygon boundary

    \pre polygon.size() >= 3
  */
  template<
  typename VertexRange,
  typename OutIterator,
  typename GeomTraits,
  typename PointMap>
  std::pair<OutIterator, bool> boundary_coordinates_2(
    const VertexRange& polygon,
    const typename GeomTraits::Point_2& query,
    OutIterator c_begin,
    const GeomTraits& traits,
    const PointMap point_map) {

    const auto result =
    internal::locate_wrt_polygon_2(polygon, query, traits, point_map);
    auto location = (*result).first;
    auto index = (*result).second;

    if (!result) {
      index = std::size_t(-1);
      location = internal::Query_point_location::UNSPECIFIED;
    }
    return internal::boundary_coordinates_2(
      polygon, query, location, index, c_begin, traits, point_map);
  }

  /*!
    \ingroup PkgBarycentricCoordinates2RefFunctions

    \brief computes 2D boundary coordinates.

    This is an overload of the function `boundary_coordinates_2`.

    \tparam VertexRange
    a model of `ConstRange` whose iterator type is `RandomAccessIterator`

    \tparam Query_2
    a model of `Kernel::Point_2`

    \tparam OutIterator
    a model of `OutputIterator` that accepts values of type `GeomTraits::FT`

    \tparam PointMap
    a model of `ReadablePropertyMap` whose key type is `VertexRange::value_type` and
    value type is `Query_2`. The default is `CGAL::Identity_property_map`.

    \param polygon
    an instance of `VertexRange` with 2D points, which form a simple polygon

    \param query
    a query point

    \param c_begin
    the beginning of the destination range with the computed coordinates

    \param point_map
    an instance of `PointMap` that maps a vertex from `polygon` to `Query_2`;
    the default initialization is provided

    \return an output iterator to the element in the destination range,
    one past the last coordinate stored + the flag indicating whether the
    query point belongs to the polygon boundary

    \pre polygon.size() >= 3
  */
  template<
  typename VertexRange,
  typename Query_2,
  typename OutIterator,
  typename PointMap = CGAL::Identity_property_map<Query_2> >
  std::pair<OutIterator, bool> boundary_coordinates_2(
    const VertexRange& polygon,
    const Query_2& query,
    OutIterator c_begin,
    const PointMap point_map = PointMap()) {

    using GeomTraits = typename Kernel_traits<Query_2>::Kernel;
    const GeomTraits traits;
    return boundary_coordinates_2(
      polygon, query, c_begin, traits, point_map);
  }

} // namespace Barycentric_coordinates
} // namespace CGAL

#endif // CGAL_BARYCENTRIC_BOUNDARY_COORDINATES_2_H
