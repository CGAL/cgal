// Copyright (c) 2021 GeometryFactory SARL (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Antonio Gomes, Dmitry Anisimov
//

#ifndef CGAL_BARYCENTRIC_COORDINATES_3_H
#define CGAL_BARYCENTRIC_COORDINATES_3_H

#include <CGAL/license/Barycentric_coordinates_3.h>

/**
* \ingroup PkgBarycentricCoordinates3Ref
* \file CGAL/Barycentric_coordinates_3.h
* A convenience header that includes all free functions and classes for
* 3D barycentric coordinates in closed form.
*/

#include <array>
#include <CGAL/Kernel_traits.h>
#include <CGAL/property_map.h>
#include <CGAL/Barycentric_coordinates_3/tetrahedron_coordinates.h>
#include <CGAL/Barycentric_coordinates_3/boundary_coordinates_3.h>
#include <CGAL/Barycentric_coordinates_3/Wachspress_coordinates_3.h>
#include <CGAL/Barycentric_coordinates_3/Discrete_harmonic_coordinates_3.h>
#include <CGAL/Barycentric_coordinates_3/Mean_value_coordinates_3.h>

namespace CGAL {
namespace Barycentric_coordinates {

/*!
  \ingroup PkgBarycentricCoordinates3RefFunctions

  \brief computes a point location from barycentric coordinates with respect to a triangle mesh.

  This function computes a point location from barycentric coordinates
  with respect to the vertices of `tmesh`

  \tparam TriangleMesh
  must be a model of the concept `FaceListGraph`

  \tparam CoordinateRange
  a range whose iterator is a model of `ForwardIterator` with value type `GeomTraits::FT`

  \tparam VertexPointMap
  a property map with `boost::graph_traits<TriangleMesh>::vertex_descriptor` as
  key type and `GeomTraits::Point_3` as value type.

  \tparam GeomTraits
  a model of `BarycentricTraits_3`, deduced from the value type of `VertexPointMap`, using `CGAL::Kernel_traits`.

  \param tmesh
  an instance of `TriangleMesh`

  \param coordinates
  barycentric coordinates of the point

  \param geom_traits
  an instance of `GeomTraits`

  \param vpm
  an instance of `VertexPointMap` that maps a vertex from `tmesh` to `GeomTraits::Point_3`

  \return point with type `boost::property_traits<VertexPointMap>::value_type`

  \pre `vertices(tmesh).size()` == `coords.size()`
*/
template<typename TriangleMesh, typename CoordinateRange, typename VertexPointMap, typename GeomTraits = typename CGAL::Kernel_traits<typename boost::property_traits<VertexPointMap>::value_type>::type>
typename boost::property_traits<VertexPointMap>::value_type apply_barycentric_coordinates(const TriangleMesh& tmesh, const CoordinateRange& coordinates, VertexPointMap vpm, GeomTraits geom_traits = GeomTraits()) {
  CGAL_precondition(vertices(tmesh).size() == coordinates.size());
  using Point = typename boost::property_traits<VertexPointMap>::value_type;
  static_assert(std::is_same_v<GeomTraits, typename Kernel_traits<Point>::Kernel>);
  using vertex_descriptor = typename boost::graph_traits<TriangleMesh>::vertex_descriptor;
  using FT = typename GeomTraits::FT;
  using Construct_cartesian_const_iterator = typename GeomTraits::Construct_cartesian_const_iterator_3;
  using Cartesian_const_iterator = typename GeomTraits::Cartesian_const_iterator_3;
  using Construct_point_3 = typename GeomTraits::Construct_point_3;

  Construct_cartesian_const_iterator construct_cci = geom_traits.construct_cartesian_const_iterator_3_object();

  std::array<FT, 3> p = { FT(0), FT(0), FT(0) };
  for (vertex_descriptor v : vertices(tmesh)) {
    const Point& pv = get(vpm, v);
    Cartesian_const_iterator cci = construct_cci(pv);
    p[0] += *cci * (*(std::begin(coordinates)));
    ++cci;
    p[1] += *cci * (*(std::begin(coordinates) + 1));
    ++cci;
    p[2] += *cci * (*(std::begin(coordinates) + 2));
    ++cci;
  }

  Construct_point_3 construct_point_3 = geom_traits.construct_point_3_object();
  return construct_point_3(p[0], p[1], p[2]);
}


/*!
  \ingroup PkgBarycentricCoordinates3RefFunctions

  \brief computes a point location from barycentric coordinates with respect to a set of points.

  This function a computes point location from barycentric coordinates
  with respect to a set of points.

  \tparam PointRange
  is a model of `ConstRange` and `RandomAccessRange` with value type `GeomTraits::Point_3`

  \tparam CoordinateRange
  a range whose iterator is a model of `ForwardIterator` with value type `GeomTraits::FT`

  \tparam GeomTraits
  a model of `BarycentricTraits_3`, deduced from the value type of `PointRange` by default.

  \param points
  a range of input points

  \param coordinates
  barycentric coordinates of the point

  \param geom_traits
  an instance of `GeomTraits`

  \return point with type `GeomTraits::Point_3`

  \pre `pts.size()` == `coords.size()`
*/
template<typename PointRange, typename CoordinateRange, typename GeomTraits = typename Kernel_traits<typename boost::range_value<PointRange>::type>::Kernel>
typename boost::range_value<PointRange>::type apply_barycentric_coordinates(const PointRange& points, const CoordinateRange& coordinates, GeomTraits geom_traits = GeomTraits()) {
  CGAL_precondition(std::distance(std::begin(points), std::end(points)) == std::distance(std::begin(coordinates), std::end(coordinates)));
  using Point = typename boost::range_value<PointRange>::type;
  static_assert(std::is_same_v<GeomTraits, typename Kernel_traits<Point>::Kernel>);
  using FT = typename GeomTraits::FT;
  using Construct_cartesian_const_iterator = typename GeomTraits::Construct_cartesian_const_iterator_3;
  using Cartesian_const_iterator = typename GeomTraits::Cartesian_const_iterator_3;
  using Construct_point_3 = typename GeomTraits::Construct_point_3;

  Construct_cartesian_const_iterator construct_cci = geom_traits.construct_cartesian_const_iterator_3_object();

  std::array<FT, 3> p = { FT(0), FT(0), FT(0) };
  for (const Point& pv : points) {
    Cartesian_const_iterator cci = construct_cci(pv);
    p[0] += *cci * (*(std::begin(coordinates)));
    ++cci;
    p[1] += *cci * (*(std::begin(coordinates) + 1));
    ++cci;
    p[2] += *cci * (*(std::begin(coordinates) + 2));
    ++cci;
  }

  Construct_point_3 construct_point_3 = geom_traits.construct_point_3_object();
  return construct_point_3(p[0], p[1], p[2]);
}
} // namespace Barycentric_coordinates
} // namespace CGAL

#endif
