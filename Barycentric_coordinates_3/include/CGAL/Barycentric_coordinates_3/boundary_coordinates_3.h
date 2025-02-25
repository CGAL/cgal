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

#ifndef CGAL_BARYCENTRIC_BOUNDARY_COORDINATES_3_H
#define CGAL_BARYCENTRIC_BOUNDARY_COORDINATES_3_H

#include <CGAL/license/Barycentric_coordinates_3.h>

#include <CGAL/Barycentric_coordinates_3/internal/utils_3.h>
#include <CGAL/Barycentric_coordinates_3/barycentric_enum_3.h>
#include <CGAL/boost/graph/property_maps.h>

namespace CGAL{
namespace Barycentric_coordinates{

  /*!
    \ingroup PkgBarycentricCoordinates3RefFunctions

    \brief computes boundary barycentric coordinates.

    This function computes boundary barycentric coordinates at a given `query` point
    with respect to the vertices of a simple `polyhedron`, that is one
    coordinate per vertex. The coordinates are stored in a destination range
    beginning at `c_begin`.

    If `query` is at the vertex, the corresponding coordinate is set to one, while
    all other coordinates are zero. If `query` is on the face, the three corresponding
    coordinates are triangle coordinates, while all other coordinates are set to zero.
    If `query` is not on the boundary, all the coordinates are set to zero.

    \tparam TriangleMesh
    must be a model of the concept `FaceListGraph`.

    \tparam GeomTraits
    a model of `BarycentricTraits_3`

    \tparam OutIterator
    a model of `OutputIterator` that accepts values of type `GeomTraits::FT`

    \tparam VertexToPointMap
    a property map with boost::graph_traits<TriangleMesh>::vertex_descriptor as
    key type and `GeomTraits::Point_3` as value type.

    \param triangle_mesh
    an instance of `TriangleMesh`, which must be a convex simplicial polyhedron

    \param query
    a query point

    \param c_begin
    the beginning of the destination range with the computed coordinates

    \param traits
    a traits class with geometric objects, predicates, and constructions;
    the default initialization is provided

    \param vertex_to_point_map
    an instance of `VertexToPointMap` that maps a vertex from `triangle_mesh` to `Point_3`;
    the default initialization is provided

    \return an output iterator to the element in the destination range,
    one past the last coordinate stored + the flag indicating whether the
    query point belongs to the polyhedron boundary

    \pre num_vertices(triangle_mesh) >= 4.
    \pre triangle_mesh is simplicial.
  */
  template<
  typename TriangleMesh,
  typename OutIterator,
  typename GeomTraits,
  typename VertexToPointMap>
  std::pair<OutIterator, bool> boundary_coordinates_3(
    const TriangleMesh& triangle_mesh,
    const typename GeomTraits::Point_3& query,
    OutIterator c_begin,
    const GeomTraits& traits,
    const VertexToPointMap vertex_to_point_map) {

    const auto edge_case = internal::locate_wrt_polyhedron(
      vertex_to_point_map, triangle_mesh, query, c_begin, traits);

    if(edge_case == internal::Edge_case::BOUNDARY)
      return {c_begin, true};
    else{
      internal::get_default(num_vertices(triangle_mesh), c_begin);
      return {c_begin, false};
    }
  }

  /*!
    \ingroup PkgBarycentricCoordinates3RefFunctions

    \brief computes boundary barycentric coordinates.

    This function computes boundary barycentric coordinates at a given `query` point
    with respect to the vertices of a simple `polyhedron`, that is one
    coordinate per vertex. The coordinates are stored in a destination range
    beginning at `c_begin`.

    If `query` is at the vertex, the corresponding coordinate is set to one, while
    all other coordinates are zero. If `query` is on the face, the three corresponding
    coordinates are triangle coordinates, while all other coordinates are set to zero.
    If `query` is not on the boundary, all the coordinates are set to zero.

    \tparam TriangleMesh
    must be a model of the concept `FaceListGraph`.

    \tparam Point_3
    a model of `Kernel::Point_3`

    \tparam OutIterator
    a model of `OutputIterator` that accepts values of type `GeomTraits::FT`

    \tparam VertexToPointMap
    a property map with boost::graph_traits<TriangleMesh>::vertex_descriptor as
    key type and Point_3 as value type. The default is `property_map_selector<TriangleMesh,
    CGAL::vertex_point_t>`.

    \param triangle_mesh
    an instance of `TriangleMesh`, which must be a convex simplicial polyhedron

    \param query
    a query point

    \param c_begin
    the beginning of the destination range with the computed coordinates

    \param vertex_to_point_map
    an instance of `VertexToPointMap` that maps a vertex from `triangle_mesh` to `Point_3`;
    the default initialization is provided

    \return an output iterator to the element in the destination range,
    one past the last coordinate stored + the flag indicating whether the
    query point belongs to the polyhedron boundary

    \pre num_vertices(triangle_mesh) >= 4.
    \pre triangle_mesh is simplicial.
  */
  template<
  typename TriangleMesh,
  typename Point_3,
  typename OutIterator,
  typename VertexToPointMap = typename boost::property_map<TriangleMesh, CGAL::vertex_point_t>::const_type>
  std::pair<OutIterator, bool> boundary_coordinates_3(
    const TriangleMesh& triangle_mesh,
    const Point_3& query,
    OutIterator c_begin,
    const VertexToPointMap vertex_to_point_map){

    using GeomTraits = typename Kernel_traits<Point_3>::Kernel;
    const GeomTraits traits;

    return boundary_coordinates_3(triangle_mesh, query, c_begin, traits, vertex_to_point_map);
  }

  template<
  typename TriangleMesh,
  typename Point_3,
  typename OutIterator>
  std::pair<OutIterator, bool> boundary_coordinates_3(
    const TriangleMesh& triangle_mesh,
    const Point_3& query,
    OutIterator c_begin){

    return boundary_coordinates_3(triangle_mesh, query, c_begin,
     get_const_property_map(CGAL::vertex_point, triangle_mesh));
  }

}
}

#endif
