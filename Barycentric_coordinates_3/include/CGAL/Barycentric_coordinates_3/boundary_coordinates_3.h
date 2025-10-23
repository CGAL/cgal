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
  with respect to the vertices of a convex simplicial polyhedron, that is one
  coordinate per vertex. The coordinates are stored in a destination range
  beginning at `oi`.

  If `query` is at the vertex, the corresponding coordinate is set to one, while
  all other coordinates are zero. If `query` is on the face, the three corresponding
  coordinates are triangle coordinates, while all other coordinates are set to zero.
  If `query` is not on the boundary, all the coordinates are set to zero.

  \tparam TriangleMesh
  must be a model of the concept `FaceListGraph`.

  \tparam GeomTraits
  a model of `BarycentricTraits_3`

  \tparam OutputIterator
  a model of `OutputIterator` that accepts values of type `GeomTraits::FT`

  \tparam VertexPointMap
  a property map with boost::graph_traits<TriangleMesh>::vertex_descriptor as
  key type and `GeomTraits::Point_3` as value type

  \param tmesh
  an instance of `TriangleMesh`

  \param query
  a query point

  \param oi
  the beginning of the destination range with the computed coordinates

  \param traits
  a traits class with geometric objects, predicates, and constructions;
  the default initialization is provided

  \param vertex_point_map
  an instance of `VertexPointMap` that maps a vertex from `tmesh` to `Point_3`;
  the default initialization is provided

  \return an output iterator to the element in the destination range,
  one past the last coordinate stored + the flag indicating whether the
  query point belongs to the polyhedron boundary

  \pre num_vertices(`tmesh`) >= 4.
  \pre is_triangle_mesh(`tmesh`).
  \pre is_closed(`tmesh`).
  \pre is_strongly_convex_3(`tmesh`).
*/
template<typename TriangleMesh,
         typename OutputIterator,
         typename GeomTraits,
         typename VertexPointMap>
std::pair<OutputIterator, bool>
boundary_coordinates_3(const TriangleMesh& tmesh,
                       const typename GeomTraits::Point_3& query,
                       OutputIterator oi,
                       const GeomTraits& traits,
                       const VertexPointMap vertex_point_map)
{
  const auto edge_case = internal::locate_wrt_polyhedron(
    vertex_point_map, tmesh, query, oi, traits);

  if(edge_case == internal::Edge_case::BOUNDARY)
    return {oi, true};
  else{
    internal::get_default(num_vertices(tmesh), oi);
    return {oi, false};
  }
}

/*!
  \ingroup PkgBarycentricCoordinates3RefFunctions

  \brief computes boundary barycentric coordinates.

  This function computes boundary barycentric coordinates at a given `query` point
  with respect to the vertices of a simple `polyhedron`, that is one
  coordinate per vertex. The coordinates are stored in a destination range
  beginning at `oi`.

  If `query` is at the vertex, the corresponding coordinate is set to one, while
  all other coordinates are zero. If `query` is on the face, the three corresponding
  coordinates are triangle coordinates, while all other coordinates are set to zero.
  If `query` is not on the boundary, all the coordinates are set to zero.

  \tparam TriangleMesh
  must be a model of the concept `FaceListGraph`.

  \tparam Point_3
  a model of `Kernel::Point_3`

  \tparam OutputIterator
  a model of `OutputIterator` that accepts values of type `GeomTraits::FT`

  \tparam VertexPointMap
  a property map with boost::graph_traits<TriangleMesh>::vertex_descriptor as
  key type and Point_3 as value type. The default is `property_map_selector<TriangleMesh,
  CGAL::vertex_point_t>`.

  \param tmesh
  an instance of `TriangleMesh`, which must be a convex simplicial polyhedron

  \param query
  a query point

  \param oi
  the beginning of the destination range with the computed coordinates

  \param vertex_point_map
  an instance of `VertexPointMap` that maps a vertex from `tmesh` to `Point_3`;
  the default initialization is provided

  \return an output iterator to the element in the destination range,
  one past the last coordinate stored + the flag indicating whether the
  query point belongs to the polyhedron boundary

  \pre is_triangle_mesh(`tmesh`)
  \pre num_vertices(`tmesh`) >= 4.
  \pre `tmesh` is simplicial.
*/
template<typename TriangleMesh,
         typename Point_3,
         typename OutputIterator,
         typename VertexPointMap = typename boost::property_map<TriangleMesh, CGAL::vertex_point_t>::const_type>
std::pair<OutputIterator, bool>
boundary_coordinates_3(const TriangleMesh& tmesh,
                       const Point_3& query,
                       OutputIterator oi,
                       const VertexPointMap vertex_point_map)
{
  using GeomTraits = typename Kernel_traits<Point_3>::Kernel;
  const GeomTraits traits;

  return boundary_coordinates_3(tmesh, query, oi, traits, vertex_point_map);
}

template<typename TriangleMesh,
         typename Point_3,
         typename OutputIterator>
std::pair<OutputIterator, bool>
boundary_coordinates_3(const TriangleMesh& tmesh,
                       const Point_3& query,
                       OutputIterator oi)
{
  return boundary_coordinates_3(tmesh, query, oi,
   get_const_property_map(CGAL::vertex_point, tmesh));
}

}
}

#endif
