// Copyright (c) 2020 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_POLYGON_MESH_PROCESSING_POLYGON_MESH_TO_POLYGON_SOUP_H
#define CGAL_POLYGON_MESH_PROCESSING_POLYGON_MESH_TO_POLYGON_SOUP_H

#include <CGAL/license/Polygon_mesh_processing/repair.h>

#include <CGAL/algorithm.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/Container_helper.h>
#include <CGAL/Dynamic_property_map.h>
#include <CGAL/property_map.h>

#include <boost/range/value_type.hpp>
#include <boost/range/reference.hpp>

namespace CGAL {
namespace Polygon_mesh_processing {

/// \ingroup PMP_repairing_grp
///
/// adds the vertices and faces of a mesh into a (possibly non-empty) polygon soup.
///
/// \tparam PolygonMesh a model of `FaceListGraph`
/// \tparam PointRange a model of the concepts `RandomAccessContainer` and
///                    `BackInsertionSequence` whose value type can be constructed from
///                    the point type of the polygon mesh
/// \tparam PolygonRange a model of the concepts `RandomAccessContainer` and `BackInsertionSequence` whose
///                      value type is itself a model of the concepts `RandomAccessContainer` and
///                      `BackInsertionSequence` whose value type is `std::size_t`
/// \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
///
/// \param mesh the mesh whose faces are being put in the polygon soup
/// \param points points making the polygons of the soup
/// \param polygons each element in the vector describes a polygon using the indices of the points in `points`
/// \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
///
/// \cgalNamedParamsBegin
///   \cgalParamNBegin{vertex_point_map}
///     \cgalParamDescription{a property map associating points to the vertices of `mesh`}
///     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
///                    as key type and `%Point_3` as value type}
///     \cgalParamDefault{`boost::get(CGAL::vertex_point, mesh)`}
///   \cgalParamNEnd
/// \cgalNamedParamsEnd
///
/// \cgalAdvancedBegin
/// `PolygonRange` can also be a model of the concepts `RandomAccessContainer` and `BackInsertionSequence`
/// whose value type is an array, but it is the user's responsability to ensure that
/// all faces have the same number of vertices, and that this number is equal to the size of the array.
/// \cgalAdvancedEnd
///
/// \sa `CGAL::Polygon_mesh_processing::orient_polygon_soup()`
/// \sa `CGAL::Polygon_mesh_processing::is_polygon_soup_a_polygon_mesh()`
/// \sa `CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh()`
///
template<typename PolygonMesh,
         typename PointRange, typename PolygonRange,
         typename NamedParameters>
void polygon_mesh_to_polygon_soup(const PolygonMesh& mesh,
                                  PointRange& points,
                                  PolygonRange& polygons,
                                  const NamedParameters& np)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor              vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor            halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor                face_descriptor;

  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::const_type      VPM;
  VPM vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_const_property_map(vertex_point, mesh));

  typedef CGAL::dynamic_vertex_property_t<std::size_t>                              Vertex_index;
  typedef typename boost::property_map<PolygonMesh, Vertex_index>::const_type       VIM;
  VIM vim = get(Vertex_index(), mesh);

  typedef typename boost::range_value<PolygonRange>::type                           Polygon;

  std::size_t index = points.size(); // so that multiple meshes can be put into the same soup
  CGAL::internal::reserve(points, points.size() + vertices(mesh).size());
  CGAL::internal::reserve(polygons, polygons.size() + faces(mesh).size());

  for(const vertex_descriptor v : vertices(mesh))
  {
    points.emplace_back(get(vpm, v));
    put(vim, v, index++);
  }

  for(const face_descriptor f : faces(mesh))
  {
    CGAL::Iterator_range<CGAL::Halfedge_around_face_iterator<PolygonMesh> > incident_halfedges =
      CGAL::halfedges_around_face(halfedge(f, mesh), mesh);

    Polygon polygon;
    CGAL::internal::resize(polygon, incident_halfedges.size());

    std::size_t pos = 0;
    for(halfedge_descriptor h : incident_halfedges)
      polygon[pos++] = get(vim, target(h, mesh));

    polygons.push_back(polygon);
  }
}

template<typename PolygonMesh, typename PointRange, typename PolygonRange>
void polygon_mesh_to_polygon_soup(const PolygonMesh& mesh,
                                  PointRange& points,
                                  PolygonRange& polygons)
{
  return polygon_mesh_to_polygon_soup(mesh, points, polygons, CGAL::parameters::all_default());
}

} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_POLYGON_MESH_TO_POLYGON_SOUP_H
