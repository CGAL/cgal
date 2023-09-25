// Copyright (c) 2023 University of Genova (Italy).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Claudio Mancinelli

#ifndef CGAL_POLYGON_MESH_PROCESSING_BSURF_LOCALLY_SHORTEST_PATH_H
#define CGAL_POLYGON_MESH_PROCESSING_BSURF_LOCALLY_SHORTEST_PATH_H

// #include <CGAL/license/Polygon_mesh_processing/bsurf.h>

#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/boost/graph/Dual.h>
#include <CGAL/Polygon_mesh_processing/locate.h>

#include <CGAL/boost/graph/dijkstra_shortest_paths.h>
#include <boost/graph/graph_traits.hpp>



namespace CGAL{
namespace Polygon_mesh_processing {

template <class TriangleMesh, class FT>
using Edge_location = std::pair< typename boost::graph_traits<TriangleMesh>::edge_descriptor, std::array<FT, 2> >;

template <typename FT, typename TriangleMesh, typename NamedParameters = parameters::Default_named_parameters>
#ifdef DOXYGEN_RUNNING
Point
construct_point(const Edge_location<TriangleMesh, FT>& loc,
#else
typename internal::Location_traits<TriangleMesh, NamedParameters>::Point
construct_point(const Edge_location<TriangleMesh,FT> & loc,
#endif
                const TriangleMesh& tm,
                const NamedParameters& np = parameters::default_values())
{
  typedef typename boost::graph_traits<TriangleMesh>::edge_descriptor            edge_descriptor;
  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type            Geom_traits;

  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::const_type  VertexPointMap;
  typedef typename boost::property_traits<VertexPointMap>::value_type            Point;
  typedef typename boost::property_traits<VertexPointMap>::reference             Point_reference;

  using parameters::get_parameter;
  using parameters::choose_parameter;

  CGAL_precondition(CGAL::is_triangle_mesh(tm));

  VertexPointMap vpm = parameters::choose_parameter(parameters::get_parameter(np, internal_np::vertex_point),
                                                    get_const_property_map(boost::vertex_point, tm));
  Geom_traits gt = choose_parameter<Geom_traits>(get_parameter(np, internal_np::geom_traits));

  edge_descriptor ed = loc.first;
  const Point_reference p0 = get(vpm, source(ed, tm));
  const Point_reference p1 = get(vpm, target(ed, tm));

  internal::Barycentric_point_constructor<Geom_traits, Point> bp_constructor;
  return bp_constructor(p0, loc.second[0], p1, loc.second[1], gt);
}

template <class FT, class TriangleMesh, class EdgeLocationRange>
void
locally_shortest_path(const Face_location<TriangleMesh, FT>& src,
                      const Face_location<TriangleMesh, FT>& tgt,
                      const TriangleMesh& tmesh,
                      EdgeLocationRange& edge_locations)
{
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;

  typename boost::property_map<TriangleMesh, CGAL::dynamic_face_property_t<face_descriptor> >::const_type
    predecessor_map =  get(CGAL::dynamic_face_property_t<face_descriptor>(), tmesh);
  typename boost::property_map<TriangleMesh, CGAL::dynamic_face_property_t<FT> >::const_type
    distance_map =  get(CGAL::dynamic_face_property_t<FT>(), tmesh);
  typename boost::property_map<TriangleMesh, CGAL::dynamic_edge_property_t<FT> >::const_type
    weight_map =  get(CGAL::dynamic_edge_property_t<FT>(), tmesh);

  Dual dual(tmesh);

  // TODO: fill the weight map

  // TODO try stopping dijkstra as soon tgt is out of the queue.
  boost::dijkstra_shortest_paths(dual, src.first,
                                 boost::distance_map(distance_map)
                                       .predecessor_map(predecessor_map)
                                       .weight_map(weight_map));

  std::vector<halfedge_descriptor> initial_path;

  auto common_halfedge=[](face_descriptor f1, face_descriptor f2, const TriangleMesh& tmesh)
  {
    halfedge_descriptor h=halfedge(f1, tmesh);
    for (int i=0; i<3;++i)
    {
      if (face(opposite(h, tmesh), tmesh)==f2)
        return h;
      h = next(h, tmesh);
    }
    CGAL_assertion(!"faces do no share a common edge");
    return halfedge_descriptor();
  };

  face_descriptor current_face = tgt.first;
  while (true)
  {
    face_descriptor prev = get(predecessor_map, current_face);
    halfedge_descriptor h=common_halfedge(current_face, prev, tmesh);
    initial_path.push_back(h);
    if (prev==src.first) break;
    current_face=prev;
  }
  std::reverse(initial_path.begin(), initial_path.end());
}


} } // CGAL::Polygon_mesh_processing

#endif
