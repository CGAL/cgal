// Copyright (c) 2018 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Sebastien Loriot

#ifndef CGAL_POLYGON_MESH_PROCESSING_MERGE_BORDER_VERTICES_H
#define CGAL_POLYGON_MESH_PROCESSING_MERGE_BORDER_VERTICES_H

#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/properties.h>
#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>
#include <CGAL/Polygon_mesh_processing/stitch_borders.h>

#include <boost/unordered_set.hpp>
#include <boost/bind.hpp>

namespace CGAL{

namespace Polygon_mesh_processing{

/// \todo document me
/// It should probably go into BGL package
template <typename PolygonMesh, typename OutputIterator>
OutputIterator
extract_boundary_cycles(PolygonMesh& pm,
                        OutputIterator out)
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;

  boost::unordered_set<halfedge_descriptor> hedge_handled;
  BOOST_FOREACH(halfedge_descriptor h, halfedges(pm))
  {
    if(is_border(h, pm) && hedge_handled.insert(h).second)
    {
      *out++=h;
      BOOST_FOREACH(halfedge_descriptor h2, halfedges_around_face(h, pm))
        hedge_handled.insert(h2);
    }
  }
  return out;
}

/// \ingroup PMP_repairing_grp
/// \todo document me
template <typename PolygonMesh, class VertexRange>
void merge_vertices(const VertexRange& vertices,
                          PolygonMesh& pm)
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;

  vertex_descriptor v_kept=*boost::begin(vertices);
  std::vector<vertex_descriptor> vertices_to_rm;

  BOOST_FOREACH(vertex_descriptor vd, vertices)
  {
    if (vd==v_kept) continue; // skip identical vertices
    if (edge(vd, v_kept, pm).second) continue; // skip null edges
    bool shall_continue=false;
    BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_target(v_kept, pm))
    {
      if (edge(vd, source(hd, pm), pm).second)
      {
        shall_continue=true;
        break;
      }
    }
    if (shall_continue) continue; // skip vertices already incident to the same vertex

    internal::update_target_vertex(halfedge(vd, pm), v_kept, pm);
    vertices_to_rm.push_back(vd);
  }

  BOOST_FOREACH(vertex_descriptor vd, vertices_to_rm)
    remove_vertex(vd, pm);
}

/// \ingroup PMP_repairing_grp
/// \todo document me
template <class PolygonMesh, class NamedParameter>
void merge_duplicated_vertices_in_boundary_cycle(typename boost::graph_traits<PolygonMesh>::halfedge_descriptor h,
                                                       PolygonMesh& pm,
                                                 const NamedParameter& np)
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
  typedef typename GetVertexPointMap<PolygonMesh, NamedParameter>::const_type Vpm;
  typedef typename boost::property_traits<Vpm>::value_type Point_3;
  Vpm vpm = choose_param(get_param(np, internal_np::vertex_point),
                         get_const_property_map(vertex_point, pm));

  // collect all the vertices of the cycle
  std::vector<vertex_descriptor> vertices;
  halfedge_descriptor start=h;
  do{
    vertices.push_back(target(h,pm));
    h=next(h, pm);
  }while(start!=h);

  // sort vertices using their point to ease the detection
  // of vertices with identical points
  CGAL::Property_map_to_unary_function<Vpm> Get_point(vpm);
  std::sort( vertices.begin(), vertices.end(),
             boost::bind(std::less<Point_3>(), boost::bind(Get_point,_1),
                                               boost::bind(Get_point, _2)) );
  std::size_t nbv=vertices.size();
  std::size_t i=1;

  std::vector< std::vector<vertex_descriptor> > identical_vertices;
  while(i!=nbv)
  {
    if (get(vpm, vertices[i]) == get(vpm, vertices[i-1]))
    {
      identical_vertices.push_back( std::vector<vertex_descriptor>() );
      identical_vertices.back().push_back(vertices[i-1]);
      identical_vertices.back().push_back(vertices[i]);
      while(++i!=nbv)
      {
        if (get(vpm, vertices[i]) == get(vpm, vertices[i-1]))
          identical_vertices.back().push_back(vertices[i]);
        else
          break;
      }
    }
    ++i;
  }
  BOOST_FOREACH(const std::vector<vertex_descriptor>& vrtcs, identical_vertices)
    merge_vertices(vrtcs, pm);
}

/// \ingroup PMP_repairing_grp
/// \todo document me
template <class PolygonMesh, class NamedParameter>
void merge_duplicated_vertices_in_boundary_cycles(      PolygonMesh& pm,
                                                  const NamedParameter& np)
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;

  std::vector<halfedge_descriptor> cycles;
  extract_boundary_cycles(pm, std::back_inserter(cycles));

  BOOST_FOREACH(halfedge_descriptor h, cycles)
    merge_duplicated_vertices_in_boundary_cycle(h, pm, np);
}

template <class PolygonMesh>
void merge_duplicated_vertices_in_boundary_cycles(PolygonMesh& pm)
{
  merge_duplicated_vertices_in_boundary_cycles(pm, parameters::all_default());
}

} } // end of CGAL::Polygon_mesh_processing

#endif //CGAL_POLYGON_MESH_PROCESSING_MERGE_BORDER_VERTICES_H
