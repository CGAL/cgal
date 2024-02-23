// Copyright (c) 2019 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Sebastien Loriot


#ifndef CGAL_POLYGON_MESH_PROCESSING_NON_MANIFOLD_FEATURE_MAP_H
#define CGAL_POLYGON_MESH_PROCESSING_NON_MANIFOLD_FEATURE_MAP_H

#include <CGAL/license/Polygon_mesh_processing/core.h>

namespace CGAL {
namespace Polygon_mesh_processing {

//TODO: right now the base name parameter mechanism will make a deep copy, we probably want to avoid that
template <class PolygonMesh>
struct Non_manifold_feature_map
{
  typedef boost::graph_traits<PolygonMesh> Graph_traits;
  typedef typename Graph_traits::vertex_descriptor vertex_descriptor;
  typedef typename Graph_traits::edge_descriptor edge_descriptor;
  typedef typename Graph_traits::halfedge_descriptor halfedge_descriptor;
  typedef dynamic_edge_property_t<std::size_t> Edge_to_id_tag;
  typedef dynamic_vertex_property_t<std::size_t> Vertex_to_id_tag;
  typedef typename boost::property_map<PolygonMesh, Edge_to_id_tag>::const_type Edge_to_nm_id;
  typedef typename boost::property_map<PolygonMesh, Vertex_to_id_tag>::const_type Vertex_to_nm_id;
  Edge_to_nm_id e_nm_id;
  Vertex_to_nm_id v_nm_id;
  std::vector< std::vector<edge_descriptor> > non_manifold_edges;
  std::vector< std::vector<vertex_descriptor> > non_manifold_vertices;

  Non_manifold_feature_map()
  {}

  template <class Vpm>
  Non_manifold_feature_map(const PolygonMesh& pm, Vpm vpm)
    : e_nm_id(get(Edge_to_id_tag(), pm))
    , v_nm_id(get(Vertex_to_id_tag(), pm))
  {
    typedef typename boost::property_traits<Vpm>::value_type Point_3;

  // detect non-manifold vertices
    std::map<Point_3, std::vector<vertex_descriptor> > vertex_map;
    for(vertex_descriptor vd : vertices(pm))
    {
      put(v_nm_id, vd, std::size_t(-1)); // init map
      vertex_map[get(vpm, vd)].push_back(vd);
    }

    for(std::pair< const Point_3,
                   std::vector<vertex_descriptor> >& p : vertex_map)
    {
      if (p.second.size()!=1)
      {
        for (vertex_descriptor vd : p.second)
          put(v_nm_id, vd, non_manifold_vertices.size());
        non_manifold_vertices.resize(non_manifold_vertices.size()+1);
        // we steal the vertor from the map
        p.second.swap(non_manifold_vertices.back());
      }
    }

  // detect non-manifold edges
    std::map< std::pair<Point_3, Point_3>, std::vector<edge_descriptor> > edge_map;
    for(edge_descriptor ed : edges(pm))
    {
      put(e_nm_id, ed, std::size_t(-1)); // init map

      halfedge_descriptor hd = halfedge(ed, pm);

      // an edge can be non-manifold only if both its vertices are non-manifold
      if ( get(v_nm_id, source(hd, pm))==std::size_t(-1) ||
           get(v_nm_id, target(hd, pm))==std::size_t(-1) ) continue;

      const Point_3& src = get(vpm, source(ed, pm));
      const Point_3& tgt = get(vpm, target(ed, pm));
      // TODO: what to do with null edges?
      if (src > tgt)
        hd = opposite(hd, pm);
      edge_map[ make_sorted_pair(src, tgt) ].push_back(edge(hd,pm));
    }

    for(std::pair< const std::pair<Point_3, Point_3>,
                   std::vector<edge_descriptor> >& p : edge_map)
    {
      if (p.second.size()!=1)
      {
        for (edge_descriptor ed : p.second)
          put(e_nm_id, ed, non_manifold_edges.size());
        non_manifold_edges.resize(non_manifold_edges.size()+1);
        // we steal the vertor from the map
        p.second.swap(non_manifold_edges.back());
      }
    }
  }

  void clear()
  {
    non_manifold_edges.clear();
    non_manifold_vertices.clear();
    e_nm_id = Edge_to_nm_id();
    v_nm_id = Vertex_to_nm_id();
  }
};

} } // end of CGAL::Polygon_mesh_processing

#endif
