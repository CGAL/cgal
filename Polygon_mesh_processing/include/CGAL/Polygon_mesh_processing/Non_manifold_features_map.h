// Copyright (c) 2019 GeometryFactory (France).
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


#ifndef CGAL_POLYGON_MESH_PROCESSING_NON_MANIFOLD_FEATURES_MAP
#define CGAL_POLYGON_MESH_PROCESSING_NON_MANIFOLD_FEATURES_MAP

namespace CGAL {
namespace Polygon_mesh_processing {

//TODO: right now the base name parameter mechanism will make a deep copy, we probably want to avoid that
template <class PolygonMesh>
struct Non_manifold_features_map
{
  typedef boost::graph_traits<PolygonMesh> Graph_traits;
  typedef typename Graph_traits::edge_descriptor edge_descriptor;
  typedef typename Graph_traits::halfedge_descriptor halfedge_descriptor;
  typedef dynamic_edge_property_t<std::size_t> Edge_to_id_tag;
  typedef typename boost::property_map<PolygonMesh, Edge_to_id_tag>::type Edge_to_nm_id;
  Edge_to_nm_id e_nm_id;
  std::vector< std::vector<edge_descriptor> > non_manifold_edges;

  Non_manifold_features_map()
  {}

  template <class Vpm>
  Non_manifold_features_map(PolygonMesh& pm, Vpm vpm)
    : e_nm_id(get(Edge_to_id_tag(), pm))
  {
    typedef typename boost::property_traits<Vpm>::value_type Point_3;

    std::map< std::pair<Point_3, Point_3>, std::vector<edge_descriptor> > edge_map;
    for(edge_descriptor ed : edges(pm))
    {
      put(e_nm_id, ed, std::size_t(-1)); // init map

      halfedge_descriptor hd = halfedge(ed, pm);
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
        std::cout << "Found a non-manifold edge, nb representatives: " << p.second.size() << "\n";
        for (edge_descriptor ed : p.second)
          put(e_nm_id, ed, non_manifold_edges.size());
        non_manifold_edges.resize(non_manifold_edges.size()+1);
        // we steal the vertor from the map
        p.second.swap(non_manifold_edges.back());
      }
    }
  }
};

} } // end of CGAL::Polygon_mesh_processing

#endif
