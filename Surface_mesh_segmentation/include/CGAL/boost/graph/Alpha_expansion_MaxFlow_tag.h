#ifndef CGAL_BOOST_GRAPH_ALPHA_EXPANSION_GRAPHCUT_MAXFLOW_H
// Copyright (c) 2014  GeometryFactory Sarl (France).
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
// Author(s)     : Ilker O. Yaz, Simon Giraudot

#define CGAL_BOOST_GRAPH_ALPHA_EXPANSION_MAXFLOW_TAG_H

#include <CGAL/license/Surface_mesh_segmentation.h>

/// \cond SKIP_IN_MANUAL

#include <CGAL/boost/graph/alpha_expansion_graphcut.h>

namespace MaxFlow
{
#include <CGAL/internal/auxiliary/graph.h>
}


namespace CGAL
{

/**
 * @brief Implements alpha-expansion graph cut algorithm.
 *
 * For underlying max-flow algorithm, it uses the MAXFLOW software implemented by Boykov & Kolmogorov.
 *  Also no pre-allocation is made.
 */
class Alpha_expansion_MaxFlow_tag
{
public:

  typedef MaxFlow::Graph::node_id Vertex_descriptor;
  
private:
  
  MaxFlow::Graph graph;

public:

  void clear_graph()
  {
    graph = MaxFlow::Graph();
  }

  Vertex_descriptor add_vertex()
  {
    return graph.add_node();
  }

  void add_tweight (Vertex_descriptor& v, double w1, double w2)
  {
    graph.add_tweights(v, w1, w2);
  }

  void init_vertices()
  {
  }

  double max_flow()
  {
    return graph.maxflow();
  }

  template <typename VertexLabelMap, typename InputVertexDescriptor>
  void update(VertexLabelMap vertex_label_map,
              const std::vector<Vertex_descriptor>& inserted_vertices,
              InputVertexDescriptor vd,
              std::size_t vertex_i,
              std::size_t alpha)
  {
    if(get(vertex_label_map, vd) != alpha
       && graph.what_segment(inserted_vertices[vertex_i]) == MaxFlow::Graph::SINK)
      put(vertex_label_map, vd, alpha);
  }
  
  void add_edge (Vertex_descriptor& v1, Vertex_descriptor& v2, double w1, double w2)
  {
    graph.add_edge(v1, v2, w1, w2);
  }
};

}//namespace CGAL

/// \endcond

#endif //CGAL_BOOST_GRAPH_ALPHA_EXPANSION_MAXFLOW_TAG_H
