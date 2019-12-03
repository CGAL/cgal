// Copyright (c) 2009-2010 INRIA Sophia-Antipolis (France).
// Copyright (c) 2014-2017 GeometryFactory Sarl (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent Rineau, Stéphane Tayeb
//
//******************************************************************************
// File Description :
//
//******************************************************************************

#ifndef CGAL_MESH_3_INTERNAL_MESH_3_HELPERS_H
#define CGAL_MESH_3_INTERNAL_MESH_3_HELPERS_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/enum.h>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/property_map.h>

#include <boost/graph/graph_traits.hpp>

#include <iostream>
#include <fstream>

namespace CGAL {

/// @cond DEVELOPERS
namespace Mesh_3 {
namespace internal {

template <typename Graph>
void dump_graph_edges(std::ostream& out, const Graph& g)
{
  typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<Graph>::edge_descriptor edge_descriptor;

  out.precision(17);
  for(edge_descriptor e : edges(g))
  {
    vertex_descriptor s = source(e, g);
    vertex_descriptor t = target(e, g);
    out << "2 " << g[s] << " " << g[t] << "\n";
  }
}

template <typename Graph>
void dump_graph_edges(const char* filename, const Graph& g)
{
  std::ofstream out(filename);
  dump_graph_edges(out, g);
}

template <typename Kernel>
struct Angle_tester
{
  template <typename vertex_descriptor, typename Graph>
  bool operator()(vertex_descriptor& v, const Graph& g) const
  {
    typedef typename boost::graph_traits<Graph>::out_edge_iterator out_edge_iterator;
    if (out_degree(v, g) != 2)
      return true;
    else
    {
      out_edge_iterator out_edge_it, out_edges_end;
      boost::tie(out_edge_it, out_edges_end) = out_edges(v, g);

      vertex_descriptor v1 = target(*out_edge_it++, g);
      vertex_descriptor v2 = target(*out_edge_it++, g);
      CGAL_assertion(out_edge_it == out_edges_end);

      const typename Kernel::Point_3& p = g[v];
      const typename Kernel::Point_3& p1 = g[v1];
      const typename Kernel::Point_3& p2 = g[v2];

      return (CGAL::angle(p1, p, p2) == CGAL::ACUTE);
    }
  }
};

template <typename Polyhedron>
struct Is_featured_edge {
  const Polyhedron* polyhedron;
  typename boost::property_map<Polyhedron, edge_is_feature_t>::type eifm;
  Is_featured_edge()
    : polyhedron(0)
  {} // required by boost::filtered_graph

  Is_featured_edge(const Polyhedron& polyhedron)
    : polyhedron(&polyhedron), eifm(get(edge_is_feature,polyhedron))
  {}

  bool operator()(typename boost::graph_traits<Polyhedron>::edge_descriptor e) const {
    return get(eifm, e);
  }
}; // end Is_featured_edge<Polyhedron>

template <typename Polyhedron>
struct Is_border_edge {
  const Polyhedron* polyhedron;
  Is_border_edge() : polyhedron(0) {} // required by boost::filtered_graph
  Is_border_edge(const Polyhedron& polyhedron) : polyhedron(&polyhedron) {}

  bool operator()(typename boost::graph_traits<Polyhedron>::edge_descriptor e) const {
    return is_border(halfedge(e, *polyhedron), *polyhedron) ||
      is_border(opposite(halfedge(e, *polyhedron), *polyhedron), *polyhedron);
  }
}; // end Is_featured_edge<Polyhedron>

template<typename Polyline_with_context,
         typename Graph>
struct Extract_polyline_with_context_visitor
{
  std::vector<Polyline_with_context>& polylines;
  const Graph& graph;

  Extract_polyline_with_context_visitor
  (const Graph& graph,
   typename std::vector<Polyline_with_context>& polylines)
    : polylines(polylines), graph(graph)
  {}

  void start_new_polyline()
  {
    polylines.push_back(Polyline_with_context());
  }

  void add_node(typename boost::graph_traits<Graph>::vertex_descriptor vd)
  {
    if(polylines.back().polyline_content.empty()) {
      polylines.back().polyline_content.push_back(graph[vd]);
    }
  }

  void add_edge(typename boost::graph_traits<Graph>::edge_descriptor ed)
  {
    typename boost::graph_traits<Graph>::vertex_descriptor
      s = source(ed, graph),
      t = target(ed, graph);
    Polyline_with_context& polyline = polylines.back();
    CGAL_assertion(!polyline.polyline_content.empty());
    if(polyline.polyline_content.back() != graph[s]) {
      polyline.polyline_content.push_back(graph[s]);
    } else if(polyline.polyline_content.back() != graph[t]) {
      // if the edge is zero-length, it is ignored
      polyline.polyline_content.push_back(graph[t]);
    }
    const typename boost::edge_bundle_type<Graph>::type &
      set_of_indices = graph[ed];
    polyline.context.adjacent_patches_ids.insert(set_of_indices.begin(),
                                                 set_of_indices.end());
  }

  void end_polyline()
  {
    // ignore degenerated polylines
    if(polylines.back().polyline_content.size() < 2)
      polylines.resize(polylines.size() - 1);
    // else {
    //   std::cerr << "Polyline with " << polylines.back().polyline_content.size()
    //             << " vertices, incident to "
    //             << polylines.back().context.adjacent_patches_ids.size()
    //             << " patches:\n ";
    //   for(auto p: polylines.back().polyline_content)
    //     std::cerr << " " << p;
    //   std::cerr << "\n";
    // }
  }
};


} // end CGAL::Mesh_3::internal
} // end CGAL::Mesh_3

/// @endcond

} // end CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_MESH_3_INTERNAL_MESH_3_HELPERS_H
