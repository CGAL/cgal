// Copyright (c) 2012-2015  GeometryFactory Sarl (France)
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
// $URL:$
// $Id:$
//
// Author(s)     : Andreas Fabri, Laurent Rineau

#ifndef CGAL_SPLIT_IN_POLYLINES_H
#define CGAL_SPLIT_IN_POLYLINES_H

#include <vector>
#include <set>
#include <iostream>
#include <fstream>
#include <boost/foreach.hpp>
#include <CGAL/number_utils.h>
#include <boost/graph/graph_traits.hpp>

namespace CGAL {
namespace internal {
namespace Mesh_3 {

template <typename Graph>
void dump_graph_edges(std::ostream& out, const Graph& g)
{
  typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<Graph>::edge_descriptor edge_descriptor;

  BOOST_FOREACH(edge_descriptor e, edges(g))
  {
    vertex_descriptor s = source(e, g);
    vertex_descriptor t = target(e, g);
    out.precision(17);
    out << "2 " << g[s] << " " << g[t] << "\n";
  }
}

template <typename Graph>
void dump_graph_edges(const char* filename, const Graph& g)
{
  std::ofstream out(filename);
  dump_graph_edges(out, g);
}

/// Splits a graph at vertices with degree higher than two.
/// The vertices are duplicated, and new incident edges created.
template <typename Graph, typename Kernel>
void split_in_polylines(Graph& G, Kernel)
{
  typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename boost::graph_traits<Graph>::vertex_iterator vertex_iterator;
  typedef typename boost::graph_traits<Graph>::out_edge_iterator out_edge_iterator;

  vertex_iterator b,e;
  boost::tie(b,e) = vertices(G);
  std::vector<vertex_descriptor> V(b,e);
  for(typename std::vector<vertex_descriptor>::iterator it = V.begin();
      it != V.end();
      ++it){
    vertex_descriptor v = *it;
    bool split = false;

    if(out_degree(v,G) > 2) {
      split = true;
    } else if(out_degree(v, G) == 2) {
      out_edge_iterator out_edge_it, out_edges_end;
      boost::tie(out_edge_it, out_edges_end) = out_edges(v, G);

      vertex_descriptor v1 = target(*out_edge_it++, G);
      vertex_descriptor v2 = target(*out_edge_it++, G);
      CGAL_assertion(out_edge_it == out_edges_end);

      const typename Kernel::Point_3& p = G[v];
      const typename Kernel::Point_3& p1 = G[v1];
      const typename Kernel::Point_3& p2 = G[v2];

      const typename Kernel::Vector_3 e1 = p1 - p;
      const typename Kernel::Vector_3 e2 = p2 - p;
      const typename Kernel::FT sc_prod = e1 * e2;
      if( sc_prod >= 0 ||   // angle < 135 degrees (3*pi/4)
          (sc_prod < 0 &&
           CGAL::square(sc_prod) < (e1 * e1) * (e2 * e2) / 2 ) )
      {
        split = true;
      }
    }

    if(split) {
      out_edge_iterator b,e;
      boost::tie(b,e) = out_edges(v,G);
      std::vector<edge_descriptor> E(b,e);
      for(unsigned int i = 1; i < E.size(); ++i){
        edge_descriptor e = E[i];
        vertex_descriptor w = target(e,G);
        remove_edge(e,G);
        vertex_descriptor vc = add_vertex(G);
        G[vc] = G[v];
        add_edge(vc,w,G);
      }
    }
  }
  CGAL_assertion_code(
  BOOST_FOREACH(vertex_descriptor v, vertices(G))
  {
    typename boost::graph_traits<Graph>::degree_size_type 
      n = out_degree(v, G);

    CGAL_assertion(n == 1 || n == 2);
  }
  BOOST_FOREACH(edge_descriptor e, edges(G))
  {
    vertex_descriptor v = target(e,G);
    vertex_descriptor w = source(e, G);
    CGAL_assertion(v != w);
    CGAL_assertion(G[v] != G[w]);
  }
                      )
}


template <typename Graph, typename Polylines_container, typename Kernel>
void split_in_polylines(Graph& G, Polylines_container& polylines, Kernel k)
{
  typedef typename Polylines_container::value_type Polyline;
  typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<Graph>::vertex_iterator vertex_iterator;
  typedef typename boost::graph_traits<Graph>::out_edge_iterator out_edge_iterator;

  // dump_graph_edges("edges.polylines.txt", G);
  split_in_polylines(G, k);
  std::set<vertex_descriptor> terminal;
  vertex_iterator b,e;
  for(boost::tie(b,e) = vertices(G); b!=e; ++b){
    if(degree(*b,G) == 1){
      terminal.insert(*b);
    }
  }

  while(! terminal.empty()){
    typename std::set<vertex_descriptor>::iterator it = terminal.begin();
    vertex_descriptor u = *it;
    terminal.erase(u);
    Polyline V;
    polylines.push_back(V);
    Polyline& polyline = polylines.back();
    polyline.push_back(G[u]);

    while(out_degree(u,G)!=0){
      CGAL_assertion(out_degree(u,G) == 1);
      out_edge_iterator b,e;
      boost::tie(b,e) = out_edges(u,G);
      vertex_descriptor v = target(*b,G);
      CGAL_assertion(G[v] != polyline.back());
      polyline.push_back(G[v]);
      remove_edge(b,G);
      u = v;
    }
    terminal.erase(u);

    if(polyline.back() == polyline.front())
    {
      CGAL_assertion(polyline.size() > 3);
      // Fake cycle. We intended that cycle to be split at polyline.front()
      // Split the line in two, arbitrary.
      std::size_t n = polyline.size() / 2;
      Polyline new_line(polyline.begin() + n,
                        polyline.end());
      polyline.resize(n+1);
      polylines.push_back(new_line);
    }
  }
  // dump_graph_edges("only-cycle-edges.polylines.txt", G);

  std::size_t nb_cycles = 0;
  // process cycles
  while(num_edges(G) != 0)
  {
    vertex_descriptor u = source(*edges(G).first, G);

    Polyline V;
    polylines.push_back(V);
    Polyline& polyline = polylines.back();
    polyline.push_back(G[u]);

    ++nb_cycles;

    CGAL_assertion_code(bool first = true);
    while(out_degree(u,G)!=0){
      CGAL_assertion(out_degree(u,G) == 1 ||
                     (first && out_degree(u, G) == 2));
      out_edge_iterator b,e;
      boost::tie(b,e) = out_edges(u,G);
      vertex_descriptor v = target(*b,G);
      polyline.push_back(G[v]);
      remove_edge(b,G);
      u = v;
      CGAL_assertion_code(first = false);
    }
  }
}

} // namespace Mesh_3
} // namespace internal
} // namespace CGAL

#endif // CGAL_SPLIT_IN_POLYLINES_H
