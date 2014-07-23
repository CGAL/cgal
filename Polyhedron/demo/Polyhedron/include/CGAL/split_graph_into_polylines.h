// Copyright (c) 2013  GeometryFactory (France). All rights reserved.
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
//
// Author(s)     : Laurent Rineau, Xiang Gao
//

#ifndef CGAL_EXTRACT_MAXIMAL_POLYLINES
#define CGAL_EXTRACT_MAXIMAL_POLYLINES

#include <boost/foreach.hpp>
#include <iostream>


namespace CGAL {

template <typename Graph, typename PointPropertyMap>
void dump_graph_edges(std::ostream& out, const Graph& g,
                      PointPropertyMap points_pmap)
{
  typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<Graph>::edge_descriptor edge_descriptor;

  BOOST_FOREACH(edge_descriptor e, edges(g))
  {
    vertex_descriptor s = source(e, g);
    vertex_descriptor t = target(e, g);
    out.precision(17);
    out << "2 " << get(points_pmap,g[s].id) << " " << get(points_pmap, g[t].id) << "\n";
  }
}

template <typename Graph, typename PointPropertyMap>
void dump_graph_edges(const char* filename, const Graph& g,
                      PointPropertyMap points_pmap)
{
  std::ofstream out(filename);
  dump_graph_edges(out, g, points_pmap);
}

struct IsTerminalDefault
{
  template <typename VertexDescriptor, typename Graph>
  bool operator ()(VertexDescriptor& , Graph& )
  {
    return false;
  }
};

/// Splits a graph at vertices with degree higher than two.
/// The vertices are duplicated, and new incident edges created.
template <typename Graph,
          typename IsTerminal>
void split_in_polylines(Graph graph,
                        IsTerminal is_terminal,
                        std::size_t /* features_id_offset */ = 0)
{
  typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename boost::graph_traits<Graph>::vertex_iterator vertex_iterator;
  typedef typename boost::graph_traits<Graph>::out_edge_iterator out_edge_iterator;

  vertex_iterator b,e;
  boost::tie(b,e) = vertices(graph);
  std::vector<vertex_descriptor> V(b,e);
  for (typename std::vector<vertex_descriptor>::iterator it = V.begin();
      it != V.end();
      ++it) {
    vertex_descriptor v = *it;
    bool split = false;

    // const typename Kernel::Point_3& p = graph[v].id;
    // std::cerr << "At point   ( " << p << " )  "
    //           << "degree=" << out_degree(v, graph) << "\n";

    if (out_degree(v, graph) > 2 || is_terminal(v, graph))
    {
      split = true;
    }

    if (split)
    {
      out_edge_iterator b, e;
      boost::tie(b, e) = out_edges(v, graph);
      std::vector<edge_descriptor> E(b, e);
      for (unsigned int i = 1; i < E.size(); ++i)
      {
        edge_descriptor e = E[i];
        vertex_descriptor w = target(e, graph);
        remove_edge(e, graph);
        vertex_descriptor vc = add_vertex(graph);
        graph[vc].id = graph[v].id;
        add_edge(vc, w, graph);
      }
    }
  }
  CGAL_assertion_code(
  BOOST_FOREACH(vertex_descriptor v, vertices(graph))
  {
    typename boost::graph_traits<Graph>::degree_size_type
      n = out_degree(v, graph);

    CGAL_assertion(n == 1 || n == 2);
  }
  BOOST_FOREACH(edge_descriptor e, edges(graph))
  {
    vertex_descriptor v = target(e, graph);
    vertex_descriptor w = source(e, graph);
    CGAL_assertion(v != w);
    CGAL_assertion(graph[v].id != graph[w].id);
  }
  ) // end of CGAL_assertion_code
}

template <typename Graph,
          typename SetOfCornerVertexIds,
          typename Visitor,
          typename IsTerminal>
std::pair<std::size_t, std::size_t>
split_in_polylines(Graph graph,
                   SetOfCornerVertexIds& corner_ids,
                   Visitor& polyline_visitor,
                   IsTerminal is_terminal,
                   std::size_t features_id_offset = 0)
{
  std::size_t curve_id = features_id_offset;

//  dump_graph_edges("edges.polylines.txt", graph, points_pmap);
  typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<Graph>::vertex_iterator vertex_iterator;
  typedef typename boost::graph_traits<Graph>::out_edge_iterator out_edge_iterator;

  split_in_polylines(graph, is_terminal, features_id_offset);
  std::set<vertex_descriptor> terminal;
  vertex_iterator b,e;
  for (boost::tie(b, e) = vertices(graph); b != e; ++b)
  {
    if (degree(*b, graph) == 1)
    {
      // std::cerr << "terminal " << graph[*b] << std::endl;
      terminal.insert(*b);
      corner_ids.insert(graph[*b].id);
    }
  }
  // std::cerr << terminal.size() << " terminal vertices\n";

  while (!terminal.empty())
  {
    typename std::set<vertex_descriptor>::iterator it = terminal.begin();
    vertex_descriptor u = *it;
    terminal.erase(u);
    polyline_visitor.onNewPolyline();
    polyline_visitor.onAddNode(graph[u].id);

    while (out_degree(u,graph) != 0)
    {
      CGAL_assertion(out_degree(u,graph) == 1);
      out_edge_iterator b, e;
      boost::tie(b, e) = out_edges(u, graph);
      vertex_descriptor v = target(*b, graph);
//      CGAL_assertion(get(points_pmap, graph[v].id) != polyline.back());
      polyline_visitor.onAddNode(graph[v].id);
      remove_edge(b, graph);
      u = v;
    }
    terminal.erase(u);

    //// Comment that piece of code. Probably Mesh_3 does not like cycles
    //// with corners, and we split them arbitrary in two polylines.
    // if(polyline.back() == polyline.front())
    // {
    //   CGAL_assertion(polyline.size() > 3);
    //   // Fake cycle. We intended that cycle to be split at polyline.front()
    //   // Split the line in two, arbitrary.
    //   std::size_t n = polyline.size() / 2;
    //   Polyline new_line(polyline.begin() + n,
    //                     polyline.end());
    //   polyline.resize(n+1);
    //   // std::cerr << "  split the line in two\n";
    //   polylines.push_back(new_line);
    // // std::cerr << "polyline with " << new_line.size() << " vertices\n";
    // }

    // std::cerr << "polyline with " << polyline.size() << " vertices\n";

    ++curve_id;
  }
  // std::cerr << polylines.size() << " polylines\n";
//  dump_graph_edges("only-cycle-edges.polylines.txt", graph, points_pmap);

  std::size_t nb_cycles = 0;
  // process cycles
  while (num_edges(graph) != 0)
  {
    vertex_descriptor u = source(*edges(graph).first, graph);

    polyline_visitor.onNewPolyline();
    polyline_visitor.onAddNode(graph[u].id);

    ++nb_cycles;

    CGAL_assertion_code(bool first = true);
    while (out_degree(u,graph) != 0)
    {
      CGAL_assertion(out_degree(u,graph) == 1 ||
                     (first && out_degree(u, graph) == 2));
      out_edge_iterator b, e;
      boost::tie(b, e) = out_edges(u, graph);
      vertex_descriptor v = target(*b, graph);
      polyline_visitor.onAddNode(graph[v].id);
      remove_edge(b, graph);
      u = v;
      CGAL_assertion_code(first = false);
    }
    // std::cerr << "cycle with " << polyline.size() - 1  << " vertices\n";
    ++curve_id;
  }
//  CGAL_assertion(polylines_of_ids.size() == (curve_id - features_id_offset));

  // std::cerr << nb_cycles << " cycles\n";
  return std::pair<std::size_t, std::size_t>(curve_id - features_id_offset,
                                             corner_ids.size());
}

} //end of namespace CGAL

#endif //CGAL_EXTRACT_MAXIMAL_POLYLINES
