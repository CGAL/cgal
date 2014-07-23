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
/// Graph must be undirected
template <typename Graph,
          typename IsTerminal>
void split_graph_into_polylines(Graph& graph,
                                IsTerminal is_terminal)
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

    if (degree(v, graph) > 2 || is_terminal(v, graph))
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

  // check all vertices are of degree 1 or 2 and that the source
  // and target of each edge are different vertices with different ids
  CGAL_assertion_code(
  BOOST_FOREACH(vertex_descriptor v, vertices(graph))
  {
    typename boost::graph_traits<Graph>::degree_size_type
      n = degree(v, graph);

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

/// Split graph into polylines delimited by node of degree 1 or more,
/// and vertices for which `is_terminal(v,graph)==true`.
/// Then the graph is visited and Visitor is called to describe the polylines
/// Graph must be undirected
template <typename Graph,
          typename Visitor,
          typename IsTerminal>
void
split_graph_into_polylines(Graph graph,
                           Visitor& polyline_visitor,
                           IsTerminal is_terminal)
{
  typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<Graph>::vertex_iterator vertex_iterator;
  typedef typename boost::graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename boost::graph_traits<Graph>::out_edge_iterator out_edge_iterator;

  // duplicate terminal vertices and vertices of degree more than 2
  split_graph_into_polylines(graph, is_terminal);

  // put polylines endpoint in a set
  std::set<vertex_descriptor> terminal;
  vertex_iterator b,e;
  for (boost::tie(b, e) = vertices(graph); b != e; ++b)
    if (degree(*b, graph) == 1) terminal.insert(*b);

  // go over all polylines and provide the description to the visitor
  while (!terminal.empty())
  {
    typename std::set<vertex_descriptor>::iterator it = terminal.begin();
    vertex_descriptor u = *it;
    terminal.erase(it);
    polyline_visitor.onNewPolyline();
    polyline_visitor.onAddNode(graph[u].id);

    while (degree(u,graph) != 0)
    {
      CGAL_assertion(degree(u,graph) == 1);
      out_edge_iterator b = out_edges(u, graph).first;
      vertex_descriptor v = target(*b, graph);
      CGAL_assertion(u!=v);
      polyline_visitor.onAddNode(graph[v].id);
      remove_edge(b, graph);
      u = v;
    }
    terminal.erase(u);
  }

  // do the same but for cycles
  while (num_edges(graph) != 0)
  {
    edge_descriptor first_edge = *edges(graph).first;
    vertex_descriptor u = source(first_edge, graph);

    polyline_visitor.onNewPolyline();
    polyline_visitor.onAddNode(graph[u].id);

    u = target(first_edge, graph);
    remove_edge(first_edge, graph);

    while (degree(u,graph) != 0)
    {
      CGAL_assertion(degree(u,graph) == 1);
      out_edge_iterator b = out_edges(u, graph).first;
      vertex_descriptor v = target(*b, graph);
      CGAL_assertion(u!=v);
      polyline_visitor.onAddNode(graph[v].id);
      remove_edge(b, graph);
      u = v;
    }
  }
}

} //end of namespace CGAL

#endif //CGAL_EXTRACT_MAXIMAL_POLYLINES
