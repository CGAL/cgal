// Copyright (c) 2013  GeometryFactory (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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

#ifndef CGAL_SPLIT_GRAPH_INTO_POLYLINES
#define CGAL_SPLIT_GRAPH_INTO_POLYLINES

#include <map> 
#include <set> 
#include <boost/foreach.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <iostream>


namespace CGAL {

struct IsTerminalDefault
{
  template <typename VertexDescriptor, typename Graph>
  bool operator ()(VertexDescriptor& , const Graph& )
  {
    return false;
  }
};

template <class Graph>
struct Dummy_visitor_for_split_graph_into_polylines
{
  void start_new_polyline(){}
  void add_node(typename boost::graph_traits<Graph>::vertex_descriptor){}
  void end_polyline(){}
};


namespace internal {

/// Splits a graph at vertices with degree higher than two and at vertices where `is_terminal  returns `true`
/// The vertices are duplicated, and new incident edges created.
/// OrigGraph must be undirected
template <typename Graph,
          typename OrigGraph,
          typename IsTerminal>
void split_graph_into_polylines(Graph& graph,
                                const OrigGraph& orig,
                                IsTerminal is_terminal)
{
  typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename boost::graph_traits<Graph>::vertex_iterator vertex_iterator;
  typedef typename boost::graph_traits<Graph>::out_edge_iterator out_edge_iterator;

  vertex_iterator b,e;
  boost::tie(b,e) = vertices(graph);
  std::vector<vertex_descriptor> V(b,e);
  BOOST_FOREACH(vertex_descriptor v, V)
  {
    if (degree(v, graph) > 2 || is_terminal(graph[v], orig))
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
            graph[vc] = graph[v];
            add_edge(vc, w, graph);
          }
      }
  }

  // check all vertices are of degree 1 or 2 and that the source
  // and target of each edge are different vertices with different ids
  CGAL_assertion_code(
                      BOOST_FOREACH(vertex_descriptor v, vertices(graph)){
                        typename boost::graph_traits<Graph>::degree_size_type
                          n = degree(v, graph);
                        CGAL_assertion( n == 0 || n == 1 || n == 2);
                      }
                      BOOST_FOREACH(edge_descriptor e, edges(graph)){
                        vertex_descriptor v = target(e, graph);
                        vertex_descriptor w = source(e, graph);
                        CGAL_assertion(v != w);
                      }
                      ) // end of CGAL_assertion_code
    }
    
} // namespace internal

  
/// Split graph into polylines delimited by vertices of degree different from 2,
/// and vertices for which `is_terminal(v,graph)==true`.
/// Then the graph is visited and Visitor is called to describe the polylines
/// Graph must be undirected
template <typename Graph,
          typename Visitor,
          typename IsTerminal>
void
split_graph_into_polylines(const Graph& graph,
                           Visitor& polyline_visitor,
                           IsTerminal is_terminal)
{
  typedef typename boost::graph_traits<Graph>::vertex_descriptor Graph_vertex_descriptor;
  typedef typename boost::graph_traits<Graph>::edge_descriptor Graph_edge_descriptor;
  
  typedef boost::adjacency_list <boost::setS, // this avoids parallel edges
                                 boost::vecS, 
                                 boost::undirectedS,
                                 Graph_vertex_descriptor > G;

  typedef typename boost::graph_traits<G>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<G>::edge_descriptor edge_descriptor;
  typedef typename boost::graph_traits<G>::out_edge_iterator out_edge_iterator;
  
  // we make a copy of the input graph
  G g;
  {
    typedef std::map<typename boost::graph_traits<Graph>::vertex_descriptor,
                     typename boost::graph_traits<G>::vertex_descriptor> V2vmap;
    V2vmap v2vmap;
    
    BOOST_FOREACH(Graph_vertex_descriptor v, vertices(graph)){
      vertex_descriptor vc = add_vertex(g);
      g[vc] = v;
      v2vmap[v] = vc; 
    }
    

    BOOST_FOREACH(Graph_edge_descriptor e, edges(graph)){
      Graph_vertex_descriptor vs = source(e,graph);
      Graph_vertex_descriptor vt = target(e,graph);
      vertex_descriptor vsc, vtc;
      if(vs == vt){
        std::cerr << "ignore self loop\n";
      }else{
        typename V2vmap::iterator it;

        if((it = v2vmap.find(vs)) == v2vmap.end()){
          vsc = add_vertex(g);
          g[vsc] = vs;
          v2vmap[vs] = vsc;
        }else{
          vsc = it->second;
        }
        if((it = v2vmap.find(vt)) == v2vmap.end()){
          vtc = add_vertex(g);
          g[vtc] = vt;
          v2vmap[vt] = vtc;
        }else{
          vtc = it->second;
        }
        add_edge(vsc,vtc,g);
      }
    }
  }  
  // duplicate terminal vertices and vertices of degree more than 2
  internal::split_graph_into_polylines(g, graph, is_terminal);
  // put polylines endpoint in a set
  std::set<vertex_descriptor> terminal;

  BOOST_FOREACH(vertex_descriptor v, vertices(g)){
    typename boost::graph_traits<Graph>::degree_size_type n = degree(v, g);
    if ( n == 1 ) terminal.insert(v);
    if ( n ==0 ){
      //isolated vertex
      polyline_visitor.start_new_polyline();
      polyline_visitor.add_node(g[v]);
      polyline_visitor.end_polyline();
    }
  }

  // go over all polylines and provide the description to the visitor
  while (!terminal.empty())
  {
    typename std::set<vertex_descriptor>::iterator it = terminal.begin();
    vertex_descriptor u = *it;
    terminal.erase(it);
    polyline_visitor.start_new_polyline();
    polyline_visitor.add_node(g[u]);
    while (degree(u,g) != 0)
    {
      CGAL_assertion(degree(u,g) == 1);
      out_edge_iterator b = out_edges(u, g).first;
      vertex_descriptor v = target(*b, g);
      CGAL_assertion(u!=v);
      polyline_visitor.add_node(g[v]);
      remove_edge(b, g);
      u = v;
    }
    terminal.erase(u);
    polyline_visitor.end_polyline();
  }

  // do the same but for cycles
  while (num_edges(g) != 0)
  {
    edge_descriptor first_edge = *edges(g).first;
    vertex_descriptor u = source(first_edge, g);

    polyline_visitor.start_new_polyline();
    polyline_visitor.add_node(g[u]);

    u = target(first_edge, g);
    remove_edge(first_edge, g);
    polyline_visitor.add_node(g[u]);

    while (degree(u,g) != 0)
    {
      CGAL_assertion(degree(u,g) == 1);
      out_edge_iterator b = out_edges(u, g).first;
      vertex_descriptor v = target(*b, g);
      CGAL_assertion(u!=v);
      polyline_visitor.add_node(g[v]);
      remove_edge(b, g);
      u = v;
    }
    polyline_visitor.end_polyline();
  }
}

/// Split graph into polylines delimited by vertices of degree different from 2.
/// Then the graph is visited and Visitor is called to describe the polylines
/// Graph must be undirected
template <typename Graph,
          typename Visitor>
void
split_graph_into_polylines(const Graph& graph,
                           Visitor& polyline_visitor)
{
  split_graph_into_polylines(graph, polyline_visitor, IsTerminalDefault());
}


} //end of namespace CGAL

#endif //CGAL_SPLIT_GRAPH_INTO_POLYLINES
