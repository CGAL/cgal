// Copyright (c) 2013  GeometryFactory (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Laurent Rineau, Xiang Gao
//

#ifndef CGAL_SPLIT_GRAPH_INTO_POLYLINES
#define CGAL_SPLIT_GRAPH_INTO_POLYLINES

#include <CGAL/disable_warnings.h>

#include <map>
#include <vector>
#include <utility>
#include <boost/graph/adjacency_list.hpp>
#include <CGAL/assertions.h>
#include <CGAL/tags.h>
#include <CGAL/Iterator_range.h>

namespace CGAL {

namespace internal{
struct IsTerminalDefault
{
  template <typename VertexDescriptor, typename Graph>
  bool operator ()(VertexDescriptor& , const Graph& )
  {
    return false;
  }
};

template<class T, class ED>
class BGL_sgip_visitor_has_add_edge
{
private:
  template<class U, U>
  class check {};

  template<class C>
  static char f(check<void(C::*)(ED), &C::add_edge>*);

  template<class C>
  static char f(check<void(C::*)(const ED&), &C::add_edge>*);

  template<class C>
  static char f(check<void(C::*)(ED&), &C::add_edge>*);

  template<class C>
  static int f(...);
public:
  static const bool value = (sizeof(f<T>(0)) == sizeof(char));
};

template <typename Visitor, typename Edge_descriptor>
void
bgl_sgip_maybe_call_visitor_add_edge_impl(Visitor&,
                                          Edge_descriptor,
                                          CGAL::Tag_false /*has_add_edge*/)
{

}

template <typename Visitor, typename Edge_descriptor>
void
bgl_sgip_maybe_call_visitor_add_edge_impl(Visitor& visitor,
                                          Edge_descriptor e,
                                          CGAL::Tag_true /*has_add_edge*/)
{
  visitor.add_edge(e);
}

template <typename Visitor, typename Edge_descriptor>
void bgl_sgip_maybe_call_visitor_add_edge(Visitor& visitor,
                                          Edge_descriptor e) {
  typedef BGL_sgip_visitor_has_add_edge<Visitor, Edge_descriptor> Has_add_edge;
  bgl_sgip_maybe_call_visitor_add_edge_impl
    ( visitor, e, CGAL::Boolean_tag<Has_add_edge::value>() );
}

template <class Graph>
struct Dummy_visitor_for_split_graph_into_polylines
{
  void start_new_polyline(){}
  void add_node(typename boost::graph_traits<Graph>::vertex_descriptor){}
  void add_edge(typename boost::graph_traits<Graph>::edge_descriptor){}
  void end_polyline(){}
};

template <typename G_copy, typename Less_on_orig_vertex_descriptors>
class Less_on_G_copy_vertex_descriptors {
  const G_copy& g_copy;
  const Less_on_orig_vertex_descriptors& less;
public:
  Less_on_G_copy_vertex_descriptors ( const G_copy& g_copy,
    const Less_on_orig_vertex_descriptors& less)
    : g_copy(g_copy), less(less) {}

  typedef typename boost::graph_traits<G_copy>::vertex_descriptor
    g_copy_vertex_descriptor;
  typedef typename boost::graph_traits<G_copy>::out_edge_iterator
    g_copy_out_edge_iterator;
  typedef typename boost::graph_traits<G_copy>::degree_size_type
    g_copy_degree_size_type;

  bool operator()(g_copy_vertex_descriptor v1,
                  g_copy_vertex_descriptor v2) const {
    if(less(g_copy[v1], g_copy[v2]))
      return true;
    else if(less(g_copy[v2], g_copy[v1]))
      return false;
    // If g_copy[v1] and g_copy[v2] are equivalent, then compare the
    // descriptors:
    if(v1 == v2) return false;
    //   - compare degrees:
    const g_copy_degree_size_type dv1 = degree(v1, g_copy);
    const g_copy_degree_size_type dv2 = degree(v2, g_copy);
    if(dv1 != dv2)
      return dv1 < dv2;
    if(dv1 == 0) return v1 < v2;
    ///  - then compare an adjacent vertex:
    g_copy_vertex_descriptor other_v1 = target(*out_edges(v1, g_copy).first,
                                               g_copy);
    g_copy_vertex_descriptor other_v2 = target(*out_edges(v2, g_copy).first,
                                               g_copy);
    return less(g_copy[other_v1], g_copy[other_v2]);
  }
}; // end class Less_on_G_copy_vertex_descriptors

// splits a graph at vertices with degree higher than two and at vertices where `is_terminal` returns `true`
// The vertices are duplicated, and new incident edges created.
// `OrigGraph` must be undirected
template <typename Graph,
          typename OrigGraph,
          typename IsTerminal>
void duplicate_terminal_vertices(Graph& graph,
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
  for(vertex_descriptor v : V)
  {
    typename boost::graph_traits<OrigGraph>::vertex_descriptor orig_v = graph[v];
    typename boost::graph_traits<Graph>::degree_size_type deg = degree(v, graph);
    if ((deg != 0 && is_terminal(orig_v, orig)) || deg > 2)
      {
        out_edge_iterator b, e;
        boost::tie(b, e) = out_edges(v, graph);
        std::vector<edge_descriptor> out_edges_of_v(b, e);
        for (unsigned int i = 1; i < out_edges_of_v.size(); ++i)
          {
            edge_descriptor e = out_edges_of_v[i];
            typename boost::graph_traits<OrigGraph>::edge_descriptor orig_e =
              graph[e];
            vertex_descriptor w = target(e, graph);
            remove_edge(e, graph);
            vertex_descriptor vc = add_vertex(graph);
            graph[vc] = orig_v;
            const std::pair<edge_descriptor, bool> pair = add_edge(vc, w, graph);
            graph[pair.first] = orig_e;
          }
        CGAL_assertion(degree(v, graph) == 1);
      }
  }

  // check all vertices are of degree 1 or 2 and that the source
  // and target of each edge are different vertices with different ids
  CGAL_assertion_code(
                      for(vertex_descriptor v : make_range(vertices(graph))){
                        typename boost::graph_traits<Graph>::degree_size_type
                          n = degree(v, graph);
                        CGAL_assertion( n == 0 || n == 1 || n == 2);
                      }
                      for(edge_descriptor e : make_range(edges(graph))){
                        vertex_descriptor v = target(e, graph);
                        vertex_descriptor w = source(e, graph);
                        CGAL_assertion(v != w);
                      }
                      ) // end of CGAL_assertion_code
} // end of duplicate_terminal_vertices

} // namespace internal

#ifndef DOXYGEN_RUNNING
template <typename Graph,
          typename Visitor,
          typename IsTerminal,
          typename LessForVertexDescriptors>
void
split_graph_into_polylines(const Graph& graph,
                           Visitor& polyline_visitor,
                           IsTerminal is_terminal,
                           LessForVertexDescriptors less);
#endif

/*!
\ingroup PkgBGLRef
splits into polylines the graph `g` at vertices of degree greater than 2
and at vertices for which `is_terminal(v,graph)==true`.
The polylines are reported using a visitor.
\tparam Graph a model of the `boost` concepts `VertexListGraph` and `EdgeListGraph`.
\tparam Visitor a class that provides:
        - <code>void start_new_polyline()</code>
          called when starting the description of a polyline.
        - <code>void add_node(typename boost::graph_traits<Graph>::%vertex_descriptor v)</code>
          called for each vertex `v` of the polyline currently described. If the polyline is closed
          this function will be called twice for the first vertex of the cycle picked (once after
          calling `start_new_polyline()` and once before the call to `end_polyline()`.
        - <code>void end_polyline()</code>
          called when the description of a polyline is finished.
\tparam IsTerminal A functor providing `bool operator()(boost::graph_traits<Graph>::%vertex_descriptor v, const Graph& g) const`
                   returning true if the vertex `v` of degree 2 is a polyline endpoint and false otherwise.

An overload without `is_terminal` is provided if no vertices but those of degree
different from 2 are polyline endpoints.

@todo Document the version with four parameters
*/
template <typename Graph,
          typename Visitor,
          typename IsTerminal>
void
split_graph_into_polylines(const Graph& graph,
                           Visitor& polyline_visitor,
                           IsTerminal is_terminal)
{
  typedef typename boost::graph_traits<Graph>::vertex_descriptor Graph_vertex_descriptor;
  std::less<Graph_vertex_descriptor> less;
  split_graph_into_polylines(graph, polyline_visitor, is_terminal, less);
}

template <typename Graph,
          typename Visitor,
          typename IsTerminal,
          typename LessForVertexDescriptors>
void
split_graph_into_polylines(const Graph& graph,
                           Visitor& polyline_visitor,
                           IsTerminal is_terminal,
                           LessForVertexDescriptors less)
{
  using boost::graph_traits;
  typedef typename graph_traits<Graph>::vertex_descriptor Graph_vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor Graph_edge_descriptor;

  typedef boost::adjacency_list <boost::setS, // this avoids parallel edges
                                 boost::vecS,
                                 boost::undirectedS,
                                 Graph_vertex_descriptor,
                                 Graph_edge_descriptor> G_copy;

  typedef typename graph_traits<G_copy>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<G_copy>::edge_descriptor edge_descriptor;
  typedef typename graph_traits<G_copy>::out_edge_iterator out_edge_iterator;

  // we make a copy of the input graph
  G_copy g_copy;
  {
    typedef std::map<typename graph_traits<Graph>::vertex_descriptor,
                     typename graph_traits<G_copy>::vertex_descriptor> V2vmap;
    V2vmap v2vmap;

    for(Graph_vertex_descriptor v : make_range(vertices(graph))){
      vertex_descriptor vc = add_vertex(g_copy);
      g_copy[vc] = v;
      v2vmap[v] = vc;
    }

    for(Graph_edge_descriptor e : make_range(edges(graph))){
      Graph_vertex_descriptor vs = source(e,graph);
      Graph_vertex_descriptor vt = target(e,graph);
      CGAL_warning_msg(vs != vt, "ignore self loops");
      if(vs != vt){
        const std::pair<edge_descriptor, bool> pair =
          add_edge(v2vmap[vs],v2vmap[vt],g_copy);
        g_copy[pair.first] = e;
      }
    }
  }
  // duplicate terminal vertices and vertices of degree more than 2
  internal::duplicate_terminal_vertices(g_copy, graph, is_terminal);

  // put polylines endpoint in a set
  typedef internal::Less_on_G_copy_vertex_descriptors<
    G_copy,
    LessForVertexDescriptors> G_copy_less;
  G_copy_less g_copy_less(g_copy, less);
  std::set<vertex_descriptor, G_copy_less> terminal(g_copy_less);

  for(vertex_descriptor v : make_range(vertices(g_copy))){
    typename graph_traits<G_copy>::degree_size_type n = degree(v, g_copy);
    if ( n == 1 ) terminal.insert(v);
    if ( n ==0 ){
      //isolated vertex
      polyline_visitor.start_new_polyline();
      polyline_visitor.add_node(g_copy[v]);
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
    polyline_visitor.add_node(g_copy[u]);
    while (degree(u,g_copy) != 0)
    {
      CGAL_assertion(degree(u,g_copy) == 1);
      out_edge_iterator b = out_edges(u, g_copy).first;
      vertex_descriptor v = target(*b, g_copy);
      CGAL_assertion(u!=v);
      polyline_visitor.add_node(g_copy[v]);
      internal::bgl_sgip_maybe_call_visitor_add_edge(polyline_visitor,
                                                     g_copy[*b]);
      if (degree(v, g_copy)==1)
        terminal.erase(v);
      remove_edge(b, g_copy);
      u = v;
    }
    polyline_visitor.end_polyline();
  }

  // do the same but for cycles
  while (num_edges(g_copy) != 0)
  {
    edge_descriptor first_edge = *edges(g_copy).first;
    vertex_descriptor u = source(first_edge, g_copy);

    polyline_visitor.start_new_polyline();
    polyline_visitor.add_node(g_copy[u]);

    u = target(first_edge, g_copy);
    polyline_visitor.add_node(g_copy[u]);
    internal::bgl_sgip_maybe_call_visitor_add_edge(polyline_visitor,
                                                   g_copy[first_edge]);
    remove_edge(first_edge, g_copy);

    while (degree(u,g_copy) != 0)
    {
      CGAL_assertion(degree(u,g_copy) == 1);
      out_edge_iterator b = out_edges(u, g_copy).first;
      vertex_descriptor v = target(*b, g_copy);
      CGAL_assertion(u!=v);
      polyline_visitor.add_node(g_copy[v]);
      internal::bgl_sgip_maybe_call_visitor_add_edge(polyline_visitor,
                                                     g_copy[*b]);
      remove_edge(b, g_copy);
      u = v;
    }
    polyline_visitor.end_polyline();
  }
}

/// \cond SKIP_IN_MANUAL
/// Split graph into polylines delimited by vertices of degree different from 2.
/// Then the graph is visited and Visitor is called to describe the polylines
/// Graph must be undirected
template <typename Graph,
          typename Visitor>
void
split_graph_into_polylines(const Graph& graph,
                           Visitor& polyline_visitor)
{
  split_graph_into_polylines(graph, polyline_visitor, internal::IsTerminalDefault());
}
/// \endcond

} //end of namespace CGAL

#include <CGAL/enable_warnings.h>

#endif //CGAL_SPLIT_GRAPH_INTO_POLYLINES
