//
//=======================================================================
// Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
//
// This file is part of the Boost Graph Library
//
// You should have received a copy of the License Agreement for the
// Boost Graph Library along with the software; see the file LICENSE.
// If not, contact Office of Research, University of Notre Dame, Notre
// Dame, IN 46556.
//
// Permission to modify the code and to distribute modified code is
// granted, provided the text of this NOTICE is retained, a notice that
// the code was modified is included with the above COPYRIGHT NOTICE and
// with the COPYRIGHT NOTICE in the LICENSE file, and that the LICENSE
// file is distributed with the modified code.
//
// LICENSOR MAKES NO REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED.
// By way of example, but not limitation, Licensor MAKES NO
// REPRESENTATIONS OR WARRANTIES OF MERCHANTABILITY OR FITNESS FOR ANY
// PARTICULAR PURPOSE OR THAT THE USE OF THE LICENSED SOFTWARE COMPONENTS
// OR DOCUMENTATION WILL NOT INFRINGE ANY PATENTS, COPYRIGHTS, TRADEMARKS
// OR OTHER RIGHTS.
//=======================================================================
//

/*
  This file implements the function

  template <class EdgeListGraph, class Size, class P, class T, class R>
  bool bellman_ford_shortest_paths(EdgeListGraph& g, Size N, 
     const bgl_named_params<P, T, R>& params)
  
 */


#ifndef BOOST_GRAPH_BELLMAN_FORD_SHORTEST_PATHS_HPP
#define BOOST_GRAPH_BELLMAN_FORD_SHORTEST_PATHS_HPP

#include <boost/config.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/relax.hpp>
#include <boost/graph/visitors.hpp>
#include <boost/graph/named_function_params.hpp>

namespace boost {

  template <class Visitor, class Graph>
  struct BellmanFordVisitorConcept {
    void constraints() {
      function_requires< CopyConstructibleConcept<Visitor> >();
      vis.examine_edge(e, g);
      vis.edge_relaxed(e, g);
      vis.edge_not_relaxed(e, g);
      vis.edge_minimized(e, g);
      vis.edge_not_minimized(e, g);
    }
    Visitor vis;
    Graph g;
    typename graph_traits<Graph>::edge_descriptor e;
  };

  template <class Visitors = null_visitor>
  class bellman_visitor {
  public:
    bellman_visitor() { }
    bellman_visitor(Visitors vis) : m_vis(vis) { }

    template <class Edge, class Graph>
    void examine_edge(Edge u, Graph& g) {
      invoke_visitors(m_vis, u, g, on_examine_edge());
    }
    template <class Edge, class Graph>
    void edge_relaxed(Edge u, Graph& g) {
      invoke_visitors(m_vis, u, g, on_edge_relaxed());      
    }
    template <class Edge, class Graph>
    void edge_not_relaxed(Edge u, Graph& g) {
      invoke_visitors(m_vis, u, g, on_edge_not_relaxed());
    }
    template <class Edge, class Graph>
    void edge_minimized(Edge u, Graph& g) {
      invoke_visitors(m_vis, u, g, on_edge_minimized());
    }
    template <class Edge, class Graph>
    void edge_not_minimized(Edge u, Graph& g) {
      invoke_visitors(m_vis, u, g, on_edge_not_minimized());
    }
  protected:
    Visitors m_vis;
  };
  template <class Visitors>
  bellman_visitor<Visitors>
  make_bellman_visitor(Visitors vis) {
    return bellman_visitor<Visitors>(vis);
  }
  typedef bellman_visitor<> default_bellman_visitor;

  template <class EdgeListGraph, class Size, class WeightMap,
            class PredecessorMap, class DistanceMap,
            class BinaryFunction, class BinaryPredicate,
            class BellmanFordVisitor>
  bool bellman_ford_shortest_paths(EdgeListGraph& g, Size N, 
                         WeightMap weight, 
                         PredecessorMap pred,
                         DistanceMap distance, 
                         BinaryFunction combine, 
                         BinaryPredicate compare,
                         BellmanFordVisitor v)
  {
    function_requires<EdgeListGraphConcept<EdgeListGraph> >();
    typedef graph_traits<EdgeListGraph> GTraits;
    typedef typename GTraits::edge_descriptor Edge;
    typedef typename GTraits::vertex_descriptor Vertex;
    function_requires<ReadWritePropertyMapConcept<DistanceMap, Vertex> >();
    function_requires<ReadablePropertyMapConcept<WeightMap, Edge> >();
    typedef typename property_traits<DistanceMap>::value_type D_value;
    typedef typename property_traits<WeightMap>::value_type W_value;

    typename GTraits::edge_iterator i, end;

    for (Size k = 0; k < N; ++k) {
      bool at_least_one_edge_relaxed = false;
      for (tie(i, end) = edges(g); i != end; ++i) {
        v.examine_edge(*i, g);
        if (relax(*i, g, weight, pred, distance, combine, compare)) {
          at_least_one_edge_relaxed = true;
          v.edge_relaxed(*i, g);
        } else
          v.edge_not_relaxed(*i, g);
      }
      if (!at_least_one_edge_relaxed)
        break;
    }

    for (tie(i, end) = edges(g); i != end; ++i)
      if (compare(combine(get(distance, source(*i, g)), 
                          get(weight, *i)),
                  get(distance, target(*i,g))))
      {
        v.edge_not_minimized(*i, g);
        return false;
      } else
        v.edge_minimized(*i, g);

    return true;
  }

  namespace detail {

    template <class EdgeListGraph, class Size, class WeightMap,
              class DistanceMap, class P, class T, class R>
    bool bellman_dispatch(EdgeListGraph& g, Size N, 
                          WeightMap weight, DistanceMap distance, 
                          const bgl_named_params<P, T, R>& params)
    {
      typedef typename property_traits<DistanceMap>::value_type D;
      bellman_visitor<> null_vis;
      dummy_property_map dummy_pred;
      return bellman_ford_shortest_paths
        (g, N, weight, 
         choose_param(get_param(params, vertex_predecessor), dummy_pred),
         distance,
         choose_param(get_param(params, distance_combine_t()),
                      closed_plus<D>()),
         choose_param(get_param(params, distance_compare_t()),
                      std::less<D>()),
         choose_param(get_param(params, graph_visitor),
                       null_vis)
         );
    }

  } // namespace detail

  template <class EdgeListGraph, class Size, class P, class T, class R>
  bool bellman_ford_shortest_paths
    (EdgeListGraph& g, Size N, 
     const bgl_named_params<P, T, R>& params)
  {                                
    return detail::bellman_dispatch
      (g, N,
       choose_const_pmap(get_param(params, edge_weight), g, edge_weight),
       choose_pmap(get_param(params, vertex_distance), g, vertex_distance),
       params);
  }

  template <class EdgeListGraph, class Size>
  bool bellman_ford_shortest_paths(EdgeListGraph& g, Size N)
  {                                
    bgl_named_params<int,int> params(0);
    return bellman_ford_shortest_paths(g, N, params);
  }

} // namespace boost

#endif // BOOST_GRAPH_BELLMAN_FORD_SHORTEST_PATHS_HPP
