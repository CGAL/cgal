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
#ifndef BOOST_GRAPH_MST_KRUSKAL_HPP
#define BOOST_GRAPH_MST_KRUSKAL_HPP

/*
 *Minimum Spanning Tree 
 *         Kruskal Algorithm
 *
 *Requirement:
 *      undirected graph
 */

#include <vector>
#include <queue>
#include <functional>

#include <boost/property_map.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <boost/graph/named_function_params.hpp>
#include <boost/pending/disjoint_sets.hpp>
#include <boost/pending/indirect_cmp.hpp>


namespace boost {

  // Kruskal's algorithm for Minimum Spanning Tree
  //
  // This is a greedy algorithm to calculate the Minimum Spanning Tree
  // for an undirected graph with weighted edges. The output will be a
  // set of edges.
  //

  namespace detail {

    template <class Graph, class OutputIterator, 
              class Rank, class Parent, class Weight>
    void
    kruskal_mst_impl(const Graph& G, 
                     OutputIterator spanning_tree_edges, 
                     Rank rank, Parent parent, Weight weight)
    {
      typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
      typedef typename graph_traits<Graph>::edge_descriptor Edge;
      function_requires<VertexListGraphConcept<Graph> >();
      function_requires<EdgeListGraphConcept<Graph> >();
      function_requires<OutputIteratorConcept<OutputIterator, Edge> >();
      function_requires<ReadWritePropertyMapConcept<Rank, Vertex> >();
      function_requires<ReadWritePropertyMapConcept<Parent, Vertex> >();
      function_requires<ReadablePropertyMapConcept<Weight, Edge> >();
      typedef typename property_traits<Weight>::value_type W_value;
      typedef typename property_traits<Rank>::value_type R_value;
      typedef typename property_traits<Parent>::value_type P_value;
      function_requires<ComparableConcept<W_value> >();
      function_requires<ConvertibleConcept<P_value, Vertex> >();
      function_requires<IntegerConcept<R_value> >();

      disjoint_sets<Rank, Parent>  dset(rank, parent);

      typename graph_traits<Graph>::vertex_iterator ui, uiend;
      for (boost::tie(ui, uiend) = vertices(G); ui != uiend; ++ui)
        dset.make_set(*ui);

      typedef indirect_cmp<Weight, std::greater<W_value> > weight_greater;
      weight_greater wl(weight);
      std::priority_queue<Edge, std::vector<Edge>, weight_greater> Q(wl);
      /*push all edge into Q*/
      typename graph_traits<Graph>::edge_iterator ei, eiend;
      for (boost::tie(ei, eiend) = edges(G); ei != eiend; ++ei) 
        Q.push(*ei);

      while (! Q.empty()) {
        Edge e = Q.top();
        Q.pop();
        Vertex u = dset.find_set(source(e, G));
        Vertex v = dset.find_set(target(e, G));
        if ( u != v ) {
          *spanning_tree_edges++ = e;
          dset.link(u, v);
        }
      }
    }

  } // namespace detail 

  // Named Parameters Variants

  template <class Graph, class OutputIterator>
  inline void 
  kruskal_minimum_spanning_tree(const Graph& g,
                                OutputIterator spanning_tree_edges)
  {
    typedef typename graph_traits<Graph>::vertices_size_type size_type;
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
    typedef typename property_map<Graph, vertex_index_t>::type index_map_t;
    typename graph_traits<Graph>::vertices_size_type
      n = num_vertices(g);
    std::vector<size_type> rank_map(n);
    std::vector<vertex_t> pred_map(n);

    detail::kruskal_mst_impl
      (g, spanning_tree_edges, 
       make_iterator_property_map(rank_map.begin(), get(vertex_index, g), rank_map[0]),
       make_iterator_property_map(pred_map.begin(), get(vertex_index, g), pred_map[0]),
       get(edge_weight, g));
  }

  template <class Graph, class OutputIterator, class P, class T, class R>
  inline void
  kruskal_minimum_spanning_tree(const Graph& g,
                                OutputIterator spanning_tree_edges, 
                                const bgl_named_params<P, T, R>& params)
  {
    typedef typename graph_traits<Graph>::vertices_size_type size_type;
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
    typename graph_traits<Graph>::vertices_size_type n;
    n = is_default_param(get_param(params, vertex_rank))
                                   ? num_vertices(g) : 1;
    std::vector<size_type> rank_map(n);
    n = is_default_param(get_param(params, vertex_predecessor))
                                   ? num_vertices(g) : 1;
    std::vector<vertex_t> pred_map(n);
    
    detail::kruskal_mst_impl
      (g, spanning_tree_edges, 
       choose_param
       (get_param(params, vertex_rank), 
        make_iterator_property_map
        (rank_map.begin(), 
         choose_pmap(get_param(params, vertex_index), g, vertex_index), rank_map[0])),
       choose_param
       (get_param(params, vertex_predecessor), 
        make_iterator_property_map
        (pred_map.begin(), 
         choose_const_pmap(get_param(params, vertex_index), g, vertex_index), 
         pred_map[0])),
       choose_const_pmap(get_param(params, edge_weight), g, edge_weight));
  }
    
} // namespace boost


#endif // BOOST_GRAPH_MST_KRUSKAL_HPP

