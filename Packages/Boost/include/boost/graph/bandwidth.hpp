// Copyright (c) Jeremy Siek 2001, Marc Wintermantel 2002
//
// Permission to use, copy, modify, distribute and sell this software
// and its documentation for any purpose is hereby granted without fee,
// provided that the above copyright notice appears in all copies and
// that both that copyright notice and this permission notice appear
// in supporting documentation.  Silicon Graphics makes no
// representations about the suitability of this software for any
// purpose.  It is provided "as is" without express or implied warranty.

#ifndef BOOST_GRAPH_BANDWIDTH_HPP
#define BOOST_GRAPH_BANDWIDTH_HPP

#include <boost/graph/graph_traits.hpp>
#include <boost/detail/numeric_traits.hpp>

namespace boost {

  template <typename Graph, typename VertexIndexMap>
  typename graph_traits<Graph>::vertices_size_type
  ith_bandwidth(typename graph_traits<Graph>::vertex_descriptor i,
                const Graph& g,
                VertexIndexMap index)
  {
    typedef typename graph_traits<Graph>::vertices_size_type size_type;
    size_type b = 0;
    typename graph_traits<Graph>::out_edge_iterator e, end;
    for (tie(e, end) = out_edges(i, g); e != end; ++e) {
      int f_i = get(index, i);
      int f_j = get(index, target(*e, g));
      using namespace std; // to call abs() unqualified
      if(f_i > f_j)
      b = std::max(b, size_type(f_i - f_j));
    }
    return b;
  }

  template <typename Graph>
  typename graph_traits<Graph>::vertices_size_type
  ith_bandwidth(typename graph_traits<Graph>::vertex_descriptor i,
                const Graph& g)
  {
    return ith_bandwidth(i, g, get(vertex_index, g));
  }

  template <typename Graph, typename VertexIndexMap>
  typename graph_traits<Graph>::vertices_size_type
  bandwidth(const Graph& g, VertexIndexMap index)
  {
    typename graph_traits<Graph>::vertices_size_type b = 0;
    typename graph_traits<Graph>::vertex_iterator i, end;
    for (tie(i, end) = vertices(g); i != end; ++i)
        b = std::max(b, ith_bandwidth(*i, g, index));
    return b;
  }

  template <typename Graph>
  typename graph_traits<Graph>::vertices_size_type
  bandwidth(const Graph& g)
  {
    return bandwidth(g, get(vertex_index, g));
  }

  template <typename Graph, typename VertexIndexMap>
  typename graph_traits<Graph>::vertices_size_type
  edgesum(const Graph& g, VertexIndexMap index_map)
  {
    typedef typename graph_traits<Graph>::vertices_size_type size_type;
    typedef typename detail::numeric_traits<size_type>::difference_type diff_t;
    size_type sum = 0;
    typename graph_traits<Graph>::edge_iterator i, end;
    for (tie(i, end) = edges(g); i != end; ++i) {
      diff_t f_u = get(index_map, source(*i, g));
      diff_t f_v = get(index_map, target(*i, g));
      using namespace std; // to call abs() unqualified
      sum += abs(f_u - f_v);
    }
    return sum;
  }
  
} // namespace boost

#endif // BOOST_GRAPH_BANDWIDTH_HPP
