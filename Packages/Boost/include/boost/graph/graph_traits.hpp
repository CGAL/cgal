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

#ifndef BOOST_GRAPH_TRAITS_HPP
#define BOOST_GRAPH_TRAITS_HPP

#include <boost/config.hpp>
#include <iterator>
#include <boost/tuple/tuple.hpp>
#include <boost/pending/ct_if.hpp>
#include <boost/iterator/iterator_categories.hpp>
#include <boost/iterator/iterator_adaptor.hpp>
#include <boost/detail/workaround.hpp>

namespace boost {

  template <typename G>
  struct graph_traits {
    typedef typename G::vertex_descriptor      vertex_descriptor;
    typedef typename G::edge_descriptor        edge_descriptor;
    typedef typename G::adjacency_iterator     adjacency_iterator;
    typedef typename G::out_edge_iterator      out_edge_iterator;
    typedef typename G::in_edge_iterator       in_edge_iterator;
    typedef typename G::vertex_iterator        vertex_iterator;
    typedef typename G::edge_iterator          edge_iterator;

    typedef typename G::directed_category      directed_category;
    typedef typename G::edge_parallel_category edge_parallel_category;
    typedef typename G::traversal_category     traversal_category;

    typedef typename G::vertices_size_type     vertices_size_type;
    typedef typename G::edges_size_type        edges_size_type;
    typedef typename G::degree_size_type       degree_size_type;

    static inline vertex_descriptor null_vertex();
  };

  template <typename G>
  inline typename graph_traits<G>::vertex_descriptor
  graph_traits<G>::null_vertex()
  {
    return G::null_vertex();
  }

  // directed_category tags
  struct directed_tag { };
  struct undirected_tag { };
  struct bidirectional_tag : public directed_tag { };

  namespace detail {
    inline bool is_directed(directed_tag) { return true; }
    inline bool is_directed(undirected_tag) { return false; }
  }

  template <typename Graph>
  bool is_directed(const Graph&) { 
    typedef typename graph_traits<Graph>::directed_category Cat;
    return detail::is_directed(Cat());
  }
  template <typename Graph>
  bool is_undirected(const Graph& g) { 
    return ! is_directed(g);
  }

  // edge_parallel_category tags
  struct allow_parallel_edge_tag {};
  struct disallow_parallel_edge_tag {};

  namespace detail {
    inline bool allows_parallel(allow_parallel_edge_tag) { return true; }
    inline bool allows_parallel(disallow_parallel_edge_tag) { return false; }
  }

  template <typename Graph>
  bool allows_parallel_edges(const Graph&) { 
    typedef typename graph_traits<Graph>::edge_parallel_category Cat;
    return detail::allows_parallel(Cat());
  }

  // traversal_category tags
  struct incidence_graph_tag { };
  struct adjacency_graph_tag { };
  struct bidirectional_graph_tag : 
    public virtual incidence_graph_tag { };
  struct vertex_list_graph_tag { };
  struct edge_list_graph_tag { };
  struct adjacency_matrix_tag { };

  //?? not the right place ?? Lee
  typedef boost::forward_traversal_tag multi_pass_input_iterator_tag;

  template <typename G>
  struct edge_property_type {
    typedef typename G::edge_property_type type;
  };
  template <typename G>
  struct vertex_property_type {
    typedef typename G::vertex_property_type type;
  };
  template <typename G>
  struct graph_property_type {
    typedef typename G::graph_property_type type;
  };

} // namespace boost

// Since pair is in namespace std, Koenig lookup will find source and
// target if they are also defined in namespace std.  This is illegal,
// but the alternative is to put source and target in the global
// namespace which causes name conflicts with other libraries (like
// SUIF).
namespace std {

  /* Some helper functions for dealing with pairs as edges */
  template <class T, class G>
  T source(pair<T,T> p, const G&) { return p.first; }

  template <class T, class G>
  T target(pair<T,T> p, const G&) { return p.second; }

}

#if defined(__GNUC__) && defined(__SGI_STL_PORT)
// For some reason g++ with STLport does not see the above definition
// of source() and target() unless we bring them into the boost
// namespace.
namespace boost {
  using std::source;
  using std::target;
}
#endif

#endif // BOOST_GRAPH_TRAITS_HPP
