//=======================================================================
// Copyright 2002 Indiana University.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
//
// This file is part of the Boost Graph Library
//
// You should have received a copy of the License Agreement for the
// Boost Graph Library along with the software; see the file LICENSE.
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

#ifndef BOOST_ADJACENCY_ITERATOR_HPP
#define BOOST_ADJACENCY_ITERATOR_HPP

#include <boost/iterator/iterator_adaptor.hpp>
#include <boost/graph/graph_traits.hpp>

namespace boost
{

  template <class Graph, class Vertex, class OutEdgeIter, class Difference>
  struct adjacency_iterator
    : iterator_adaptor<
          adjacency_iterator<Graph,Vertex,OutEdgeIter,Difference>
        , OutEdgeIter
        , Vertex
        , use_default
        , Vertex
        , Difference
      >
  {
      typedef iterator_adaptor<
          adjacency_iterator<Graph,Vertex,OutEdgeIter,Difference>
        , OutEdgeIter
        , Vertex
        , use_default
        , Vertex
        , Difference
      > super_t;
      
      inline adjacency_iterator() {}
      inline adjacency_iterator(OutEdgeIter const& i, const Graph* g) : super_t(i), m_g(g) { }

      inline Vertex
      dereference() const
        { return target(*this->base(), *m_g); }

      const Graph* m_g;
  };

  template <class Graph,
            class Vertex = typename graph_traits<Graph>::vertex_descriptor,
            class OutEdgeIter=typename graph_traits<Graph>::out_edge_iterator>
  class adjacency_iterator_generator
  {
    typedef typename boost::detail::iterator_traits<OutEdgeIter>
      ::difference_type difference_type;
  public:
      typedef adjacency_iterator<Graph,Vertex,OutEdgeIter,difference_type> type;
  };

  template <class Graph, class Vertex, class InEdgeIter, class Difference>
  struct inv_adjacency_iterator
    : iterator_adaptor<
          inv_adjacency_iterator<Graph,Vertex,InEdgeIter,Difference>
        , InEdgeIter
        , Vertex
        , use_default
        , Vertex
        , Difference
      >
    {
      typedef iterator_adaptor<
                  inv_adjacency_iterator<Graph,Vertex,InEdgeIter,Difference>
                , InEdgeIter
                , Vertex
                , use_default
                , Vertex
                , Difference
              > super_t;

      inline inv_adjacency_iterator() { }
      inline inv_adjacency_iterator(InEdgeIter const& i, const Graph* g) : super_t(i), m_g(g) { }

      inline Vertex
      dereference() const
        { return source(*this->base(), *m_g); }

      const Graph* m_g;
    };

  template <class Graph,
            class Vertex = typename graph_traits<Graph>::vertex_descriptor,
            class InEdgeIter = typename graph_traits<Graph>::in_edge_iterator>
  class inv_adjacency_iterator_generator {
    typedef typename boost::detail::iterator_traits<InEdgeIter>
      ::difference_type difference_type;
  public:
      typedef inv_adjacency_iterator<Graph, Vertex, InEdgeIter, difference_type> type;
  };

} // namespace boost

#endif // BOOST_DETAIL_ADJACENCY_ITERATOR_HPP
