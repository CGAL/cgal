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
#ifndef BOOST_GRAPH_DETAIL_ADJ_LIST_EDGE_ITERATOR_HPP
#define BOOST_GRAPH_DETAIL_ADJ_LIST_EDGE_ITERATOR_HPP

#include <iterator>
#include <utility>
#include <boost/detail/workaround.hpp>

#if BOOST_WORKAROUND( __IBMCPP__, <= 600 )
#  define BOOST_GRAPH_NO_OPTIONAL
#endif

#ifdef BOOST_GRAPH_NO_OPTIONAL
#  define BOOST_GRAPH_MEMBER .
#else 
#  define BOOST_GRAPH_MEMBER ->
#  include <boost/optional.hpp>
#endif // ndef BOOST_GRAPH_NO_OPTIONAL

namespace boost {

  namespace detail {

    template <class VertexIterator, class OutEdgeIterator, class Graph>
    class adj_list_edge_iterator {
      typedef adj_list_edge_iterator self;
    public:
      typedef std::forward_iterator_tag iterator_category;
      typedef typename OutEdgeIterator::value_type value_type;
      typedef typename OutEdgeIterator::reference  reference;
      typedef typename OutEdgeIterator::pointer    pointer;
      typedef typename OutEdgeIterator::difference_type difference_type;
      typedef difference_type distance_type;

      inline adj_list_edge_iterator() {}

      inline adj_list_edge_iterator(const self& x) 
      : vBegin(x.vBegin), vCurr(x.vCurr), vEnd(x.vEnd),
        edges(x.edges), m_g(x.m_g) { }

      template <class G>
      inline adj_list_edge_iterator(VertexIterator b, 
                                    VertexIterator c,
                                    VertexIterator e,
                                    const G& g) 
        : vBegin(b), vCurr(c), vEnd(e), m_g(&g) {
        if ( vCurr != vEnd ) {
          while ( vCurr != vEnd && out_degree(*vCurr, *m_g) == 0 )
            ++vCurr;
          if ( vCurr != vEnd )
            edges = out_edges(*vCurr, *m_g);
        }
      }

      /*Note:
        In the directed graph cases, it is fine. 
        For undirected graphs, one edge go through twice.
      */
      inline self& operator++() {
        ++edges BOOST_GRAPH_MEMBER first;
        if (edges BOOST_GRAPH_MEMBER first == edges BOOST_GRAPH_MEMBER second) 
        {
          ++vCurr;
          while ( vCurr != vEnd && out_degree(*vCurr, *m_g) == 0 )
            ++vCurr;
          if ( vCurr != vEnd )
            edges = out_edges(*vCurr, *m_g);
        }
        return *this;
      }
      inline self operator++(int) {
        self tmp = *this;
        ++(*this);
        return tmp;
      }
      inline value_type operator*() const 
      { return *edges BOOST_GRAPH_MEMBER first; } 
      inline bool operator==(const self& x) const {
        return vCurr == x.vCurr 
          && (vCurr == vEnd 
              || edges BOOST_GRAPH_MEMBER first == x.edges BOOST_GRAPH_MEMBER first);
      }
      inline bool operator!=(const self& x) const {
        return vCurr != x.vCurr 
          || (vCurr != vEnd 
              && edges BOOST_GRAPH_MEMBER first != x.edges BOOST_GRAPH_MEMBER first);
      }
    protected:
      VertexIterator vBegin;
      VertexIterator vCurr;
      VertexIterator vEnd;

#ifdef BOOST_GRAPH_NO_OPTIONAL
      std::pair<OutEdgeIterator, OutEdgeIterator> edges;
#else
      boost::optional<std::pair<OutEdgeIterator, OutEdgeIterator> >
        edges;
#endif // ndef BOOST_GRAPH_NO_OPTIONAL
      const Graph* m_g;
    };

  } // namespace detail

}

#undef BOOST_GRAPH_MEMBER

#endif // BOOST_GRAPH_DETAIL_ADJ_LIST_EDGE_ITERATOR_HPP
