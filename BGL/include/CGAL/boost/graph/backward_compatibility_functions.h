// Copyright (c) 2013  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_BOOST_GRAPH_BACKWARD_COMPATIBILITY_FUNCTIONS_H
#define CGAL_BOOST_GRAPH_BACKWARD_COMPATIBILITY_FUNCTIONS_H

namespace CGAL {

  template<class Graph>
  typename boost::graph_traits<Graph>::edge_descriptor
  opposite_edge(typename boost::graph_traits<Graph>::edge_descriptor e
                , const Graph& g)
  {
    typename boost::graph_traits<Graph>::halfedge_descriptor h = halfedge(e, g);
    return edge(opposite(h,g), g);
  }


  template<class Graph>
  typename boost::graph_traits<Graph>::edge_descriptor
  next_edge(typename boost::graph_traits<Graph>::edge_descriptor e
            , const Graph& g)
  {
    typename boost::graph_traits<Graph>::halfedge_descriptor h = halfedge(e, g);
    return edge(next(h, g), g);
}


  template<class Graph>
  typename boost::graph_traits<Graph>::edge_descriptor
  next_edge_cw(typename boost::graph_traits<Graph>::edge_descriptor e
               , const Graph& g)
  {
    typename boost::graph_traits<Graph>::halfedge_descriptor h = halfedge(e, g);
    return edge(opposite(next(h, g), g), g);
  }


  template<class Graph>
  typename boost::graph_traits<Graph>::edge_descriptor
  next_edge_ccw(typename boost::graph_traits<Graph>::edge_descriptor e
                , const Graph& g)
  {
    typename boost::graph_traits<Graph>::halfedge_descriptor h = halfedge(e, g);

    return edge(prev(opposite(h, g), g), g);
  }

  template <class Graph>
  struct halfedge_graph_traits;

   template<class Graph>
   std::pair<typename CGAL::halfedge_graph_traits<Graph>::undirected_edge_iterator, typename CGAL::halfedge_graph_traits<Graph>::undirected_edge_iterator>
   undirected_edges(const Graph& g)
   {
     return edges(g);
   }

}  //end of namespace CGAL

#endif //CGAL_BOOST_GRAPH_BACKWARD_COMPATIBILITY_FUNCTIONS_H
