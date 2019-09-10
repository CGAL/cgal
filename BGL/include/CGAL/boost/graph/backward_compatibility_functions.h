// Copyright (c) 2013  GeometryFactory (France).  All rights reserved.
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
// SPDX-License-Identifier: LGPL-3.0+
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
