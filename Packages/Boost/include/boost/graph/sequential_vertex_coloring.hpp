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
#ifndef BOOST_GRAPH_SEQUENTIAL_VERTEX_COLORING_HPP
#define BOOST_GRAPH_SEQUENTIAL_VERTEX_COLORING_HPP

#include <vector>
#include <boost/graph/graph_traits.hpp>

/* This algorithm is to find coloring of a graph

   Algorithm: 
   Let G = (V,E) be a graph with vertices (somehow) ordered v_1, v_2, ...,
   v_n. For k = 1, 2, ..., n the sequential algorithm assigns v_k to the
   smallest possible color. 

   Reference:

   Thomas F. Coleman and Jorge J. More, Estimation of sparse Jacobian
   matrices and graph coloring problems. J. Numer. Anal. V20, P187-209, 1983

   v_k is stored as o[k] here. 

   The color of the vertex v will be stored in color[v].
   i.e., vertex v belongs to coloring color[v] */

namespace boost {
  template <class VertexListGraph, class OrderPA, class ColorMap>
  typename graph_traits<VertexListGraph>::size_type
  sequential_vertex_coloring(const VertexListGraph& G, OrderPA order, 
                             ColorMap color)
  {
    using graph_traits;
    using boost::tie;
    typedef graph_traits<VertexListGraph> GraphTraits;
    typedef typename GraphTraits::vertex_descriptor Vertex;
    typedef typename GraphTraits::size_type size_type;
    
    size_type max_color = 0;
    const size_type V = num_vertices(G);

    // We need to keep track of which colors are used by
    // adjacent vertices. We do this by marking the colors
    // that are used. The mark array contains the mark
    // for each color. The length of mark is the
    // number of vertices since the maximum possible number of colors
    // is the number of vertices.
    std::vector<size_type> mark(V, numeric_limits_max(max_color));
    
    //Initialize colors 
    typename GraphTraits::vertex_iterator v, vend;
    for (tie(v, vend) = vertices(G); v != vend; ++v)
      put(color, *v, V-1);
    
    //Determine the color for every vertex one by one
    for ( size_type i = 0; i < V; i++) {
      Vertex current = get(order,i);
      typename GraphTraits::adjacency_iterator v, vend;
      
      //Mark the colors of vertices adjacent to current.
      //i can be the value for marking since i increases successively
      for (tie(v,vend) = adjacent_vertices(current, G); v != vend; ++v)
        mark[get(color,*v)] = i; 
      
      //Next step is to assign the smallest un-marked color
      //to the current vertex.
      size_type j = 0;

      //Scan through all useable colors, find the smallest possible
      //color which is not used by neighbors.  Note that if mark[j]
      //is equal to i, color j is used by one of the current vertex's
      //neighbors.
      while ( j < max_color && mark[j] == i ) 
        ++j;
      
      if ( j == max_color )  //All colors are used up. Add one more color
        ++max_color;

      //At this point, j is the smallest possible color
      put(color, current, j);  //Save the color of vertex current
      
    }
    
    return max_color;
  }
}

#endif
