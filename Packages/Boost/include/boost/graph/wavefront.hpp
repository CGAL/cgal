//
//=======================================================================
// Copyright 2002 Marc Wintermantel (wintermantel@imes.mavt.ethz.ch)
// ETH Zurich, Center of Structure Technologies (www.imes.ethz.ch/st)
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

#ifndef BOOST_GRAPH_WAVEFRONT_HPP
#define BOOST_GRAPH_WAVEFRONT_HPP

#include <boost/graph/graph_traits.hpp>
#include <boost/detail/numeric_traits.hpp>
#include <boost/graph/bandwidth.hpp>
#include <cmath>
#include <vector>

namespace boost {

  template <typename Graph, typename VertexIndexMap>
  typename graph_traits<Graph>::vertices_size_type
  ith_wavefront(typename graph_traits<Graph>::vertex_descriptor i,
                const Graph& g,
                VertexIndexMap index)
  {
    typename graph_traits<Graph>::vertex_descriptor v, w;
    typename graph_traits<Graph>::vertices_size_type b = 1;
    typename graph_traits<Graph>::out_edge_iterator edge_it2, edge_it2_end; 
    typename graph_traits<Graph>::vertices_size_type index_i = index[i];
    std::vector<bool> rows_active(num_vertices(g), false);

    rows_active[index_i] = true;
      
      typename graph_traits<Graph>::vertex_iterator ui, ui_end;
      for (tie(ui, ui_end) = vertices(g); ui != ui_end; ++ui)
      {
        v = *ui;
          if(index[v] <= index_i)
            {
              for (tie(edge_it2, edge_it2_end) = out_edges(v, g); edge_it2 != edge_it2_end; ++edge_it2)
              {
                w = target(*edge_it2, g);
                if( (index[w] >= index_i) && (!rows_active[index[w]]) )
                  {
                    b++;
                    rows_active[index[w]] = true;
                  }
              }
            }
      }
 
    return b;
  }


  template <typename Graph>
  typename graph_traits<Graph>::vertices_size_type
  ith_wavefront(typename graph_traits<Graph>::vertex_descriptor i,
                const Graph& g)
  {
    return ith_wavefront(i, g, get(vertex_index, g));
  }


  template <typename Graph, typename VertexIndexMap>
  typename graph_traits<Graph>::vertices_size_type
  max_wavefront(const Graph& g, VertexIndexMap index)
  {
    typename graph_traits<Graph>::vertices_size_type b = 0;
    typename graph_traits<Graph>::vertex_iterator i, end;
    for (tie(i, end) = vertices(g); i != end; ++i)
      b = std::max(b, ith_wavefront(*i, g, index));
    return b;
  }

  template <typename Graph>
  typename graph_traits<Graph>::vertices_size_type
  max_wavefront(const Graph& g)
  {
    return max_wavefront(g, get(vertex_index, g));
  }


  template <typename Graph, typename VertexIndexMap>
  double
  aver_wavefront(const Graph& g, VertexIndexMap index)
  {
    double b = 0;
    typename graph_traits<Graph>::vertex_iterator i, end;
    for (tie(i, end) = vertices(g); i != end; ++i)
      b += ith_wavefront(*i, g, index);

    b /= num_vertices(g);
    return b;
  }

  template <typename Graph>
  double
  aver_wavefront(const Graph& g)
  {
    return aver_wavefront(g, get(vertex_index, g));
  }


  template <typename Graph, typename VertexIndexMap>
  double
  rms_wavefront(const Graph& g, VertexIndexMap index)
  {
    double b = 0;
    typename graph_traits<Graph>::vertex_iterator i, end;
    for (tie(i, end) = vertices(g); i != end; ++i)
      b += std::pow(double ( ith_wavefront(*i, g, index) ), 2.0);

    b /= num_vertices(g);

    return std::sqrt(b);
  }

  template <typename Graph>
  double
  rms_wavefront(const Graph& g)
  {
    return rms_wavefront(g, get(vertex_index, g));
  }
 
  
} // namespace boost

#endif // BOOST_GRAPH_WAVEFRONT_HPP
