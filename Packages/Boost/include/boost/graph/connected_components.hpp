//
//=======================================================================
// Copyright 1997-2001 University of Notre Dame.
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
#ifndef BOOST_GRAPH_CONNECTED_COMPONENTS_HPP
#define BOOST_GRAPH_CONNECTED_COMPONENTS_HPP

#include <boost/config.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_concepts.hpp>

#include <boost/static_assert.hpp>

namespace boost {

  namespace detail {

    // This visitor is used both in the connected_components algorithm
    // and in the kosaraju strong components algorithm during the
    // second DFS traversal.
    template <class ComponentsMap>
    class components_recorder : public dfs_visitor<>
    {
      typedef typename property_traits<ComponentsMap>::value_type comp_type;
    public:
      components_recorder(ComponentsMap c, 
                          comp_type& c_count)
        : m_component(c), m_count(c_count) {}

      template <class Vertex, class Graph>
      void start_vertex(Vertex, Graph&) {
        if (m_count == std::numeric_limits<comp_type>::max())
          m_count = 0; // start counting components at zero
        else
          ++m_count;
      }
      template <class Vertex, class Graph>
      void discover_vertex(Vertex u, Graph&) {
        put(m_component, u, m_count);
      }
    protected:
      ComponentsMap m_component;
      comp_type& m_count;
    };

  } // namespace detail

  // This function computes the connected components of an undirected
  // graph using a single application of depth first search.

  template <class Graph, class ComponentMap, class P, class T, class R>
  inline typename property_traits<ComponentMap>::value_type
  connected_components(const Graph& g, ComponentMap c, 
                       const bgl_named_params<P, T, R>& params)
  {
    typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
    function_requires< WritablePropertyMapConcept<ComponentMap, Vertex> >();
    typedef typename boost::graph_traits<Graph>::directed_category directed;
    BOOST_STATIC_ASSERT((boost::is_same<directed, undirected_tag>::value));

    typedef typename property_traits<ComponentMap>::value_type comp_type;
    // c_count initialized to "nil" (with nil represented by max())
    comp_type c_count(std::numeric_limits<comp_type>::max());
    detail::components_recorder<ComponentMap> vis(c, c_count);
    depth_first_search(g, params.visitor(vis));
    return c_count + 1;
  }

  template <class Graph, class ComponentMap>
  inline typename property_traits<ComponentMap>::value_type
  connected_components(const Graph& g, ComponentMap c)
  {
    typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
    function_requires< WritablePropertyMapConcept<ComponentMap, Vertex> >();
    typedef typename boost::graph_traits<Graph>::directed_category directed;
    BOOST_STATIC_ASSERT((boost::is_same<directed, undirected_tag>::value));

    typedef typename property_traits<ComponentMap>::value_type comp_type;
    // c_count initialized to "nil" (with nil represented by max())
    comp_type c_count(std::numeric_limits<comp_type>::max());
    detail::components_recorder<ComponentMap> vis(c, c_count);
    depth_first_search(g, visitor(vis));
    return c_count + 1;
  }

  
} // namespace boost


#endif // BOOST_GRAPH_CONNECTED_COMPONENTS_HPP
