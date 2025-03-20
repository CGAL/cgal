// Copyright (c) 2025  GeometryFactory (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Jane Tournois, Andreas Fabri
//

#ifndef CGAL_BOOST_GRAPH_DIJKSTRA_SHORTEST_PATH_H
#define CGAL_BOOST_GRAPH_DIJKSTRA_SHORTEST_PATH_H

#include <boost/graph/dijkstra_shortest_paths.hpp>

#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <boost/property_map/property_map.hpp>
#include <CGAL/boost/graph/properties.h>

#include <list>
#include <unordered_map>

namespace CGAL {
namespace internal {

  /// An exception used while catching a throw that stops Dijkstra's algorithm
  /// once the shortest path to a target has been found.
  class Dijkstra_end_exception : public std::exception
  {
    const char* what() const throw ()
    {
      return "Dijkstra shortest path: reached the target vertex.";
    }
  };

  /// Visitor to stop Dijkstra's algorithm once the given target turns 'BLACK',
  /// that is when the target has been examined through all its incident edges and
  /// the shortest path is thus known.
  template<typename Graph, typename VertexEdgeMap>
  class Stop_at_target_Dijkstra_visitor : public boost::default_dijkstra_visitor
  {
    using vertex_descriptor = typename boost::graph_traits<Graph>::vertex_descriptor;
    using edge_descriptor = typename boost::graph_traits<Graph>::edge_descriptor;

  public:
    vertex_descriptor destination_vd;
    VertexEdgeMap& relaxed_edges;

    Stop_at_target_Dijkstra_visitor(vertex_descriptor destination_vd,
                                    VertexEdgeMap& relaxed_edges)
      : destination_vd(destination_vd), relaxed_edges(relaxed_edges)
    {}


    void edge_relaxed(const edge_descriptor& e, const Graph& g) const
    {
      relaxed_edges[target(e, g)] = e;
    }

    void finish_vertex(const vertex_descriptor& vd, const Graph& /* g*/) const
    {
      if (vd == destination_vd)
        throw Dijkstra_end_exception();
    }
  };
} // namespace internal

/*!
* \ingroup PkgBGLTraversal
* computes the shortest path between two vertices in a graph `g`, where the vertices must belong to the same connected component of `g`.
*
* @tparam Graph a model of the concept `HalfedgeListGraph`
* @tparam OutputIterator an output iterator with value type `boost::graph_traits<Graph>::%halfedge_descriptor`
* @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
*
* @param vs source vertex
* @param vt target vertex
* @param g the graph
* @param halfedge_sequence_oit the output iterator holding the output sequence
* of halfedges that form the shortest path from `vs` to `vt` on `g`
* @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*
*  \cgalNamedParamsBegin
*   \cgalParamNBegin{edge_weight_map}
*    \cgalParamDescription{a property map associating to each edge in the graph its weight or ``length''.
*                          The weights must all be non-negative.}
*    \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<PolygonMesh>::%edge_descriptor`
*                   as key type and a value type which as specified in the named parameter `distance_map`of the function
                    <A href="https://www.boost.org/doc/libs/release/libs/graph/doc/dijkstra_shortest_paths.html">`boost::graph::dijkstra_shortest_paths()`</A>,
                    with any model of `RingNumberType` fulfilling the requirements. }
*    \cgalParamDefault{`get(boost::edge_weight, mesh)`}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{vertex_index_map}
*     \cgalParamDescription{a property map associating to each vertex of `g` a unique index between `0` and `num_vertices(g) - 1`}
*     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
*                    as key type and `std::size_t` as value type}
*     \cgalParamDefault{an automatically indexed internal map}
*   \cgalParamNEnd
*  \cgalNamedParamsEnd
*/
template<typename Graph,
         typename OutputIterator,
         typename NamedParameters = parameters::Default_named_parameters>
OutputIterator dijkstra_shortest_path(
  const typename boost::graph_traits<Graph>::vertex_descriptor vs,//source
  const typename boost::graph_traits<Graph>::vertex_descriptor vt,//target
  const Graph& g,
  OutputIterator halfedge_sequence_oit,
  const NamedParameters& np = parameters::default_values())
{
  using vertex_descriptor = typename boost::graph_traits<Graph>::vertex_descriptor;
  using halfedge_descriptor = typename boost::graph_traits<Graph>::halfedge_descriptor;
  using edge_descriptor = typename boost::graph_traits<Graph>::edge_descriptor;

  using Pred_umap = std::unordered_map<vertex_descriptor, vertex_descriptor>;
  using Pred_pmap = boost::associative_property_map<Pred_umap>;

  using parameters::get_parameter;
  using parameters::choose_parameter;

  const auto w_map = choose_parameter(get_parameter(np, internal_np::edge_weight),
                                      get(boost::edge_weight, g));
  const auto vim = get_initialized_vertex_index_map(g, np);

  Pred_umap predecessor;
  Pred_pmap pred_pmap(predecessor);

  using VEMap = std::unordered_map<vertex_descriptor, edge_descriptor>;
  VEMap relaxed_edges_map;
  internal::Stop_at_target_Dijkstra_visitor<Graph, VEMap> vis(vt, relaxed_edges_map);
  try
  {
    boost::dijkstra_shortest_paths(g, vs,
      boost::predecessor_map(pred_pmap)
      .visitor(vis)
      .weight_map(w_map)
      .vertex_index_map(vim));
  }
  catch (const internal::Dijkstra_end_exception& ){}

  std::list<halfedge_descriptor> path;
  vertex_descriptor t = vt;
  do {
    path.push_front(halfedge(relaxed_edges_map[t],g));
    t = get(pred_pmap, t);
  }while (t != vs);
  for(auto he : path){
    *halfedge_sequence_oit++ = he;
  }
  return halfedge_sequence_oit;
}

} // namespace CGAL


#endif //CGAL_BOOST_GRAPH_DIJKSTRA_SHORTEST_PATH_H
