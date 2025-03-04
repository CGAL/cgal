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

#ifndef CGAL_BOOST_GRAPH_SHORTEST_PATH_H
#define CGAL_BOOST_GRAPH_SHORTEST_PATH_H

#include <boost/graph/dijkstra_shortest_paths.hpp>

#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <boost/property_map/property_map.hpp>
#include <CGAL/boost/graph/properties.h>

#include <vector>
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
  class Stop_at_target_Dijkstra_visitor : boost::default_dijkstra_visitor
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

    void initialize_vertex(const vertex_descriptor& /*s*/, const Graph& /*g*/) const {}
    void examine_vertex(const vertex_descriptor& /*s*/, const Graph& /*g*/) const {}
    void examine_edge(const edge_descriptor& /*e*/, const Graph& /*g*/) const {}
    void edge_relaxed(const edge_descriptor& e, const Graph& g) const
    {
      relaxed_edges[target(e, g)] = e;
    }
    void discover_vertex(const vertex_descriptor& /*s*/, const Graph& /*g*/) const {}
    void edge_not_relaxed(const edge_descriptor& /*e*/, const Graph& /*g*/) const {}
    void finish_vertex(const vertex_descriptor& vd, const Graph& /* g*/) const
    {
      if (vd == destination_vd)
        throw Dijkstra_end_exception();
    }
  };
} // namespace internal

/*!
* \ingroup PkgBGLTraversal
* Computes the shortest path between two vertices in a graph `g`
*
*@tparam Graph a model of the concept `HalfedgeListGraph`
* @param vs source vertex
* @param vt target vertex
* @param g the graph
* @param halfedge_sequence_oit the output iterator holding the output sequence
* of halfedges that form the shortest path from `vs` to `vt` on `g`
* @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
* ** edge_weight_map : @todo deal with input mesh with no internal pmap.
*                       default in boost is `get(boost::edge_weight, mesh)`
* ** vertex_index_map : @todo deal with input mesh with no internal pmap.
*
*/
template<typename Graph,
         typename OutputIterator,
         typename NamedParameters = parameters::Default_named_parameters>
OutputIterator shortest_path_between_two_vertices(
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
      .weight_map(w_map));
  }
  catch (const internal::Dijkstra_end_exception& ){}

  // Walk back from target to source and collect vertices along the way
  struct vertex_on_path
  {
    vertex_descriptor vertex;
    bool is_constrained;
  };

  std::vector<vertex_descriptor> constrained_vertices = { vs };
  vertex_descriptor t = vt;
  std::vector<vertex_on_path> path;
  do
  {
    const bool is_new_vertex = (constrained_vertices.end()
      == std::find(constrained_vertices.begin(), constrained_vertices.end(), t));

    vertex_on_path vop;
    vop.vertex = t;
    vop.is_constrained = !is_new_vertex;
    path.push_back(vop);

    t = get(pred_pmap, t);
  }
  while (t != vs);

  // Add the last vertex
  vertex_on_path vop;
  vop.vertex = constrained_vertices.back();
  vop.is_constrained = true;
  path.push_back(vop);

  // Display path
  for (auto path_it = path.begin(); path_it != path.end() - 1; ++path_it)
  {
    const auto map_it = vis.relaxed_edges.find(path_it->vertex);
    if (map_it != vis.relaxed_edges.end())
      *halfedge_sequence_oit++ = halfedge(map_it->second, g);
  }
  return halfedge_sequence_oit;
}

} // namespace CGAL


#endif //CGAL_BOOST_GRAPH_SHORTEST_PATH_H
