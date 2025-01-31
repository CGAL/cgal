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
  template<typename Mesh>
  class Stop_at_target_Dijkstra_visitor : boost::default_dijkstra_visitor
  {
    using vertex_descriptor = typename boost::graph_traits<Mesh>::vertex_descriptor;
    using edge_descriptor = typename boost::graph_traits<Mesh>::edge_descriptor;

  public:
    vertex_descriptor destination_vd;

    void initialize_vertex(const vertex_descriptor& /*s*/, const Mesh& /*mesh*/) const {}
    void examine_vertex(const vertex_descriptor& /*s*/, const Mesh& /*mesh*/) const {}
    void examine_edge(const edge_descriptor& /*e*/, const Mesh& /*mesh*/) const {}
    void edge_relaxed(const edge_descriptor& /*e*/, const Mesh& /*mesh*/) const {}
    void discover_vertex(const vertex_descriptor& /*s*/, const Mesh& /*mesh*/) const {}
    void edge_not_relaxed(const edge_descriptor& /*e*/, const Mesh& /*mesh*/) const {}
    void finish_vertex(const vertex_descriptor& vd, const Mesh& /* mesh*/) const
    {
      if (vd == destination_vd)
        throw Dijkstra_end_exception();
    }
  };
} // namespace internal

/*!
*@tparam Mesh a model of the concept `HalfedgeListGraph`
* @param vs source vertex
* @param vt target vertex
* @param mesh the mesh
* @param halfedge_sequence the sequence of halfedges that form the shortest path on `mesh`
*/
template<typename Mesh,
         typename OutputIterator,
         typename NamedParameters = parameters::Default_named_parameters>
OutputIterator shortest_path_between_two_vertices(
  const typename boost::graph_traits<Mesh>::vertex_descriptor vs,//source
  const typename boost::graph_traits<Mesh>::vertex_descriptor vt,//target
  const Mesh& mesh,
  OutputIterator halfedge_sequence_oit,
  const NamedParameters& np = parameters::default_values())
{
  using vertex_descriptor = typename boost::graph_traits<Mesh>::vertex_descriptor;
  using halfedge_descriptor = typename boost::graph_traits<Mesh>::halfedge_descriptor;

  using Pred_umap = std::unordered_map<vertex_descriptor, vertex_descriptor>;
  using Pred_pmap = boost::associative_property_map<Pred_umap>;

  Pred_umap predecessor;
  Pred_pmap pred_pmap(predecessor);

  internal::Stop_at_target_Dijkstra_visitor<Mesh> vis;
  vis.destination_vd = vt;
  try
  {
    boost::dijkstra_shortest_paths(mesh, vs,
      boost::predecessor_map(pred_pmap).visitor(vis));
  }
  catch (const std::exception& e){}

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
    const std::pair<halfedge_descriptor, bool>
      h = halfedge((path_it + 1)->vertex, path_it->vertex, mesh);
    if (h.second)
      *halfedge_sequence_oit++ = h.first;
  }
  return halfedge_sequence_oit;
}

} // namespace CGAL


#endif //CGAL_BOOST_GRAPH_SHORTEST_PATH_H
