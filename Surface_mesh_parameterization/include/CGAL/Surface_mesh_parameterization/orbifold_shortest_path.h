// Copyright (c) 2016  GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_SURFACE_MESH_PARAMETERIZATION_SHORTEST_PATH_H
#define CGAL_SURFACE_MESH_PARAMETERIZATION_SHORTEST_PATH_H

#include <CGAL/license/Surface_mesh_parameterization.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/assertions.h>

#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/unordered_map.hpp>

#include <exception>
#include <list>
#include <utility>
#include <iostream>
#include <fstream>

namespace CGAL {

namespace Surface_mesh_parameterization {

namespace internal {

class Dijkstra_end_exception : public std::exception
{
  const char* what() const throw ()
  {
#ifdef CGAL_SMP_ORBIFOLD_DEBUG
    return "Dijkstra: reached the target vertex";
#else
    return "";
#endif
  }
};

template<typename TriangleMesh, typename Seam_container>
void output_shortest_paths_to_selection_file(const TriangleMesh& mesh,
                                             const Seam_container& seams,
                                             std::ofstream& os)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::edge_descriptor   edge_descriptor;

  boost::unordered_map<vertex_descriptor, int> index_map;

  int counter = 0;
  for(vertex_descriptor vd : vertices(mesh)) {
    index_map[vd] = counter++;
  }

  os << std::endl /* vertices */ << std::endl /* faces */;

  for(edge_descriptor ed : seams) {
    // could be made more efficient...
    os << index_map[source(ed, mesh)] << " " << index_map[target(ed, mesh)] << " ";
  }

  os << std::endl;
}

// Visitor to stop Dijkstra when a target turns 'BLACK' (the point has been examined
// through all its edges)
template<typename TriangleMesh>
class Stop_at_target_Dijkstra_visitor : boost::default_dijkstra_visitor
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor  vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::edge_descriptor    edge_descriptor;

  vertex_descriptor destination_vd;

public:
  Stop_at_target_Dijkstra_visitor(vertex_descriptor destination_vd)
    : destination_vd(destination_vd)
  { }

  void initialize_vertex(const vertex_descriptor& /*s*/, const TriangleMesh& /*mesh*/) const { }
  void examine_vertex(const vertex_descriptor& /*s*/, const TriangleMesh& /*mesh*/) const { }
  void examine_edge(const edge_descriptor& /*e*/, const TriangleMesh& /*mesh*/) const { }
  void edge_relaxed(const edge_descriptor& /*e*/, const TriangleMesh& /*mesh*/) const { }
  void discover_vertex(const vertex_descriptor& /*s*/, const TriangleMesh& /*mesh*/) const { }
  void edge_not_relaxed(const edge_descriptor& /*e*/, const TriangleMesh& /*mesh*/) const { }
  void finish_vertex(const vertex_descriptor &vd, const TriangleMesh& /* mesh*/) const
  {
    if(vd == destination_vd)
      throw Dijkstra_end_exception();
  }
};

} // namespace internal

/// \ingroup PkgSurfaceMeshParameterizationOrbifoldHelperFunctions
///
/// computes the shortest path between `source` and `target` over `mesh`, using
/// <a href="https://www.boost.org/doc/libs/release/libs/graph/doc/dijkstra_shortest_paths.html">
/// boost::dijkstra_shortest_paths()</a>.
///
/// \tparam TriangleMesh A triangle mesh, model of `FaceListGraph` and `HalfedgeListGraph`.
/// \tparam EdgeOutputIterator A model of `OutputIterator` with value type
///                            `boost::graph_traits<TriangleMesh>::%edge_descriptor`.
///
/// \param mesh the triangular mesh to be parameterized
/// \param source, target the extremities of the path to be computed
/// \param oi the output iterator
///
/// \pre `source` and `target` are vertices of `mesh`.
/// \pre `source != target`
template<typename TriangleMesh, typename EdgeOutputIterator>
void compute_shortest_paths_between_two_cones(const TriangleMesh& mesh,
                                              typename boost::graph_traits<TriangleMesh>::vertex_descriptor source,
                                              typename boost::graph_traits<TriangleMesh>::vertex_descriptor target,
                                              EdgeOutputIterator oi)
{
  if(source == target) {
    std::cerr << "Warning: the source and target are identical in 'shortest_path' " << std::endl;
    return;
  }

  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor  vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::edge_descriptor    edge_descriptor;

  typedef internal::Stop_at_target_Dijkstra_visitor<TriangleMesh>        Stop_visitor;

  typedef boost::unordered_map<vertex_descriptor, vertex_descriptor>     Pred_umap;
  typedef boost::associative_property_map<Pred_umap>                     Pred_pmap;

  Pred_umap predecessor;
  Pred_pmap pred_pmap(predecessor);

  Stop_visitor vis(target);

  try {
    boost::dijkstra_shortest_paths(mesh, source, boost::predecessor_map(pred_pmap).visitor(vis));
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }

  // Draw the path from target to source and collect the edges along the way
  vertex_descriptor s, t = target;
  do {
    s = get(pred_pmap, t);
    std::pair<edge_descriptor, bool> e = edge(s, t, mesh);
    CGAL_assertion(e.second);
    if(e.second) {
      *oi++ = e.first;
    }
    t = s;
  } while (s != source);
}

/// \ingroup PkgSurfaceMeshParameterizationOrbifoldHelperFunctions
///
/// Given a range `[first; beyond[` of cones (described as vertex descriptors),
/// compute the shortest path for all pairs of consecutive entries in the range
/// and add them to the container `seams`.
///
/// \tparam TriangleMesh A triangle mesh, model of `FaceListGraph` and `HalfedgeListGraph`.
/// \tparam InputConesForwardIterator A model of `ForwardIterator` with value type
///                                   `boost::graph_traits<TriangleMesh>::%vertex_descriptor`.
/// \tparam SeamContainer A model of `SequenceContainer` with value type
///                       `boost::graph_traits<TriangleMesh>::%edge_descriptor`.
///
/// \param mesh the triangular mesh on which paths are computed
/// \param first, beyond a range of cones
/// \param seams a container that will store the paths, as a sequence of edges of the mesh.
///
/// \pre `std::distance(first,beyond) > 1`
template<typename TriangleMesh, typename InputConesForwardIterator, typename SeamContainer>
void compute_shortest_paths_between_cones(const TriangleMesh& mesh,
                                          InputConesForwardIterator first, InputConesForwardIterator beyond,
                                          SeamContainer& seams)
{
  CGAL_precondition(std::distance(first, beyond) == 3 || std::distance(first, beyond) == 4);
  InputConesForwardIterator last = --beyond;
  for(; first!=last; ++first) {
    InputConesForwardIterator next = first;
    ++next;
    compute_shortest_paths_between_two_cones(mesh, *first, *next, std::back_inserter(seams));
  }

  std::ofstream out("shortest_path.selection.txt");
#ifdef CGAL_SMP_ORBIFOLD_DEBUG
  internal::output_shortest_paths_to_selection_file(mesh, seams, out);
#endif
}

} // namespace Surface_mesh_parameterization

} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_SURFACE_MESH_PARAMETERIZATION_SHORTEST_PATH_H
