// Copyright (c) 2016  GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_SURFACE_MESH_PARAMETERIZATION_INTERNAL_SHORTEST_PATH_H
#define CGAL_SURFACE_MESH_PARAMETERIZATION_INTERNAL_SHORTEST_PATH_H

#include <CGAL/license/Surface_mesh_parameterization.h>

#include <CGAL/assertions.h>

#include <boost/foreach.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/unordered_map.hpp>

#include <exception>
#include <list>
#include <utility>

namespace CGAL {

namespace Surface_mesh_parameterization {

namespace internal {

class Dijkstra_end_exception : public std::exception
{
  const char* what() const throw ()
  {
    return "Dijkstra: reached the target vertex";
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
  BOOST_FOREACH(vertex_descriptor vd, vertices(mesh)) {
    index_map[vd] = counter++;
  }

  os << std::endl /* vertices */ << std::endl /* faces */;

  BOOST_FOREACH(edge_descriptor ed, seams) {
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

template<typename TriangleMesh, typename OutputIterator>
void compute_shortest_paths_between_two_cones(const TriangleMesh& mesh,
           typename boost::graph_traits<TriangleMesh>::vertex_descriptor source,
           typename boost::graph_traits<TriangleMesh>::vertex_descriptor target,
           OutputIterator oi)
{
  if(source == target) {
    std::cout << "Warning: the source and target are identical in 'shortest_path' " << std::endl;
    return;
  }

  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor  vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::edge_descriptor    edge_descriptor;

  typedef Stop_at_target_Dijkstra_visitor<TriangleMesh>                  Stop_visitor;

  typedef boost::unordered_map<vertex_descriptor, vertex_descriptor>     Pred_umap;
  typedef boost::associative_property_map<Pred_umap>                     Pred_pmap;

  Pred_umap predecessor;
  Pred_pmap pred_pmap(predecessor);

  Stop_visitor vis(target);

  try {
    boost::dijkstra_shortest_paths(mesh, source, boost::predecessor_map(pred_pmap).visitor(vis));
  } catch (const std::exception& e) {
    std::cout << e.what() << std::endl;
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

template<typename TriangleMesh, typename Cones_vector, typename Seam_container>
void compute_shortest_paths_between_cones(const TriangleMesh& mesh,
                                          const Cones_vector& cones,
                                          Seam_container& seams)
{
  CGAL_precondition(cones.size() == 3 || cones.size() == 4);
  for(std::size_t i=0; i<cones.size() - 1; ++i) {
    compute_shortest_paths_between_two_cones(mesh, cones[i], cones[i+1],
        std::back_inserter(seams));
  }

  std::ofstream out("shortest_path.selection.txt");
  output_shortest_paths_to_selection_file(mesh, seams, out);
}

} // namespace internal

} // namespace Surface_mesh_parameterization

} // namespace CGAL

#endif // CGAL_SURFACE_MESH_PARAMETERIZATION_INTERNAL_SHORTEST_PATH_H
