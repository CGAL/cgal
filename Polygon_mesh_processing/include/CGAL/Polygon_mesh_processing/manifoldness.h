// Copyright (c) 2015-2019 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Sebastien Loriot,
//                 Mael Rouxel-Labbé
//
#ifndef CGAL_POLYGON_MESH_PROCESSING_MANIFOLDNESS_H
#define CGAL_POLYGON_MESH_PROCESSING_MANIFOLDNESS_H

#include <CGAL/license/Polygon_mesh_processing/repair.h>

#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>

#include <CGAL/algorithm.h>
#include <CGAL/assertions.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Dynamic_property_map.h>
#include <CGAL/iterator.h>
#include <CGAL/property_map.h>

#include <iterator>
#include <utility>

namespace CGAL {
namespace Polygon_mesh_processing {

/// \ingroup PMP_repairing_grp
/// returns whether a vertex of a polygon mesh is non-manifold.
///
/// @tparam PolygonMesh a model of `HalfedgeListGraph`
///
/// @param v a vertex of `pm`
/// @param pm a triangle mesh containing `v`
///
/// \warning This function has linear runtime with respect to the size of the mesh.
///
/// \sa `duplicate_non_manifold_vertices()`
///
/// \return `true` if the vertex is non-manifold, `false` otherwise.
template <typename PolygonMesh>
bool is_non_manifold_vertex(typename boost::graph_traits<PolygonMesh>::vertex_descriptor v,
                            const PolygonMesh& pm)
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;

  typedef CGAL::dynamic_halfedge_property_t<bool>                                       Halfedge_property_tag;
  typedef typename boost::property_map<PolygonMesh, Halfedge_property_tag>::const_type  Visited_halfedge_map;

  // Dynamic pmaps do not have default initialization values (yet)
  Visited_halfedge_map visited_halfedges = get(Halfedge_property_tag(), pm);
  for(halfedge_descriptor h : halfedges(pm))
    put(visited_halfedges, h, false);

  std::size_t incident_null_faces_counter = 0;
  for(halfedge_descriptor h : halfedges_around_target(v, pm))
  {
    put(visited_halfedges, h, true);
    if(CGAL::is_border(h, pm))
      ++incident_null_faces_counter;
  }

  if(incident_null_faces_counter > 1)
  {
    // The vertex is the sole connection between two connected components --> non-manifold
    return true;
  }

  for(halfedge_descriptor h : halfedges(pm))
  {
    if(v == target(h, pm))
    {
      // Haven't seen that halfedge yet ==> more than one umbrella incident to 'v' ==> non-manifold
      if(!get(visited_halfedges, h))
        return true;
    }
  }

  return false;
}

/// \ingroup PMP_repairing_grp
/// collects the non-manifold vertices (if any) present in the mesh. A non-manifold vertex `v` is returned
/// via one incident halfedge `h` such that `target(h, pm) = v` for all the umbrellas that `v` appears in
/// (an <i>umbrella</i> being the set of faces incident to all the halfedges reachable by walking around `v`
/// using `hnext = prev(opposite(h, pm), pm)`, starting from `h`).
///
/// @tparam PolygonMesh a model of `HalfedgeListGraph`
/// @tparam OutputIterator a model of `OutputIterator` holding objects of type
///                         `boost::graph_traits<PolygonMesh>::%halfedge_descriptor`
///
/// @param pm a triangle mesh
/// @param out the output iterator that collects halfedges incident to `v`
///
/// \sa `is_non_manifold_vertex()`
/// \sa `duplicate_non_manifold_vertices()`
///
/// \return the output iterator.
template <typename PolygonMesh, typename OutputIterator>
OutputIterator non_manifold_vertices(const PolygonMesh& pm,
                                     OutputIterator out)
{
  // Non-manifoldness can appear either:
  // - if 'pm' is pinched at a vertex. While traversing the incoming halfedges at this vertex,
  //   we will meet strictly more than one border halfedge.
  // - if there are multiple umbrellas around a vertex. In that case, we will find a non-visited
  //   halfedge that has for target a vertex that is already visited.

  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor                  vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor                halfedge_descriptor;

  typedef CGAL::dynamic_vertex_property_t<halfedge_descriptor>                          Vertex_halfedge_tag;
  typedef typename boost::property_map<PolygonMesh, Vertex_halfedge_tag>::const_type    Visited_vertex_map;
  typedef CGAL::dynamic_halfedge_property_t<bool>                                       Halfedge_property_tag;
  typedef typename boost::property_map<PolygonMesh, Halfedge_property_tag>::const_type  Visited_halfedge_map;

  std::set<vertex_descriptor> known_nm_vertices;
  Visited_vertex_map visited_vertices = get(Vertex_halfedge_tag(), pm);
  Visited_halfedge_map visited_halfedges = get(Halfedge_property_tag(), pm);

  halfedge_descriptor null_h = boost::graph_traits<PolygonMesh>::null_halfedge();

  // Dynamic pmaps do not have default initialization values (yet).
  for(vertex_descriptor v : vertices(pm))
    put(visited_vertices, v, null_h);

  for(halfedge_descriptor h : halfedges(pm))
    put(visited_halfedges, h, false);

  for(halfedge_descriptor h : halfedges(pm))
  {
    // If 'h' is not visited yet, we walk around the target of 'h' and mark these
    // halfedges as visited. Thus, if we are here and the target is already marked as visited,
    // it means that the vertex is non manifold.
    if(!get(visited_halfedges, h))
    {
      put(visited_halfedges, h, true);

      bool is_non_manifold = false;
      const vertex_descriptor v = target(h, pm);

      if(get(visited_vertices, v) != null_h) // already seen this vertex, but not from this star
      {
        is_non_manifold = true;

        // If this is the second time we visit that vertex and the first star was manifold, we have
        // never reported the first star, but must do so now.
        if(known_nm_vertices.count(v) == 0)
          *out++ = get(visited_vertices, v); // that's a halfedge of the first star we've seen 'v' in
      }
      else
      {
        // First time we meet this vertex, only remember the halfedge we found the vertex with.
        put(visited_vertices, v, h);
      }

      // While walking the star of this halfedge, if we meet a border halfedge more than once,
      // it means the mesh is pinched and we are also in the case of a non-manifold situation.
      halfedge_descriptor ih = h, done = ih;
      int border_counter = 0;
      do
      {
        put(visited_halfedges, ih, true);
        if(is_border(ih, pm))
          ++border_counter;

        ih = prev(opposite(ih, pm), pm);
      }
      while(ih != done);

      if(border_counter > 1)
        is_non_manifold = true;

      if(is_non_manifold)
      {
        *out++ = h;
        known_nm_vertices.insert(v);
      }
    }
  }

  return out;
}

// pretty much the same as 'PMP::non_manifold_vertices()', but consider the geometry
// instead of the combinatorics ({combinatorial non-manifold vertices} C {geometrical non-manifold vertices})
template <typename PolygonMesh, typename NMPContainer, typename NMVM, typename CGAL_PMP_NP_TEMPLATE_PARAMETERS>
std::size_t geometrically_non_manifold_vertices(const PolygonMesh& pmesh,
                                                NMPContainer& nm_points, // m[vertex] = {halfedges}
                                                NMVM nm_marks,
                                                const CGAL_PMP_NP_CLASS& np)
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor                   vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor                 halfedge_descriptor;

  typedef typename GetGeomTraits<PolygonMesh, CGAL_PMP_NP_CLASS>::type                   Geom_traits;
//  Geom_traits gt = choose_parameter<Geom_traits>(get_parameter(np, internal_np::geom_traits));

  typedef typename GetVertexPointMap<PolygonMesh, CGAL_PMP_NP_CLASS>::const_type         VertexPointMap;
  VertexPointMap vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                                        get_const_property_map(vertex_point, pmesh));

  const bool verbose = choose_parameter(get_parameter(np, internal_np::verbose), false);

  typedef typename boost::property_traits<VertexPointMap>::reference                    Point_ref;
  typedef typename boost::property_traits<VertexPointMap>::value_type                   Point;

  typedef CGAL::dynamic_halfedge_property_t<bool>                                       Halfedge_property_tag;
  typedef typename boost::property_map<PolygonMesh, Halfedge_property_tag>::const_type  Visited_halfedge_map;

  CGAL_precondition(nm_points.empty());

  Visited_halfedge_map visited_halfedges = get(Halfedge_property_tag(), pmesh);
  std::unordered_map<Point, halfedge_descriptor> visited_points;

  for(halfedge_descriptor h : halfedges(pmesh))
  {
    // If 'h' is not visited yet, we walk around the target of 'h' and mark these
    // halfedges as visited. Thus, if we are here and the target is already marked as visited,
    // it means that the vertex is non manifold.
    if(get(visited_halfedges, h))
      continue;

    put(visited_halfedges, h, true);

    bool is_non_manifold_due_to_multiple_umbrellas = false;
    const vertex_descriptor v = target(h, pmesh);
    const Point_ref p = get(vpm, v);

    const auto visited_itb = visited_points.emplace(p, h);
    if(!visited_itb.second) // already seen this point, but not from this star
    {
      is_non_manifold_due_to_multiple_umbrellas = true;
      put(nm_marks, v, true);

      if(verbose)
        std::cout << p << std::endl;

      // if this is the second time we visit that vertex and the first star was manifold, we have
      // not marked the vertex as non-manifold from the first star
      const auto nm_itb = nm_points.emplace(p, std::set<halfedge_descriptor>{h});
      if(nm_itb.second) // successful insertion
      {
        const halfedge_descriptor h_from_another_star = visited_itb.first->second;
        nm_itb.first->second.insert(h_from_another_star);
        put(nm_marks, target(h_from_another_star, pmesh), true);
      }
      else
      {
        nm_itb.first->second.insert(h);
      }
    }

    // While walking the star of this halfedge, if we meet a border halfedge more than once,
    // it means the mesh is pinched and we are also in the case of a non-manifold situation
    halfedge_descriptor ih = h, done = ih;
    int border_counter = 0;
    do
    {
      put(visited_halfedges, ih, true);
      if(is_border(ih, pmesh))
        ++border_counter;

      ih = prev(opposite(ih, pmesh), pmesh);
    }
    while(ih != done);

    if(border_counter > 1 && !is_non_manifold_due_to_multiple_umbrellas)
    {
      put(nm_marks, v, true);
      nm_points[p].insert(h); // might or might not have been an empty vector before

      if(verbose)
        std::cout << "pinched star centered at " << p << std::endl;
    }
  }

  return nm_points.size();
}

template <typename PolygonMesh, typename NMVM, typename CGAL_PMP_NP_TEMPLATE_PARAMETERS>
std::size_t geometrically_non_manifold_vertices(const PolygonMesh& pmesh,
                                                NMVM& nm_points,
                                                const CGAL_PMP_NP_CLASS& np)
{
  typedef CGAL::dynamic_vertex_property_t<bool>                               Mark;
  typedef typename boost::property_map<PolygonMesh, Mark>::const_type         Marked_vertices;
  Marked_vertices nm_marks = get(Mark(), pmesh);
  for(auto v : vertices(pmesh))
    put(nm_marks, v, false);

  return geometrically_non_manifold_vertices(pmesh, nm_points, nm_marks, np);
}

template <typename PolygonMesh, typename NMVM>
std::size_t geometrically_non_manifold_vertices(const PolygonMesh& pmesh,
                                                NMVM& nm_points)
{
  return geometrically_non_manifold_vertices(pmesh, nm_points, parameters::all_default());
}

template <typename PolygonMesh, typename CGAL_PMP_NP_TEMPLATE_PARAMETERS>
std::size_t geometrically_non_manifold_vertices(const PolygonMesh& pmesh,
                                                const CGAL_PMP_NP_CLASS& np)
{
  typedef typename property_map_selector<PolygonMesh, boost::vertex_point_t>::const_type VPM;
  typedef typename boost::property_traits<VPM>::value_type                               Point;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor                 halfedge_descriptor;

  std::unordered_map<Point, std::set<halfedge_descriptor> > nm_points;
  return geometrically_non_manifold_vertices(pmesh, nm_points, np);
}

template <typename PolygonMesh>
std::size_t geometrically_non_manifold_vertices(const PolygonMesh& pmesh)
{
  return geometrically_non_manifold_vertices(pmesh, parameters::all_default());
}

} // namespace Polygon_mesh_processing
} // namespace CGAL

// for backward compatibility, and the functions above are required 'repair_manifoldness.h'
#include <CGAL/Polygon_mesh_processing/repair_manifoldness.h>

#endif // CGAL_POLYGON_MESH_PROCESSING_MANIFOLDNESS_H
