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
/// checks whether a vertex of a polygon mesh is non-manifold.
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
/// via one incident halfedge `h` such that `target(h, pm) = v` for all the umbrellas that `v` apppears in
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

  typedef CGAL::dynamic_vertex_property_t<bool>                                         Vertex_bool_tag;
  typedef typename boost::property_map<PolygonMesh, Vertex_bool_tag>::const_type        Known_manifold_vertex_map;
  typedef CGAL::dynamic_vertex_property_t<halfedge_descriptor>                          Vertex_halfedge_tag;
  typedef typename boost::property_map<PolygonMesh, Vertex_halfedge_tag>::const_type    Visited_vertex_map;
  typedef CGAL::dynamic_halfedge_property_t<bool>                                       Halfedge_property_tag;
  typedef typename boost::property_map<PolygonMesh, Halfedge_property_tag>::const_type  Visited_halfedge_map;

  Known_manifold_vertex_map known_nm_vertices = get(Vertex_bool_tag(), pm);
  Visited_vertex_map visited_vertices = get(Vertex_halfedge_tag(), pm);
  Visited_halfedge_map visited_halfedges = get(Halfedge_property_tag(), pm);

  halfedge_descriptor null_h = boost::graph_traits<PolygonMesh>::null_halfedge();

  // Dynamic pmaps do not have default initialization values (yet)
  for(vertex_descriptor v : vertices(pm))
  {
    put(known_nm_vertices, v, false);
    put(visited_vertices, v, null_h);
  }
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

      vertex_descriptor v = target(h, pm);

      if(get(visited_vertices, v) != null_h) // already seen this vertex, but not from this star
      {
        is_non_manifold = true;

        // if this is the second time we visit that vertex and the first star was manifold, we have
        // never reported the first star, but we must now
        if(!get(known_nm_vertices, v))
          *out++ = get(visited_vertices, v); // that's a halfedge of the first star we've seen 'v' in
      }
      else
      {
        // first time we meet this vertex, just mark it so, and keep the halfedge we found the vertex with in memory
        put(visited_vertices, v, h);
      }

      // While walking the star of this halfedge, if we meet a border halfedge more than once,
      // it means the mesh is pinched and we are also in the case of a non-manifold situation
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
        put(known_nm_vertices, v, true);
      }
    }
  }

  return out;
}

} // namespace Polygon_mesh_processing
} // namespace CGAL

// for backward compatibility, and the functions above are required 'repair_manifoldness.h'
#include <CGAL/Polygon_mesh_processing/repair_manifoldness.h>

#endif // CGAL_POLYGON_MESH_PROCESSING_MANIFOLDNESS_H
