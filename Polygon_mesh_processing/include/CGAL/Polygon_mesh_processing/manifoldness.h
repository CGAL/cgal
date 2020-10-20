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
#include <map>
#include <utility>
#include <vector>

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

namespace internal {

template <typename G>
struct Vertex_collector
{
  typedef typename boost::graph_traits<G>::vertex_descriptor      vertex_descriptor;

  bool has_old_vertex(const vertex_descriptor v) const { return collections.count(v) != 0; }
  void tag_old_vertex(const vertex_descriptor v)
  {
    CGAL_precondition(!has_old_vertex(v));
    collections[v];
  }

  void collect_vertices(vertex_descriptor v1, vertex_descriptor v2)
  {
    std::vector<vertex_descriptor>& verts = collections[v1];
    if(verts.empty())
      verts.push_back(v1);
    verts.push_back(v2);
  }

  template<typename OutputIterator>
  void dump(OutputIterator out)
  {
    typedef std::pair<const vertex_descriptor, std::vector<vertex_descriptor> > Pair_type;
    for(const Pair_type& p : collections)
      *out++ = p.second;
  }

  void dump(Emptyset_iterator) { }

  std::map<vertex_descriptor, std::vector<vertex_descriptor> > collections;
};

template <typename PolygonMesh, typename VPM, typename ConstraintMap>
typename boost::graph_traits<PolygonMesh>::vertex_descriptor
create_new_vertex_for_sector(typename boost::graph_traits<PolygonMesh>::halfedge_descriptor sector_begin_h,
                             typename boost::graph_traits<PolygonMesh>::halfedge_descriptor sector_last_h,
                             PolygonMesh& pm,
                             const VPM& vpm,
                             const ConstraintMap& cmap)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor    vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor  halfedge_descriptor;

  vertex_descriptor old_vd = target(sector_begin_h, pm);
  vertex_descriptor new_vd = add_vertex(pm);
  put(vpm, new_vd, get(vpm, old_vd));

  put(cmap, new_vd, true);

  set_halfedge(new_vd, sector_begin_h, pm);
  halfedge_descriptor h = sector_begin_h;
  do
  {
    set_target(h, new_vd, pm);

    if(h == sector_last_h)
      break;
    else
      h = prev(opposite(h, pm), pm);
  }
  while(h != sector_begin_h); // for safety
  CGAL_assertion(h != sector_begin_h);

  return new_vd;
}

template <typename PolygonMesh, typename NamedParameters>
std::size_t make_umbrella_manifold(typename boost::graph_traits<PolygonMesh>::halfedge_descriptor h,
                                   PolygonMesh& pm,
                                   internal::Vertex_collector<PolygonMesh>& dmap,
                                   const NamedParameters& np)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor    vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor  halfedge_descriptor;

  using parameters::get_parameter;
  using parameters::choose_parameter;

  typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::type VertexPointMap;
  VertexPointMap vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                                        get_property_map(vertex_point, pm));

  typedef typename internal_np::Lookup_named_param_def<internal_np::vertex_is_constrained_t,
                                                       NamedParameters,
                                                       Static_boolean_property_map<vertex_descriptor, false> // default (no constraint pmap)
                                                       >::type                  VerticesMap;
  VerticesMap cmap = choose_parameter(get_parameter(np, internal_np::vertex_is_constrained),
                                      Static_boolean_property_map<vertex_descriptor, false>());

  std::size_t nb_new_vertices = 0;

  vertex_descriptor old_v = target(h, pm);
  put(cmap, old_v, true); // store the duplicates

  // count the number of borders
  int border_counter = 0;
  halfedge_descriptor ih = h, done = ih, border_h = h;
  do
  {
    if(is_border(ih, pm))
    {
      border_h = ih;
      ++border_counter;
    }

    ih = prev(opposite(ih, pm), pm);
  }
  while(ih != done);

  bool is_non_manifold_within_umbrella = (border_counter > 1);
  if(!is_non_manifold_within_umbrella)
  {
    const bool first_time_meeting_v = !dmap.has_old_vertex(old_v);
    if(first_time_meeting_v)
    {
      // The star is manifold, so if it is the first time we have met that vertex,
      // there is nothing to do, we just keep the same vertex.
      set_halfedge(old_v, h, pm); // to ensure halfedge(old_v, pm) stays valid
      dmap.tag_old_vertex(old_v); // so that we know we have met old_v already, next time, we'll have to duplicate
    }
    else
    {
      // This is not the canonical star associated to 'v'.
      // Create a new vertex, and move the whole star to that new vertex
      halfedge_descriptor last_h = opposite(next(h, pm), pm);
      vertex_descriptor new_v = create_new_vertex_for_sector(h, last_h, pm, vpm, cmap);
      dmap.collect_vertices(old_v, new_v);
      nb_new_vertices = 1;
    }
  }
  // if there is more than one sector, look at each sector and split them away from the main one
  else
  {
    // the first manifold sector, described by two halfedges
    halfedge_descriptor sector_start_h = border_h;
    CGAL_assertion(is_border(border_h, pm));

    bool should_stop = false;
    bool is_main_sector = true;
    do
    {
      CGAL_assertion(is_border(sector_start_h, pm));

      // collect the sector and split it away if it must be
      halfedge_descriptor sector_last_h = sector_start_h;
      do
      {
        halfedge_descriptor next_h = prev(opposite(sector_last_h, pm), pm);

        if(is_border(next_h, pm))
          break;

        sector_last_h = next_h;
      }
      while(sector_last_h != sector_start_h);
      CGAL_assertion(!is_border(sector_last_h, pm));
      CGAL_assertion(sector_last_h != sector_start_h);

      halfedge_descriptor next_start_h = prev(opposite(sector_last_h, pm), pm);

      // there are multiple CCs incident to this particular vertex, and we should create a new vertex
      // if it's not the first umbrella around 'old_v' or not the first sector, but only not if it's
      // both the first umbrella and first sector.
      bool must_create_new_vertex = (!is_main_sector || dmap.has_old_vertex(old_v));

      // In any case, we must set up the next pointer correctly
      set_next(sector_start_h, opposite(sector_last_h, pm), pm);

      if(must_create_new_vertex)
      {
        vertex_descriptor new_v = create_new_vertex_for_sector(sector_start_h, sector_last_h, pm, vpm, cmap);
        dmap.collect_vertices(old_v, new_v);
        ++nb_new_vertices;
      }
      else
      {
        // We are in the first sector and first star, ensure that halfedge(old_v, pm) stays valid
        set_halfedge(old_v, sector_start_h, pm);
      }

      is_main_sector = false;
      sector_start_h = next_start_h;
      should_stop = (sector_start_h == border_h);
    }
    while(!should_stop);
  }

  return nb_new_vertices;
}

} // end namespace internal

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

/// \ingroup PMP_repairing_grp
/// duplicates all the non-manifold vertices of the input mesh.
///
/// @tparam PolygonMesh a model of `HalfedgeListGraph` and `MutableHalfedgeGraph`
/// @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
///
/// @param pm the surface mesh to be repaired
/// @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
///
/// \cgalNamedParamsBegin
///   \cgalParamNBegin{vertex_point_map}
///     \cgalParamDescription{a property map associating points to the vertices of `pm`}
///     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
///                    as key type and `%Point_3` as value type}
///     \cgalParamDefault{`boost::get(CGAL::vertex_point, pm)`}
///     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
///                     must be available in `PolygonMesh`.}
///   \cgalParamNEnd
///
///   \cgalParamNBegin{vertex_is_constrained_map}
///     \cgalParamDescription{a property map containing the constrained-or-not status of each vertex of `pm`.}
///     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
///                    as key type and `bool` as value type. It must be default constructible.}
///     \cgalParamDefault{a default property map where no vertex is constrained}
///     \cgalParamExtra{`put(vcm, v, true)` will be called for each duplicated
///                     vertices, as well as the original non-manifold vertex in the input mesh.}
///   \cgalParamNEnd
///
///   \cgalParamNBegin{output_iterator}
///     \cgalParamDescription{an output iterator to collect the duplicated vertices}
///     \cgalParamType{a model of `OutputIterator` with value type `std::vector<vertex_descriptor>`}
///     \cgalParamDefault{unused}
///     \cgalParamExtra{The first vertex of each vector is a non-manifold vertex of the input mesh,
///                     followed by the new vertices that were created to fix the given non-manifold configuration.}
///   \cgalParamNEnd
/// \cgalNamedParamsEnd
///
/// \return the number of vertices created.
template <typename PolygonMesh, typename NamedParameters>
std::size_t duplicate_non_manifold_vertices(PolygonMesh& pm,
                                            const NamedParameters& np)
{
  using parameters::get_parameter;
  using parameters::choose_parameter;

  typedef boost::graph_traits<PolygonMesh>                            GT;
  typedef typename GT::halfedge_descriptor                            halfedge_descriptor;

  typedef typename internal_np::Lookup_named_param_def<internal_np::output_iterator_t,
                                                       NamedParameters,
                                                       Emptyset_iterator>::type Output_iterator;

  Output_iterator out = choose_parameter(get_parameter(np, internal_np::output_iterator),
                                         Emptyset_iterator());

  std::vector<halfedge_descriptor> non_manifold_cones;
  non_manifold_vertices(pm, std::back_inserter(non_manifold_cones));

  internal::Vertex_collector<PolygonMesh> dmap;
  std::size_t nb_new_vertices = 0;
  if(!non_manifold_cones.empty())
  {
    for(halfedge_descriptor h : non_manifold_cones)
      nb_new_vertices += internal::make_umbrella_manifold(h, pm, dmap, np);

    dmap.dump(out);
  }

  return nb_new_vertices;
}

template <class PolygonMesh>
std::size_t duplicate_non_manifold_vertices(PolygonMesh& pm)
{
  return duplicate_non_manifold_vertices(pm, parameters::all_default());
}

} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_MANIFOLDNESS_H
