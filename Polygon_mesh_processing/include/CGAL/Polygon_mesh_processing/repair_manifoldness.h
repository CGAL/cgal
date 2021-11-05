// Copyright (c) 2019-2020 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé,
//                 Sébastien Loriot

#ifndef CGAL_POLYGON_MESH_PROCESSING_REPAIR_MANIFOLDNESS_H
#define CGAL_POLYGON_MESH_PROCESSING_REPAIR_MANIFOLDNESS_H

#include <CGAL/license/Polygon_mesh_processing/repair.h>

#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>
#include <CGAL/Polygon_mesh_processing/internal/Repair/helper.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/clip.h>
#include <CGAL/Polygon_mesh_processing/clip_self_intersecting.h>
#include <CGAL/Polygon_mesh_processing/detect_features.h>
#include <CGAL/Polygon_mesh_processing/manifoldness.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Canvas/BGL_canvas.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_halfedge_graph_segment_primitive.h>
#include <CGAL/assertions.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/selection.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/Heat_method_3/Surface_mesh_geodesic_distances_3.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/IO/facets_in_complex_2_to_triangle_mesh.h>
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/number_utils.h>
#include <CGAL/Origin.h>
#include <CGAL/property_map.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/utility.h>

#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <stack>
#include <utility>
#include <vector>

namespace CGAL {
namespace Polygon_mesh_processing {

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/// Combinatorial treatment

namespace internal {

template <typename G>
struct Vertex_collector
{
  typedef typename boost::graph_traits<G>::vertex_descriptor                   vertex_descriptor;
  typedef std::map<vertex_descriptor, std::vector<vertex_descriptor> >         Collection;

  const Collection& collection() const { return _coll; }

  bool has_old_vertex(const vertex_descriptor v) const { return _coll.count(v) != 0; }
  void tag_old_vertex(const vertex_descriptor v)
  {
    CGAL_precondition(!has_old_vertex(v));
    _coll[v];
  }

  void collect_vertices(vertex_descriptor v1, vertex_descriptor v2)
  {
    std::vector<vertex_descriptor>& verts = _coll[v1];
    if(verts.empty())
      verts.push_back(v1);
    verts.push_back(v2);
  }

  template<typename OutputIterator>
  void dump(OutputIterator out)
  {
    typedef std::pair<const vertex_descriptor, std::vector<vertex_descriptor> > Pair_type;
    for(const Pair_type& p : _coll)
      *out++ = p.second;
  }

  void dump(Emptyset_iterator) { }

  Collection _coll;
};

// Replaces the current target vertex of a fan of faces, described by two halfedges, with a new vertex
template <typename PolygonMesh, typename VPM, typename ConstraintMap>
typename boost::graph_traits<PolygonMesh>::vertex_descriptor
create_new_vertex_for_sector(typename boost::graph_traits<PolygonMesh>::halfedge_descriptor sector_begin_h,
                             typename boost::graph_traits<PolygonMesh>::halfedge_descriptor sector_last_h,
                             PolygonMesh& pmesh,
                             const VPM& vpm,
                             const ConstraintMap& cmap)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor    vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor  halfedge_descriptor;

  vertex_descriptor old_vd = target(sector_begin_h, pmesh);
  vertex_descriptor new_vd = add_vertex(pmesh);
  put(vpm, new_vd, get(vpm, old_vd));

  put(cmap, new_vd, true);

  set_halfedge(new_vd, sector_begin_h, pmesh);
  halfedge_descriptor h = sector_begin_h;
  do
  {
    set_target(h, new_vd, pmesh);

    if(h == sector_last_h)
      break;
    else
      h = prev(opposite(h, pmesh), pmesh);
  }
  while(h != sector_begin_h); // for safety only
  CGAL_postcondition(h != sector_begin_h);

  return new_vd;
}

template <typename PolygonMesh, typename VPM, typename CMAP>
std::size_t make_umbrella_manifold(typename boost::graph_traits<PolygonMesh>::halfedge_descriptor h,
                                   PolygonMesh& pmesh,
                                   internal::Vertex_collector<PolygonMesh>& dmap,
                                   VPM vpm,
                                   CMAP cmap)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor    vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor  halfedge_descriptor;

  std::size_t nb_new_vertices = 0;

  vertex_descriptor old_v = target(h, pmesh);
  put(cmap, old_v, true); // store the duplicates

  // count the number of borders
  int border_counter = 0;
  halfedge_descriptor ih = h, done = ih, border_h = h;
  do
  {
    if(is_border(ih, pmesh))
    {
      border_h = ih;
      ++border_counter;
    }

    ih = prev(opposite(ih, pmesh), pmesh);
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
      set_halfedge(old_v, h, pmesh); // to ensure halfedge(old_v, pmesh) stays valid
      dmap.tag_old_vertex(old_v); // so that we know we have met old_v already, next time, we'll have to duplicate
    }
    else
    {
      // This is not the canonical star associated to 'v'.
      // Create a new vertex, and move the whole star to that new vertex
      halfedge_descriptor last_h = opposite(next(h, pmesh), pmesh);
      vertex_descriptor new_v = create_new_vertex_for_sector(h, last_h, pmesh, vpm, cmap);
      dmap.collect_vertices(old_v, new_v);
      nb_new_vertices = 1;
    }
  }
  // if there is more than one sector, look at each sector and split them away from the main one
  else
  {
    // the first manifold sector, described by two halfedges
    halfedge_descriptor sector_start_h = border_h;
    CGAL_assertion(is_border(border_h, pmesh));

    bool should_stop = false;
    bool is_main_sector = true;
    do
    {
      CGAL_assertion(is_border(sector_start_h, pmesh));

      // collect the sector and split it away if it must be
      halfedge_descriptor sector_last_h = sector_start_h;
      do
      {
        halfedge_descriptor next_h = prev(opposite(sector_last_h, pmesh), pmesh);

        if(is_border(next_h, pmesh))
          break;

        sector_last_h = next_h;
      }
      while(sector_last_h != sector_start_h);
      CGAL_assertion(!is_border(sector_last_h, pmesh));
      CGAL_assertion(sector_last_h != sector_start_h);

      halfedge_descriptor next_start_h = prev(opposite(sector_last_h, pmesh), pmesh);

      // there are multiple CCs incident to this particular vertex, and we should create a new vertex
      // if it's not the first umbrella around 'old_v' or not the first sector, but only not if it's
      // both the first umbrella and first sector.
      bool must_create_new_vertex = (!is_main_sector || dmap.has_old_vertex(old_v));

      // In any case, we must set up the next pointer correctly
      set_next(sector_start_h, opposite(sector_last_h, pmesh), pmesh);

      if(must_create_new_vertex)
      {
        vertex_descriptor new_v = create_new_vertex_for_sector(sector_start_h, sector_last_h, pmesh, vpm, cmap);
        dmap.collect_vertices(old_v, new_v);
        ++nb_new_vertices;
      }
      else
      {
        // We are in the first sector and first star, ensure that halfedge(old_v, pmesh) stays valid
        set_halfedge(old_v, sector_start_h, pmesh);
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
/// duplicates all the non-manifold vertices of the input mesh.
///
/// @tparam PolygonMesh a model of `HalfedgeListGraph` and `MutableHalfedgeGraph`
/// @tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
///
/// @param pm the surface mesh to be repaired
/// @param np optional \ref pmp_namedparameters "Named Parameters" described below
///
/// \cgalNamedParamsBegin
///    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `pmesh`.
///       The type of this map is model of `ReadWritePropertyMap`.
///       If this parameter is omitted, an internal property map for
///       `CGAL::vertex_point_t` should be available in `PolygonMesh`
///    \cgalParamEnd
///   \cgalParamBegin{vertex_is_constrained_map} a writable property map with `vertex_descriptor`
///     as key and `bool` as `value_type`. `put(pmap, v, true)` will be called for each duplicated
///     vertices, as well as the original non-manifold vertex in the input mesh.
///  \cgalParamEnd
///   \cgalParamBegin{output_iterator} a model of `OutputIterator` with value type
///      `std::vector<vertex_descriptor>`. The first vertex of each vector is a non-manifold vertex
///       of the input mesh, followed by the new vertices that were created to fix this precise
///       non-manifold configuration.
///  \cgalParamEnd
/// \cgalNamedParamsEnd
///
/// \return the number of vertices created.
template <typename PolygonMesh, typename NamedParameters>
std::size_t duplicate_non_manifold_vertices(PolygonMesh& pmesh,
                                            const NamedParameters& np)
{
  using parameters::get_parameter;
  using parameters::choose_parameter;

  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor         vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor       halfedge_descriptor;

  typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::type VertexPointMap;
  VertexPointMap vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                                        get_property_map(vertex_point, pmesh));

  typedef typename internal_np::Lookup_named_param_def<
                                  internal_np::vertex_is_constrained_t,
                                  NamedParameters,
                                  Static_boolean_property_map<vertex_descriptor, false> // default (no constraint pmap)
                                  >::type                                      VerticesMap;
  VerticesMap cmap = choose_parameter(get_parameter(np, internal_np::vertex_is_constrained),
                                      Static_boolean_property_map<vertex_descriptor, false>());

  typedef typename internal_np::Lookup_named_param_def<
                                  internal_np::output_iterator_t,
                                  NamedParameters,
                                  Emptyset_iterator>::type                     Output_iterator;
  Output_iterator out = choose_parameter(get_parameter(np, internal_np::output_iterator),
                                         Emptyset_iterator());

  std::vector<halfedge_descriptor> non_manifold_umbrellas;
  non_manifold_vertices(pmesh, std::back_inserter(non_manifold_umbrellas));

  internal::Vertex_collector<PolygonMesh> dmap;
  std::size_t nb_new_vertices = 0;
  if(!non_manifold_umbrellas.empty())
  {
    for(halfedge_descriptor h : non_manifold_umbrellas)
      nb_new_vertices += internal::make_umbrella_manifold(h, pmesh, dmap, vpm, cmap);

    dmap.dump(out);
  }

  return nb_new_vertices;
}

template <class PolygonMesh>
std::size_t duplicate_non_manifold_vertices(PolygonMesh& pmesh)
{
  return duplicate_non_manifold_vertices(pmesh, parameters::all_default());
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/// Geometrical treatment

// main:
// @todo VPM for heat method / heat method on a FFG

enum NM_TREATMENT
{
  SEPARATE = 0,
  CLIP,
  MERGE
};

namespace internal {

template <class PolygonMesh, class NamedParameters>
void preprocess_mesh_for_treatment_of_non_manifold_vertices(PolygonMesh& pmesh,
                                                            const bool join_when_looking_from_outside,
                                                            const NamedParameters& np)
{
  namespace pred = CGAL::Polygon_mesh_processing::Corefinement;

  // extract types from NPs
  typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type               Traits;
  typedef typename CGAL::GetInitializedVertexIndexMap<PolygonMesh, NamedParameters>::type VertexIndexMap;
  typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::type           VPM;

  typedef typename Traits::Point_3                                                 Point_3;

  typedef boost::graph_traits<PolygonMesh>                                         Graph_traits;
  typedef typename Graph_traits::vertex_descriptor                                 vertex_descriptor;
  typedef typename Graph_traits::edge_descriptor                                   edge_descriptor;
  typedef typename Graph_traits::halfedge_descriptor                               halfedge_descriptor;
  typedef typename Graph_traits::face_descriptor                                   face_descriptor;

  typedef std::size_t                                                              PID;

  using parameters::choose_parameter;
  using parameters::get_parameter;

  VertexIndexMap vim = CGAL::get_initialized_vertex_index_map(pmesh, np);
  VPM vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_property_map(vertex_point, pmesh));

  std::unordered_map<Point_3, PID> point_ids; // Point_3 to unique-id
  std::vector<Point_3> points; // All the points: id -> Point_3
  std::vector<PID> pids(num_vertices(pmesh)); // maps each vertex to its id (through vim)

  // identifies identical vertices
  points.reserve(num_vertices(pmesh));
  for(vertex_descriptor v : vertices(pmesh))
  {
    auto insert_res = point_ids.emplace(get(vpm, v), points.size());
    if(insert_res.second)
      points.push_back(get(vpm, v));

    pids[get(vim, v)] = insert_res.first->second;
  }

  std::cout << "points.size() " << points.size() << "\n";
  std::cout << "vertices.size() " << vertices(pmesh).size() << "\n";

  // collect duplicated edges
  std::vector<std::unordered_map<std::size_t, std::vector<edge_descriptor> > > edge_map(points.size());
  for(edge_descriptor ed : edges(pmesh))
  {
    std::pair<PID,PID> e_pids = CGAL::make_sorted_pair(pids[get(vim, source(ed, pmesh))],
                                                       pids[get(vim, target(ed, pmesh))]);
    edge_map[e_pids.first][e_pids.second].push_back(ed);
  }

  std::vector<std::vector<halfedge_descriptor> > sorted_hedges; // use a list?
  for(PID pid0 = 0; pid0<points.size(); ++pid0)
  {
    for(const auto& pid1_and_edges : edge_map[pid0])
    {
      if(pid1_and_edges.second.size() > 1)
      {
        std::cout << "processing a non-manifold edge...\n";

        for(edge_descriptor ed : pid1_and_edges.second)
        {
          // @todo remove the limitation or return without doing nothing
          CGAL_assertion(!is_border(ed, pmesh));
          if(is_border(ed, pmesh))
            return;
        }

        PID pid1 = pid1_and_edges.first;
        sorted_hedges.resize(sorted_hedges.size()+1);

        for(edge_descriptor ed : pid1_and_edges.second)
        {
          halfedge_descriptor h = halfedge(ed, pmesh);
          if(!is_border(h, pmesh))
            sorted_hedges.back().push_back(h);

          h = opposite(h, pmesh);
          if(!is_border(h, pmesh))
            sorted_hedges.back().push_back(h);
        }

        const Point_3& ref = get(vpm, target(next(sorted_hedges.back()[0], pmesh), pmesh));
        auto less = [&ref, &pmesh, &points, vpm, pid0, pid1](halfedge_descriptor h1, halfedge_descriptor h2)
        {
          return pred::sorted_around_edge<Traits>(points[pid0], points[pid1], ref,
                                                  get(vpm, target(next(h1,pmesh), pmesh)),
                                                  get(vpm, target(next(h2, pmesh), pmesh)));
        };

        std::sort(sorted_hedges.back().begin()+1, sorted_hedges.back().end(), less);

        // we have the faces of sorted_hedges.back()[0](f0) and sorted_hedges.back()[1](f1)
        // that are consecutive along the edge defined by(points[pid0],points[pid1]),
        // and in counterclockwise order when looking from points[pids0].
        // To know if the volume between f0 and f1 is inside or outside we need to look at
        // the order of the vertices in f0 and we enter the interior of the volume
        // if in f0 id(src(sorted_hedges.back()[0]) == pid0

        const bool sweeping_interior = (pids[get(vim,source(sorted_hedges.back()[0], pmesh))] == pid0);
        if(sweeping_interior == join_when_looking_from_outside)
          std::rotate(sorted_hedges.back().begin(), sorted_hedges.back().begin()+1, sorted_hedges.back().end());

        // @todo: if we want to avoid the rotate we should change sorted_hedges.back()[0] to be
      }
    }
  }

  // collect vertices per pid for edges that we are going to update
  std::set<PID> pids_to_update;
  std::vector< std::set<vertex_descriptor> > vertices_per_pid(points.size());
  for(const std::vector<halfedge_descriptor>& hedges : sorted_hedges)
  {
    for(halfedge_descriptor h : hedges)
    {
      PID src_id = pids[get(vim, source(h, pmesh))];
      PID tgt_id = pids[get(vim, target(h, pmesh))];

      vertices_per_pid[ src_id].insert(source(h, pmesh));
      vertices_per_pid[ tgt_id ].insert(target(h, pmesh));

      pids_to_update.insert(src_id);
      pids_to_update.insert(tgt_id);
    }
  }

  // pick one vertex for all vertex with the same point
  for(PID pid : pids_to_update)
  {
    std::cout << "PID vertices:";
    CGAL_assertion(!vertices_per_pid[pid].empty());
    vertex_descriptor ref = *std::prev(vertices_per_pid[pid].end());
    for(vertex_descriptor v : vertices_per_pid[pid])
    {
      std::cout << " " << v;
      halfedge_descriptor h = halfedge(v, pmesh);
      if(target(h, pmesh) != ref)
      {
        for(halfedge_descriptor h2up : halfedges_around_target(h, pmesh))
          set_target(h2up, ref, pmesh);
      }
      set_halfedge(v, Graph_traits::null_halfedge(), pmesh); // to indicate that the vertex has not been processed
    }
    vertices_per_pid[pid].erase(std::prev(vertices_per_pid[pid].end()));

    std::cout << " -- ref = " << ref << "\n";
  }

  // update halfedges
  auto replace_halfedge = [&pmesh](halfedge_descriptor h_new, halfedge_descriptor h_old)
  {
    halfedge_descriptor prv = prev(h_old, pmesh);
    halfedge_descriptor nxt = next(h_old, pmesh);
    face_descriptor f = face(h_old, pmesh);
    vertex_descriptor tgt = target(h_old, pmesh);

    set_face(h_new, f, pmesh);
    set_next(h_new, nxt, pmesh);
    set_next(prv, h_new, pmesh);
    set_halfedge(f, h_new, pmesh);
    set_target(h_new, tgt, pmesh);
    set_face(h_old, Graph_traits::null_face(), pmesh);
  };

//   std::size_t kk=0;
  std::vector<std::vector<halfedge_descriptor> > hedges_to_postprocess;
  hedges_to_postprocess.reserve(sorted_hedges.size());
  for(const std::vector<halfedge_descriptor>& hedges : sorted_hedges)
  {
    hedges_to_postprocess.resize(hedges_to_postprocess.size()+1);

    CGAL_assertion(hedges.size()%2 == 0);
    std::set<edge_descriptor> edges_to_remove;
    for(std::size_t i=0; i<hedges.size(); i+=2)
    {
      halfedge_descriptor h1 = hedges[i], h2 = hedges[i+1];

//      std::ofstream debug("results/wedge_"+std::to_string(kk++)+".polylines.txt");
//      debug << "3 " << pmesh.point(source(h1,pmesh)) << " " << pmesh.point(target(next(h1,pmesh),pmesh)) << " " << pmesh.point(target(h1,pmesh)) << "\n";
//      debug << "3 " << pmesh.point(source(h2,pmesh)) << " " << pmesh.point(target(next(h2,pmesh),pmesh)) << " " << pmesh.point(target(h2,pmesh)) << "\n";

      if(opposite(h1, pmesh) != h2)
      {
        // we create a new halfedge
        halfedge_descriptor h_new = halfedge(add_edge(pmesh), pmesh);
        replace_halfedge(h_new, hedges[i]);
        replace_halfedge(opposite(h_new, pmesh), hedges[i+1]);
        edges_to_remove.insert(edge(h1, pmesh));
        edges_to_remove.insert(edge(h2, pmesh));
        hedges_to_postprocess.back().push_back(h_new);
        hedges_to_postprocess.back().push_back(opposite(h_new, pmesh));
      }
      else
      {
        hedges_to_postprocess.back().push_back(h1);
        hedges_to_postprocess.back().push_back(h2);
      }
    }

    for(edge_descriptor ed : edges_to_remove)
      remove_edge(ed, pmesh);
  }

  // set vertex per umbrella
  for(const std::vector<halfedge_descriptor>& hedges : hedges_to_postprocess)
  {
    for(halfedge_descriptor h : hedges)
    {
      if(halfedge(target(h, pmesh), pmesh) == Graph_traits::null_halfedge())
      {
        // vertex needs to be updated
        PID pid = pids[ get(vim, target(h, pmesh))];

        if(!vertices_per_pid[pid].empty())
        {
          set_target(h, *std::prev(vertices_per_pid[pid].end()), pmesh);
          vertices_per_pid[pid].erase(std::prev(vertices_per_pid[pid].end()));
        }

        set_halfedge(target(h, pmesh), h, pmesh);

        for(halfedge_descriptor hat : halfedges_around_target(h, pmesh))
          set_target(hat, target(h, pmesh), pmesh);
      }
    }
  }

  // @todo: find a better way to erase ref vertices
  std::vector<vertex_descriptor> vertices_to_remove;
  for(vertex_descriptor v : vertices(pmesh))
    if(halfedge(v, pmesh) == Graph_traits::null_halfedge())
      vertices_to_remove.push_back(v);

  std::cout <<"removing " << vertices_to_remove.size() << " vertices\n";

  for(vertex_descriptor v : vertices_to_remove)
  {
    std::cout << "removing " << v << "\n";
    remove_vertex(v, pmesh);
  }

  // post-processing for edges that need to be split as still non-manifold:
  // look for edges with the same vertices as endpoint
  // @todo: This will not work for boundary hedges
  for(const std::vector<halfedge_descriptor>& hedges : hedges_to_postprocess)
  {
    std::map<std::pair<vertex_descriptor,vertex_descriptor>, std::vector<halfedge_descriptor> > l_edge_map;
    for(std::size_t k=0; k<hedges.size(); k+=2)
      l_edge_map[CGAL::make_sorted_pair(source(hedges[k], pmesh),
                                        target(hedges[k], pmesh))].push_back(hedges[k]);

    for(const auto& vrts_and_hedges : l_edge_map)
    {
      if(vrts_and_hedges.second.size() > 1)
      {
        Point_3 midpt = midpoint(get(vpm, source(vrts_and_hedges.second.front(), pmesh)),
                                 get(vpm, target(vrts_and_hedges.second.front(), pmesh)));
        for(halfedge_descriptor h_to_split : vrts_and_hedges.second)
        {
          halfedge_descriptor hnew = Euler::split_edge(h_to_split, pmesh);
          put(vpm, target(hnew, pmesh), midpt);
          Euler::split_face(hnew, next(next(hnew, pmesh), pmesh), pmesh);
          hnew = opposite(h_to_split, pmesh);
          Euler::split_face(hnew, next(next(hnew, pmesh), pmesh), pmesh);
        }
      }
    }
  }
}

template <typename NMEdgeRange, typename PolygonMesh, typename VPM, typename GeomTraits>
void sample_non_manifold_edges(const NMEdgeRange& nm_edges,
                               const typename GeomTraits::FT r,
                               PolygonMesh& pmesh,
                               VPM vpm,
                               const GeomTraits& gt)
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor       halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::edge_descriptor           edge_descriptor;

  typedef typename boost::property_traits<VPM>::value_type                     Point;
  typedef typename boost::property_traits<VPM>::reference                      Point_ref;

  typedef typename GeomTraits::FT                                              FT;
  typedef typename GeomTraits::Vector_3                                        Vector;

  CGAL_precondition(! CGAL_NTS is_zero(r));

  for(const edge_descriptor e : nm_edges)
  {
    halfedge_descriptor h = halfedge(e, pmesh);

    // to get a canonical direction and the same subdivision
    if(get(vpm, source(h, pmesh)) > get(vpm, target(h, pmesh)))
      h = opposite(h, pmesh);

    const Point_ref sp = get(vpm, source(h, pmesh));
    const Point_ref tp = get(vpm, target(h, pmesh));

    // Need the balls of radius `r` to cover the halfedge
    const FT length = CGAL_NTS approximate_sqrt(gt.compute_squared_distance_3_object()(sp, tp));
    const int interval_n = std::ceil(CGAL::to_double(length / r));
    int points_n = interval_n + 1;

    if(points_n <= 2)
      continue;

    points_n -= 2;

    Vector d{sp, tp};
    d = gt.construct_scaled_vector_3_object()(d, 1 / FT(interval_n));

#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
    std::cout << "Sample edge " << sp << " " << tp << " with " << points_n << " extra points" << std::endl;
#endif

    for(int i=1; i<=points_n; ++i)
    {
      const Point new_pos = gt.construct_translated_point_3_object()(sp, i * d);
      halfedge_descriptor new_h = Euler::split_edge(h, pmesh);
      put(vpm, target(new_h, pmesh), new_pos);
    }
  }

  Polygon_mesh_processing::triangulate_faces(faces(pmesh), pmesh); // @speed
}

template <typename DistanceMap, typename VertexSet,
          typename PolygonMesh, typename VPM, typename GT>
void compute_distances_with_Campen(DistanceMap distances,
                                   const VertexSet& sources,
                                   const PolygonMesh& pmesh,
                                   const VPM vpm,
                                   const GT& gt)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor         vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor       halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::edge_descriptor           edge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor           face_descriptor;

  typedef typename GT::FT                                                      FT;
  typedef typename boost::property_traits<VPM>::value_type                     Point;
  typedef typename boost::property_traits<VPM>::reference                      Point_ref;
  typedef typename GT::Vector_3                                                Vector;

  using Metric_field = CGAL::Canvas::Metric_fields::Euclidean_metric_field<GT>;
  using Canvas = CGAL::Canvas::BGL_canvas<PolygonMesh, VPM, Metric_field, GT>;

  Metric_field* metric_field = new Metric_field(1.0, 1.0, 1.0);
  Canvas canvas(pmesh, vpm, metric_field, gt);

  canvas.initialize();
  canvas.initialize_seeds(sources);
  canvas.paint(distances);
}

template <typename VertexSet,
          typename PolygonMesh, typename VPM, typename GeomTraits>
void split_too_long_incident_edges(const VertexSet& sources,
                                   const typename GeomTraits::FT radius,
                                   PolygonMesh& pmesh,
                                   VPM vpm,
                                   const GeomTraits& gt)
{
  using FT = typename GeomTraits::FT;
  using Point_3 = typename GeomTraits::Point_3;

  using vertex_descriptor = typename boost::graph_traits<PolygonMesh>::vertex_descriptor;
  using halfedge_descriptor = typename boost::graph_traits<PolygonMesh>::halfedge_descriptor;
  using edge_descriptor = typename boost::graph_traits<PolygonMesh>::edge_descriptor;

  const FT sq_radius = CGAL::square(radius);

  std::set<edge_descriptor> edges_to_treat;
  for(const vertex_descriptor v : sources)
  {
    for(const halfedge_descriptor h : halfedges_around_target(halfedge(v, pmesh), pmesh))
    {
      const vertex_descriptor vs = source(h, pmesh);
      if(sources.find(vs) == sources.end())
        continue;

      const FT sq_length = gt.compute_squared_distance_3_object()(get(vpm, vs), get(vpm, v));
      if(sq_length > sq_radius)
        edges_to_treat.insert(edge(h, pmesh));
    }
  }

  for(edge_descriptor e : edges_to_treat)
  {
    const Point_3 new_pos = gt.construct_midpoint_3_object()(get(vpm, source(e, pmesh)),
                                                             get(vpm, target(e, pmesh)));
    halfedge_descriptor new_h = CGAL::Euler::split_edge(halfedge(e, pmesh), pmesh);
    put(vpm, target(new_h, pmesh), new_pos);
  }

  Polygon_mesh_processing::triangulate_faces(faces(pmesh), pmesh); // @speed
}

template <typename HalfedgeSet, typename EdgeSet, typename FaceSet,
          typename PolygonMesh, typename VPM, typename GeomTraits>
void geodesic_refine_and_select_faces(const HalfedgeSet& nm_vertices,
                                      const EdgeSet& nm_edges,
                                      const typename GeomTraits::FT radius,
                                      PolygonMesh& pmesh,
                                      VPM vpm,
                                      const GeomTraits& gt,
                                      FaceSet& selected_faces)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor         vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor       halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::edge_descriptor           edge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor           face_descriptor;

  typedef typename GeomTraits::FT                                              FT;
  typedef typename boost::property_traits<VPM>::value_type                     Point;
  typedef typename boost::property_traits<VPM>::reference                      Point_ref;
  typedef typename GeomTraits::Vector_3                                        Vector;

  typedef CGAL::dynamic_vertex_property_t<FT>                                  Distance_tag;
  typedef typename boost::property_map<PolygonMesh, Distance_tag>::type        Vertex_distance_map;
  Vertex_distance_map distances = get(Distance_tag(), pmesh);

#ifdef CGAL_PMP_REPAIR_NM_USE_HEAT_METHOD
  Heat_method_3::Surface_mesh_geodesic_distances_3<PolygonMesh, CGAL::Heat_method_3::Intrinsic_Delaunay, VPM> smgd(pmesh, vpm);
#else
  std::set<vertex_descriptor> sources;
#endif

#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
  auto vim = CGAL::get_initialized_vertex_index_map(pmesh);
  std::set<std::size_t> unique_source_ids;

  std::ofstream out_s("results/HM_sources.xyz");
  out_s.precision(17);
#endif

  for(const halfedge_descriptor h : nm_vertices)
  {
#ifdef CGAL_PMP_REPAIR_NM_USE_HEAT_METHOD
    smgd.add_source(target(h, pmesh));
#else
    sources.insert(target(h, pmesh));
#endif

#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
    unique_source_ids.insert(get(vim, target(h, pmesh)));
    out_s << get(vpm, target(h, pmesh)) << "\n";
#endif
  }

  // @todo most extremities of nm_edges are already in nm_vertices
  for(const edge_descriptor e : nm_edges)
  {
#ifdef CGAL_PMP_REPAIR_NM_USE_HEAT_METHOD
    smgd.add_source(source(e, pmesh));
    smgd.add_source(target(e, pmesh));
#else
    sources.insert(source(e, pmesh));
    sources.insert(target(e, pmesh));
#endif

#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
    unique_source_ids.insert(get(vim, source(e, pmesh)));
    unique_source_ids.insert(get(vim, target(e, pmesh)));
    out_s << get(vpm, source(e, pmesh)) << "\n";
    out_s << get(vpm, target(e, pmesh)) << "\n";
#endif
  }

#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
  std::ofstream out_sel("results/HM_sources.selection.txt");
  for(std::size_t i : unique_source_ids)
    out_sel << i << " ";
  out_sel << std::endl;
#endif

  // Avoid big faces being iso-refined due to to their vertices all being sources
#ifdef CGAL_PMP_REPAIR_NM_USE_HEAT_METHOD
  split_too_long_incident_edges(smgd.sources(), radius, pmesh, vpm, gt);
#else
  split_too_long_incident_edges(sources, radius, pmesh, vpm, gt);
#endif

#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
  CGAL::IO::write_polygon_mesh("results/post_long_edge_split.off", pmesh, CGAL::parameters::stream_precision(17));
#endif

#ifdef CGAL_PMP_REPAIR_NM_USE_HEAT_METHOD
  smgd.estimate_geodesic_distances(distances);
#else
  compute_distances_with_Campen(distances, sources, pmesh, vpm, gt);
#endif

  for(const vertex_descriptor& v : vertices(pmesh))
    std::cout << "Vertex " << v << " pos " << pmesh.point(v) << " Distance: " << get(distances, v) << std::endl;

  // Corefine the mesh with a piecewise(face)-linear approximation of the union of the geodesic circles
  std::set<face_descriptor> faces_to_consider;
  std::set<halfedge_descriptor> edges_to_split;

  std::stack<halfedge_descriptor> halfedges_to_consider;

  for(halfedge_descriptor h : nm_vertices)
    halfedges_to_consider.push(h);

  for(edge_descriptor e : nm_edges)
  {
    halfedges_to_consider.push(halfedge(e, pmesh));
    halfedges_to_consider.push(opposite(halfedge(e, pmesh), pmesh));
  }

  typedef CGAL::dynamic_edge_property_t<bool>                                  Considered_tag;
  typedef typename boost::property_map<PolygonMesh, Considered_tag>::type      Considered_edge_map;
  Considered_edge_map considered_edges = get(Considered_tag(), pmesh);

  while(!halfedges_to_consider.empty())
  {
    const halfedge_descriptor curr_h = halfedges_to_consider.top();
    halfedges_to_consider.pop();

    const edge_descriptor curr_e = edge(curr_h, pmesh);
    if(get(considered_edges, curr_e))
      continue;

    put(considered_edges, curr_e, true);

    const bool is_s_in = (get(distances, source(curr_e, pmesh)) <= radius);
    const bool is_t_in = (get(distances, target(curr_e, pmesh)) <= radius);

    if(is_s_in && is_t_in)
    {
      if(!is_border(curr_h, pmesh))
        faces_to_consider.insert(face(curr_h, pmesh));
      if(!is_border(opposite(curr_h, pmesh), pmesh))
        faces_to_consider.insert(face(opposite(curr_h, pmesh), pmesh));
    }

    if(is_s_in != is_t_in)
      edges_to_split.insert(curr_h);

    for(const halfedge_descriptor adj_h : CGAL::halfedges_around_source(curr_h, pmesh))
    {
      if(get(considered_edges, edge(adj_h, pmesh))) // already visited
        continue;

      // want the source of the new halfedges to be equal to 'source(curr_h, pmesh)'
      halfedges_to_consider.push(opposite(adj_h, pmesh));
    }
  }

#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
  std::cout << edges_to_split.size() << " edges to split" << std::endl;
#endif

#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
  CGAL::IO::write_polygon_mesh("results/pre_geodesy.off", pmesh, parameters::stream_precision(17));
#endif

  // Don't insert points too close to the isovalue, it can create small edges,
  // which isotropic_remeshing has trouble dealing with
  const FT eps = 0.05 * radius;

  // Actual split
  for(halfedge_descriptor h : edges_to_split)
  {
    const vertex_descriptor vs = source(h, pmesh);
    const vertex_descriptor vt = target(h, pmesh);
    const Point_ref spt = get(vpm, vs);
    const Point_ref tpt = get(vpm, vt);
    Vector tsv = gt.construct_vector_3_object()(tpt, spt);

    const FT dist_at_vs = get(distances, vs);
    const FT dist_at_vt = get(distances, vt);
    if(std::abs(radius - dist_at_vs) < eps || std::abs(radius - dist_at_vt) < eps) // nothing to do
      continue;

    Point new_p;
    if(dist_at_vs < dist_at_vt)
    {
      CGAL_assertion(dist_at_vs < radius && radius <= dist_at_vt);
      const FT lambda = (radius - dist_at_vs) / (dist_at_vt - dist_at_vs);
      new_p = gt.construct_translated_point_3_object()(
                spt, gt.construct_scaled_vector_3_object()(tsv, - lambda));
    }
    else
    {
      CGAL_assertion(dist_at_vt < radius && radius <= dist_at_vs);
      const FT lambda = (radius - dist_at_vt) / (dist_at_vs - dist_at_vt);
      new_p = gt.construct_translated_point_3_object()(
                tpt, gt.construct_scaled_vector_3_object()(tsv, lambda));
    }

    // avoid inexact constructions creating degenerate faces, which is problematic later during fairing
    // --
    if(new_p == spt || new_p == tpt)
      continue;

    if(!is_border(h, pmesh))
    {
      const Point_ref opt = get(vpm, target(next(h, pmesh), pmesh));
      if(gt.collinear_3_object()(opt, spt, new_p) || gt.collinear_3_object()(opt, tpt, new_p))
        continue;
    }

    if(!is_border(opposite(h, pmesh), pmesh))
    {
      const Point_ref opt = get(vpm, target(next(opposite(h, pmesh), pmesh), pmesh));
      if(gt.collinear_3_object()(opt, spt, new_p) || gt.collinear_3_object()(opt, tpt, new_p))
        continue;
    }
    // --

    halfedge_descriptor new_h = Euler::split_edge_and_incident_faces(h, pmesh);
    put(vpm, target(new_h, pmesh), new_p);
    put(distances, target(new_h, pmesh), radius);

    if(!is_border(h, pmesh))
    {
      faces_to_consider.insert(face(h, pmesh));
      faces_to_consider.insert(face(new_h, pmesh));
    }

    if(!is_border(opposite(h, pmesh), pmesh))
    {
      faces_to_consider.insert(face(opposite(h, pmesh), pmesh));
      faces_to_consider.insert(face(opposite(new_h, pmesh), pmesh));
    }
  }

#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
  std::cout << faces_to_consider.size() << " faces to consider" << std::endl;

  static int fi = 0;
  std::stringstream oss;
  oss << "results/faces_to_consider_" << fi++ << ".off" << std::ends;
  dump_cc(faces_to_consider, pmesh, oss.str().c_str());
#endif

  for(const face_descriptor f : faces_to_consider)
  {
    bool is_face_in = true;
    for(halfedge_descriptor h : CGAL::halfedges_around_face(halfedge(f, pmesh), pmesh))
    {
      if(get(distances, target(h, pmesh)) > radius)
      {
        is_face_in = false;
        break;
      }
    }

    if(is_face_in)
      selected_faces.insert(f);
  }

#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
  CGAL::IO::write_polygon_mesh("results/post_geodesy.off", pmesh, parameters::stream_precision(17));
#endif
}

template <typename HalfedgeSet, typename EdgeSet, typename FaceSet,
          typename PolygonMesh, typename VPM, typename GeomTraits>
void select_treatment_faces(const HalfedgeSet& nm_vertices,
                            const EdgeSet& nm_edges,
                            const typename GeomTraits::FT radius,
                            PolygonMesh& pmesh,
                            VPM vpm,
                            const GeomTraits& gt,
                            FaceSet& selected_faces)
{
  typedef Boolean_property_map<FaceSet>                                        Face_selection_map;
  Face_selection_map sf_pm(selected_faces);

  // 3*radius because fair() needs a continuous sampling, so we will remesh a larger area,
  // and then only keep a subset of the selected faces
  geodesic_refine_and_select_faces(nm_vertices, nm_edges, 3 * radius, pmesh, vpm, gt, selected_faces);

  // Regularize and sanitize for both settings
  regularize_face_selection_borders(pmesh, sf_pm, 0.5 /*weight*/,
                                    CGAL::parameters::prevent_unselection(true)
                                                     .vertex_point_map(vpm)
                                                     .geom_traits(gt));
#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
  std::cout << selected_faces.size() << " selected faces after regularization" << std::endl;
#endif

  // Sanitize to ensure no non-manifold vertices on the border
  expand_face_selection_for_removal(selected_faces, pmesh, sf_pm);

#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
  std::cout << selected_faces.size() << " selected faces" << std::endl;
  dump_cc(selected_faces, pmesh, "results/face_selection.off");
#endif
}

template <typename VertexSet, typename EdgeSet, typename FaceSet, typename FT,
          typename PolygonMesh, typename VPM, typename GT>
void remesh_selection(VertexSet& nm_vertices,
                      const EdgeSet& nm_edges,
                      const FaceSet& faces,
                      const FT radius,
                      PolygonMesh& pmesh,
                      VPM vpm,
                      const GT& gt)
{
  typedef typename boost::graph_traits<PolygonMesh>::edge_descriptor           edge_descriptor;

  typedef Boolean_property_map<VertexSet>                                      CVPM;
  CVPM cvm(nm_vertices);

  typedef Boolean_property_map<EdgeSet>                                        CEPM;
  EdgeSet features_edges;
  CEPM cem(features_edges);
  Polygon_mesh_processing::detect_sharp_edges(pmesh, 60, cem,
                                              CGAL::parameters::vertex_point_map(vpm)
                                                               .geom_traits(gt));

  // Prevent non-manifold edges from being remeshed
  for(edge_descriptor e : nm_edges)
    put(cem, e, true);

  // Not protecting constraints because they are either:
  // - nm edges and in that case they have been sampled and should not be affected
  // - sharp edges living in the work area, and in that case they should be split/merged
  CGAL::Polygon_mesh_processing::isotropic_remeshing(faces, 0.5 * radius, pmesh,
                                                     CGAL::parameters::vertex_point_map(vpm)
                                                                      .geom_traits(gt)
                                                                      .number_of_iterations(5)
                                                                      .number_of_relaxation_steps(5)
                                                                      .do_project(false)
                                                                      .protect_constraints(false)
                                                                      .collapse_constraints(false)
                                                                      .vertex_is_constrained_map(cvm)
                                                                      .edge_is_constrained_map(cem));
}

template <typename HalfedgeSet, typename EdgeSet, typename FaceSet, typename FT,
          typename PolygonMesh, typename VPM, typename GeomTraits>
void treat_with_fairing(const HalfedgeSet& nm_vertices,
                        const EdgeSet& nm_edges,
                        FaceSet& selected_faces,
                        const FT radius,
                        PolygonMesh& pmesh,
                        VPM vpm,
                        const GeomTraits& gt)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor         vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor       halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::edge_descriptor           edge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor           face_descriptor;

  typedef std::set<vertex_descriptor>                                          Vertex_set;
  typedef Boolean_property_map<Vertex_set>                                     Vertex_selection_map;

  Vertex_set selected_vertices;
  Vertex_selection_map sv_pm(selected_vertices);

  for(const halfedge_descriptor h : nm_vertices)
    selected_vertices.insert(target(h, pmesh));

  for(const edge_descriptor e : nm_edges)
  {
    selected_vertices.insert(source(e, pmesh));
    selected_vertices.insert(target(e, pmesh));
  }

  // Remesh to improve the stability and quality of fairing
  internal::remesh_selection(selected_vertices, nm_edges, selected_faces, radius, pmesh, vpm, gt);

#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
  CGAL::IO::write_polygon_mesh("results/pre_fairing-remeshed.off", pmesh, parameters::stream_precision(17));
#endif

#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
  std::ofstream out_sv_base("results/selected_vertices_base.xyz");
  out_sv_base.precision(17);
  for(const vertex_descriptor v : selected_vertices)
    out_sv_base << pmesh.point(v) << "\n";
  out_sv_base.close();
#endif

  // Get the vertex selection for fairing
  // @fixme? border vertices
  expand_vertex_selection(selected_vertices, pmesh, 2 /*iterations*/, sv_pm,
                          std::inserter(selected_vertices, selected_vertices.end()));

#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
  std::ofstream out("results/selected_vertices.xyz");
  out.precision(17);
  for(const vertex_descriptor v : selected_vertices)
    out << pmesh.point(v) << "\n";
  out.close();

  // clean selection, for the garbage-less mesh 'pre_fairing-remeshed.off'
  std::ofstream out_sel("results/fairing_vertices.selection.txt");
  int i = 0;
  for(vertex_descriptor v : vertices(pmesh))
  {
    if(selected_vertices.count(v) != 0)
      out_sel << i << " ";

    ++i;
  }

  out_sel << std::endl << std::endl << std::endl;
  out_sel.close();
#endif

  bool res = CGAL::Polygon_mesh_processing::fair(pmesh, selected_vertices,
                                                 CGAL::parameters::fairing_continuity(2)
                                                                  .vertex_point_map(vpm)
                                                                  .geom_traits(gt));

  if(!res)
    std::cerr << "Warning: fairing failed" << std::endl;

#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
  CGAL::IO::write_polygon_mesh("results/faired.off", pmesh, CGAL::parameters::stream_precision(17));
#endif
}

} // namespace internal

namespace experimental {

template <typename PolygonMesh, typename NamedParameters>
std::size_t repair_non_manifoldness(PolygonMesh& pmesh,
                                    const NM_TREATMENT treatment,
                                    const NamedParameters& np)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor         vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor       halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::edge_descriptor           edge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor           face_descriptor;

  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::type       VertexPointMap;
  VertexPointMap vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                                        get_property_map(vertex_point, pmesh));

  typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type           Geom_traits;
  Geom_traits gt = choose_parameter<Geom_traits>(get_parameter(np, internal_np::geom_traits));

  typedef typename Geom_traits::FT                                             FT;
  typedef typename boost::property_traits<VertexPointMap>::value_type          Point;

  const FT radius = 0.1; // @todo np (+ default value?)

  CGAL_precondition(! CGAL_NTS is_zero(radius));

  // Ensure that halfedge(v, pmesh) is canonical and describes a unique umbrella (+ no pinched stars)
  duplicate_non_manifold_vertices(pmesh, np);

  if(treatment == SEPARATE)
    internal::preprocess_mesh_for_treatment_of_non_manifold_vertices(pmesh, false /*join_when_looking_from_outside*/, np);
  else if(treatment == MERGE)
    internal::preprocess_mesh_for_treatment_of_non_manifold_vertices(pmesh, true /*join_when_looking_from_outside*/, np);

  // @fixme is this actually needed?
  PMP::stitch_borders(pmesh);

#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
  CGAL::IO::write_polygon_mesh("results/preprocessed_input.off", pmesh, parameters::stream_precision(17));
#endif

  std::set<edge_descriptor> nm_edges;
  non_manifold_edges(pmesh, std::inserter(nm_edges, nm_edges.end()), np);

  // If the radius is small compare to edge size, we need to ensure the distance field is correct
  //
  // Not a problem in the separation strategy since pushing away extremities pushes away
  // the whole edge.
  //
  // !! Don't change this '0.5' without changing the other one (in iso_remesh)
  internal::sample_non_manifold_edges(nm_edges, 0.5 * radius, pmesh, vpm, gt);

#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
  CGAL::IO::write_polygon_mesh("results/post_sampling.off", pmesh, parameters::stream_precision(17));
#endif

  // @speed avoid this by just updating the set during sampling
  // ---
  nm_edges.clear();
  non_manifold_edges(pmesh, std::inserter(nm_edges, nm_edges.end()), np);
  // ---

  // @todo there is a lot of redundancy between non-manifold vertices and non-manifold edges
  std::set<halfedge_descriptor> nm_vertices;
  geometrically_non_manifold_vertices(pmesh, std::inserter(nm_vertices, nm_vertices.end()), np);

#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
  std::ofstream out_nmv("results/nm_vertices.xyz");
  out_nmv.precision(17);

  std::cout << nm_vertices.size() << " vertice(s) to treat:" << std::endl;
  for(const halfedge_descriptor h : nm_vertices)
  {
    std::cout << "NM vertex at pos (" << get(vpm, target(h, pmesh)) << ") "
              << "canonical halfedge: " << h << std::endl;
    out_nmv << pmesh.point(target(h, pmesh)) << "\n";
  }
  out_nmv.close();

  std::ofstream out_nme("results/nm_edges.polylines.txt");
  out_nme.precision(17);

  std::cout << nm_edges.size() << " edge(s) to treat:" << std::endl;
  for(const edge_descriptor e : nm_edges)
  {
    std::cout << "NM edge (" << get(vpm, source(e, pmesh)) << ") - ("
                             << get(vpm, target(e, pmesh)) << ")" << std::endl;
    out_nme << "2 " << get(vpm, source(e, pmesh)) << " " << get(vpm, target(e, pmesh)) << "\n";
  }
  out_nme.close();
#endif

  // @todo is that even meaningful
  const std::size_t initial_n = nm_vertices.size() + nm_edges.size();
  if(initial_n == 0)
    return initial_n;

  // Construct the zone(s) that will be modified
  std::set<face_descriptor> selected_faces;
  internal::select_treatment_faces(nm_vertices, nm_edges, radius, pmesh, vpm, gt, selected_faces);

#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
  internal::dump_cc(selected_faces, pmesh, "results/selected_faces.off");
  CGAL::IO::write_polygon_mesh("results/pre_fairing.off", pmesh, parameters::stream_precision(17));
#endif

  // @speed avoid recomputing it (this is because geodesic refinement messes up the halfedges)
  nm_vertices.clear();
  geometrically_non_manifold_vertices(pmesh, std::inserter(nm_vertices, nm_vertices.end()), np);

  // Update the mesh to remove non-manifoldness
  internal::treat_with_fairing(nm_vertices, nm_edges, selected_faces, radius, pmesh, vpm, gt);

#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
  IO::write_polygon_mesh("results/post_fairing.off", pmesh, CGAL::parameters::stream_precision(17));
#endif

  CGAL_postcondition(is_valid_polygon_mesh(pmesh, true));

  CGAL_postcondition_code(nm_edges.clear();)
  CGAL_postcondition_code(non_manifold_edges(pmesh, std::inserter(nm_edges, nm_edges.end()), np);)
  CGAL_postcondition(nm_edges.empty());

  CGAL_postcondition_code(nm_vertices.clear();)
  CGAL_postcondition_code(geometrically_non_manifold_vertices(pmesh, std::inserter(nm_vertices, nm_vertices.end()), np);)
  CGAL_postcondition(nm_vertices.empty());

  return initial_n;
}

template <typename PolygonMesh>
std::size_t repair_non_manifoldness(PolygonMesh& pmesh,
                                    const NM_TREATMENT treatment)
{
  return repair_non_manifoldness(pmesh, treatment, parameters::all_default());
}

template <typename PolygonMesh>
std::size_t repair_non_manifoldness(PolygonMesh& pmesh)
{
  return repair_non_manifoldness(pmesh, SEPARATE, parameters::all_default());
}

} // namespace experimental
} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_REPAIR_MANIFOLDNESS_H
