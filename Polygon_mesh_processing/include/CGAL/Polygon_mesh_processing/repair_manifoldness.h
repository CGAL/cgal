// Copyright (c) 2019-2020 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

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
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/manifoldness.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>

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

#include <boost/shared_ptr.hpp>

#include <iostream>
#include <iterator>
#include <fstream>
#include <map>
#include <stack>
#include <vector>
#include <utility>

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
                                  >::type                            VerticesMap;
  VerticesMap cmap = choose_parameter(get_parameter(np, internal_np::vertex_is_constrained),
                                      Static_boolean_property_map<vertex_descriptor, false>());

  typedef typename internal_np::Lookup_named_param_def<
                                  internal_np::output_iterator_t,
                                  NamedParameters,
                                  Emptyset_iterator>::type           Output_iterator;
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
//
enum NM_TREATMENT
{
  SEPARATE = 0,
  CLIP,
  MERGE
};

namespace internal {

template <typename PolygonMesh, typename VPM, typename GeomTraits>
void sample_non_manifold_edges(const typename GeomTraits::FT r,
                               PolygonMesh& pmesh,
                               VPM vpm,
                               const GeomTraits& gt)
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor       halfedge_descriptor;

  typedef typename boost::property_traits<VPM>::value_type                     Point;
  typedef typename boost::property_traits<VPM>::reference                      Point_ref;

  typedef typename GeomTraits::FT                                              FT;
  typedef typename GeomTraits::Vector_3                                        Vector;

  CGAL_precondition(! CGAL_NTS is_zero(r));

  std::map<std::pair<Point, Point>, std::vector<halfedge_descriptor> > nm_edges;

  for(const halfedge_descriptor h : halfedges(pmesh))
  {
    Point_ref sp = get(vpm, source(h, pmesh));
    Point_ref tp = get(vpm, target(h, pmesh));

    if(sp > tp)
      continue;

    auto is_insert_successful = nm_edges.emplace(std::make_pair(sp, tp),
                                                 std::initializer_list<halfedge_descriptor>{h});

    if(!is_insert_successful.second)
      is_insert_successful.first->second.push_back(h);
  }

  for(const auto& e : nm_edges)
  {
    if(e.second.size() == 1) // manifold edge, nothing to do
      continue;

    const Point& sp = e.first.first;
    const Point& tp = e.first.second;

    // Need the balls of radius `r` to cover the halfedge
    const FT length = CGAL_NTS approximate_sqrt(gt.compute_squared_distance_3_object()(sp, tp));
    const int interval_n = std::ceil(CGAL::to_double(length / r));
    int points_n = interval_n + 1;

    if(points_n <= 2)
      continue;

    points_n -= 2;

    Vector d{sp, tp};
    d = gt.construct_scaled_vector_3_object()(d, 1 / FT(interval_n));

    std::cout << "length | interval_n | d " << length << " " << interval_n << " " << CGAL::approximate_sqrt(d.squared_length()) << std::endl;
    std::cout << interval_n * CGAL::approximate_sqrt(d.squared_length()) << std::endl;

#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
    std::cout << "Sample edge " << sp << " " << tp << " with " << points_n << " extra points" << std::endl;
#endif

    for(int i=1; i<=points_n; ++i)
    {
      const Point new_pos = gt.construct_translated_point_3_object()(sp, i * d);
      for(halfedge_descriptor h : e.second)
      {
        halfedge_descriptor new_h = Euler::split_edge_and_incident_faces(h, pmesh);
        put(vpm, target(new_h, pmesh), new_pos);
      }
    }
  }

#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
  CGAL::write_polygon_mesh("results/post_sampling.off", pmesh, parameters::stream_precision(17));
#endif
}

template <typename NMPointContainer, typename SelectedFaceContainer,
          typename PolygonMesh, typename VPM, typename GeomTraits>
void geodesic_refine_and_select_faces(const NMPointContainer& nm_points,
                                      const typename GeomTraits::FT radius,
                                      PolygonMesh& pmesh,
                                      VPM vpm,
                                      const GeomTraits& gt,
                                      SelectedFaceContainer& selected_faces)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor         vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor       halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::edge_descriptor           edge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor           face_descriptor;

  typedef typename GeomTraits::FT                                              FT;
  typedef typename boost::property_traits<VPM>::value_type                     Point;
  typedef typename boost::property_traits<VPM>::reference                      Point_ref;
  typedef typename GeomTraits::Vector_3                                        Vector;

  Heat_method_3::Surface_mesh_geodesic_distances_3<PolygonMesh, CGAL::Heat_method_3::Intrinsic_Delaunay, VPM> smgd(pmesh, vpm);

  std::cout << "sources" << std::endl;
  for(const halfedge_descriptor h : nm_points)
  {
    std::cout << get(vpm, target(h, pmesh)) << std::endl;
    smgd.add_source(target(h, pmesh));
  }

  typedef CGAL::dynamic_vertex_property_t<FT>                                  Distance_tag;
  typedef typename boost::property_map<PolygonMesh, Distance_tag>::type        Vertex_distance_map;
  Vertex_distance_map distances = get(Distance_tag(), pmesh);

  smgd.estimate_geodesic_distances(distances);

  // Corefine the mesh with a piecewise(face)-linear approximation of the union of the geodesic circles
  std::set<face_descriptor> faces_to_consider;
  std::vector<halfedge_descriptor> edges_to_split;

  std::stack<halfedge_descriptor> halfedges_to_consider;

  for(halfedge_descriptor h : nm_points)
    halfedges_to_consider.push(h);

  typedef CGAL::dynamic_edge_property_t<bool>                                  Considered_tag;
  typedef typename boost::property_map<PolygonMesh, Considered_tag>::type      Considered_edge_map;
  Considered_edge_map considered_edges = get(Considered_tag(), pmesh);

  while(!halfedges_to_consider.empty())
  {
    const halfedge_descriptor curr_h = halfedges_to_consider.top();
    halfedges_to_consider.pop();

    const edge_descriptor curr_e = edge(curr_h, pmesh);
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
      edges_to_split.push_back(curr_h);

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
    if(dist_at_vs == radius || dist_at_vt == radius) // nothing to do
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
}

template <typename NMPointContainer, typename SelectedFaceContainer,
          typename PolygonMesh, typename VPM, typename GeomTraits>
void select_faces_with_expansion(const NMPointContainer& nm_points,
                                 const int expand_selection_k,
                                 PolygonMesh& pmesh,
                                 VPM vpm,
                                 const GeomTraits& gt,
                                 SelectedFaceContainer& selected_faces,
                                 Boolean_property_map<SelectedFaceContainer>& sf_pm)
{
  select_incident_faces(nm_points, pmesh, std::inserter(selected_faces, selected_faces.end()));

#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
  std::cout << selected_faces.size() << " faces in non-manifold umbrellas" << std::endl;
#endif

  expand_face_selection(selected_faces, pmesh, expand_selection_k,
                        sf_pm, std::inserter(selected_faces, selected_faces.end()));
#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
  std::cout << selected_faces.size() << " selected faces after expansion" << std::endl;
#endif
}

template <typename NMPointContainer, typename SelectedFaceContainer,
          typename PolygonMesh, typename VPM, typename GeomTraits>
void select_treatment_faces(const NMPointContainer& nm_points,
                            const NM_TREATMENT treatment,
                            const int expand_selection_k,
                            const bool regularize_selection,
                            const typename GeomTraits::FT radius,
                            PolygonMesh& pmesh,
                            VPM vpm,
                            const GeomTraits& gt,
                            SelectedFaceContainer& selected_faces)
{
  typedef Boolean_property_map<SelectedFaceContainer>                          Face_selection_map;
  Face_selection_map sf_pm(selected_faces);

  if(! CGAL_NTS is_zero(radius))
    geodesic_refine_and_select_faces(nm_points, radius, pmesh, vpm, gt, selected_faces);
  else
    select_faces_with_expansion(nm_points, expand_selection_k, pmesh, vpm, gt, selected_faces, sf_pm);

  // Regularize and sanitize for both settings
  if(regularize_selection)
  {
    regularize_face_selection_borders(pmesh, sf_pm, 0.5,
                                      CGAL::parameters::prevent_unselection(true)
                                                       .vertex_point_map(vpm)
                                                       .geom_traits(gt));
#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
    std::cout << selected_faces.size() << " selected faces after regularization" << std::endl;
#endif
  }

  // Sanitize to ensure no non-manifold vertices on the border
  expand_face_selection_for_removal(selected_faces, pmesh, sf_pm);

#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
  dump_cc(selected_faces, pmesh, "results/face_selection.off");
  std::cout << selected_faces.size() << " selected faces in final selection" << std::endl;
#endif
}

// @todo something nicer?
template <typename PolygonMesh, typename VPM, typename GeomTraits>
typename boost::property_traits<VPM>::value_type
construct_adjusted_position(const typename boost::graph_traits<PolygonMesh>::vertex_descriptor v,
                            const PolygonMesh& pmesh,
                            VPM vpm,
                            const GeomTraits& gt)
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor       halfedge_descriptor;

  typedef typename GeomTraits::Vector_3                                        Vector_3;
  typedef typename boost::property_traits<VPM>::reference                      Point_ref;

  const typename boost::graph_traits<PolygonMesh>::degree_size_type d = degree(v, pmesh);
  CGAL_assertion(d > 1);

  const Point_ref pt = get(vpm, v);

  Vector_3 move{CGAL::NULL_VECTOR};
  for(auto hn : CGAL::halfedges_around_target(v, pmesh))
  {
    if(is_border(hn, pmesh))
      continue;

    const Point_ref opt1 = get(vpm, source(prev(hn, pmesh), pmesh));
    const Point_ref opt2 = get(vpm, source(hn, pmesh));
    const auto mpt = gt.construct_midpoint_3_object()(opt1, opt2);

    // - "/2" to get a point within the adjacent face
    // - (d-1) is the number of faces in the sector
    move = gt.construct_sum_of_vectors_3_object()(
             move, gt.construct_scaled_vector_3_object()(Vector_3{pt, mpt}, 0.5 / (d - 1)));
  }

  return gt.construct_translated_point_3_object()(pt, move);
}

template <typename SelectedFaceContainer, typename PolygonMesh, typename VPM, typename GeomTraits>
void treat_with_separation(const SelectedFaceContainer& selected_faces,
                           PolygonMesh& pmesh,
                           VPM vpm,
                           const GeomTraits& gt,
                           const int smoothing_iterations)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor         vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor       halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor           face_descriptor;

  typedef typename boost::property_traits<VPM>::value_type                     Point;

  // map rather than dynamic pmap as it should be a small set
  std::vector<std::pair<vertex_descriptor, Point> > updated_positions;

  for(face_descriptor f : selected_faces)
  {
    for(halfedge_descriptor h : CGAL::halfedges_around_face(halfedge(f, pmesh), pmesh))
    {
      const vertex_descriptor v = target(h, pmesh);
      bool is_unconstrained_vertex = true;
      for(face_descriptor inc_f : CGAL::faces_around_target(h, pmesh))
      {
        // a vertex not in a selected face must be constrained
        if(inc_f != boost::graph_traits<PolygonMesh>::null_face() && selected_faces.count(inc_f) == 0)
        {
          is_unconstrained_vertex = false;
          break;
        }
      }

      if(is_unconstrained_vertex)
        updated_positions.emplace_back(v, CGAL::ORIGIN);
    }
  }

  for(int i=0; i<smoothing_iterations; ++i)
  {
    for(auto& e : updated_positions)
      e.second = construct_adjusted_position(e.first, pmesh, vpm, gt);

    for(const auto& e : updated_positions)
      put(vpm, e.first, e.second);
  }
}

template <typename FFG, typename PolygonMesh, typename VPM, typename GT, typename NewFacesContainer>
bool fill_holes(const FFG& ffg,
                PolygonMesh& pmesh,
                VPM vpm,
                const GT& gt,
                NewFacesContainer& new_faces)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor         vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor       halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor           face_descriptor;

  typedef typename boost::property_traits<VPM>::value_type                     Point;

  bool all_went_well = true;

  std::set<halfedge_descriptor> visited_halfedges;
  std::vector<halfedge_descriptor> borders_to_fill_reps;

  for(halfedge_descriptor h : halfedges(ffg))
  {
    // seek a halfedge on an unvisited border, on the boundary of the selection zone,
    // but not on the border of the mesh
    if(!CGAL::is_border(h, ffg) || CGAL::is_border(h, pmesh) || !visited_halfedges.insert(h).second)
      continue;

#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
    std::cout << "border representative " << h << std::endl;
#endif

    // keep a single representative by boundary cycle
    bool is_boundary_cycle_eligible_for_filling = true;
    halfedge_descriptor done = h;
    do
    {
      // Post face deletion, all is left are internal halfedges
      if(CGAL::is_border(h, pmesh))
      {
        is_boundary_cycle_eligible_for_filling = false;
        // could break here, but continue to mark all the cycle's halfedges as visited
      }

      visited_halfedges.insert(h);
      h = next(h, ffg);
    }
    while(h != done);

    if(is_boundary_cycle_eligible_for_filling)
      borders_to_fill_reps.push_back(h);
  }

  // Dig
  for(face_descriptor f : faces(ffg))
    Euler::remove_face(halfedge(f, pmesh), pmesh);

  // Fill where it makes sense
  for(halfedge_descriptor bh : borders_to_fill_reps)
  {
    CGAL_assertion(!CGAL::is_border(bh, pmesh) && CGAL::is_border(opposite(bh, pmesh), pmesh));

    // @todo two pass with delaunay&whatnot?
    std::vector<vertex_descriptor> unused_new_vertices;
    auto res = triangulate_refine_and_fair_hole(pmesh, opposite(bh, pmesh),
                                                std::back_inserter(new_faces),
                                                std::back_inserter(unused_new_vertices),
                                                parameters::vertex_point_map(vpm)
                                                           .geom_traits(gt));

    if(!std::get<0>(res))
    {
#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
      std::cout << "Warning: Failed hole filling" << std::endl;
#endif
      all_went_well = false;
    }
  }

  CGAL_postcondition(is_valid_polygon_mesh(pmesh));

  return all_went_well;
}

template <typename SelectedFaceContainer, typename PolygonMesh, typename VPM, typename GeomTraits>
void treat_with_clipping(const SelectedFaceContainer& selected_faces,
                         PolygonMesh& pmesh,
                         VPM vpm,
                         const GeomTraits& gt)
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor       halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::edge_descriptor           edge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor           face_descriptor;

  // split the selected faces into CCs
  typedef CGAL::dynamic_edge_property_t<bool>                                  Constraint_tag;
  typedef typename boost::property_map<PolygonMesh, Constraint_tag>::type      Constrained_edges;
  Constrained_edges cepm = get(Constraint_tag(), pmesh);

  for(const edge_descriptor e : edges(pmesh))
  {
    const halfedge_descriptor h = halfedge(e, pmesh);
    const face_descriptor f = face(h, pmesh);
    const bool is_f_selected = (!is_border(h, pmesh) && selected_faces.count(f));
    const face_descriptor fn = face(opposite(h, pmesh), pmesh);
    const bool is_fn_selected = (!is_border(opposite(h, pmesh), pmesh) && selected_faces.count(fn));
    const bool is_constrained = (is_f_selected != is_fn_selected);
    put(cepm, e, is_constrained);
  }

  typedef typename boost::graph_traits<PolygonMesh>::faces_size_type           faces_size_type;
  typedef CGAL::dynamic_face_property_t<faces_size_type>                       Face_property_tag;
  typedef typename boost::property_map<PolygonMesh, Face_property_tag>::type   Patch_ids_map;
  Patch_ids_map patch_ids_map = get(Face_property_tag(), pmesh);

  const std::size_t ccn = Polygon_mesh_processing::connected_components(
                            pmesh, patch_ids_map, parameters::edge_is_constrained_map(cepm));

  // Deal with each CC independently
  bool all_went_well = true;
  for(std::size_t i=0; i<ccn; ++i)
  {
    // FFG is easy on the eyes, but it's not the cheapest (it's still iterating over the underlying graph)
    CGAL::Face_filtered_graph<PolygonMesh> ffg(pmesh, i, patch_ids_map);

    if(is_empty(ffg) || selected_faces.count(*(faces(ffg).begin())) == 0)
      continue; // not an interesting cc

    std::vector<face_descriptor> new_faces;
    if(!fill_holes(ffg, pmesh, vpm, gt, new_faces))
      all_went_well = false;

    // Probably not necessary, but just to be safe in case the ID map has weird default initialization
    for(const face_descriptor f : new_faces)
      put(patch_ids_map, f, i);
  }

  return all_went_well;
}

template <typename SelectedFaceContainer, typename PolygonMesh, typename VPM, typename GeomTraits>
void treat_with_merge(const SelectedFaceContainer& /*selected_faces*/,
                      PolygonMesh& /*pmesh*/,
                      VPM /*vpm*/,
                      const GeomTraits& /*gt*/)
{
  CGAL_assertion(false); // @todo ?
}

} // namespace internal

template <typename PolygonMesh, typename NamedParameters>
std::size_t repair_non_manifoldness(PolygonMesh& pmesh,
                                    const NM_TREATMENT treatment,
                                    const NamedParameters& np)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor         vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor       halfedge_descriptor;
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

  // @todo nps
  const int expand_selection_k = 2;
  const bool regularize_selection = true;
  const FT radius = 0.1;
  const int smoothing_iterations = expand_selection_k;

  // Ensure that halfedge(v, pmesh) is canonical and unique umbrella (+ no pinched stars)
  duplicate_non_manifold_vertices(pmesh, np);

  if(treatment == CLIP && ! CGAL_NTS is_zero(radius))
  {
    // If the radius is small compare to edge size, we need to ensure a whole nm edge
    // will be clipped away, and not just the area around the extremities.
    //
    // Not a problem in the separation strategy since pushing away extremities pushes away
    // the whole edge.
    internal::sample_non_manifold_edges(radius, pmesh, vpm, gt);
  }

  // Collect the non-manifold vertices
  std::set<halfedge_descriptor> nm_points;
  geometrically_non_manifold_vertices(pmesh, std::inserter(nm_points, nm_points.end()), np);
  const std::size_t initial_n = nm_points.size();

#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
  std::cout << nm_points.size() << " case(s) to treat:" << std::endl;
  for(const halfedge_descriptor h : nm_points)
  {
    std::cout << "NM vertex at pos (" << get(vpm, target(h, pmesh)) << ") "
              << "canonical halfedge: " << h << std::endl;
  }
#endif

  // Construct the zone(s) that will be modified
  std::set<face_descriptor> selected_faces;
  internal::select_treatment_faces(nm_points, treatment, expand_selection_k, expand_selection_k, radius,
                                   pmesh, vpm, gt, selected_faces);

#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
  internal::dump_cc(selected_faces, pmesh, "results/selected_faces.off");
  CGAL::write_polygon_mesh("results/pre_treatment.off", pmesh, parameters::stream_precision(17));
#endif

  // Update the mesh to remove non-manifoldness
  if(treatment == SEPARATE)
    internal::treat_with_separation(selected_faces, pmesh, vpm, gt, smoothing_iterations);
  else if(treatment == CLIP)
    internal::treat_with_clipping(selected_faces, pmesh, vpm, gt);
  else
    internal::treat_with_merge(selected_faces, pmesh, vpm, gt);

#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
  write_polygon_mesh("post_treatment.off", pmesh, CGAL::parameters::stream_precision(17));
#endif

  CGAL_postcondition(is_valid_polygon_mesh(pmesh, true));
  CGAL_postcondition_code(nm_points.clear();)
  CGAL_postcondition_code(geometrically_non_manifold_vertices(pmesh, std::inserter(nm_points, nm_points.end()), np);)
  CGAL_postcondition(nm_points.empty());

  return initial_n;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TriangleMesh, class GeomTraits>
class Offset_function
{
  typedef AABB_halfedge_graph_segment_primitive<TriangleMesh> Primitive;
  typedef AABB_traits<GeomTraits, Primitive> Traits;
  typedef AABB_tree<Traits> Tree;

public:
  template <typename EdgeRange>
  Offset_function(TriangleMesh& tm,
                  const EdgeRange& edges,
                  double offset_distance)
    : m_tree_ptr(new Tree(std::begin(edges), std::end(edges), tm) ),
      m_offset_distance(offset_distance)
  {}

  double operator()(const typename GeomTraits::Point_3& p) const
  {
    typename GeomTraits::Point_3 closest_point = m_tree_ptr->closest_point(p);
    double distance = sqrt(squared_distance(p, closest_point));

    return m_offset_distance - distance;
  }

private:
  boost::shared_ptr<Tree> m_tree_ptr;
  double m_offset_distance;
};

template <typename PolygonMesh, typename NamedParameters>
void repair_non_manifold_edges_with_clipping(PolygonMesh& pmesh,
                                             const NamedParameters& np)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor    vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor  halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::edge_descriptor      edge_descriptor;

  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::type       VertexPointMap;
  VertexPointMap vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                                        get_property_map(vertex_point, pmesh));

  typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type           Geom_traits;
  Geom_traits gt = choose_parameter<Geom_traits>(get_parameter(np, internal_np::geom_traits));

  // default triangulation for Surface_mesher
  typedef typename CGAL::Surface_mesher::Surface_mesh_default_triangulation_3_generator<Geom_traits>::Type Tr;

  // c2t3
  typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;
  typedef typename Geom_traits::Sphere_3 Sphere;
  typedef typename Geom_traits::Point_3 Point;
  typedef typename Geom_traits::FT FT;

  typedef typename boost::property_traits<VertexPointMap>::reference           Point_ref;

  typedef Offset_function<PolygonMesh, Geom_traits> Offset_function;
  typedef CGAL::Implicit_surface_3<Geom_traits, Offset_function> Surface_3;

  double offset_distance = 0.1;
  double angle_bound = 20;
  double radius_bound = 0.1;
  double distance_bound = 0.01;

  CGAL::Bbox_3 bbox = CGAL::Polygon_mesh_processing::bbox_3(pmesh);

  Point center((bbox.xmax() + bbox.xmin())/2,
               (bbox.ymax() + bbox.ymin())/2,
               (bbox.zmax() + bbox.zmin())/2);
  double sqrad = 0.6 * std::sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                                 CGAL::square(bbox.ymax() - bbox.ymin()) +
                                 CGAL::square(bbox.zmax() - bbox.zmin()))
                + offset_distance;
  sqrad = CGAL::square(sqrad);

  std::cout << "Offset distance = " << offset_distance << "\n";
  std::cout << "Bounding sphere center = " << center << "\n";
  std::cout << "Bounding sphere squared radius = " << sqrad << "\n";
  std::cout << "Angular bound " << angle_bound << "\n";
  std::cout << "Radius bound " << radius_bound << "\n";
  std::cout << "Distance bound " << distance_bound << "\n";

  std::map<std::pair<Point, Point>, std::vector<halfedge_descriptor> > nm_edges;

  for(const halfedge_descriptor h : halfedges(pmesh))
  {
    Point_ref sp = get(vpm, source(h, pmesh));
    Point_ref tp = get(vpm, target(h, pmesh));

    if(sp > tp)
      continue;

    auto is_insert_successful = nm_edges.emplace(std::make_pair(sp, tp),
                                                 std::initializer_list<halfedge_descriptor>{h});

    if(!is_insert_successful.second)
      is_insert_successful.first->second.push_back(h);
  }

  std::cout << nm_edges.size() << " nm edges" << std::endl;

  std::set<edge_descriptor> nmes;
  for(const auto& e : nm_edges)
  {
    if(e.second.size() == 1)
      continue;
    for(const halfedge_descriptor h : e.second)
      nmes.insert(edge(h, pmesh));
  }

  Offset_function offset_function(pmesh, nmes, offset_distance);

  Tr tr;
  C2t3 c2t3(tr);

  // defining the surface
  Surface_3 surface(offset_function, Sphere(center, sqrad)); // bounding sphere

  // defining meshing criteria
  CGAL::Surface_mesh_default_criteria_3<Tr> criteria(angle_bound, radius_bound, distance_bound);

  // meshing surface
  std::cout << nmes.size() << " nm edges" << std::endl;
  std::cout << "Make..." << std::endl;
  CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Manifold_tag());
  std::cout << "Final number of points: " << tr.number_of_vertices() << "\n";

  // write to file
  std::string output_name = "offset";
  output_name.resize(output_name.size()-4); // strip .off
  std::stringstream sstr;
  sstr << "results/" << output_name << "_OD" << offset_distance
                                    << "_AB" << angle_bound
                                    << "_RB" << radius_bound
                                    << "_DB" << distance_bound
                                    << ".off";

  std::cout << "Writing result in " << sstr.str() << "\n";

  std::ofstream output(sstr.str().c_str());
  CGAL::output_surface_facets_to_off(output, c2t3);

  PolygonMesh nm_edge_offset;
  CGAL::facets_in_complex_2_to_triangle_mesh(c2t3, nm_edge_offset);

  if(is_outward_oriented(nm_edge_offset))
    reverse_face_orientations(nm_edge_offset);

  //extend the bbox a bit to avoid border cases
  const double xd = (std::max)(1., 0.01 * (bbox.xmax() - bbox.xmin()));
  const double yd = (std::max)(1., 0.01 * (bbox.ymax() - bbox.ymin()));
  const double zd = (std::max)(1., 0.01 * (bbox.zmax() - bbox.zmin()));
  bbox = CGAL::Bbox_3(bbox.xmin()-xd, bbox.ymin()-yd, bbox.zmin()-zd,
                      bbox.xmax()+xd, bbox.ymax()+yd, bbox.zmax()+zd);

  typename Geom_traits::Iso_cuboid_3 ic(bbox);
  PolygonMesh bbox_mesh;
  make_hexahedron(ic[0], ic[1], ic[2], ic[3], ic[4], ic[5], ic[6], ic[7], bbox_mesh);
  triangulate_faces(bbox_mesh);

  copy_face_graph(bbox_mesh, nm_edge_offset);
  write_polygon_mesh("results/clipper.off", nm_edge_offset, parameters::stream_precision(17));

  generic_clip(pmesh, nm_edge_offset);

  write_polygon_mesh("results/clipped.off", pmesh, CGAL::parameters::stream_precision(17));
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

} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_REPAIR_MANIFOLDNESS_H
