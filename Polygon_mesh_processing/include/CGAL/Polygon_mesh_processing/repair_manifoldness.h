// Copyright (c) 2020 GeometryFactory (France).
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
#include <CGAL/Polygon_mesh_processing/manifoldness.h>
#include <CGAL/Polygon_mesh_processing/remove_self_intersections.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>

#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/utility.h>

#include <iterator>
#include <fstream>
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
  while(h != sector_begin_h); // for safety
  CGAL_assertion(h != sector_begin_h);

  return new_vd;
}

template <typename PolygonMesh, typename NamedParameters>
std::size_t make_umbrella_manifold(typename boost::graph_traits<PolygonMesh>::halfedge_descriptor h,
                                   PolygonMesh& pmesh,
                                   internal::Vertex_collector<PolygonMesh>& dmap,
                                   const NamedParameters& np)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor    vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor  halfedge_descriptor;

  using parameters::get_parameter;
  using parameters::choose_parameter;

  typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::type VertexPointMap;
  VertexPointMap vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                                        get_property_map(vertex_point, pmesh));

  typedef typename internal_np::Lookup_named_param_def<internal_np::vertex_is_constrained_t,
                                                       NamedParameters,
                                                       Constant_property_map<vertex_descriptor, bool> // default (no constraint pmap)
                                                       >::type                  VerticesMap;
  VerticesMap cmap = choose_parameter(get_parameter(np, internal_np::vertex_is_constrained),
                                      Constant_property_map<vertex_descriptor, bool>(false));

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

  typedef boost::graph_traits<PolygonMesh>                            GT;
  typedef typename GT::halfedge_descriptor                            halfedge_descriptor;

  typedef typename internal_np::Lookup_named_param_def<internal_np::output_iterator_t,
                                                       NamedParameters,
                                                       Emptyset_iterator>::type Output_iterator;

  Output_iterator out = choose_parameter(get_parameter(np, internal_np::output_iterator),
                                         Emptyset_iterator());

  std::vector<halfedge_descriptor> non_manifold_cones;
  non_manifold_vertices(pmesh, std::back_inserter(non_manifold_cones));

  internal::Vertex_collector<PolygonMesh> dmap;
  std::size_t nb_new_vertices = 0;
  if(!non_manifold_cones.empty())
  {
    for(halfedge_descriptor h : non_manifold_cones)
      nb_new_vertices += internal::make_umbrella_manifold(h, pmesh, dmap, np);

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

namespace internal {

template <typename PolygonMesh>
bool is_star_without_border(const typename boost::graph_traits<PolygonMesh>::halfedge_descriptor h,
                            const PolygonMesh& g)
{
  CGAL::Halfedge_around_target_iterator<PolygonMesh> havib, havie;
  for(std::tie(havib, havie) = CGAL::halfedges_around_target(h, g); havib != havie; ++havib)
  {
    if(is_border(*havib, g))
      return false;
  }
  return true;
}

template <typename UmbrellaContainer,
          typename PolygonMesh,
          typename VPM,
          typename GeomTraits>
void separate_umbrellas(const UmbrellaContainer& umbrellas,
                        PolygonMesh& pmesh,
                        const VPM vpm,
                        const GeomTraits& /*gt*/)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor        vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor      halfedge_descriptor;

#define CGAL_PMP_TREAT_NM_VERTICES_COMBINATORIAL_SEPARATE_WITH_MOVE
#ifdef CGAL_PMP_TREAT_NM_VERTICES_COMBINATORIAL_SEPARATE_WITH_MOVE
  // @todo less naive since this can create self-intersections
  // compute normal + move into normal direction?

  for(const halfedge_descriptor h : umbrellas)
  {
    vertex_descriptor v = target(h, pmesh);

#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
  std::cout << "Move " << v << " (" << get(vpm, v) << ")" << std::endl;
#endif

    typename GeomTraits::Point_3 c = CGAL::centroid(pmesh.point(v),
                                                    pmesh.point(source(h, pmesh)),
                                                    pmesh.point(target(next(h, pmesh), pmesh)));

    pmesh.point(v) = c + 0.99 * (pmesh.point(v) - c);
  }
#else
  // separate by removing the incident star + hole filling if possible
  // @todo
#endif
}

template <typename HalfedgeContainer_A, typename HalfedgeContainer_B,
          typename VPM_A, typename VPM_B,
          typename TriangleMesh,
          typename OutputIterator,
          typename GeomTraits>
bool two_borders_hole_fill(const HalfedgeContainer_A& bhv_A,
                           const TriangleMesh& pmesh_A,
                           const VPM_A vpm_A,
                           const HalfedgeContainer_B& bhv_B,
                           const TriangleMesh& pmesh_B,
                           const VPM_B vpm_B,
                           OutputIterator out,
                           const GeomTraits& gt)
{
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor       halfedge_descriptor;

  typedef typename boost::property_traits<VPM_A>::value_type                    Point;
  typedef typename boost::property_traits<VPM_A>::reference                     Point_ref_A;
  typedef typename boost::property_traits<VPM_B>::reference                     Point_ref_B;

  typedef typename HalfedgeContainer_A::const_iterator                          HCit_A;
  typedef typename HalfedgeContainer_B::const_iterator                          HCit_B;

  typedef CGAL::Triple<int, int, int>                                           Face_indices;

  CGAL_precondition(bhv_A.size() != 0);
  CGAL_precondition(bhv_B.size() != 0);

  HCit_A canon_A;
  HCit_B canon_B;

  double best_score = std::numeric_limits<double>::max(); // might not be the best idea...

  // @todo avoid the O(n^2) complexity (but the borders are likely small, so...)
  for(HCit_A it_A=bhv_A.begin(), A_end=bhv_A.end(); it_A!=A_end; ++it_A)
  {
    const Point_ref_A sap = get(vpm_A, source(*it_A, pmesh_A));
    const Point_ref_A tap = get(vpm_A, target(*it_A, pmesh_A));

    for(HCit_B it_B=bhv_B.begin(), B_end=bhv_B.end(); it_B!=B_end; ++it_B)
    {
      const Point_ref_B sbp = get(vpm_B, source(*it_B, pmesh_B));
      const Point_ref_B tbp = get(vpm_B, target(*it_B, pmesh_B));

      const double score = gt.compute_squared_distance_3_object()(sap, tbp)
                         + gt.compute_squared_distance_3_object()(tap, sbp);
      if(score < best_score)
      {
        best_score = score;
        canon_A = it_A;
        canon_B = it_B;
      }
    }
  }

  // polyline
  std::vector<Point> hole_points, third_points;

  // _____   _________   ________
  //      | | <------ | |
  //      | | canon_A | |
  //      | |         | |
  //      | | canon_B | |
  //      | | ------> | |
  // -----   --------   -------

  const Point_ref_A sap = get(vpm_A, source(*canon_A, pmesh_A));
  const Point_ref_A tap = get(vpm_A, target(*canon_A, pmesh_A));
  const Point_ref_B sbp = get(vpm_B, source(*canon_B, pmesh_B));
  const Point_ref_B tbp = get(vpm_B, target(*canon_B, pmesh_B));

#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
  std::cout << "best score:" << std::endl;
  std::cout << sap << std::endl;
  std::cout << tap << std::endl;
  std::cout << sbp << std::endl;
  std::cout << tbp << std::endl;
#endif

  // Walk A's border
  HCit_A it_A = canon_A, last_A = std::prev(bhv_A.end());
  do
  {
    hole_points.push_back(get(vpm_A, source(*it_A, pmesh_A)));
    third_points.push_back(get(vpm_A, target(next(opposite(*it_A, pmesh_A), pmesh_A), pmesh_A)));
    it_A = (it_A == last_A) ? bhv_A.begin() : std::next(it_A);
  }
  while(it_A != canon_A);

  // vertical
  hole_points.push_back(sap);
  third_points.push_back(CGAL::midpoint(sbp, tap));

  // Walk B's border
  HCit_B it_B = canon_B, last_B = std::prev(bhv_B.end());
  do
  {
    hole_points.push_back(get(vpm_B, source(*it_B, pmesh_B)));
    third_points.push_back(get(vpm_B, target(next(opposite(*it_B, pmesh_B), pmesh_B), pmesh_B)));
    it_B = (it_B == last_B) ? bhv_B.begin() : std::next(it_B);
  }
  while(it_B != canon_B);

  hole_points.push_back(sbp);
  third_points.push_back(CGAL::midpoint(tbp, sap));

  CGAL_assertion(hole_points.size() == third_points.size());

  std::vector<Face_indices> patch;
  triangulate_hole_polyline(hole_points, third_points, std::back_inserter(patch));

  if(patch.empty())
  {
#ifndef CGAL_HOLE_FILLING_DO_NOT_USE_DT3
#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
    std::cout << "Failed to fill a hole using Delaunay search space.\n";
#endif

    triangulate_hole_polyline(hole_points, third_points, std::back_inserter(patch),
                              parameters::use_delaunay_triangulation(false));
#endif // CGAL_HOLE_FILLING_DO_NOT_USE_DT3
    if(patch.empty())
    {
#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
      std::cout << "Failed to fill a hole using the whole search space.\n";
#endif
      return false;
    }
  }

#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
  std::ofstream pout("results/patch.off");
  pout << std::setprecision(17);
  pout << 3 * patch.size() << " " << patch.size() << " 0\n";

  for(const Face_indices& f : patch)
  {
    pout << hole_points[f.first] << "\n";
    pout << hole_points[f.second] << "\n";
    pout << hole_points[f.third] << "\n";
  }

  for(std::size_t i=0, ps=patch.size(); i<ps; ++i)
    pout << "3 " << 3*i << " " << 3*i+1 << " " << 3*i+2 << "\n";
#endif

  for(const Face_indices& face : patch)
  {
    *out++ = std::initializer_list<Point>{ hole_points[face.first],
                                           hole_points[face.second],
                                           hole_points[face.third] };
  }

  return true;
}

template <typename PolygonMesh,
          typename VPM,
          typename GeomTraits>
bool merge_umbrellas(const typename boost::graph_traits<PolygonMesh>::halfedge_descriptor h1,
                     const typename boost::graph_traits<PolygonMesh>::halfedge_descriptor h2,
                     PolygonMesh& pmesh,
                     const VPM vpm,
                     const GeomTraits& gt)
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor      halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor          face_descriptor;

  typedef typename GeomTraits::Point_3                                        Point;

#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
  std::cout << "Merging umbrellas around " << target(h1, pmesh) << " (" << h1 << " and " << h2 << ")" << std::endl;
#endif

  CGAL_precondition(target(h1, pmesh) == target(h2, pmesh));

#define CGAL_PMP_TREAT_NM_VERTICES_REMOVE_STARS
#ifdef CGAL_PMP_TREAT_NM_VERTICES_REMOVE_STARS
  // remove incident stars

  std::vector<halfedge_descriptor> bhv_1, bhv_2;
  std::set<face_descriptor> faces_to_delete;

  CGAL::Halfedge_around_target_iterator<PolygonMesh> havib, havie;
  for(std::tie(havib, havie) = CGAL::halfedges_around_target(h1, pmesh); havib != havie; ++havib)
  {
    bhv_1.push_back(*havib);
    faces_to_delete.insert(face(*havib, pmesh));
  }

  for(std::tie(havib, havie) = CGAL::halfedges_around_target(h2, pmesh); havib != havie; ++havib)
  {
    bhv_2.push_back(*havib);
    faces_to_delete.insert(face(*havib, pmesh));
  }
#else
  // clip with a sphere?

#endif

  // make sure that the holes are topological disks
  std::vector<std::vector<Point> > point_patch;
  point_patch.reserve(2 * bhv_1.size());
  if(!two_borders_hole_fill(bhv_1, pmesh, vpm, bhv_2, pmesh, vpm, std::back_inserter(point_patch), gt))
    return false;

  if(!check_patch_sanity<PolygonMesh>(point_patch))
  {
#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
    std::cout << "this patch is INSANE!" << std::endl;
#endif
    return false;
  }

#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
  std::cout << "Replacing " << faces_to_delete.size() << " faces with " << point_patch.size() << " new faces" << std::endl;
#endif

  replace_faces_with_patch(faces_to_delete, point_patch, pmesh, vpm);

  return true;
}

} // namespace internal

enum NM_TREATMENT
{
  SEPARATE = 0,
  MERGE
};

template <typename PolygonMesh, typename NamedParameters>
void treat_non_manifold_vertices(PolygonMesh& pmesh,
                                 const NM_TREATMENT treatment,
                                 const NamedParameters& np)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor        vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor      halfedge_descriptor;

  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::type      VertexPointMap;
  VertexPointMap vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                                        get_property_map(vertex_point, pmesh));

  typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type          Geom_traits;
  Geom_traits gt = choose_parameter<Geom_traits>(get_parameter(np, internal_np::geom_traits));

#define CGAL_PMP_TREAT_NM_VERTICES_COMBINATORIAL_NM_ONLY // @todo test the geometric version
#ifdef CGAL_PMP_TREAT_NM_VERTICES_COMBINATORIAL_NM_ONLY
  // not worth getting a dynamic pmap: there are usually few non-manifold vertices
  std::unordered_map<vertex_descriptor, std::vector<halfedge_descriptor> > non_manifold_umbrellas;

  std::vector<halfedge_descriptor> non_manifold_cones;
  non_manifold_vertices(pmesh, std::back_inserter(non_manifold_cones));

  for(const halfedge_descriptor h : non_manifold_cones)
    non_manifold_umbrellas[target(h, pmesh)].push_back(h);
#else
  typedef typename Geom_traits::Point_3                                       Point_3;

  // a little brutal for now, we can do something better akin to non_manifold_vertices()
  std::map<Point_3, std::vector<halfedge_descriptor> > non_manifold_umbrellas;

  typedef CGAL::dynamic_halfedge_property_t<bool>                                       Halfedge_property_tag;
  typedef typename boost::property_map<PolygonMesh, Halfedge_property_tag>::const_type  Visited_halfedge_map;

  Visited_halfedge_map visited_halfedges = get(Halfedge_property_tag(), pmesh);
  for(const halfedge_descriptor h : halfedges(pmesh))
    put(visited_halfedges, h, false);

  for(const halfedge_descriptor h : halfedges(pmesh))
  {
    if(is_border(h, pmesh) || get(visited_halfedges, h))
      continue;

    non_manifold_umbrellas[get(vpm, target(h, pmesh))].push_back(h);

    // move ccw
    halfedge_descriptor ih = h, done = ih;
    do
    {
      put(visited_halfedges, ih, true);
      ih = prev(opposite(ih, pmesh), pmesh);
      if(is_border(ih, pmesh))
        break;
    }
    while(ih != done);

    // move cw
    if(ih == done) // no
    {
      ih = h;
      do
      {
        put(visited_halfedges, ih, true);
        ih = opposite(next(ih, pmesh), pmesh);
        if(is_border(ih, pmesh))
          break;
      }
      while(ih != done);
    }
  }
#endif

  if(non_manifold_umbrellas.empty())
    return;

  for(const auto& e : non_manifold_umbrellas)
  {
    const std::vector<halfedge_descriptor>& cones = e.second;

    CGAL_assertion(cones.size() >= 2);
    CGAL_assertion(std::set<halfedge_descriptor>(cones.begin(), cones.end()).size() == cones.size());

    if(treatment == SEPARATE)
      internal::separate_umbrellas(cones, pmesh, vpm, gt);

    if(cones.size() > 2 ||
       !internal::is_star_without_border(cones.front(), pmesh) ||
       !internal::is_star_without_border(cones.back(), pmesh))
    {
#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
      std::cout << "MERGING treatment requested, but impossible" << std::endl;
#endif
      internal::separate_umbrellas(cones, pmesh, vpm, gt);
    }

    internal::merge_umbrellas(cones.front(), cones.back(), pmesh, vpm, gt);
  }

  CGAL_postcondition_code(non_manifold_cones.clear();)
  CGAL_postcondition_code(non_manifold_vertices(pmesh, std::back_inserter(non_manifold_cones));)
  CGAL_postcondition(non_manifold_cones.empty());
}

template <typename PolygonMesh>
void treat_non_manifold_vertices(PolygonMesh& pmesh,
                                 const NM_TREATMENT treatment)
{
  return treat_non_manifold_vertices(pmesh, treatment, parameters::all_default());
}

template <typename PolygonMesh>
void treat_non_manifold_vertices(PolygonMesh& pmesh)
{
  return treat_non_manifold_vertices(pmesh, MERGE, parameters::all_default());
}

} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_REPAIR_MANIFOLDNESS_H
