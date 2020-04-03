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
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/manifoldness.h>
#include <CGAL/Polygon_mesh_processing/repair_self_intersections.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>

#include <CGAL/boost/graph/Euler_operations.h>
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

  std::vector<halfedge_descriptor> non_manifold_umbrellas;
  non_manifold_vertices(pmesh, std::back_inserter(non_manifold_umbrellas));

  internal::Vertex_collector<PolygonMesh> dmap;
  std::size_t nb_new_vertices = 0;
  if(!non_manifold_umbrellas.empty())
  {
    for(halfedge_descriptor h : non_manifold_umbrellas)
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

enum NM_TREATMENT
{
  SEPARATE = 0,
  MERGE
};

namespace internal {

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Merging treatment

template <typename PointRangeContainer, typename PolygonMesh, typename PVM, typename VPM>
bool add_faces(const PointRangeContainer& pr_range,
               PolygonMesh& pmesh,
               PVM& point_to_vs, // std::map<Point, vertex_descriptor>
               const VPM vpm) // vertex to point
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor        vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor          face_descriptor;

  typedef typename boost::property_traits<VPM>::value_type                    Point;

  const vertex_descriptor null_v = boost::graph_traits<PolygonMesh>::null_vertex();

  for(const auto& pr : pr_range)
  {
    std::vector<vertex_descriptor> vf;
    for(const Point& p : pr)
    {
      bool success;
      typename std::map<Point, vertex_descriptor>::iterator it;
      std::tie(it, success) = point_to_vs.insert(std::make_pair(p, null_v));
      vertex_descriptor& v = it->second;

      if(success) // first time we meet that point, means it`s an interior point and we need to make a new vertex
      {
        v = add_vertex(pmesh);
        put(vpm, v, p);
      }

      vf.push_back(v);
    }

    face_descriptor f = Euler::add_face(vf, pmesh);
    if(f == boost::graph_traits<PolygonMesh>::null_face())
    {
      CGAL_assertion(false);
      return false;
    }
  }

  return true;
}

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

// @todo relax that once pinched star merging is done
template <typename HalfedgeContainer, typename PolygonMesh>
bool are_umbrellas_mergeable(const HalfedgeContainer& umbrellas,
                             const PolygonMesh& pmesh)
{
  return !(umbrellas.size() > 2 ||
           !internal::is_star_without_border(umbrellas.front(), pmesh) ||
           !internal::is_star_without_border(umbrellas.back(), pmesh));
}

template <typename PolygonMesh>
typename boost::graph_traits<PolygonMesh>::halfedge_descriptor
split_edge_and_triangulate_incident_faces(typename boost::graph_traits<PolygonMesh>::halfedge_descriptor h,
                                          PolygonMesh& pmesh)
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor      halfedge_descriptor;

  halfedge_descriptor res = Euler::split_edge(h, pmesh);
  Euler::split_face(res, next(h, pmesh), pmesh);

  if(!is_border(edge(h, pmesh), pmesh))
  {
    halfedge_descriptor opp_h = opposite(h, pmesh);
    Euler::split_face(opp_h, next(next(opp_h, pmesh), pmesh), pmesh);
  }

  return res;
}

// Merging combinatorially changes the star, so we don't want to have two adjacent nm vertices
template <typename NMVM, typename PolygonMesh, typename VPM, typename GeomTraits>
void enforce_non_manifold_vertex_separation(NMVM nm_marks,
                                            PolygonMesh& pmesh,
                                            VPM vpm,
                                            const GeomTraits& gt)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor        vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor      halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::edge_descriptor          edge_descriptor;

  typedef typename boost::property_traits<VPM>::reference                     Point_ref;
  typedef typename boost::property_traits<VPM>::value_type                    Point;

  std::set<edge_descriptor> edges_to_split;
  for(const edge_descriptor e : edges(pmesh))
  {
    if(get(nm_marks, source(e, pmesh)) && get(nm_marks, target(e, pmesh)))
    {
#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
      std::cout << "edge that needs to be combinatorially separated: " << e
                << " vertices: " << source(e, pmesh) << " " << target(e, pmesh) << std::endl;
#endif
      edges_to_split.insert(e);
    }
  }

  for(const edge_descriptor e : edges_to_split)
  {
    halfedge_descriptor h = halfedge(e, pmesh);
    if(is_border(h, pmesh))
      h = opposite(h, pmesh);

    halfedge_descriptor new_h = split_edge_and_triangulate_incident_faces(h, pmesh);

    const vertex_descriptor vs = source(h, pmesh);
    const vertex_descriptor vt = target(h, pmesh);
    const Point_ref sp = get(vpm, vs);
    const Point_ref tp = get(vpm, vt);
    const Point mp = gt.construct_midpoint_3_object()(sp, tp);

    put(vpm, target(new_h, pmesh), mp);
  }
}

// @todo that's a rough version of something smarter to be still done
template <typename PolygonMesh, typename VPM, typename GeomTraits>
typename boost::graph_traits<PolygonMesh>::halfedge_descriptor
dig_hole(const typename boost::graph_traits<PolygonMesh>::halfedge_descriptor h,
         const typename GeomTraits::FT radius,
         PolygonMesh& pmesh,
         VPM vpm,
         const GeomTraits& gt)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor        vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor      halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor          face_descriptor;

  typedef typename boost::property_traits<VPM>::reference                     Point_ref;
  typedef typename boost::property_traits<VPM>::value_type                    Point;

  typedef typename GeomTraits::Vector_3                                       Vector;

  const vertex_descriptor vt = target(h, pmesh);
  const Point_ref tpt = get(vpm, vt);

  halfedge_descriptor ih = h, done = h;
  do
  {
    CGAL_assertion(target(ih, pmesh) == target(h, pmesh));

    const halfedge_descriptor next_ih = prev(opposite(ih, pmesh), pmesh);
    const vertex_descriptor vs = source(ih, pmesh);

    halfedge_descriptor ih_to_split = ih;
    if(is_border(ih_to_split, pmesh))
      ih_to_split = opposite(ih_to_split, pmesh);

    const halfedge_descriptor new_ih = split_edge_and_triangulate_incident_faces(ih_to_split, pmesh);

    // note that below uses 'ih', which always points to nm point
    const Point_ref spt = get(vpm, vs);
    Vector tsv(tpt, spt);
    internal::normalize(tsv, gt);

    const Point new_pt = tpt + radius * tsv;
    put(vpm, target(new_ih, pmesh), new_pt);

    ih = next_ih;
  }
  while(ih != done);

  halfedge_descriptor anchor_h = h;
  if(is_border(h, pmesh))
    anchor_h = opposite(next(opposite(h, pmesh), pmesh), pmesh);
  else
    anchor_h = opposite(prev(h, pmesh), pmesh);

  // now delete the faces
  Iterator_range<Face_around_target_iterator<PolygonMesh> > fat = faces_around_target(h, pmesh);
  std::vector<face_descriptor> faces_to_delete(std::begin(fat), std::end(fat));

  for(face_descriptor f : faces_to_delete)
    Euler::remove_face(halfedge(f, pmesh), pmesh);

  CGAL_postcondition(is_border(opposite(anchor_h, pmesh), pmesh));
  return opposite(anchor_h, pmesh);
}

// compute the faces to remove
template <typename UmbrellaContainer, typename PolygonMesh, typename VPM, typename GeomTraits>
std::vector<typename boost::graph_traits<PolygonMesh>::halfedge_descriptor>
dig_holes(UmbrellaContainer& umbrellas,
          const typename GeomTraits::FT radius,
          PolygonMesh& pmesh,
          VPM vpm,
          const GeomTraits& gt)
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor      halfedge_descriptor;

  std::vector<halfedge_descriptor> holes;
  for(const halfedge_descriptor h : umbrellas)
  {
    // @tmp gotta compute, say 1/3rd of the smallest incident edge, something like that
    halfedge_descriptor hole_h = dig_hole(h, radius, pmesh, vpm, gt);
    holes.push_back(hole_h);
  }

  std::ofstream out_pr("results/post_ref.off");
  out_pr << std::setprecision(17) << pmesh;
  out_pr.close();

  return holes;
}

// Currently performed with a hack: transform the borders into a single border
// by manually adding a patch between the closest edges
template <typename HalfedgeContainer_A, typename HalfedgeContainer_B,
          typename VPM_A, typename VPM_B,
          typename PolygonMesh,
          typename OutputIterator,
          typename GeomTraits>
bool two_borders_hole_fill(const HalfedgeContainer_A& bhv_A,
                           const PolygonMesh& pmesh_A,
                           const VPM_A vpm_A,
                           const HalfedgeContainer_B& bhv_B,
                           const PolygonMesh& pmesh_B,
                           const VPM_B vpm_B,
                           OutputIterator out,
                           const GeomTraits& gt)
{
  typedef typename boost::property_traits<VPM_A>::value_type                  Point;
  typedef typename boost::property_traits<VPM_A>::reference                   Point_ref_A;
  typedef typename boost::property_traits<VPM_B>::reference                   Point_ref_B;

  typedef typename HalfedgeContainer_A::const_iterator                        HCit_A;
  typedef typename HalfedgeContainer_B::const_iterator                        HCit_B;

  typedef CGAL::Triple<int, int, int>                                         Face_indices;

  CGAL_precondition(bhv_A.size() >= 3);
  CGAL_precondition(bhv_B.size() >= 3);

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
  std::cout << "sap: " << sap << std::endl;
  std::cout << "tap: " << tap << std::endl;
  std::cout << "sbp: " << sbp << std::endl;
  std::cout << "tbp: " << tbp << std::endl;
#endif

  std::ofstream border_out("results/hole_border.polylines.txt");
  border_out.precision(17);

  // Only a single ID is needed, the rest is 0, hole.size()-1, and quad_id + 1
  std::size_t quad_id = static_cast<std::size_t>(-1);

  // Walk A's border
  HCit_A last_A = std::prev(bhv_A.end());
  HCit_A it_A = (canon_A == last_A) ? bhv_A.begin() : std::next(canon_A);
  do
  {
    if(!hole_points.empty())
      border_out << "2 " << hole_points.back() << " " << get(vpm_A, source(*it_A, pmesh_A)) << std::endl;

    hole_points.push_back(get(vpm_A, source(*it_A, pmesh_A)));
    std::cout << "point on hole: " << hole_points.back() << std::endl;

    third_points.push_back(get(vpm_A, target(next(opposite(*it_A, pmesh_A), pmesh_A), pmesh_A)));
    it_A = (it_A == last_A) ? bhv_A.begin() : std::next(it_A);
  }
  while(it_A != canon_A);

  // vertical
  border_out << "2 " << hole_points.back() << " " << sap << std::endl;

  quad_id = hole_points.size();
  hole_points.push_back(sap);
  third_points.push_back(CGAL::midpoint(sbp, tap));

  // Walk B's border
  HCit_B last_B = std::prev(bhv_B.end());
  HCit_B it_B = (canon_B == last_B) ? bhv_B.begin() : std::next(canon_B);
  do
  {
    border_out << "2 " << hole_points.back() << " " << get(vpm_B, source(*it_B, pmesh_B)) << std::endl;
    hole_points.push_back(get(vpm_B, source(*it_B, pmesh_B)));
    third_points.push_back(get(vpm_B, target(next(opposite(*it_B, pmesh_B), pmesh_B), pmesh_B)));
    it_B = (it_B == last_B) ? bhv_B.begin() : std::next(it_B);
  }
  while(it_B != canon_B);

  border_out << "2 " << hole_points.back() << " " << sbp << std::endl;
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

  // add the missing quad
  patch.emplace_back(quad_id, 0, quad_id+1);
  patch.emplace_back(quad_id+1, 0, hole_points.size() - 1);

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

// @todo shouldn't the orientation always be correct by construction?
// Do this with combinatorics if it turns out it might not be
template <typename Patch,
          typename PolygonMesh,
          typename VPM,
          typename GeomTraits>
bool fix_patch_orientation(Patch& point_patch,
                           const typename boost::graph_traits<PolygonMesh>::halfedge_descriptor h1,
                           const typename boost::graph_traits<PolygonMesh>::halfedge_descriptor h2,
                           const PolygonMesh& pmesh,
                           const VPM vpm,
                           const GeomTraits& gt)
{
  typedef typename boost::property_traits<VPM>::reference                       Point_ref;
  typedef typename boost::property_traits<VPM>::value_type                      Point;

  bool ok_orientation_1 = false, ok_orientation_2 = false;

  const Point_ref h1sp = get(vpm, source(h1, pmesh));
  const Point_ref h1tp = get(vpm, target(h1, pmesh));
  const Point_ref h2sp = get(vpm, source(h2, pmesh));
  const Point_ref h2tp = get(vpm, target(h2, pmesh));

#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
  std::cout << "h1sp: " << h1sp << std::endl;
  std::cout << "h1tp: " << h1tp << std::endl;
  std::cout << "h2sp: " << h2sp << std::endl;
  std::cout << "h2tp: " << h2tp << std::endl;
#endif

  for(const auto& face : point_patch)
  {
    for(int i=0; i<3; ++i)
    {
      const Point& p1 = face[i];
      const Point& p2 = face[(i+1)%3];

      if(gt.equal_3_object()(p1, h1sp) && gt.equal_3_object()(p2, h1tp))
        ok_orientation_1 = true;
      if(gt.equal_3_object()(p1, h2sp) && gt.equal_3_object()(p2, h2tp))
        ok_orientation_2 = true;
    }
  }

  std::cout << "orientations: " << ok_orientation_1 << " " << ok_orientation_2 << std::endl;

  if(ok_orientation_1 != ok_orientation_2)
    return false;

  if(!ok_orientation_1)
  {
    for(auto& face : point_patch)
      std::swap(face[0], face[1]);
  }

  return true;
}

template <typename PolygonMesh,
          typename VPM,
          typename GeomTraits>
bool merge_holes(const typename boost::graph_traits<PolygonMesh>::halfedge_descriptor bh1,
                 const typename boost::graph_traits<PolygonMesh>::halfedge_descriptor bh2,
                 PolygonMesh& pmesh,
                 const VPM vpm,
                 const GeomTraits& gt)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor        vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor      halfedge_descriptor;

  typedef typename GeomTraits::Point_3                                        Point;

  CGAL_precondition(is_border(bh1, pmesh) && is_border(bh2, pmesh));

  const Iterator_range<Halfedge_around_face_iterator<PolygonMesh> > haf1 = CGAL::halfedges_around_face(bh1, pmesh);
  std::vector<halfedge_descriptor> bhv_1(std::begin(haf1), std::end(haf1));

  const Iterator_range<Halfedge_around_face_iterator<PolygonMesh> > haf2 = CGAL::halfedges_around_face(bh2, pmesh);
  std::vector<halfedge_descriptor> bhv_2(std::begin(haf2), std::end(haf2));

  // keep a correspondence map for when we will add the patch
  std::map<Point, vertex_descriptor> points_to_vs;
  for(const halfedge_descriptor h : bhv_1)
    points_to_vs[get(vpm, target(h, pmesh))] = target(h, pmesh);

  for(const halfedge_descriptor h : bhv_2)
    points_to_vs[get(vpm, target(h, pmesh))] = target(h, pmesh);

#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
  std::cout << "Merge holes (" << bh1 << " and " << bh2 << ")" << std::endl;
  std::cout << "holes of size: " << bhv_1.size() << " && " << bhv_2.size() << std::endl;
#endif

  // make sure that the holes are topological disks
  std::vector<std::vector<Point> > point_patch;
  point_patch.reserve(2 * bhv_1.size());

  bool success = two_borders_hole_fill(bhv_1, pmesh, vpm, bhv_2, pmesh, vpm,
                                       std::back_inserter(point_patch), gt);
  if(!success)
    return false;

#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
  dump_tentative_hole(point_patch, "results/tentative_patch.off");
  std::cout << point_patch.size() << " new faces" << std::endl;
#endif

  if(!check_patch_sanity<PolygonMesh>(point_patch))
  {
#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
    std::cout << "this patch is INSANE!" << std::endl;
#endif
    return false;
  }

  // @todo not useful?
  success = fix_patch_orientation(point_patch, bhv_1.front(), bhv_2.front(), pmesh, vpm, gt);
  if(!success)
  {
#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
    std::cout << "Incompatible orientations" << std::endl;
#endif
    return false; // can't find an orientation compatible with both borders
  }

  return add_faces(point_patch, pmesh, points_to_vs, vpm);
}

template <typename PolygonMesh, typename VPM, typename GeomTraits>
bool fill_hole(const typename boost::graph_traits<PolygonMesh>::halfedge_descriptor h,
               PolygonMesh& pmesh,
               VPM vpm,
               const GeomTraits& gt)
{
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor                    face_descriptor;

  CGAL_precondition(is_border(h, pmesh));

  std::vector<face_descriptor> patch_faces;
  triangulate_hole(pmesh, h, std::back_inserter(patch_faces),
                   parameters::vertex_point_map(vpm).geom_traits(gt));

  if(patch_faces.empty())
  {
#ifndef CGAL_HOLE_FILLING_DO_NOT_USE_DT3
#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
    std::cout << "Failed to fill a hole using Delaunay search space.\n";
#endif

    triangulate_hole(pmesh, h, std::back_inserter(patch_faces),
                     parameters::use_delaunay_triangulation(false).vertex_point_map(vpm).geom_traits(gt));

#endif // CGAL_HOLE_FILLING_DO_NOT_USE_DT3
    if(patch_faces.empty())
    {
#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
      std::cout << "Failed to fill a hole using the whole search space.\n";
#endif
      return false;
    }
  }

  return true;
}

template <typename HoleHalfedgeContainer, typename PolygonMesh, typename VPM, typename GeomTraits>
bool fill_holes(const HoleHalfedgeContainer& holes,
                PolygonMesh& pmesh,
                VPM vpm,
                const GeomTraits& gt)
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor                halfedge_descriptor;

  for(const halfedge_descriptor h : holes)
    if(!fill_hole(h, pmesh, vpm, gt))
      return false;

  return true;
}

template <typename UmbrellaContainer, typename PolygonMesh, typename VPM, typename GeomTraits>
bool treat_umbrellas(UmbrellaContainer& umbrellas,
                     const NM_TREATMENT treatment,
                     const typename GeomTraits::FT radius,
                     PolygonMesh& pmesh,
                     VPM vpm,
                     const GeomTraits& gt)
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor                halfedge_descriptor;

  // @todo
  // treat_pinched_vertices(umbrellas, radius, treatment, pmesh, vpm, gt);

  std::vector<halfedge_descriptor> holes;
  holes = dig_holes(umbrellas, radius, pmesh, vpm, gt);

  if(treatment == SEPARATE)
    return fill_holes(holes, pmesh, vpm, gt);

  if(holes.size() != 2)
  {
#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
    std::cout << "merging treatment requested, but configuration makes it impossible" << std::endl;
#endif
    return fill_holes(umbrellas, pmesh, vpm, gt);
  }
  else
  {
    const bool success = merge_holes(holes.front(), holes.back(), pmesh, vpm, gt);
    if(!success)
    {
#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
      std::cout << "merging failed, falling back to separating strategy" << std::endl;
#endif
      return fill_holes(holes, pmesh, vpm, gt);
    }

    return true;
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Collect

// pretty much the same as 'PMP::non_manifold_vertices()', but consider the geometry instead of the combinatorics
// ({combinatorial non-manifold vertices} <= {geometrical non-manifold vertices}) so that creates a few changes
template <typename NMPContainer, typename NMVM,
          typename PolygonMesh, typename VPM, typename GeomTraits>
void geometrically_non_manifold_vertices(NMPContainer& nm_points, // m[vertex] = {halfedges}
                                         NMVM nm_marks,
                                         const PolygonMesh& pmesh,
                                         const VPM vpm,
                                         const GeomTraits& /*gt*/)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor                  vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor                halfedge_descriptor;

  typedef typename boost::property_traits<VPM>::reference                               Point_ref;
  typedef typename boost::property_traits<VPM>::value_type                              Point;

  typedef CGAL::dynamic_halfedge_property_t<bool>                                       Halfedge_property_tag;
  typedef typename boost::property_map<PolygonMesh, Halfedge_property_tag>::const_type  Visited_halfedge_map;

  CGAL_precondition(nm_points.empty());

  Visited_halfedge_map visited_halfedges = get(Halfedge_property_tag(), pmesh);

  // @todo here again a could be map with a vertex as key
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

      // if this is the second time we visit that vertex and the first star was manifold, we have
      // not marked the vertex as non-manifold from the first star
      const auto nm_itb = nm_points.emplace(p, std::vector<halfedge_descriptor>{h});
      if(nm_itb.second) // successful insertion
      {
        const halfedge_descriptor h_from_another_star = visited_itb.first->second;
        nm_itb.first->second.push_back(h_from_another_star);
        put(nm_marks, target(h_from_another_star, pmesh), true);
      }
      else
      {
        nm_itb.first->second.push_back(h);
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
      nm_points[p].push_back(h); // might or might not have been an empty vector before
    }
  }
}

} // namespace internal

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

  typedef typename Geom_traits::FT                                            FT;
  typedef typename Geom_traits::Point_3                                       Point;

  const FT radius = 0.25; // @todo automatic or np

  // Collect the non-manifold vertices
  typedef CGAL::dynamic_vertex_property_t<bool>                               Mark;
  typedef typename boost::property_map<PolygonMesh, Mark>::type               Marked_vertices;
  Marked_vertices nm_marks = get(Mark(), pmesh);
  for(vertex_descriptor v : vertices(pmesh))
    put(nm_marks, v, false);

  typedef std::vector<halfedge_descriptor>                                    Cones;
  std::unordered_map<Point, Cones> nm_points;

  // @todo could be made lighter with a vertex as key, using a custom equal
  internal::geometrically_non_manifold_vertices(nm_points, nm_marks, pmesh, vpm, gt);

#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
  std::cout << nm_points.size() << " cases to treat" << std::endl;
#endif

  // If the treatment is merging, there will be combinatorics change and the cone halfedges
  // might become invalid. To ensure that they do stay invalid, we must ensure that no two
  // non-manifold vertices share an edge
  if(treatment == MERGE)
  {
    internal::enforce_non_manifold_vertex_separation(nm_marks, pmesh, vpm, gt);

    std::ofstream("results/comb_separated.off") << std::setprecision(17) << pmesh;
  }

  for(const auto& e : nm_points)
  {
    const std::vector<halfedge_descriptor>& umbrellas = e.second;
    CGAL_assertion(umbrellas.size() >= 2);
    CGAL_assertion(std::set<halfedge_descriptor>(umbrellas.begin(), umbrellas.end()).size() == umbrellas.size());

#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
    std::cout << "NM vertex with " << umbrellas.size() << " incident umbrellas:";
    for(const halfedge_descriptor h : umbrellas)
      std::cout << " " << h;
    std::cout << std::endl;
#endif

    internal::treat_umbrellas(umbrellas, treatment, radius, pmesh, vpm, gt);
    std::cin.get();
  }

  CGAL_postcondition_code(nm_points.clear();)
  CGAL_postcondition_code(internal::geometrically_non_manifold_vertices(nm_points, nm_marks, pmesh, vpm, gt);)
  CGAL_postcondition(nm_points.empty());
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
