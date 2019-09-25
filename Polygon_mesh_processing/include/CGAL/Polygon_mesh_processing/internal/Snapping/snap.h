// Copyright (c) 2018, 2019 GeometryFactory (France).
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
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Sebastien Loriot
//                 Mael Rouxel-Labbé

#ifndef CGAL_POLYGON_MESH_PROCESSING_INTERNAL_SNAP_H
#define CGAL_POLYGON_MESH_PROCESSING_INTERNAL_SNAP_H

#include <CGAL/license/Polygon_mesh_processing/repair.h>

#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>
#include <CGAL/Polygon_mesh_processing/border.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_halfedge_graph_segment_primitive.h>
#include <CGAL/assertions.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/box_intersection_d.h>
#include <CGAL/circulator.h>
#include <CGAL/Dynamic_property_map.h>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/number_utils.h>

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/member.hpp>

#ifdef CGAL_LINKED_WITH_TBB
# include <tbb/concurrent_vector.h>
# include <tbb/parallel_for.h>
#endif

#include <functional>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <type_traits>
#include <utility>
#include <unordered_set>
#include <vector>

namespace CGAL {
namespace Polygon_mesh_processing {
namespace internal {

// Assigns at each vertex the 'tolerance' value as tolerance, but bounded by a percentage of the length of its shortest incident edge
template <typename HalfedgeRange,
          typename ToleranceMap,
          typename PolygonMesh,
          typename SourceNamedParameters>
void assign_tolerance_with_local_edge_length_bound(const HalfedgeRange& hrange,
                                                   ToleranceMap& tol_pmap,
                                                   const typename GetGeomTraits<PolygonMesh, SourceNamedParameters>::type::FT tolerance,
                                                   PolygonMesh& mesh,
                                                   const SourceNamedParameters& snp)
{
  using parameters::get_parameter;
  using parameters::choose_parameter;

  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor                vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor              halfedge_descriptor;

  typedef typename GetVertexPointMap<PolygonMesh, SourceNamedParameters>::type        SVPM;
  typedef typename GetGeomTraits<PolygonMesh, SourceNamedParameters>::type            GT;
  typedef typename GT::FT                                                             FT;

  GT gt = choose_parameter(get_parameter(snp, internal_np::geom_traits), GT());
  SVPM svpm = choose_parameter(get_parameter(snp, internal_np::vertex_point),
                               get_property_map(vertex_point, mesh));

  for(halfedge_descriptor hd : hrange)
  {
    const vertex_descriptor vd = target(hd, mesh);
    CGAL::Halfedge_around_target_iterator<PolygonMesh> hit, hend;
    boost::tie(hit, hend) = CGAL::halfedges_around_target(vd, mesh);
    CGAL_assertion(hit != hend);

    FT sq_length = gt.compute_squared_distance_3_object()(get(svpm, source(*hit, mesh)),
                                                          get(svpm, target(*hit, mesh)));
    FT min_sq_dist = sq_length;
    ++hit;

    for(; hit!=hend; ++hit)
    {
      sq_length = gt.compute_squared_distance_3_object()(get(svpm, source(*hit, mesh)),
                                                         get(svpm, target(*hit, mesh)));

      if(sq_length < min_sq_dist)
        min_sq_dist = sq_length;
    }

#ifdef CGAL_PMP_SNAP_DEBUG
    std::cout << "tolerance at vd: " /*<< vd */ << " [" << get(svpm, vd) << "]: min of "
              << 0.9 * CGAL::approximate_sqrt(min_sq_dist) << " AND " << tolerance << std::endl;
#endif
    put(tol_pmap, vd, CGAL::min<FT>(0.9 * CGAL::approximate_sqrt(min_sq_dist), tolerance));
  }
}

template <typename HalfedgeRange,
          typename ToleranceMap,
          typename PolygonMesh>
void assign_tolerance_with_local_edge_length_bound(const HalfedgeRange& hrange,
                                                   ToleranceMap& tol_pmap,
                                                   const typename GetGeomTraits<PolygonMesh>::type::FT tolerance,
                                                   PolygonMesh& mesh)
{
  return assign_tolerance_with_local_edge_length_bound(hrange, tol_pmap, tolerance, mesh, CGAL::parameters::all_default());
}

template <typename PolygonMesh, typename GeomTraits>
struct Snapping_pair
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor       vertex_descriptor;
  typedef typename GeomTraits::FT                                            FT;

  Snapping_pair(const vertex_descriptor vs_, const vertex_descriptor vt_, const FT sq_dist_)
    : vs(vs_), vt(vt_), sq_dist(sq_dist_)
  { }

  vertex_descriptor vs;
  vertex_descriptor vt;
  FT sq_dist;
};

template <typename PolygonMesh, typename GeomTraits,
          typename DistanceMultiIndexContainer,
          typename SVPM, typename TVPM,
          typename ToleranceMap,
          typename Box>
struct Vertex_proximity_report
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor       vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::edge_descriptor         edge_descriptor;

  typedef typename GeomTraits::FT                                            FT;
  typedef typename boost::property_traits<SVPM>::value_type                  Point;

  Vertex_proximity_report(DistanceMultiIndexContainer& snapping_pairs,
                          const SVPM& svpm, const PolygonMesh& smesh,
                          const TVPM& tvpm, const PolygonMesh& tmesh,
                          const ToleranceMap& tol_pmap,
                          const GeomTraits& gt)
    :
      m_snapping_pairs(snapping_pairs),
      is_same_mesh((&smesh == &tmesh)),
      svpm(svpm), tvpm(tvpm),
      tol_pmap(tol_pmap),
      gt(gt)
  { }

  void operator()(const Box& a, const Box& b)
  {
    vertex_descriptor va = a.info();
    vertex_descriptor vb = b.info();

    if(is_same_mesh && va == vb)
      return;

    const Point sp = get(svpm, va);
    const Point tp = get(tvpm, vb);
    const FT tol = get(tol_pmap, va);

    // Don't reject a '0' distance, it still needs to lock the points in place
    const FT sq_dist = gt.compute_squared_distance_3_object()(sp, tp);
    CGAL::Comparison_result res = CGAL::compare(sq_dist, tol * tol);

#ifdef CGAL_PMP_SNAP_DEBUG_PP
    std::cout << "distance between " /*<< va*/ << " [" << sp << "] and "
                                     /*<< vb*/ << " [" << tp << "]: " << sq_dist
                                     << " (bound: " << tol*tol << ") larger? " << (res == CGAL::LARGER)
                                     << std::endl;
#endif

    if(res == CGAL::LARGER)
      return;

    m_snapping_pairs.insert(Snapping_pair<PolygonMesh, GeomTraits>(va, vb, sq_dist));
  }

private:
  DistanceMultiIndexContainer& m_snapping_pairs;

  const bool is_same_mesh;
  const SVPM& svpm;
  const TVPM& tvpm;
  const ToleranceMap& tol_pmap;
  const GeomTraits& gt;
};

} // namespace internal

namespace experimental {

// This is the function if you know what you're doing with the ranges
//
// \ingroup PMP_repairing_grp
//
// Attempts to snap the vertices in `source_hrange` onto the vertices in `target_hrange`.
// A vertex of the source range is only snapped to a vertex of the target mesh
// if its distance to the target mesh vertex is smaller than a user-chosen bound.
// If a source vertex can be snapped onto multiple vertices of the target
// range, the closest target vertex is used.
// If multiple vertices within the source range could be snapped onto the same target vertex,
// the source vertex closest to the target vertex will be snapped. Other source vertices will try to snap
// to free target vertices until they find a valid target, or until there is no more target
// to snap to (within the tolerance).
//
// @warning This function does not give any guarantee on the conformity between the source and target meshes after the snapping.
// @warning This function does not merge vertices together.
//
// @tparam PolygonMesh a model of `FaceListGraph` and `MutableFaceGraph`
// @tparam SourceHalfedgeRange a model of `Range` with value type `boost::graph_traits<PolygonMesh>::%halfedge_descriptor`
// @tparam TargetHalfedgeRange a model of `Range` with value type `boost::graph_traits<PolygonMesh>::%halfedge_descriptor`
// @tparam ToleranceMap a model of `ReadablePropertyMap` with key type `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
//                      and value type `GetGeomTraits<PolygonMesh, SourceNamedParameters>::type::FT`
// @tparam SourceNamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
// @tparam TargetNamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
//
// @param source_hrange a range of vertices of the source mesh whose positions can be changed.
//                      the vertices must be border vertices of `smesh`.
// @param smesh the source mesh whose border vertices might be moved
// @param target_hrange a range of vertices of the target mesh which are potential new positions
//                      for the vertices in the source range
// @param tmesh the target mesh to which the vertices in `target_hrange` belong
// @param tol_pmap a tolerance map associating to each vertex of the source range a tolerance value:
//               potential projection targets are sought in a sphere centered at the vertex and
//               whose radius is the tolerance value.
// @param snp optional \ref pmp_namedparameters "Named Parameters" related to the source mesh,
//            amongst those described below:
//
// \cgalNamedParamsBegin
//    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of the source mesh.
//                                      The type of this map is model of `ReadWritePropertyMap`.
//                                      If this parameter is omitted, an internal property map for
//                                      `CGAL::vertex_point_t` must be available in `PolygonMesh`
//    \cgalParamEnd
//    \cgalParamBegin{geom_traits} a geometric traits class instance.
//       The traits class must provide the nested types `Point_3` and `Vector_3`,
//       and the nested functors :
//         - `Construct_bbox_3` to construct a bounding box of a point,
//         - `Compute_squared_distance_3` to compute the distance between two points,
//
//       and, for each functor `Foo`, a function `Foo foo_object()`
//   \cgalParamEnd
// \cgalNamedParamsEnd
//
// @param tnp optional \ref pmp_namedparameters "Named Parameters" related to the target mesh,
//            amongst those described below:
//
// \cgalNamedParamsBegin
//    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of the target mesh.
//                                      The type of this map is model of `ReadablePropertyMap`.
//                                      If this parameter is omitted, an internal property map for
//                                      `CGAL::vertex_point_t` must be available in `PolygonMesh`
//    \cgalParamEnd
// \cgalNamedParamsEnd
//
// @return the number of snapped vertices
//
// @sa `merge_duplicated_vertices_in_boundary_cycles()`
//
template <typename PolygonMesh,
          typename SourceHalfedgeRange, typename TargetHalfedgeRange,
          typename ToleranceMap,
          typename SourceNamedParameters, typename TargetNamedParameters>
std::size_t snap_vertex_range_onto_vertex_range(const SourceHalfedgeRange& source_hrange,
                                                PolygonMesh& smesh,
                                                const TargetHalfedgeRange& target_hrange,
                                                const PolygonMesh& tmesh,
                                                const ToleranceMap& tol_pmap,
                                                const SourceNamedParameters& snp,
                                                const TargetNamedParameters& tnp)
{
  using parameters::get_parameter;
  using parameters::choose_parameter;

  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor      vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor    halfedge_descriptor;

  typedef CGAL::Box_intersection_d::Box_with_info_d<double, 3, vertex_descriptor>     Box;

  typedef typename GetVertexPointMap<PolygonMesh, SourceNamedParameters>::type        SVPM;
  typedef typename GetVertexPointMap<PolygonMesh, TargetNamedParameters>::const_type  TVPM;
  typedef typename boost::property_traits<TVPM>::value_type                           Point;

  typedef typename GetGeomTraits<PolygonMesh, SourceNamedParameters>::type            GT;
  typedef typename GT::FT                                                             FT;

  if(is_empty_range(source_hrange.begin(), source_hrange.end()) ||
     is_empty_range(target_hrange.begin(), target_hrange.end()))
    return 0;

  CGAL_static_assertion((std::is_same<Point, typename GT::Point_3>::value));

  GT gt = choose_parameter(get_parameter(snp, internal_np::geom_traits), GT());

  SVPM svpm = choose_parameter(get_parameter(snp, internal_np::vertex_point),
                               get_property_map(vertex_point, smesh));
  TVPM tvpm = choose_parameter(get_parameter(tnp, internal_np::vertex_point),
                               get_const_property_map(vertex_point, tmesh));

#ifdef CGAL_PMP_SNAP_DEBUG
  std::cout << "Snapping vertices to vertices. Range sizes: "
            << std::distance(source_hrange.begin(), source_hrange.end()) << " and "
            << std::distance(target_hrange.begin(), target_hrange.end()) << std::endl;

  if(&smesh == &tmesh)
    std::cout << "same mesh!" << std::endl;
#endif

  // Try to snap vertices
  std::vector<Box> boxes;
  std::unordered_set<vertex_descriptor> unique_vertices;
  for(halfedge_descriptor hd : source_hrange)
  {
    const vertex_descriptor vd = target(hd, smesh);
    if(!unique_vertices.insert(vd).second)
      continue; // if 'vd' appears multiple times on the border, move it only once

    // only making the box a little larger, but the final tolerance is not changed
    const double eps = 1.01 * CGAL::to_double(get(tol_pmap, vd));

    const Bbox_3 pb = gt.construct_bbox_3_object()(get(svpm, vd));
    const Bbox_3 b(pb.xmin() - eps, pb.ymin() - eps, pb.zmin() - eps,
                   pb.xmax() + eps, pb.ymax() + eps, pb.zmax() + eps);
    boxes.push_back(Box(b, vd));
  }

  std::vector<Box> target_boxes;
  for(halfedge_descriptor hd : target_hrange)
  {
    const vertex_descriptor vd = target(hd, tmesh);
    const Point& p = get(tvpm, vd);
    target_boxes.push_back(Box(gt.construct_bbox_3_object()(p), vd));
  }

  // Use a multi index to sort easily by sources, targets, AND distance.
  // Then, look up the distances in increasing order, and snap whenever the source and the target
  // have both not been snapped yet.
  typedef internal::Snapping_pair<PolygonMesh, GT>                                    Snapping_pair;
  typedef boost::multi_index::multi_index_container<
    Snapping_pair,
    boost::multi_index::indexed_by<
      boost::multi_index::ordered_non_unique<
        BOOST_MULTI_INDEX_MEMBER(Snapping_pair, vertex_descriptor, vs)>,
      boost::multi_index::ordered_non_unique<
        BOOST_MULTI_INDEX_MEMBER(Snapping_pair, vertex_descriptor, vt)>,
      boost::multi_index::ordered_non_unique<
        BOOST_MULTI_INDEX_MEMBER(Snapping_pair, FT, sq_dist)>
    >
  >                                                                                   Snapping_pair_container;

  typedef internal::Vertex_proximity_report<PolygonMesh, GT, Snapping_pair_container,
                                            SVPM, TVPM, ToleranceMap, Box>            Reporter;

  Snapping_pair_container snapping_pairs;
  Reporter vpr(snapping_pairs, svpm, smesh, tvpm, tmesh, tol_pmap, gt);

  // Shenanigans to pass a reference as callback (which is copied by value by 'box_intersection_d')
  std::function<void(const Box&, const Box&)> callback(std::ref(vpr));

  CGAL::box_intersection_d(boxes.begin(), boxes.end(),
                           target_boxes.begin(), target_boxes.end(),
                           callback);

#ifdef CGAL_PMP_SNAP_DEBUG
    std::cout << snapping_pairs.size() << " snappable pair(s)!" << std::endl;
#endif

  if(snapping_pairs.empty())
    return 0;

  // Sorted views of the container
  typedef typename Snapping_pair_container::template nth_index<0>::type   Container_by_source;
  typedef typename Snapping_pair_container::template nth_index<1>::type   Container_by_target;
  typedef typename Snapping_pair_container::template nth_index<2>::type   Container_by_distance;

  Container_by_source& container_by_source = snapping_pairs.template get<0>();
  Container_by_target& container_by_target = snapping_pairs.template get<1>();
  Container_by_distance& container_by_dist = snapping_pairs.template get<2>();

  std::size_t counter = 0;

  CGAL_assertion_code(FT prev = -1;)
  while(!container_by_dist.empty())
  {
    const Snapping_pair& sp = *(container_by_dist.begin());
    const vertex_descriptor vs = sp.vs;
    const vertex_descriptor vt = sp.vt;
    CGAL_assertion(sp.sq_dist >= prev);
    CGAL_assertion_code(prev = sp.sq_dist;)

#ifdef CGAL_PMP_SNAP_DEBUG_PP
    std::cout << "Snapping " /*<< vs*/ << " (" << get(svpm, vs) << ") "
              << " to " /*<< vt*/ << " (" << get(tvpm, vt) << ") at dist: " << sp.sq_dist << std::endl;
#endif

    // Collect all the source vertices projecting onto that target vertex
   ++counter;
    put(svpm, vs, get(tvpm, vt));

    // 'vs' and 'vt' cannot be used anymore, remove them from the container
    container_by_source.erase(vs);
    container_by_target.erase(vt);
  }

#ifdef CGAL_PMP_SNAP_DEBUG
    std::cout << "Snapped " << counter << " pair(s)!" << std::endl;
#endif

  return counter;
}

template <typename PolygonMesh, typename SourceHalfedgeRange, typename TargetHalfedgeRange,
          typename SourceNamedParameters, typename TargetNamedParameters>
std::size_t snap_vertex_range_onto_vertex_range(const SourceHalfedgeRange& source_hrange,
                                                PolygonMesh& smesh,
                                                const TargetHalfedgeRange& target_hrange,
                                                const PolygonMesh& tmesh,
                                                const SourceNamedParameters& snp,
                                                const TargetNamedParameters& tnp)
{
  typedef typename GetGeomTraits<PolygonMesh, SourceNamedParameters>::type            GT;
  typedef typename GT::FT                                                             FT;
  typedef CGAL::dynamic_vertex_property_t<FT>                                         Vertex_property_tag;
  typedef typename boost::property_map<PolygonMesh, Vertex_property_tag>::type        Tolerance_map;

  Tolerance_map tol_pmap = get(Vertex_property_tag(), smesh);
  const FT tol_mx(std::numeric_limits<double>::max());
  internal::assign_tolerance_with_local_edge_length_bound(source_hrange, tol_pmap, tol_mx, smesh, snp);

  return snap_vertex_range_onto_vertex_range(source_hrange, smesh, target_hrange, tmesh, tol_pmap, snp, tnp);
}

template <typename PolygonMesh, typename SourceHalfedgeRange, typename TargetHalfedgeRange, typename ToleranceMap>
std::size_t snap_vertex_range_onto_vertex_range(const SourceHalfedgeRange& source_hrange,
                                                PolygonMesh& smesh,
                                                const TargetHalfedgeRange& target_hrange,
                                                const PolygonMesh& tmesh,
                                                const ToleranceMap& tol_pmap)
{
  return snap_vertex_range_onto_vertex_range(source_hrange, smesh, target_hrange, tmesh, tol_pmap,
                                             CGAL::parameters::all_default(),
                                             CGAL::parameters::all_default());
}

template <typename PolygonMesh, typename SourceHalfedgeRange, typename TargetHalfedgeRange>
std::size_t snap_vertex_range_onto_vertex_range(const SourceHalfedgeRange& source_hrange,
                                                PolygonMesh& smesh,
                                                const TargetHalfedgeRange& target_hrange,
                                                const PolygonMesh& tmesh)
{
  return snap_vertex_range_onto_vertex_range(source_hrange, smesh, target_hrange, tmesh,
                                             CGAL::parameters::all_default(),
                                             CGAL::parameters::all_default());
}

template <typename PolygonMesh>
std::size_t snap_vertex_range_onto_vertex_range(PolygonMesh& smesh, const PolygonMesh& tmesh)
{
  return snap_vertex_range_onto_vertex_range(halfedges(smesh), smesh, halfedges(tmesh), tmesh);
}

template <typename PolygonMesh, typename HalfedgeRange, typename ToleranceMap>
std::size_t snap_border_vertices_onto_vertex_range(PolygonMesh& smesh,
                                                   const HalfedgeRange& target_hrange,
                                                   const PolygonMesh& tmesh,
                                                   const ToleranceMap& tol_pmap)
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor    halfedge_descriptor;

  std::vector<halfedge_descriptor> border_vertices;
  border_halfedges(smesh, std::back_inserter(border_vertices));

  return snap_vertex_range_onto_vertex_range(border_vertices, smesh, target_hrange, tmesh, tol_pmap);
}

template <typename PolygonMesh, typename HalfedgeRange>
std::size_t snap_border_vertices_onto_vertex_range(PolygonMesh& smesh,
                                                   const HalfedgeRange& target_hrange,
                                                   const PolygonMesh& tmesh)
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor    halfedge_descriptor;

  std::vector<halfedge_descriptor> border_vertices;
  border_halfedges(smesh, std::back_inserter(border_vertices));

  return snap_vertex_range_onto_vertex_range(border_vertices, smesh, target_hrange, tmesh);
}

template <typename PolygonMesh, typename ToleranceMap>
std::size_t snap_border_vertices_onto_border_vertices(PolygonMesh& smesh,
                                                      const PolygonMesh& tmesh,
                                                      const ToleranceMap& tol_pmap)
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor    halfedge_descriptor;

  std::vector<halfedge_descriptor> sborder_vertices;
  border_halfedges(smesh, std::back_inserter(sborder_vertices));

  std::vector<halfedge_descriptor> tborder_vertices;
  border_halfedges(tmesh, std::back_inserter(tborder_vertices));

  return snap_vertex_range_onto_vertex_range(sborder_vertices, smesh, tborder_vertices, tmesh, tol_pmap);
}

// \ingroup PMP_repairing_grp
//
// Attempts to snap the border vertices of the source mesh onto the vertices of the target mesh.
//
// A vertex of the source range is only snapped to a vertex of the target mesh
// if its distance to the target mesh vertex is smaller than a user-chosen bound.
// If any source vertex can be snapped onto multiple vertices of the target
// range, the closest one is chosen.
// If multiple vertices within the source range could be snapped onto the same target vertex,
// the source vertex closest to the target vertex will be snapped. Other source vertices will try to snap
// to free target vertices until they find a valid target, or until there is no more target
// to snap to (within the tolerance).
//
// @tparam PolygonMesh a model of `FaceListGraph` and `MutableFaceGraph`
// @tparam ToleranceMap a model of `ReadablePropertyMap` with key type `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
//                      and value type the number type associated with the traits of the mesh.
//
// @param smesh the source mesh whose border vertices might be moved
// @param tmesh the target mesh whose vertices are potential projection targets
// @param tol_pmap a tolerance map associating to each vertex of the source range a tolerance value:
//               potential projection targets are sought in a sphere centered at the vertex and
//               whose radius is the tolerance value.
//
// @pre `smesh` and `tmesh` are different meshes
//
// \return the number of snapped vertices
//
template <typename PolygonMesh, typename ToleranceMap>
std::size_t snap_border_vertices_onto_vertex_range(PolygonMesh& smesh,
                                                   const PolygonMesh& tmesh,
                                                   const ToleranceMap& tol_pmap)
{
  return snap_border_vertices_onto_vertex_range(smesh, halfedges(tmesh), tmesh, tol_pmap);
}

template <typename PolygonMesh>
std::size_t snap_border_vertices_onto_vertex_range(PolygonMesh& smesh, const PolygonMesh& tmesh)
{
  return snap_border_vertices_onto_vertex_range(smesh, halfedges(tmesh), tmesh);
}

} // namespace experimental

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

namespace internal {

// Adapted from <CGAL/internal/AABB_tree/AABB_traversal_traits.h>
template <typename AABBTraits>
class Projection_traits
{
  typedef typename AABBTraits::FT                               FT;
  typedef typename AABBTraits::Point_3                          Point;
  typedef typename AABBTraits::Primitive                        Primitive;
  typedef typename AABBTraits::Bounding_box                     Bounding_box;
  typedef typename AABBTraits::Primitive::Id                    Primitive_id;
  typedef typename AABBTraits::Point_and_primitive_id           Point_and_primitive_id;
  typedef typename AABBTraits::Object_and_primitive_id          Object_and_primitive_id;

  typedef CGAL::AABB_node<AABBTraits>                           Node;

public:
  explicit Projection_traits(const AABBTraits& traits,
                             const FT sq_tolerance)
    : m_traits(traits),
      m_closest_point_initialized(false),
      m_sq_dist(sq_tolerance)
  {}

  bool go_further() const { return true; }

  void intersection(const Point& query, const Primitive& primitive)
  {
    // skip the primitive if one of its endpoints is the query
    const typename Primitive::Datum& s = primitive.datum(m_traits.shared_data());
    if(m_traits.equal_3_object()(s[0], query) || m_traits.equal_3_object()(s[1], query))
      return;

    if(!m_closest_point_initialized)
    {
      m_closest_point_initialized = true;
      m_closest_primitive = primitive.id();
      m_closest_point = primitive.reference_point(m_traits.shared_data());
      m_sq_dist = m_traits.squared_distance_object()(query, m_closest_point);
    }

    Point new_closest_point = m_traits.closest_point_object()(query, primitive, m_closest_point);
    if(!m_traits.equal_3_object()(new_closest_point, m_closest_point))
    {
      m_closest_primitive = primitive.id();
      m_closest_point = new_closest_point;
      m_sq_dist = m_traits.squared_distance_object()(query, m_closest_point);
    }
  }

  bool do_intersect(const Point& query, const Node& node) const
  {
    return m_traits.compare_distance_object()(query, node.bbox(), m_sq_dist) == CGAL::SMALLER;
  }

  const Point& closest_point() const { return m_closest_point; }
  typename Primitive::Id closest_primitive_id() const { return m_closest_primitive; }
  bool closest_point_initialized() const { return m_closest_point_initialized; }

private:
  const AABBTraits& m_traits;
  bool m_closest_point_initialized;
  Point m_closest_point;
  FT m_sq_dist;

  typename Primitive::Id m_closest_primitive;
};

// The UniqueVertex is a pair of a container of vertex_descriptor and FT, representing
// vertices with the same geometric position and their associated snapping tolerance
// (defined as the minimum of the tolerances of the vertices of the bunch)
template <typename UniqueVertex, typename TriangleMesh,
          typename EdgeToSplitMap, typename AABBTree, typename VPM, typename GT>
void find_splittable_edge(const UniqueVertex& unique_vertex,
                          EdgeToSplitMap& edges_to_split,
                          const AABBTree* aabb_tree_ptr,
                          const TriangleMesh& /* pms */,
                          const VPM& vpms,
                          const TriangleMesh& pmt,
                          const VPM& vpmt,
                          const GT& gt)
{
  CGAL_USE(vpmt);

  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor   halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::edge_descriptor       edge_descriptor;

  typedef typename GT::FT                                                   FT;
  typedef typename boost::property_traits<VPM>::value_type                  Point;

  typedef typename AABBTree::AABB_traits                                    AABB_traits;

  CGAL_assertion(!unique_vertex.first.empty());

  const Point& query = get(vpms, unique_vertex.first.front()); // by construction the whole range has the same position
  const FT sq_tolerance = CGAL::square(unique_vertex.second);

  // use the source halfedge as hint
  internal::Projection_traits<AABB_traits> traversal_traits(aabb_tree_ptr->traits(), sq_tolerance);
  aabb_tree_ptr->traversal(query, traversal_traits);

  if(!traversal_traits.closest_point_initialized())
    return;

  const Point& closest_p = traversal_traits.closest_point();

  // The filtering in the AABB tree checks the dist query <-> node bbox, which might be smaller than
  // the actual distance between the query <-> closest point
  const FT sq_dist_to_closest = gt.compute_squared_distance_3_object()(query, closest_p);
  bool is_close_enough = (sq_dist_to_closest <= sq_tolerance);

#ifdef CGAL_PMP_SNAP_DEBUG_PP
  std::cout << "  Closest point: (" << closest_p << ")" << std::endl
            << "  at distance " << gt.compute_squared_distance_3_object()(query, closest_p)
            << " with squared tolerance: " << sq_tolerance
            << " && close enough? " << is_close_enough << ")" << std::endl;
#endif

  if(!is_close_enough)
    return;

#ifdef CGAL_PMP_SNAP_DEBUG_PP
  std::cout << "\t and is thus beneath tolerance" << std::endl;
#endif

  edge_descriptor closest_e = traversal_traits.closest_primitive_id();
  CGAL_assertion(get(vpmt, source(closest_e, pmt)) != query &&
                 get(vpmt, target(closest_e, pmt)) != query);

  halfedge_descriptor closest_h = halfedge(closest_e, pmt);
  CGAL_assertion(is_border(edge(closest_h, pmt), pmt));

  if(!is_border(closest_h, pmt))
    closest_h = opposite(closest_h, pmt);

#ifdef CGAL_PMP_SNAP_DEBUG_PP
  std::cout << "splitting position: " << query << std::endl;
#endif

  // a map because we need to know if the same halfedge is split multiple times
#ifdef CGAL_LINKED_WITH_TBB
  typename EdgeToSplitMap::accessor acc;
  edges_to_split.insert(acc, closest_h);
  acc->second.emplace_back(&(unique_vertex.first), closest_p);
#else
  edges_to_split[closest_h].emplace_back(&(unique_vertex.first), closest_p);
#endif
}

#ifdef CGAL_LINKED_WITH_TBB
template <typename PointWithToleranceContainer,
          typename TriangleMesh, typename EdgeToSplitMap,
          typename AABBTree, typename VPM, typename GT>
struct Find_splittable_edge_for_parallel_for
{
  Find_splittable_edge_for_parallel_for(const PointWithToleranceContainer& points_with_tolerance,
                                        EdgeToSplitMap& edges_to_split,
                                        const AABBTree* aabb_tree_ptr,
                                        const TriangleMesh& pms,
                                        const VPM& vpms,
                                        const TriangleMesh& pmt,
                                        const VPM& vpmt,
                                        const GT& gt)
    :
      m_points_with_tolerance(points_with_tolerance),
      m_edges_to_split(edges_to_split), m_aabb_tree_ptr(aabb_tree_ptr),
      m_pms(pms), m_vpms(vpms), m_pmt(pmt), m_vpmt(vpmt), m_gt(gt)
  { }

  void operator()(const tbb::blocked_range<size_t>& r) const
  {
    for(std::size_t i=r.begin(); i!=r.end(); ++i)
    {
      find_splittable_edge(m_points_with_tolerance[i], m_edges_to_split,
                           m_aabb_tree_ptr,
                           m_pms, m_vpms, m_pmt, m_vpmt,
                           m_gt);
    }
  }

private:
  const PointWithToleranceContainer& m_points_with_tolerance;
  EdgeToSplitMap& m_edges_to_split;
  const AABBTree* m_aabb_tree_ptr;
  const TriangleMesh& m_pms;
  const VPM m_vpms;
  const TriangleMesh& m_pmt;
  const VPM m_vpmt;
  const GT& m_gt;
};
#endif

} // namespace internal

namespace experimental {

// Try to merge vertices in 'hrange' (vd = target(hd, pm)) onto the edges in the 'erange',
// non-conformingly
template <typename HalfedgeRange, typename TriangleMesh,
          typename ToleranceMap, typename NamedParameters>
std::size_t snap_vertex_range_onto_vertex_range_non_conforming(const HalfedgeRange& source_hrange,
                                                               TriangleMesh& pms,
                                                               const HalfedgeRange& target_hrange,
                                                               TriangleMesh& pmt,
                                                               const ToleranceMap& tol_pmap,
                                                               const NamedParameters& nps,
                                                               const NamedParameters& npt)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor     vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor   halfedge_descriptor;

  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::type   VPM;
  typedef typename boost::property_traits<VPM>::value_type                  Point;
  typedef typename boost::property_traits<VPM>::reference                   Point_ref;

  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type       GT;
  typedef typename GT::FT                                                   FT;

  typedef CGAL::AABB_halfedge_graph_segment_primitive<TriangleMesh, VPM>    Primitive;
  typedef CGAL::AABB_traits<GT, Primitive>                                  AABB_Traits;
  typedef CGAL::AABB_tree<AABB_Traits>                                      AABB_tree;

  using parameters::get_parameter;
  using parameters::choose_parameter;

  // Need to have a triangle mesh to ensure that the refinement point of a border edge 'e' has visibility
  // of the third p oint of the face 'f' incident to 'e'
  CGAL_precondition(is_triangle_mesh(pmt));

  std::size_t snapped_n = 0;

  VPM vpms = choose_parameter(get_parameter(nps, internal_np::vertex_point), get_property_map(vertex_point, pms));
  VPM vpmt = choose_parameter(get_parameter(npt, internal_np::vertex_point), get_property_map(vertex_point, pmt));

  const GT gt = choose_parameter(get_parameter(nps, internal_np::geom_traits), GT());

  // start by snapping vertices together to simplify things
  snapped_n = snap_vertex_range_onto_vertex_range(source_hrange, pms, target_hrange, pmt, tol_pmap, nps, npt);

#ifdef CGAL_PMP_SNAP_DEBUG
  std::cout << "Snapping vertices to edges. Range sizes: "
            << std::distance(source_hrange.begin(), source_hrange.end()) << " and "
            << std::distance(target_hrange.begin(), target_hrange.end()) << std::endl;
#endif

  // Collect border points that can be projected onto a border edge
  //
  // If multiple vertices on the source border have the same geometric position, they are moved
  // all at once and we use the smallest tolerance from all source vertices with that position

  // Pair of:
  // - vertices with the same geometric position
  // - the smallest tolerance for the corresponding geometric position
  typedef std::vector<vertex_descriptor>                                      Vertex_container;
  typedef std::pair<Vertex_container, FT>                                     Unique_vertex;
  typedef std::vector<Unique_vertex>                                          Unique_vertices_with_tolerance;

  // Map[unique geometric position, iterator of unique vertex]
  // This is not simply a map[unique_vertices, tolerance] sorted on the geometric position
  // because we need to use a vector in the tbb::parallel_for()
  typedef typename Unique_vertices_with_tolerance::iterator                   Unique_vertex_iterator;
  typedef std::map<Point, Unique_vertex_iterator>                             Unique_positions_with_iterator;

#ifdef CGAL_PMP_SNAP_DEBUG
  std::cout << "Gather unique points in source range..." << std::endl;
#endif

  Unique_positions_with_iterator unique_positions;

  Unique_vertices_with_tolerance unique_vertices_with_tolerance;
  unique_vertices_with_tolerance.reserve(source_hrange.size()); // ensures that iterators stay valid

  for(halfedge_descriptor hd : source_hrange)
  {
    // Skip the source vertex if its two incident halfedges are geometrically identical (it means that
    // the two halfedges are already stitchable and we don't want this common vertex to be used
    // to split a halfedge somewhere else)
    if(get(vpms, source(hd, pms)) == get(vpms, target(next(hd, pms), pms)))
      continue;

    const vertex_descriptor v = target(hd, pms);
    const Point_ref query = get(vpms, v);
    const FT tolerance = get(tol_pmap, v);

#ifdef CGAL_PMP_SNAP_DEBUG_PP
      std::cout << "Query: " << v << " (" << query << ")" << std::endl;
#endif

    // The iterator will be set after insertion
    std::pair<typename Unique_positions_with_iterator::iterator, bool> is_insert_successful =
      unique_positions.insert(std::make_pair(query, Unique_vertex_iterator()));

    if(is_insert_successful.second) // successful insertion <=> never met this point
    {
      std::vector<vertex_descriptor> vv {{ v }};
      unique_vertices_with_tolerance.emplace_back(vv, tolerance);
      is_insert_successful.first->second = std::prev(unique_vertices_with_tolerance.end());
    }
    else // point was already met
    {
      CGAL_assertion(is_insert_successful.first->second != Unique_vertex_iterator());
      Unique_vertex& uv = *(is_insert_successful.first->second);

      std::vector<vertex_descriptor>& vv = uv.first;
      CGAL_assertion(std::find(vv.begin(), vv.end(), v) == vv.end());
      vv.push_back(v);

      uv.second = (std::min)(uv.second, tolerance);
    }
  }

  // Since we're inserting primitives one by one, we can't pass this shared data in the constructor of the tree
  AABB_Traits aabb_traits;
  aabb_traits.set_shared_data(pmt, vpmt);
  AABB_tree aabb_tree(aabb_traits);

  for(halfedge_descriptor hd : target_hrange)
  {
    CGAL_precondition(is_border(edge(hd, pmt), pmt));
    aabb_tree.insert(Primitive(edge(hd, pmt), pmt, vpmt));
  }

  // Now, check which edges of the target range ought to be split by source vertices
#ifdef CGAL_PMP_SNAP_DEBUG_PP
  std::cout << "Collect edges to split w/ " << unique_vertices_with_tolerance.size() << " unique vertices" << std::endl;
#endif

  typedef std::pair<const Vertex_container*, Point>                               Vertices_with_new_position;
  typedef std::vector<Vertices_with_new_position>                                 VNP_container;

#ifdef CGAL_LINKED_WITH_TBB
 #ifdef CGAL_PMP_SNAP_DEBUG
  std::cout << "Parallel find splittable edges!" << std::endl;
 #endif

  typedef tbb::concurrent_hash_map<halfedge_descriptor, VNP_container>            Concurrent_edge_to_split_container;
  typedef internal::Find_splittable_edge_for_parallel_for<Unique_vertices_with_tolerance,
                                                          TriangleMesh,
                                                          Concurrent_edge_to_split_container,
                                                          AABB_tree, VPM, GT>     Functor;

  Concurrent_edge_to_split_container edges_to_split;
  Functor f(unique_vertices_with_tolerance, edges_to_split, &aabb_tree, pms, vpms, pmt, vpmt, gt);
  tbb::parallel_for(tbb::blocked_range<std::size_t>(0, unique_vertices_with_tolerance.size()), f);
#else // CGAL_LINKED_WITH_TBB
 #ifdef CGAL_PMP_SNAP_DEBUG
  std::cout << "Sequential find splittable edges!" << std::endl;
 #endif

  std::map<halfedge_descriptor, VNP_container> edges_to_split; // @todo hash ?
  for(const Unique_vertex& uv : unique_vertices_with_tolerance)
    internal::find_splittable_edge(uv, edges_to_split, &aabb_tree, pms, vpms, pmt, vpmt, gt);
#endif // CGAL_LINKED_WITH_TBB

#ifdef CGAL_PMP_SNAP_DEBUG
  std::cout << "split " << edges_to_split.size() << " edges" << std::endl;
#endif

  typedef std::pair<const halfedge_descriptor, VNP_container>                     Splitters;
  for(Splitters& splitter : edges_to_split)
  {
    const halfedge_descriptor h = splitter.first;
    CGAL_assertion(is_border(h, pmt));

    VNP_container& splitters = splitter.second;

    if(splitters.size() > 1)
    {
#ifdef CGAL_PMP_SNAP_DEBUG_PP
      std::cout << " MUST SORT ON BORDER!" << std::endl;
#endif
      /// \todo Sorting the projected positions is too simple (for example, this doesn't work
      ///       for a zigzaging polyline). Rather, points should be sorted following
      ///       the order along the matching border. Note that this requires identifying
      ///       matching polylines because a sequence of zigzaging points might not all project
      ///       onto the same halfedge.

      const Point_ref hsp = get(vpmt, source(h, pmt));
      std::sort(splitters.begin(), splitters.end(),
                [&](const Vertices_with_new_position& l, const Vertices_with_new_position& r) -> bool
                {
                  return gt.less_distance_to_point_3_object()(hsp, l.second, r.second);
                });
    }

    halfedge_descriptor h_to_split = h;
    for(const Vertices_with_new_position& vnp : splitters)
    {
      // vnp.first is an iterator in the container of std::pair<vertices with same position, tolerance>
      const Vertex_container& vertices_to_move = *(vnp.first);
      const Point& p = vnp.second;

      CGAL_assertion(!vertices_to_move.empty());

#ifdef CGAL_PMP_SNAP_DEBUG_PP
      std::cout << "Splitting " << edge(h_to_split, pmt) << " |||| "
                << " Vs " << source(h_to_split, pmt) << " (" << pmt.point(source(h_to_split, pmt)) << ")"
                << " --- Vt " << target(h_to_split, pmt) << " (" << pmt.point(target(h_to_split, pmt)) << ")" << std::endl;
      std::cout << "  with pos " << p << std::endl;
      std::cout << vertices_to_move.size() << " vertices to move" << std::endl;
#endif

      halfedge_descriptor res = CGAL::Euler::split_edge(h_to_split, pmt);
      ++snapped_n;

      h_to_split = next(res, pmt); // inserting new points ordered from the source to the target of (the initial) 'h'

      // Update positions of the the source mesh and set the position of the new point in the target mesh
      put(vpmt, target(res, pmt), p);
      for(const vertex_descriptor v : vertices_to_move)
        put(vpms, v, p);

      // Now need to triangulate the quad created by the edge split
      const halfedge_descriptor face_split_h1 = opposite(h_to_split, pmt);
      const halfedge_descriptor face_split_h2 = prev(prev(face_split_h1, pmt), pmt);

      CGAL::Euler::split_face(face_split_h1, face_split_h2, pmt);
    }
  }

  return snapped_n;
}

template <typename HalfedgeRange, typename TriangleMesh, typename ToleranceMap>
std::size_t snap_vertex_range_onto_vertex_range_non_conforming(const HalfedgeRange& source_hrange,
                                                               TriangleMesh& pms,
                                                               const HalfedgeRange& target_hrange,
                                                               TriangleMesh& pmt,
                                                               const ToleranceMap& tol_pmap)
{
  return snap_vertex_range_onto_vertex_range_non_conforming(source_hrange, pms, target_hrange, pmt,
                                                            tol_pmap, parameters::all_default(),
                                                            parameters::all_default());
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

// Attempt to snap the vertices of the border of 'pm1' onto border edges of 'pm2'.

template <typename TriangleMesh>
std::size_t snap_border_vertices_non_conforming(TriangleMesh& pm1,
                                                TriangleMesh& pm2)
{
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor        halfedge_descriptor;

  typedef typename GetGeomTraits<TriangleMesh>::type                             GT;
  typedef typename GT::FT                                                        FT;

  typedef CGAL::dynamic_vertex_property_t<FT>                                    Vertex_property_tag;
  typedef typename boost::property_map<TriangleMesh, Vertex_property_tag>::type  Tolerance_map;

  std::vector<halfedge_descriptor> border_vertices1;
  border_halfedges(pm1, std::back_inserter(border_vertices1));

  std::vector<halfedge_descriptor> border_vertices2;
  if(&pm1 == &pm2)
    border_vertices2 = border_vertices1;
  else
    border_halfedges(pm2, std::back_inserter(border_vertices2));

  Tolerance_map tol_pmap = get(Vertex_property_tag(), pm1);
  const FT tol_mx(std::numeric_limits<double>::max());
  internal::assign_tolerance_with_local_edge_length_bound(border_vertices1, tol_pmap, tol_mx, pm1);

  return snap_vertex_range_onto_vertex_range_non_conforming(border_vertices1, pm1,
                                                            border_vertices2, pm2, tol_pmap,
                                                            CGAL::parameters::all_default(),
                                                            CGAL::parameters::all_default());
}

template <typename TriangleMesh, typename ToleranceMap, typename NamedParameters>
std::size_t snap_border_vertices_non_conforming(TriangleMesh& pm1,
                                                TriangleMesh& pm2,
                                                const ToleranceMap& tol_pmap,
                                                const NamedParameters& np1,
                                                const NamedParameters& np2)
{
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor        halfedge_descriptor;

  std::vector<halfedge_descriptor> border_vertices1;
  border_halfedges(pm1, std::back_inserter(border_vertices1));

  std::vector<halfedge_descriptor> border_vertices2;
  border_halfedges(pm2, std::back_inserter(border_vertices2));

  return snap_vertex_range_onto_vertex_range_non_conforming(border_vertices1, pm1,
                                                            border_vertices2, pm2,
                                                            tol_pmap, np1, np2);
}

template <typename TriangleMesh, typename ToleranceMap>
std::size_t snap_border_vertices_non_conforming(TriangleMesh& pm1,
                                                TriangleMesh& pm2,
                                                const ToleranceMap& tol_pmap)
{
  return snap_border_vertices_non_conforming(pm1, pm2, tol_pmap,
                                             CGAL::parameters::all_default(), CGAL::parameters::all_default());
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

// Same as above, but with a single mesh

template <typename TriangleMesh>
std::size_t snap_border_vertices_non_conforming(TriangleMesh& pm)
{
  return snap_border_vertices_non_conforming(pm, pm);
}

template <typename TriangleMesh, typename ToleranceMap, typename NamedParameters>
std::size_t snap_border_vertices_non_conforming(TriangleMesh& pm,
                                                const ToleranceMap& tol_pmap,
                                                const NamedParameters& np)
{
  return snap_border_vertices_non_conforming(pm, pm, tol_pmap, np, np);
}

template <typename TriangleMesh, typename ToleranceMap>
std::size_t snap_border_vertices_non_conforming(TriangleMesh& pm,
                                                const ToleranceMap& tol_pmap)
{
  return snap_border_vertices_non_conforming(pm, pm, tol_pmap,
                                             CGAL::parameters::all_default(), CGAL::parameters::all_default());
}

} // end namespace experimental

} // end namespace Polygon_mesh_processing

} // end namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_INTERNAL_SNAP_H
