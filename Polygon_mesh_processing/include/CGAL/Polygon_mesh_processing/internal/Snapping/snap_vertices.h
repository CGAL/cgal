// Copyright (c) 2018, 2019 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_POLYGON_MESH_PROCESSING_SNAPPING_SNAP_VERTICES_H
#define CGAL_POLYGON_MESH_PROCESSING_SNAPPING_SNAP_VERTICES_H

#include <CGAL/license/Polygon_mesh_processing/repair.h>

#ifdef CGAL_PMP_SNAP_DEBUG_PP
 #ifndef CGAL_PMP_SNAP_DEBUG
  #define CGAL_PMP_SNAP_DEBUG
 #endif
#endif

#include <CGAL/Polygon_mesh_processing/internal/Snapping/helper.h>
#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>
#include <CGAL/Polygon_mesh_processing/border.h>

#include <CGAL/assertions.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/selection.h>
#include <CGAL/box_intersection_d.h>
#include <CGAL/circulator.h>
#include <CGAL/Dynamic_property_map.h>
#include <CGAL/iterator.h>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/number_utils.h>
#include <CGAL/utility.h>
#include <CGAL/Real_timer.h>
#include <CGAL/tags.h>

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/member.hpp>

#ifdef CGAL_LINKED_WITH_TBB
# include <tbb/blocked_range.h>
# include <tbb/concurrent_vector.h>
# include <tbb/parallel_for.h>
#endif

#include <functional>
#ifdef CGAL_PMP_SNAP_DEBUG_OUTPUT
#include <fstream>
#endif
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <type_traits>
#include <unordered_set>
#include <utility>
#include <vector>

namespace CGAL {
namespace Polygon_mesh_processing {
namespace internal {

// Data entry of the multi-index container
template <typename PolygonMesh, typename GeomTraits>
struct Snapping_pair
{
  typedef typename GeomTraits::FT                                                     FT;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor              halfedge_descriptor;
  typedef std::vector<halfedge_descriptor>                                            Vertex_container;
  typedef std::pair<Vertex_container, FT>                                             Unique_vertex;
  typedef const Unique_vertex*                                                        Unique_vertex_ptr;

  Snapping_pair(Unique_vertex_ptr uv_a_, Unique_vertex_ptr uv_b_, const FT sq_dist_)
    : uv_a(uv_a_), uv_b(uv_b_), sq_dist(sq_dist_)
  { }

  Unique_vertex_ptr uv_a;
  Unique_vertex_ptr uv_b;
  FT sq_dist;
};

// Functor that just forwards the pair of the two intersecting boxes
template <class OutputIterator>
struct Intersecting_boxes_pairs_report
{
  Intersecting_boxes_pairs_report(OutputIterator it) :  m_iterator(it) { }

  template <class Box>
  void operator()(const Box* b, const Box* c) const {
    *m_iterator++ = std::make_pair(b->info(), c->info());
  }

  mutable OutputIterator m_iterator;
};

template <typename SnappingPairContainer,
          typename PolygonMesh, typename GeomTraits,
          typename VPM_A, typename VPM_B,
          typename ToleranceMap_A, typename ToleranceMap_B,
          typename VertexPatchMap_A, typename VertexPatchMap_B,
          typename Box
#ifdef CGAL_LINKED_WITH_TBB
          , typename UniqueVertexPairContainer = void
          , typename ToKeepContainer = void
#endif
          >
struct Vertex_proximity_report
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor                vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor              halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::edge_descriptor                  edge_descriptor;

  typedef typename GeomTraits::FT                                                     FT;

  typedef std::vector<halfedge_descriptor>                                            Vertex_container;
  typedef std::pair<Vertex_container, FT>                                             Unique_vertex;
  typedef const Unique_vertex*                                                        Unique_vertex_ptr;
  typedef std::pair<Unique_vertex_ptr, Unique_vertex_ptr>                             Unique_vertex_pair;

  Vertex_proximity_report(SnappingPairContainer& snapping_pairs,
                          const PolygonMesh& tm_A,
                          ToleranceMap_A tolerance_map_A,
                          VertexPatchMap_A vertex_patch_map_A,
                          const VPM_A& vpm_A,
                          const PolygonMesh& tm_B,
                          ToleranceMap_B tolerance_map_B,
                          VertexPatchMap_B vertex_patch_map_B,
                          const VPM_B& vpm_B,
                          const GeomTraits& gt
#ifdef CGAL_LINKED_WITH_TBB
                          , const UniqueVertexPairContainer* uv_pairs = nullptr
                          , ToKeepContainer* to_keep = nullptr
#endif
                          )
    :
      m_tm_A(tm_A),
      m_tolerance_map_A(tolerance_map_A),
      m_vertex_patch_map_A(vertex_patch_map_A),
      m_vpm_A(vpm_A),
      m_tm_B(tm_B),
      m_tolerance_map_B(tolerance_map_B),
      m_vertex_patch_map_B(vertex_patch_map_B),
      m_vpm_B(vpm_B),
      m_gt(gt),
      m_snapping_pairs(snapping_pairs)
#ifdef CGAL_LINKED_WITH_TBB
    , m_uv_pairs(uv_pairs)
    , m_to_keep(to_keep)
#endif
  { }

  bool are_equal_vertices(Unique_vertex_ptr uv_a, Unique_vertex_ptr uv_b) const
  {
    if(&m_tm_A != &m_tm_B) // must be the same mesh
      return false;

    if(uv_a->first.size() != uv_b->first.size())
      return false;

    CGAL_assertion(!uv_a->first.empty() && !uv_b->first.empty());
    const halfedge_descriptor ha = uv_a->first.front();
    return (std::find(uv_b->first.begin(), uv_b->first.end(), ha) != uv_b->first.end());
  }

  std::pair<FT, bool> do_keep(Unique_vertex_ptr uv_a, Unique_vertex_ptr uv_b) const
  {
    const Vertex_container& vs_a = uv_a->first;
    const Vertex_container& vs_b = uv_b->first;
    const vertex_descriptor va = target(vs_a.front(), m_tm_A);
    const vertex_descriptor vb = target(vs_b.front(), m_tm_B);

    // Check for patch compatibility
    if(get(m_vertex_patch_map_A, va) != get(m_vertex_patch_map_B, vb))
      return std::make_pair(-1, false); // ignore

    // Reject if same mesh and already grouped together
    if(are_equal_vertices(uv_a, uv_b))
      return std::make_pair(-1, false);

    const FT tol_a = uv_a->second;
    const FT tol_b = uv_b->second;
    CGAL_assertion(tol_a >= FT(0));
    CGAL_assertion(tol_b >= FT(0));

    const FT upper_bound_squared = CGAL::square(0.5*(tol_a + tol_b));
    CGAL::Comparison_result res =
      m_gt.compare_squared_distance_3_object()(get(m_vpm_A, va), get(m_vpm_B, vb), upper_bound_squared);

#ifdef CGAL_PMP_SNAP_DEBUG_PP
    std::cout << "squared distance between "
              << va << " [" << get(m_vpm_A, va) << "] and "
              << vb << " [" << get(m_vpm_B, vb) << "]: " << sq_dist
              << " (tol a/b: " << tol_a << " " << tol_b << ") larger? " << (res == CGAL::LARGER)
              << std::endl;
#endif

    if(res == CGAL::LARGER)
      return std::make_pair(-1, false); // ignore
    else
    {
      const FT sq_dist = (min)(m_gt.compute_squared_distance_3_object()(get(m_vpm_A, va), get(m_vpm_B, vb)),
                               upper_bound_squared);
      return std::make_pair(sq_dist, true); // keep
    }
  }

  // that's the sequential version
  void operator()(const Box* a, const Box* b)
  {
    Unique_vertex_ptr uv_a = a->info(), uv_b = b->info();
    const std::pair<FT, bool> res = do_keep(uv_a, uv_b);
    if(res.second)
      m_snapping_pairs.insert(Snapping_pair<PolygonMesh, GeomTraits>(uv_a, uv_b, res.first));
  }

#ifdef CGAL_LINKED_WITH_TBB
  // that's the parallel version
  void operator()(const tbb::blocked_range<std::size_t>& r) const
  {
    CGAL_assertion(m_uv_pairs != nullptr);
    CGAL_assertion(m_to_keep != nullptr);

    for(std::size_t ri = r.begin(); ri != r.end(); ++ri)
    {
      const Unique_vertex_pair& hp = (*m_uv_pairs)[ri];
      m_to_keep->operator[](ri) = do_keep(hp.first, hp.second);
    }
  }
#endif

private:
  const PolygonMesh& m_tm_A;
  ToleranceMap_A m_tolerance_map_A;
  VertexPatchMap_A m_vertex_patch_map_A;
  VPM_A m_vpm_A;
  const PolygonMesh& m_tm_B;
  ToleranceMap_B m_tolerance_map_B;
  VertexPatchMap_B m_vertex_patch_map_B;
  VPM_B m_vpm_B;
  const GeomTraits& m_gt;

  SnappingPairContainer& m_snapping_pairs; // will only be filled here in the sequential setting

#ifdef CGAL_LINKED_WITH_TBB
  const UniqueVertexPairContainer* m_uv_pairs;
  ToKeepContainer* m_to_keep;
#endif
};

template <typename ConcurrencyTag = CGAL::Sequential_tag,
          typename HalfedgeRange_A, typename HalfedgeRange_B, typename PolygonMesh,
          typename ToleranceMap_A, typename ToleranceMap_B,
          typename VertexPatchMap_A, typename VertexPatchMap_B,
          typename LockableVerticesOutputIterator, typename LockableHalfedgesOutputIterator,
          typename NamedParameters_A, typename NamedParameters_B>
std::size_t snap_vertices_two_way(const HalfedgeRange_A& halfedge_range_A,
                                  PolygonMesh& tm_A,
                                  ToleranceMap_A tolerance_map_A,
                                  VertexPatchMap_A vertex_patch_map_A,
                                  const HalfedgeRange_B& halfedge_range_B,
                                  PolygonMesh& tm_B,
                                  ToleranceMap_B tolerance_map_B,
                                  VertexPatchMap_B vertex_patch_map_B,
                                  LockableVerticesOutputIterator lockable_vps_out,
                                  LockableHalfedgesOutputIterator lockable_ha_out,
                                  LockableHalfedgesOutputIterator lockable_hb_out,
                                  const bool is_second_mesh_fixed,
                                  const NamedParameters_A& np_A,
                                  const NamedParameters_B& np_B)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor                vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor              halfedge_descriptor;

  typedef typename GetVertexPointMap<PolygonMesh, NamedParameters_A>::type            VPM_A;
  typedef typename GetVertexPointMap<PolygonMesh, NamedParameters_B>::type            VPM_B;

  typedef typename GetGeomTraits<PolygonMesh, NamedParameters_A>::type                GT;
  typedef typename GT::FT                                                             FT;
  typedef typename boost::property_traits<VPM_B>::value_type                          Point;

  typedef std::vector<halfedge_descriptor>                                            Vertex_container;
  typedef std::pair<Vertex_container, FT>                                             Unique_vertex;
  typedef const Unique_vertex*                                                        Unique_vertex_ptr;
  typedef std::map<std::pair<Point, std::size_t>, Unique_vertex>                      Unique_positions;

  typedef Box_intersection_d::ID_FROM_BOX_ADDRESS                                             Box_policy;
  typedef CGAL::Box_intersection_d::Box_with_info_d<double, 3, Unique_vertex_ptr, Box_policy> Box;

  using parameters::get_parameter;
  using parameters::choose_parameter;

  CGAL_static_assertion((std::is_same<Point, typename GT::Point_3>::value));

  GT gt = choose_parameter<GT>(get_parameter(np_A, internal_np::geom_traits));
  VPM_A vpm_A = choose_parameter(get_parameter(np_A, internal_np::vertex_point),
                                 get_property_map(vertex_point, tm_A));
  VPM_B vpm_B = choose_parameter(get_parameter(np_B, internal_np::vertex_point),
                                 get_property_map(vertex_point, tm_B));

#ifdef CGAL_PMP_SNAP_DEBUG
  std::cout << "Finding snappables vertices. Range sizes: "
            << std::distance(halfedge_range_A.begin(), halfedge_range_A.end()) << " and "
            << std::distance(halfedge_range_B.begin(), halfedge_range_B.end()) << std::endl;
#endif

  if(is_empty_range(halfedge_range_A.begin(), halfedge_range_A.end()) ||
     is_empty_range(halfedge_range_B.begin(), halfedge_range_B.end()))
    return 0;

  // Vertex-Vertex snapping is performed as follows:
  // - Identify points which are already equal and group them together so that they are moved together
  // - Create a single box for these points

  std::vector<Box> boxes_A;
  boxes_A.reserve(halfedge_range_A.size());
  std::vector<Box> boxes_B;
  boxes_B.reserve(halfedge_range_B.size());

  Unique_positions unique_positions_A;
  for(halfedge_descriptor h : halfedge_range_A)
  {
    const vertex_descriptor v = target(h, tm_A);
    const FT tolerance = get(tolerance_map_A, v);

#ifdef CGAL_PMP_SNAP_DEBUG_PP
    std::cout << "Pos (A): " << v << " (" << get(vpm_A, v) << ")" << std::endl;
#endif

    Vertex_container nvc {{ h }};
    std::pair<typename Unique_positions::iterator, bool> is_insert_successful =
      unique_positions_A.insert(std::make_pair(
                                  std::make_pair(get(vpm_A, v), get(vertex_patch_map_A, v)), // point and patch id
                                  std::make_pair(nvc, tolerance)));

    if(!is_insert_successful.second) // point was already met
    {
      Unique_vertex& uv = is_insert_successful.first->second;
      Vertex_container& vc = uv.first; // second is the tolerance
      CGAL_assertion(std::find(vc.begin(), vc.end(), h) == vc.end());
      vc.push_back(h);
      uv.second = (std::min)(uv.second, tolerance);
    }
  }

  // same for tm_B (@todo avoid all that for self snapping + use self_intersection_d)
  Unique_positions unique_positions_B;
  for(halfedge_descriptor h : halfedge_range_B)
  {
    const vertex_descriptor v = target(h, tm_B);
    const FT tolerance = get(tolerance_map_B, v);

#ifdef CGAL_PMP_SNAP_DEBUG_PP
    std::cout << "Pos (B): " << v << " (" << get(vpm_B, v) << ")" << std::endl;
#endif

    Vertex_container nvc {{ h }};
    std::pair<typename Unique_positions::iterator, bool> is_insert_successful =
      unique_positions_B.insert(std::make_pair(
                                  std::make_pair(get(vpm_B, v), get(vertex_patch_map_B, v)), // point and patch id
                                  std::make_pair(nvc, tolerance)));

    if(!is_insert_successful.second) // point was already met
    {
      Unique_vertex& uv = is_insert_successful.first->second;
      Vertex_container& vc = uv.first; // second is the tolerance
      CGAL_assertion(std::find(vc.begin(), vc.end(), h) == vc.end());
      vc.push_back(h);
      uv.second = (std::min)(uv.second, tolerance);
    }
  }

  // Actually build the boxes now
  for(const auto& p : unique_positions_A)
  {
#ifdef CGAL_PMP_SNAP_DEBUG_PP
    std::cout << "Unique_vertex (A), pos: " << p.first.first << " vertices:";
    for(const halfedge_descriptor h : p.second.first)
      std::cout << " " << target(h, tm_A);
    std::cout << std::endl;
#endif

    const Unique_vertex& ev = p.second;
    CGAL_assertion(!ev.first.empty());

    // this only makes the box a little larger to ease intersection computations,
    // the final tolerance is not changed
    const double eps = 1.01 * CGAL::to_double(ev.second);
    const Bbox_3 pb = gt.construct_bbox_3_object()(p.first.first);
    const Bbox_3 b(pb.xmin() - eps, pb.ymin() - eps, pb.zmin() - eps,
                   pb.xmax() + eps, pb.ymax() + eps, pb.zmax() + eps);
    boxes_A.push_back(Box(b, &ev));
  }

  for(const auto& p : unique_positions_B)
  {
#ifdef CGAL_PMP_SNAP_DEBUG_PP
    std::cout << "Unique_vertex (B), pos: " << p.first.first << " vertices:";
    for(const halfedge_descriptor h : p.second.first)
      std::cout << " " << target(h, tm_B);
    std::cout << std::endl;
#endif

    const Unique_vertex& ev = p.second;
    CGAL_assertion(!ev.first.empty());

    const double eps = 1.01 * CGAL::to_double(ev.second);
    const Bbox_3 pb = gt.construct_bbox_3_object()(p.first.first);
    const Bbox_3 b(pb.xmin() - eps, pb.ymin() - eps, pb.zmin() - eps,
                   pb.xmax() + eps, pb.ymax() + eps, pb.zmax() + eps);
    boxes_B.push_back(Box(b, &ev));
  }

  // @fixme bench and don't use ptrs if not useful
  std::vector<const Box*> boxes_A_ptr;
  boxes_A_ptr.reserve(boxes_A.size());
  for(const Box& b : boxes_A)
    boxes_A_ptr.push_back(&b);

  std::vector<const Box*> boxes_B_ptr;
  boxes_B_ptr.reserve(boxes_B.size());
  for(const Box& b : boxes_B)
    boxes_B_ptr.push_back(&b);

  // Use a multi index to sort easily by sources, targets, AND distance.
  // Then, look up the distances in increasing order, and snap whenever the source and the target
  // have both not been snapped yet.
  typedef internal::Snapping_pair<PolygonMesh, GT>                                    Snapping_pair;
  typedef boost::multi_index::multi_index_container<
    Snapping_pair,
    boost::multi_index::indexed_by<
      boost::multi_index::ordered_non_unique<
        BOOST_MULTI_INDEX_MEMBER(Snapping_pair, Unique_vertex_ptr, uv_a)>,
      boost::multi_index::ordered_non_unique<
        BOOST_MULTI_INDEX_MEMBER(Snapping_pair, Unique_vertex_ptr, uv_b)>,
      boost::multi_index::ordered_non_unique<
        BOOST_MULTI_INDEX_MEMBER(Snapping_pair, FT, sq_dist)>
    >
  >                                                                                   Snapping_pair_container;

  Snapping_pair_container snapping_pairs;

#if !defined(CGAL_LINKED_WITH_TBB)
  CGAL_static_assertion_msg (!(std::is_convertible<ConcurrencyTag, Parallel_tag>::value),
                             "Parallel_tag is enabled but TBB is unavailable.");
#else
  if(std::is_convertible<ConcurrencyTag, Parallel_tag>::value)
  {
    typedef std::pair<Unique_vertex_ptr, Unique_vertex_ptr>                           Unique_vertex_pair;
    typedef tbb::concurrent_vector<Unique_vertex_pair>                                Unique_vertex_pairs;
    typedef std::back_insert_iterator<Unique_vertex_pairs>                            UVP_output_iterator;

    Unique_vertex_pairs uv_pairs;
    Intersecting_boxes_pairs_report<UVP_output_iterator> callback(std::back_inserter(uv_pairs));

    CGAL::Real_timer timer;
    timer.start();

    // Grab the boxes that are interesecting but don't do any extra filtering (in parallel)
    CGAL::box_intersection_d<CGAL::Parallel_tag>(boxes_A_ptr.begin(), boxes_A_ptr.end(),
                                                 boxes_B_ptr.begin(), boxes_B_ptr.end(),
                                                 callback);

    std::cout << "time for box_d: " << timer.time() << std::endl;

    // Actually filter the range of intersecting boxes now (in parallel)
    typedef std::vector<std::pair<FT, bool> >                                         Filters;
    typedef Vertex_proximity_report<Snapping_pair_container,
                                    PolygonMesh, GT, VPM_A, VPM_B,
                                    ToleranceMap_A, ToleranceMap_B,
                                    VertexPatchMap_A, VertexPatchMap_B,
                                    Box, Unique_vertex_pairs, Filters>                Reporter;

    Filters to_keep(uv_pairs.size());
    Reporter proximity_filterer(snapping_pairs, tm_A, tolerance_map_A, vertex_patch_map_A, vpm_A,
                                                tm_B, tolerance_map_B, vertex_patch_map_B, vpm_B,
                                gt, &uv_pairs, &to_keep);
    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, uv_pairs.size()), proximity_filterer);

    // Now fill the multi index, sequentially
    for(std::size_t i=0, uvps = uv_pairs.size(); i<uvps; ++i)
      if(to_keep[i].second)
        snapping_pairs.insert(Snapping_pair(uv_pairs[i].first, uv_pairs[i].second, to_keep[i].first));
  }
  else
#endif
  {
    typedef Vertex_proximity_report<Snapping_pair_container, PolygonMesh, GT,
                                    VPM_A, VPM_B, ToleranceMap_A, ToleranceMap_B,
                                    VertexPatchMap_A, VertexPatchMap_B, Box>          Reporter;

    Reporter vpr(snapping_pairs, tm_A, tolerance_map_A, vertex_patch_map_A, vpm_A,
                                 tm_B, tolerance_map_B, vertex_patch_map_B, vpm_B, gt);

    // Shenanigans to pass a reference as callback (which is copied by value by 'box_intersection_d')
    std::function<void(const Box*, const Box*)> callback(std::ref(vpr));

    CGAL::box_intersection_d<CGAL::Sequential_tag>(boxes_A_ptr.begin(), boxes_A_ptr.end(),
                                                   boxes_B_ptr.begin(), boxes_B_ptr.end(),
                                                   callback);
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  /// Done collecting; start matching
  //////////////////////////////////////////////////////////////////////////////////////////////////

  if(snapping_pairs.empty())
    return 0;

  typedef std::pair<Unique_vertex_ptr, Unique_vertex_ptr>                             Unique_vertex_pair;
  std::vector<Unique_vertex_pair> snappable_vertices_pairs;

  // Sorted views of the container
  typedef typename Snapping_pair_container::template nth_index<0>::type               Container_by_source;
  typedef typename Snapping_pair_container::template nth_index<1>::type               Container_by_target;
  typedef typename Snapping_pair_container::template nth_index<2>::type               Container_by_distance;

  Container_by_source& container_by_source = snapping_pairs.template get<0>();
  Container_by_target& container_by_target = snapping_pairs.template get<1>();
  Container_by_distance& container_by_dist = snapping_pairs.template get<2>();

  // @todo this could be a best-match algorithm rather than a greedy algorithm,
  // but is the increased complexity -and thus runtime- worth it?
  CGAL_assertion_code(FT prev = -1;)
  while(!container_by_dist.empty())
  {
    const Snapping_pair& sp = *(container_by_dist.begin());
    Unique_vertex_ptr uv_a = sp.uv_a;
    Unique_vertex_ptr uv_b = sp.uv_b;
    CGAL_assertion(sp.sq_dist >= prev);
    CGAL_assertion_code(prev = sp.sq_dist;)

#ifdef CGAL_PMP_SNAP_DEBUG_PP
    const Point& pa = get(vpm_A, target(uv_a->first.front(), tm_A));
    const Point& pb = get(vpm_B, target(uv_b->first.front(), tm_B));
    std::cout << "Snapping (" << pa << ") to (" << pb << ") at dist: " << sp.sq_dist << std::endl;
    std::cout << "#verts A: " << uv_a->first.size() << " #verts B: " << uv_b->first.size() << std::endl;
#endif

    snappable_vertices_pairs.emplace_back(uv_a, uv_b);

    // 'va' and 'vb' cannot be used anymore, remove them from the container
    container_by_source.erase(uv_a);
    container_by_target.erase(uv_b);

    if((&tm_A == &tm_B))
    {
      container_by_source.erase(uv_b);
      container_by_target.erase(uv_a);
    }
  }

#ifdef CGAL_PMP_SNAP_DEBUG
  std::cout << snappable_vertices_pairs.size() << " snappable pair(s)" << std::endl;
#endif

  //////////////////////////////////////////////////////////////////////////////////////////////////
  /// Done matching; start snapping
  //////////////////////////////////////////////////////////////////////////////////////////////////

  std::size_t counter = 0;

#ifdef CGAL_PMP_SNAP_DEBUG_OUTPUT
  std::ofstream out_edges("results/snappable.polylines.txt");
  out_edges.precision(17);
#endif

  for(const Unique_vertex_pair& uvp : snappable_vertices_pairs)
  {
    Unique_vertex_ptr uv_a = uvp.first;
    Unique_vertex_ptr uv_b = uvp.second;
    const Vertex_container& vs_a = uv_a->first;
    const Vertex_container& vs_b = uv_b->first;
    const vertex_descriptor va = target(vs_a.front(), tm_A);
    const vertex_descriptor vb = target(vs_b.front(), tm_B);

#ifdef CGAL_PMP_SNAP_DEBUG_OUTPUT
    out_edges << "2 " << tm_A.point(va) << " " << tm_B.point(vb) << std::endl;
#endif

    if(!gt.equal_3_object()(get(vpm_A, va), get(vpm_B, vb)))
    {
      if(is_second_mesh_fixed)
      {
        for(const halfedge_descriptor ha : vs_a)
          put(vpm_A, target(ha, tm_A), get(vpm_B, vb));
      }
      else
      {
        // Pick a point that is on the segment [va; vb], with a ratio based on the respective tolerances
        const FT tol_s = uv_a->second;
        const FT tol_t = uv_b->second;
        CGAL_assertion(tol_s != FT(0) || tol_t != FT(0));

        const FT lambda = tol_t / (tol_s + tol_t);
        const Point new_p = get(vpm_A, va) + lambda * (get(vpm_B, vb) - get(vpm_A, va));
#ifdef CGAL_PMP_SNAP_DEBUG_PP
        std::cout << "new position of " << va << " " << vb << " --> " << new_p << std::endl;
#endif

        for(const halfedge_descriptor ha : vs_a)
          put(vpm_A, target(ha, tm_A), new_p);

        for(const halfedge_descriptor hb : vs_b)
          put(vpm_B, target(hb, tm_B), new_p);
      }

      ++counter;
    }
  }

#ifdef CGAL_PMP_SNAP_DEBUG_OUTPUT
  out_edges.close();
#endif

#ifdef CGAL_PMP_SNAP_DEBUG
  std::cout << "Snapped " << counter << " pair(s)!" << std::endl;
#endif

  //////////////////////////////////////////////////////////////////////////////////////////////////
  /// Done snapping; start analyzing
  //////////////////////////////////////////////////////////////////////////////////////////////////

  // Below is used in non-conformal snapping
  // @todo could avoid doing it if not required
  //
  // Now that vertex-vertex snapping has been performed, look around to see if we can already
  // lock some vertices and halfedges...
  //
  // #1 : If a pair of edges (1 incident to va, the other to vb) are already fully matching,
  //      we don't want to project anything more onto any of those 2 edges
  // #2 : If pairs on either side of the two matching vertices are compatible (not necessary fully matching),
  //      then the two vertices should be locked as no better match can be obtained
  // #3 : If a pair of incident edges are not fully matching, but still have compatible directions
  //      (i.e. collinear and opposite directions), then we don't want to project onto the shorter
  //      of the two
  for(const Unique_vertex_pair& uvp : snappable_vertices_pairs)
  {
    Unique_vertex_ptr uv_a = uvp.first;
    Unique_vertex_ptr uv_b = uvp.second;
    const Vertex_container& vs_a = uv_a->first;
    const Vertex_container& vs_b = uv_b->first;

    // Quadratic, but all halfedges in vs_a and vs_b point to the same point. There shouldn't be many.
    // @fixme this assumes compatible orientation...
    for(const halfedge_descriptor ha : vs_a)
    {
      for(const halfedge_descriptor hb : vs_b)
      {
        if(!is_border(ha, tm_A) || !is_border(hb, tm_B))
          continue;

        const vertex_descriptor va = target(ha, tm_A);
        const vertex_descriptor vb = target(hb, tm_B);

        // The two folloing halfedges might not be in the range to snap, but it doesn't matter
        const halfedge_descriptor nha = next(ha, tm_A);
        const halfedge_descriptor nhb = next(hb, tm_B);
        const bool is_stitchable_left = gt.equal_3_object()(get(vpm_A, source(ha, tm_A)),
                                                            get(vpm_B, target(nhb, tm_B)));
        const bool is_stitchable_right = gt.equal_3_object()(get(vpm_A, target(nha, tm_A)),
                                                             get(vpm_B, source(hb, tm_B)));
        bool is_collinear_left = false;
        bool is_collinear_right = false;

        if(is_stitchable_left)
        {
          *lockable_ha_out++ = ha; // #1
          *lockable_hb_out++ = nhb; // #1
          is_collinear_left = true;
        }
        else
        {
          const vertex_descriptor lva = source(ha, tm_A);
          const vertex_descriptor lvb = target(nhb, tm_B);

          is_collinear_left = is_collinear_with_tolerance(get(vpm_A, va), get(vpm_A, lva), get(vpm_B, lvb), gt);
          if(is_collinear_left) // but not stitchable
          {
            if(gt.less_distance_to_point_3_object()(get(vpm_A, va), get(vpm_A, lva), get(vpm_B, lvb)))
              *lockable_ha_out++ = ha; // #2
            else
              *lockable_hb_out++ = nhb; // #2
          }
        }

        if(is_stitchable_right)
        {
          *lockable_ha_out++ = nha; // #1
          *lockable_hb_out++ = hb; // #1
          is_collinear_right = true;
        }
        else
        {
          const vertex_descriptor rva = target(nha, tm_A);
          const vertex_descriptor rvb = source(hb, tm_B);

          is_collinear_right = is_collinear_with_tolerance(get(vpm_A, va), get(vpm_A, rva), get(vpm_B, rvb), gt);
          if(is_collinear_right) // but not stitchable
          {
            if(gt.less_distance_to_point_3_object()(get(vpm_A, va), get(vpm_A, rva), get(vpm_B, rvb)))
              *lockable_ha_out++ = nha; // #2
            else
              *lockable_hb_out++ = hb; // #2
          }
        }

        if(is_collinear_left && is_collinear_right) // #3
        {
#ifdef CGAL_PMP_SNAP_DEBUG
          std::cout << va << " (" << get(vpm_A, va) << ") and " << vb << " (" << get(vpm_B, vb) << ") are locked vertices" << std::endl;
#endif
          *lockable_vps_out++ = std::make_pair(va, vb);
        }
      }
    }
  }

  return counter;
}

// Convenience overload for snap_borders
template <typename ConcurrencyTag = CGAL::Sequential_tag,
          typename HalfedgeRange_A, typename HalfedgeRange_B, typename PolygonMesh,
          typename ToleranceMap_A, typename ToleranceMap_B,
          typename LockableVerticesOutputIterator, typename LockableHalfedgesOutputIterator,
          typename NamedParameters_A, typename NamedParameters_B>
std::size_t snap_vertices_two_way(const HalfedgeRange_A& halfedge_range_A,
                                  PolygonMesh& tm_A,
                                  ToleranceMap_A tolerance_map_A,
                                  const HalfedgeRange_B& halfedge_range_B,
                                  PolygonMesh& tm_B,
                                  ToleranceMap_B tolerance_map_B,
                                  LockableVerticesOutputIterator lockable_vps_out,
                                  LockableHalfedgesOutputIterator lockable_ha_out,
                                  LockableHalfedgesOutputIterator lockable_hb_out,
                                  const NamedParameters_A& np_A,
                                  const NamedParameters_B& np_B)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor             vertex_descriptor;

  using CGAL::parameters::choose_parameter;
  using CGAL::parameters::get_parameter;

  return snap_vertices_two_way(halfedge_range_A, tm_A, tolerance_map_A,
                               choose_parameter(get_parameter(np_A, internal_np::vertex_incident_patches),
                                                Constant_property_map<vertex_descriptor, std::size_t>(-1)),
                               halfedge_range_B, tm_B, tolerance_map_B,
                               choose_parameter(get_parameter(np_B, internal_np::vertex_incident_patches),
                                                Constant_property_map<vertex_descriptor, std::size_t>(-1)),
                               lockable_vps_out, lockable_ha_out, lockable_hb_out,
                               choose_parameter(get_parameter(np_B, internal_np::do_lock_mesh), false),
                               np_A, np_B);
}

} // namespace internal

namespace experimental {

// \ingroup PMP_repairing_grp
//
// Attempts to snap the vertices in `halfedge_range_A` and `halfedge_range_B`.
// A vertex from the first range and a vertex from the second range are only snapped
// if the distance between both vertices is smaller than the sum of their user-prescribed tolerance.
// All such pairs are collected and processed in a greedy order: the two vertices are moved
// to a common vertex (more on this below), and both vertices are locked and cannot be used
// in any other matches.
//
// \warning This function does not give any guarantee on the conformity between the two meshes after the snapping.
// \warning This function does not merge vertices or the meshes, it is purely geometric.
//
// \tparam PolygonMesh a model of `FaceListGraph` and `MutableFaceGraph`
// \tparam HalfedgeRange_A a model of `Range` with value type `boost::graph_traits<PolygonMesh>::%halfedge_descriptor`
// \tparam HalfedgeRange_B a model of `Range` with value type `boost::graph_traits<PolygonMesh>::%halfedge_descriptor`
// \tparam ToleranceMap_A a model of `ReadablePropertyMap` with key type `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
//                        and value type `GetGeomTraits<PolygonMesh, NamedParameters_A>::type::FT`
// \tparam ToleranceMap_B a model of `ReadablePropertyMap` with key type `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
//                        and value type `GetGeomTraits<PolygonMesh, NamedParameters_A>::type::FT`
// \tparam NamedParameters_A a sequence of \ref bgl_namedparameters "Named Parameters"
// \tparam NamedParameters_B a sequence of \ref bgl_namedparameters "Named Parameters"
//
// \param halfedge_range_A a range of halfedges of the first mesh defining a set of vertices (as targets of the halfeges)
// \param tm_A the first mesh to which the vertices in `halfedge_range_A` belong
// \param tolerance_map_A a tolerance map associating to each vertex of the first range a tolerance value
// \param halfedge_range_B a range of vertices of the second mesh defining a set of vertices (as targets of the halfeges)
// \param tolerance_map_B a tolerance map associating to each vertex of the second range a tolerance value
// \param tm_B the target mesh to which the vertices in `halfedge_range_B` belong
// \param np_A an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
//
// \cgalNamedParamsBegin
//   \cgalParamBegin{vertex_point_map}
//     the property map with the points associated to the vertices of the source mesh.
//     The type of this map is model of `ReadWritePropertyMap`. If this parameter is omitted,
//     an internal property map for `CGAL::vertex_point_t` must be available in `PolygonMesh`.
//   \cgalParamEnd
//   \cgalParamBegin{geom_traits}
//     a geometric traits class instance, must be a model of `Kernel`
//   \cgalParamEnd
// \cgalNamedParamsEnd
//
// \param np_B an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
//
// \cgalNamedParamsBegin
//   \cgalParamBegin{vertex_point_map}
//     the property map with the points associated to the vertices of the target mesh.
//     The type of this map is model of `ReadablePropertyMap`. If this parameter is omitted,
//     an internal property map for `CGAL::vertex_point_t` must be available in `PolygonMesh`.
//   \cgalParamEnd
// \cgalNamedParamsEnd
//
// \return the number of snapped vertex pairs
//
template <typename HalfedgeRange_A, typename HalfedgeRange_B, typename PolygonMesh,
          typename ToleranceMap_A, typename ToleranceMap_B,
          typename NamedParameters_A, typename NamedParameters_B>
std::size_t snap_vertices(const HalfedgeRange_A& halfedge_range_A,
                          PolygonMesh& tm_A,
                          ToleranceMap_A tolerance_map_A,
                          const HalfedgeRange_B& halfedge_range_B,
                          PolygonMesh& tm_B,
                          ToleranceMap_B tolerance_map_B,
                          const NamedParameters_A& np_A,
                          const NamedParameters_B& np_B)
{
  CGAL::Emptyset_iterator unused_output_iterator;

  return internal::snap_vertices_two_way(halfedge_range_A, tm_A, tolerance_map_A,
                                         halfedge_range_B, tm_B, tolerance_map_B,
                                         unused_output_iterator, unused_output_iterator,
                                         unused_output_iterator, np_A, np_B);
}

template <typename HalfedgeRange_A, typename HalfedgeRange_B, typename PolygonMesh,
          typename ToleranceMap_A, typename ToleranceMap_B>
std::size_t snap_vertices(const HalfedgeRange_A& halfedge_range_A,
                          PolygonMesh& tm_A,
                          ToleranceMap_A tolerance_map_A,
                          const HalfedgeRange_B& halfedge_range_B,
                          PolygonMesh& tm_B,
                          ToleranceMap_B tolerance_map_B)
{
  return snap_vertices(halfedge_range_A, tm_A, tolerance_map_A, halfedge_range_B, tm_B, tolerance_map_B,
                       CGAL::parameters::all_default(), CGAL::parameters::all_default());
}

template <typename HalfedgeRange_A, typename HalfedgeRange_B, typename PolygonMesh,
          typename T_A, typename Tag_A, typename Base_A,
          typename T_B, typename Tag_B, typename Base_B>
std::size_t snap_vertices(const HalfedgeRange_A& halfedge_range_A,
                          PolygonMesh& tm_A,
                          const HalfedgeRange_B& halfedge_range_B,
                          PolygonMesh& tm_B,
                          const CGAL::Named_function_parameters<T_A, Tag_A, Base_A>& np_A,
                          const CGAL::Named_function_parameters<T_B, Tag_B, Base_B>& np_B)
{
  typedef CGAL::Named_function_parameters<T_A, Tag_A, Base_A>                         NamedParameters_A;
  typedef typename GetGeomTraits<PolygonMesh, NamedParameters_A>::type                GT;
  typedef typename GT::FT                                                             FT;
  typedef CGAL::dynamic_vertex_property_t<FT>                                         Vertex_property_tag;
  typedef typename boost::property_map<PolygonMesh, Vertex_property_tag>::type        Tolerance_map;

  const FT max_tol((std::numeric_limits<double>::max)());

  Tolerance_map tolerance_map_A = get(Vertex_property_tag(), tm_A);
  internal::assign_tolerance_with_local_edge_length_bound(halfedge_range_A, tolerance_map_A, max_tol, tm_A, np_A);
  Tolerance_map tolerance_map_B = get(Vertex_property_tag(), tm_B);
  internal::assign_tolerance_with_local_edge_length_bound(halfedge_range_B, tolerance_map_B, max_tol, tm_B, np_B);

  return snap_vertices(halfedge_range_A, tm_A, tolerance_map_A, halfedge_range_B, tm_B, tolerance_map_B, np_A, np_B);
}

template <typename HalfedgeRange_A, typename HalfedgeRange_B, typename PolygonMesh>
std::size_t snap_vertices(const HalfedgeRange_A& halfedge_range_A,
                          PolygonMesh& tm_A,
                          const HalfedgeRange_B& halfedge_range_B,
                          PolygonMesh& tm_B)
{
  return snap_vertices(halfedge_range_A, tm_A, halfedge_range_B, tm_B,
                       parameters::all_default(), parameters::all_default());
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// Border convenience overloads
///////////////////////////////////////////////////////////////////////////////////////////////////

template <typename PolygonMesh, typename ToleranceMap_A, typename ToleranceMap_B>
std::size_t snap_border_vertices(PolygonMesh& tm_A,
                                 ToleranceMap_A tolerance_map_A,
                                 PolygonMesh& tm_B,
                                 ToleranceMap_B tolerance_map_B)
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor    halfedge_descriptor;

  std::vector<halfedge_descriptor> border_A;
  border_halfedges(tm_A, std::back_inserter(border_A));
  std::vector<halfedge_descriptor> border_B;
  border_halfedges(tm_B, std::back_inserter(border_B));

  return snap_vertices(border_A, tm_A, tolerance_map_A,
                       border_B, tm_B, tolerance_map_B);
}

template <typename PolygonMesh>
std::size_t snap_border_vertices(PolygonMesh& tm_A,
                                 PolygonMesh& tm_B)
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor    halfedge_descriptor;

  std::vector<halfedge_descriptor> border_vertices_A;
  border_halfedges(tm_A, std::back_inserter(border_vertices_A));

  std::vector<halfedge_descriptor> border_vertices_B;
  border_halfedges(tm_B, std::back_inserter(border_vertices_B));

  return snap_vertices(border_vertices_A, tm_A, border_vertices_B, tm_B);
}

template <typename PolygonMesh, typename ToleranceMap>
std::size_t snap_border_vertices(PolygonMesh& tm, ToleranceMap tolerance_map)
{
  return snap_border_vertices(tm, tolerance_map, tm, tolerance_map);
}

template <typename PolygonMesh>
std::size_t snap_border_vertices(PolygonMesh& tm)
{
  return snap_border_vertices(tm, tm);
}

} // namespace experimental
} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_SNAPPING_SNAP_VERTICES_H
