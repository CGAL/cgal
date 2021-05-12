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
//                 Sebastien Loriot

#ifndef CGAL_POLYGON_MESH_PROCESSING_SNAPPING_SNAP_H
#define CGAL_POLYGON_MESH_PROCESSING_SNAPPING_SNAP_H

#include <CGAL/license/Polygon_mesh_processing/repair.h>

#ifdef CGAL_PMP_SNAP_DEBUG_PP
 #ifndef CGAL_PMP_SNAP_DEBUG
  #define CGAL_PMP_SNAP_DEBUG
 #endif
#endif

#include <CGAL/Polygon_mesh_processing/internal/Snapping/helper.h>
#include <CGAL/Polygon_mesh_processing/internal/Snapping/snap_vertices.h>

#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_halfedge_graph_segment_primitive.h>
#include <CGAL/assertions.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/selection.h>
#include <CGAL/circulator.h>
#include <CGAL/Dynamic_property_map.h>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/number_utils.h>
#include <CGAL/Real_timer.h>
#include <CGAL/utility.h>

#ifdef CGAL_EIGEN3_ENABLED
#include <Eigen/Dense>
#endif

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/concurrent_vector.h>
#include <tbb/parallel_for.h>
#endif

#include <functional>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <type_traits>
#include <unordered_set>
#include <utility>
#include <vector>

#ifdef CGAL_PMP_SNAP_DEBUG_OUTPUT
#include <fstream>
#endif

namespace CGAL {
namespace Polygon_mesh_processing {
namespace internal {

template <typename HalfedgeRange, typename TriangleMesh,
          typename ToleranceMap, typename NamedParameters>
void simplify_range(HalfedgeRange& halfedge_range,
                    TriangleMesh& tm,
                    ToleranceMap tolerance_map,
                    const NamedParameters& np)
{
  // @todo the simplification below is too naive and will fail to treat complicated zigzaging
  // with long edges.
  // A real simplification would create a fuzzy area around the border (say the dilation
  // of the border polyline by a sphere of radius that continuously interpolates the tolerance
  // at each vertex). The border is then simplified to use the least amount of edges possible
  // in that fuzzy zone. Need to look into cartography papers, probably.

  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor           vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor         halfedge_descriptor;

  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type             GT;
  typedef typename GT::FT                                                         FT;

  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::type         VPM;
  typedef typename boost::property_traits<VPM>::reference                         Point_ref;

  using parameters::get_parameter;
  using parameters::choose_parameter;

  const GT gt = choose_parameter<GT>(get_parameter(np, internal_np::geom_traits));
  VPM vpm = choose_parameter(get_parameter(np, internal_np::vertex_point), get_property_map(vertex_point, tm));

  typedef CGAL::dynamic_halfedge_property_t<bool>                                 Halfedge_bool_tag;
  typedef typename boost::property_map<TriangleMesh, Halfedge_bool_tag>::type     Range_halfedges;

  Range_halfedges range_halfedges = get(Halfedge_bool_tag(), tm);
  for(halfedge_descriptor h : halfedge_range)
    put(range_halfedges, h, true);

  CGAL_postcondition_code(const std::size_t initial_n = halfedge_range.size();)

  std::set<halfedge_descriptor> edges_to_test(halfedge_range.begin(), halfedge_range.end());

  int collapsed_n = 0;
  while(!edges_to_test.empty())
  {
    const halfedge_descriptor h = *(edges_to_test.begin());
    edges_to_test.erase(edges_to_test.begin());

    const vertex_descriptor vs = source(h, tm);
    const vertex_descriptor vt = target(h, tm);
    const Point_ref ps = get(vpm, vs);
    const Point_ref pt = get(vpm, vt);

    // @fixme what if the source vertex is not to be snapped? Tolerance cannot be obtained...
    // and where should the post-collapse vertex be since we can't move the source vertex...
    // --> simply don't collapse?
    const FT min_tol = (std::min)(get(tolerance_map, vs), get(tolerance_map, vt));
    const FT max_tol = (std::max)(get(tolerance_map, vs), get(tolerance_map, vt));

    if(gt.compare_squared_distance_3_object()(ps,pt,CGAL::square(max_tol))==SMALLER)
    {
      const halfedge_descriptor prev_h = prev(h, tm);
      const halfedge_descriptor next_h = next(h, tm);

      // check that the border has at least 4 edges not to create degenerate volumes
      if(border_size(h, tm) >= 4)
      {
        const FT h_sq_length = gt.compute_squared_distance_3_object()(ps, pt);
        vertex_descriptor v = Euler::collapse_edge(edge(h, tm), tm);

        put(vpm, v, gt.construct_midpoint_3_object()(ps, pt));
        put(tolerance_map, v, min_tol + 0.5 * CGAL::approximate_sqrt(h_sq_length));

        if(get(range_halfedges, prev_h))
          edges_to_test.insert(prev_h);
        if(get(range_halfedges, next_h))
          edges_to_test.insert(next_h);

        ++collapsed_n;
      }
    }
  }

#ifdef CGAL_PMP_SNAP_DEBUG
  std::cout << "collapsed " << collapsed_n << " edges" << std::endl;
#endif

  std::vector<halfedge_descriptor> new_range;
  new_range.reserve(halfedge_range.size());

  for(halfedge_descriptor h : halfedges(tm))
    if(get(range_halfedges, h))
      new_range.push_back(h);

  halfedge_range = HalfedgeRange(new_range.begin(), new_range.end());

  CGAL_postcondition(halfedge_range.size() == initial_n - collapsed_n);
}

// Adapted from <CGAL/internal/AABB_tree/AABB_traversal_traits.h>
template <typename TriangleMesh,
          typename VPMS, typename VPMT, typename FacePatchMap,
          typename AABBTraits>
class Projection_traits
{
  typedef typename AABBTraits::FT                                          FT;
  typedef typename AABBTraits::Point_3                                     Point;
  typedef typename AABBTraits::Primitive                                   Primitive;
  typedef typename AABBTraits::Bounding_box                                Bounding_box;
  typedef typename AABBTraits::Primitive::Id                               Primitive_id;
  typedef typename AABBTraits::Point_and_primitive_id                      Point_and_primitive_id;
  typedef typename AABBTraits::Object_and_primitive_id                     Object_and_primitive_id;

  typedef CGAL::AABB_node<AABBTraits>                                      Node;

  typedef typename AABBTraits::Geom_traits                                 Geom_traits;
  typedef typename Geom_traits::Vector_3                                   Vector;
  typedef typename Geom_traits::Plane_3                                    Plane;

  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor  halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::edge_descriptor      edge_descriptor;

#ifdef CGAL_PMP_SNAP_USE_ANISOTROPIC_DISTANCE
  void build_metric()
  {
    Vector hv(get(m_vpm_S, source(m_h, m_tm_S)), get(m_vpm_S, target(m_h, m_tm_S)));
    Vector n = CGAL::Polygon_mesh_processing::compute_face_normal(face(opposite(m_h, m_tm_S), m_tm_S), m_tm_S,
                                                                  CGAL::parameters::geom_traits(m_gt)
                                                                                   .vertex_point_map(m_vpm_S));
    Vector pn = m_gt.construct_cross_product_vector_3_object()(hv, n);
    Plane pl(get(m_vpm_S, target(m_h, m_tm_S)), pn);

    Vector b1 = pl.base1();
    Vector b2 = pl.base2();

    internal::normalize(b1, m_gt);
    internal::normalize(b2, m_gt);
    internal::normalize(pn, m_gt);

    Eigen::Matrix3d eigen_m;
    eigen_m(0,0) = b1.x(); eigen_m(0,1) = b2.x(); eigen_m(0,2) = pn.x();
    eigen_m(1,0) = b1.y(); eigen_m(1,1) = b2.y(); eigen_m(1,2) = pn.y();
    eigen_m(2,0) = b1.z(); eigen_m(2,1) = b2.z(); eigen_m(2,2) = pn.z();

    Eigen::Matrix3d eigen_diag = Eigen::Matrix3d::Zero();
    eigen_diag(0,0) = 1;
    eigen_diag(1,1) = 1;
    eigen_diag(2,2) = 100; // we scale by 1/sqrt(lambda_0) in the direction of n

    Eigen::Matrix3d eigen_mtransp = eigen_m.transpose();
    m_metric = eigen_m * eigen_diag * eigen_mtransp;
  }
#endif

  FT squared_anisotropic_distance(const Point& p, const Point& q) const
  {
#ifdef CGAL_PMP_SNAP_USE_ANISOTROPIC_DISTANCE
    Vector pq(p, q);
    Eigen::Vector3d ev;
    ev(0) = pq.x();
    ev(1) = pq.y();
    ev(2) = pq.z();

    return ev.transpose() * m_metric * ev;
#else
    return CGAL::squared_distance(p, q);
#endif
  }

  // The idea behind this function is to give priority to projections that are parallel to the edge 'e0'
  // whom the query is the target of. This is done by considering an (anisotropic) distance that
  // increases the distance if we are not projecting along the vector orthogonal to both the edge 'e0'
  // and to the normal of the face, and also by checking the scalar product between the edges.
  // The latter is because the distance to the projection point does not entirely reflect
  // the alignment of the edges.
  //
  // This is WIP and likely to evolve...
  bool is_better_than_current_best(const FT sq_dist,
                                   const FT scalar_product)
  {
#ifdef CGAL_PMP_SNAP_DEBUG_PP
    std::cout << "is_better_than_current_best()" << std::endl;
    std::cout << "sq_dist: " << sq_dist << std::endl;
    std::cout << "m_sq_adist: " << m_sq_adist << std::endl;
    std::cout << "scalar products: " << scalar_product << " " << m_sp_with_closest_edge << std::endl;
#endif

    // Automatically accept much closer tentative targets; automatically reject much farther tentative targets
    const FT lambda = 4; // comparing squared distances
    if(lambda * sq_dist < m_sq_adist)
      return true;
    if(lambda * m_sq_adist < sq_dist)
      return false;

    return (scalar_product > m_sp_with_closest_edge);
  }

public:
  explicit Projection_traits(const halfedge_descriptor h,
                             const std::size_t patch_id,
                             const FT sq_tolerance,
                             const TriangleMesh& tm_S,
                             const VPMS vpm_S,
                             const TriangleMesh& tm_T,
                             const FacePatchMap face_patch_map_T,
                             const VPMT vpm_T,
                             const AABBTraits& aabb_traits)
    :
      m_h(h),
      m_patch_id(patch_id),
      m_sq_tol(sq_tolerance),
      m_tm_S(tm_S), m_vpm_S(vpm_S),
      m_tm_T(tm_T), m_face_patch_map_T(face_patch_map_T), m_vpm_T(vpm_T),
      m_is_same_mesh((&tm_S == &tm_T)),
      m_continue(true),
      m_closest_point_initialized(false),
      m_traits(aabb_traits),
      m_gt(Geom_traits()) // blame AABB's traits management
  {
    // not fully convinced that it is still useful after adding scalar product shenanigans
// #define CGAL_PMP_SNAP_USE_ANISOTROPIC_DISTANCE
#ifdef CGAL_PMP_SNAP_USE_ANISOTROPIC_DISTANCE
    build_metric();
#endif

    Vector hv(get(vpm_S, source(h, tm_S)), get(vpm_S, target(h, tm_S)));
    m_direction = hv;
    CGAL::Polygon_mesh_processing::internal::normalize(m_direction, m_gt);
  }

  bool go_further() const { return m_continue; }

  void intersection(const Point& query, const Primitive& primitive)
  {
#ifdef CGAL_PMP_SNAP_DEBUG_PP
    std::cout << "~~~~ intersection with primitive: " << primitive.id() << std::endl;
    std::cout << get(m_vpm_T, source(primitive.id(), m_tm_T)) << std::endl;
    std::cout << get(m_vpm_T, target(primitive.id(), m_tm_T)) << std::endl;
#endif

    halfedge_descriptor h = halfedge(primitive.id(), m_tm_T);
    if(is_border(h, m_tm_T))
      h = opposite(h, m_tm_T);

#ifdef CGAL_PMP_SNAP_DEBUG_PP
      std::cout << "patches: " << get(m_face_patch_map_T, face(h, m_tm_T)) << " " << m_patch_id << std::endl;
#endif

    if(get(m_face_patch_map_T, face(h, m_tm_T)) != m_patch_id)
      return;

    const typename Primitive::Datum& s = primitive.datum(m_traits.shared_data());
    if(m_traits.equal_3_object()(s[0], query) || m_traits.equal_3_object()(s[1], query))
    {
      // If we are NOT using the same mesh and the query is (geometrically) equal to one extremity
      // of the target edge, we don't want to move the source point away from the target point
      // (because we cannot move the target point).
      if(!m_is_same_mesh)
      {
#ifdef CGAL_PMP_SNAP_DEBUG_PP
        std::cout << "This vertex is stuck because it is equal to a vertex on the target mesh" << std::endl;
#endif
        m_closest_point_initialized = false;
        m_continue = false;
      }

      // skip the primitive if one of its endpoints is the query
      return;
    }

    // We are searching for a point on the target border that is close to target(h).
    // To try and avoid foldings, we penalize potential matches that are roughly in the direction of 'h'
    // by using an anisotropic distance that is the Euclidean on the plane orthogonal to the direction,
    // and much larger when traveling in the rough direction of 'h'.
    //
    // Note that we apply this to points that fall within the tolerance ball, so if there is
    // only a single candidate that is in a direction we don't like, we still take it.
    //
    // @todo should the construct_projected_point_3 be anisotropic too?
    const Point new_closest_point = m_gt.construct_projected_point_3_object()(
      CGAL::internal::Primitive_helper<AABBTraits>::get_datum(primitive, m_traits), query);

#ifdef CGAL_PMP_SNAP_DEBUG_PP
    std::cout << "closest point to primitive: " << new_closest_point << std::endl;
#endif

    const FT sq_ad_to_tentative_closest_pt = squared_anisotropic_distance(query, new_closest_point);

    Vector tentative_dir(get(m_vpm_T, source(primitive.id(), m_tm_T)),
                         get(m_vpm_T, target(primitive.id(), m_tm_T)));
    CGAL::Polygon_mesh_processing::internal::normalize(tentative_dir, m_gt);
    const FT sp_with_tentative = CGAL::abs(tentative_dir * m_direction);

    if(!m_closest_point_initialized ||
       is_better_than_current_best(sq_ad_to_tentative_closest_pt, sp_with_tentative))
    {
#ifdef CGAL_PMP_SNAP_DEBUG_PP
      std::cout << "is better!" << std::endl;
#endif

      m_closest_point_initialized = true;

      m_closest_primitive = primitive.id();
      m_closest_point = new_closest_point;
      m_sq_adist = sq_ad_to_tentative_closest_pt;
      m_sp_with_closest_edge = sp_with_tentative;
    }
  }

  bool do_intersect(const Point& query, const Node& node) const
  {
    // This should NOT be the anisotropic distance, as we want to test all targets within the tolerance
    return m_traits.compare_distance_object()(query, node.bbox(), m_sq_tol) == CGAL::SMALLER;
  }

  const Point& closest_point() const { return m_closest_point; }
  typename Primitive::Id closest_primitive_id() const { return m_closest_primitive; }
  bool closest_point_initialized() const { return m_closest_point_initialized; }
  FT scalar_product_with_best() const { return m_sp_with_closest_edge; }

private:
  halfedge_descriptor m_h;
  const std::size_t m_patch_id;
  const FT m_sq_tol;
  Vector m_direction;
#ifdef CGAL_PMP_SNAP_USE_ANISOTROPIC_DISTANCE
  Eigen::Matrix3d m_metric;
#endif

  const TriangleMesh& m_tm_S;
  VPMS m_vpm_S;
  const TriangleMesh& m_tm_T;
  FacePatchMap m_face_patch_map_T;
  VPMT m_vpm_T;
  const bool m_is_same_mesh;

  bool m_continue;
  bool m_closest_point_initialized;
  typename Primitive::Id m_closest_primitive;
  Point m_closest_point;
  FT m_sq_adist;
  FT m_sp_with_closest_edge; // scalar product between m_h and the current best edge candidate

  const AABBTraits& m_traits;
  Geom_traits m_gt;
};

template <typename ConcurrencyTag>
struct Edges_to_split_map_inserter // Parallel
{
#ifdef CGAL_LINKED_WITH_TBB
  template <typename EdgeToSplitMap, typename HalfedgeDescriptor, typename Point>
  void operator()(EdgeToSplitMap& edges_to_split,
                  const HalfedgeDescriptor closest_h,
                  const HalfedgeDescriptor h,
                  const Point& closest_p)
  {
    typename EdgeToSplitMap::accessor acc;
    edges_to_split.insert(acc, closest_h);
    acc->second.emplace_back(h, closest_p);
  }
#endif
};

template <>
struct Edges_to_split_map_inserter<CGAL::Sequential_tag>
{
  template <typename EdgeToSplitMap, typename HalfedgeDescriptor, typename Point>
  void operator()(EdgeToSplitMap& edges_to_split,
                  const HalfedgeDescriptor closest_h,
                  const HalfedgeDescriptor h,
                  const Point& closest_p)
  {
    edges_to_split[closest_h].emplace_back(h, closest_p);
  }
};

// The UniqueVertex is a pair of a container of vertex_descriptor and FT, representing
// vertices with the same geometric position and their associated snapping tolerance
// (defined as the minimum of the tolerances of the vertices of the bunch)
template <typename ConcurrencyTag = CGAL::Sequential_tag,
          typename VertexWithTolerance, typename TriangleMesh,
          typename EdgeToSplitMap, typename AABBTree,
          typename VertexPatchMap_S, typename FacePatchMap_T,
          typename VPMS, typename VPMT, typename GT>
void find_splittable_edge(const VertexWithTolerance& vertex_with_tolerance,
                          EdgeToSplitMap& edges_to_split,
                          const AABBTree* aabb_tree_ptr,
                          const TriangleMesh& tm_S,
                          VertexPatchMap_S vertex_patch_map_S,
                          VPMS vpm_S,
                          const TriangleMesh& tm_T,
                          FacePatchMap_T face_patch_map_T,
                          VPMT vpm_T,
                          const GT& gt)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor           vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor         halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::edge_descriptor             edge_descriptor;

  typedef typename GT::FT                                                         FT;
  typedef typename boost::property_traits<VPMS>::value_type                       Point;
  typedef typename boost::property_traits<VPMS>::reference                        Point_ref;

  typedef typename AABBTree::AABB_traits                                          AABB_traits;

  typedef internal::Projection_traits<TriangleMesh, VPMS, VPMT, FacePatchMap_T, AABB_traits> Projection_traits;

  // by construction the whole range has the same position
  const halfedge_descriptor h = vertex_with_tolerance.first;
  const vertex_descriptor v = target(h, tm_S);
  const Point_ref query = get(vpm_S, v);
  const FT sq_tolerance = CGAL::square(vertex_with_tolerance.second);
  const std::size_t patch_id = get(vertex_patch_map_S, v);

#ifdef CGAL_PMP_SNAP_DEBUG_PP
  std::cout << "--------------------------- Query: " << v << " (" << query << ")" << std::endl;
#endif

  Projection_traits traversal_traits(h, patch_id, sq_tolerance, tm_S, vpm_S,
                                     tm_T, face_patch_map_T, vpm_T, aabb_tree_ptr->traits());
  aabb_tree_ptr->traversal(query, traversal_traits);

  if(!traversal_traits.closest_point_initialized())
  {
#ifdef CGAL_PMP_SNAP_DEBUG_PP
  std::cout << "Couldn't find any single projection point" << std::endl;
#endif
    return;
  }

  const Point& closest_p = traversal_traits.closest_point();

  // The filtering in the AABB tree checks the dist query <-> node bbox, which might be smaller than
  // the actual distance between the query <-> closest point
  edge_descriptor closest_e = traversal_traits.closest_primitive_id();
  bool is_close_enough =
    gt.compare_squared_distance_3_object()(query,
                                           gt.construct_segment_3_object()(get(vpm_T, source(closest_e, tm_T)),
                                                                           get(vpm_T, target(closest_e, tm_T))),
                                           sq_tolerance) != LARGER;

#ifdef CGAL_PMP_SNAP_DEBUG_PP
  std::cout << "  Closest point: (" << closest_p << ")" << std::endl
            << "  on edge: " << traversal_traits.closest_primitive_id()
            << "  at sq distance " << gt.compute_squared_distance_3_object()(query, closest_p)
            << " with squared tolerance: " << sq_tolerance
            << " && close enough? " << is_close_enough << std::endl;
#endif

  if(!is_close_enough)
    return;


  CGAL_assertion(get(vpm_T, source(closest_e, tm_T)) != query &&
                 get(vpm_T, target(closest_e, tm_T)) != query);

  halfedge_descriptor closest_h = halfedge(closest_e, tm_T);
  CGAL_assertion(is_border(edge(closest_h, tm_T), tm_T));

  if(!is_border(closest_h, tm_T))
    closest_h = opposite(closest_h, tm_T);

  // Using a map because we need to know if the same halfedge is split multiple times
  Edges_to_split_map_inserter<ConcurrencyTag>()(edges_to_split, closest_h, vertex_with_tolerance.first, closest_p);
}

#ifdef CGAL_LINKED_WITH_TBB
template <typename PointWithToleranceContainer,
          typename TriangleMesh, typename EdgeToSplitMap, typename AABBTree,
          typename VertexPatchMap_S, typename FacePatchMap_T,
          typename VPMS, typename VPMT, typename GT>
struct Find_splittable_edge_for_parallel_for
{
  Find_splittable_edge_for_parallel_for(const PointWithToleranceContainer& points_with_tolerance,
                                        EdgeToSplitMap& edges_to_split,
                                        const AABBTree* aabb_tree_ptr,
                                        const TriangleMesh& tm_S,
                                        const VertexPatchMap_S vertex_patch_map_S,
                                        const VPMS vpm_S,
                                        const TriangleMesh& tm_T,
                                        const FacePatchMap_T face_patch_map_T,
                                        const VPMT vpm_T,
                                        const GT& gt)
    :
      m_points_with_tolerance(points_with_tolerance),
      m_edges_to_split(edges_to_split), m_aabb_tree_ptr(aabb_tree_ptr),
      m_tm_S(tm_S), m_vertex_patch_map_S(vertex_patch_map_S), m_vpm_S(vpm_S),
      m_tm_T(tm_T), m_face_patch_map_T(face_patch_map_T), m_vpm_T(vpm_T),
      m_gt(gt)
  { }

  void operator()(const tbb::blocked_range<size_t>& r) const
  {
    for(std::size_t i=r.begin(); i!=r.end(); ++i)
    {
      find_splittable_edge<CGAL::Parallel_tag>(m_points_with_tolerance[i], m_edges_to_split, m_aabb_tree_ptr,
                                               m_tm_S, m_vertex_patch_map_S, m_vpm_S,
                                               m_tm_T, m_face_patch_map_T, m_vpm_T, m_gt);
    }
  }

private:
  const PointWithToleranceContainer& m_points_with_tolerance;
  EdgeToSplitMap& m_edges_to_split;
  const AABBTree* m_aabb_tree_ptr;
  const TriangleMesh& m_tm_S;
  const VertexPatchMap_S m_vertex_patch_map_S;
  const VPMS m_vpm_S;
  const TriangleMesh& m_tm_T;
  const FacePatchMap_T m_face_patch_map_T;
  const VPMT m_vpm_T;
  const GT& m_gt;
};
#endif

template <typename EdgesToSplitContainer,
          typename TriangleMesh, typename GeomTraits,
          typename VPMS, typename VPMT>
std::size_t split_edges(EdgesToSplitContainer& edges_to_split,
                        TriangleMesh& tm_S,
                        VPMS vpm_S,
                        TriangleMesh& tm_T,
                        VPMT vpm_T,
                        const GeomTraits& gt,
                        const bool is_source_mesh_fixed) // when snapping is B --> A and the mesh B is fixed
{
#ifdef CGAL_PMP_SNAP_DEBUG
  std::cout << "split " << edges_to_split.size() << " edges" << std::endl;
#endif

  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor           vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor         halfedge_descriptor;

  typedef typename boost::property_traits<VPMS>::value_type                       Point;
  typedef typename boost::property_traits<VPMT>::reference                        Point_ref;
  typedef typename GeomTraits::Vector_3                                           Vector;

  typedef std::pair<halfedge_descriptor, Point>                                   Vertex_with_new_position;
  typedef std::vector<Vertex_with_new_position>                                   Vertices_with_new_position;
  typedef std::pair<const halfedge_descriptor, Vertices_with_new_position>        Edge_and_splitters;

  std::size_t snapped_n = 0;

  CGAL::Real_timer timer;
  timer.start();

  for(Edge_and_splitters& es : edges_to_split)
  {
    halfedge_descriptor h_to_split = es.first;
    CGAL_assertion(is_border(h_to_split, tm_T));

    Vertices_with_new_position& splitters = es.second;

    if(splitters.size() > 1)
    {
#ifdef CGAL_PMP_SNAP_DEBUG_PP
      std::cout << " _______ Multiple splitting points on the same halfedge" << std::endl;
#endif

      const Point_ref hsp = get(vpm_T, source(h_to_split, tm_T));
      std::sort(splitters.begin(), splitters.end(),
                [&](const Vertex_with_new_position& l, const Vertex_with_new_position& r) -> bool
                {
                  return gt.less_distance_to_point_3_object()(hsp, l.second, r.second);
                });
    }

    // Inserting new points ordered from the source to the target of (the initial) 'h'
    bool first_split = true;
    Point previous_split_position = get(vpm_S, *(vertices(tm_S).begin())); // dummy value to avoid "used uninitialized" warning
    for(const Vertex_with_new_position& vnp : splitters)
    {
      const halfedge_descriptor splitter_h = vnp.first;
      const vertex_descriptor splitter_v = target(splitter_h, tm_S);
      const Point new_position = is_source_mesh_fixed ? get(vpm_S, splitter_v) : vnp.second;

      bool do_split = true;

      // Some splits can create degenerate faces, avoid that
      if((new_position == get(vpm_T, target(h_to_split, tm_T))) ||
         (new_position == get(vpm_T, source(h_to_split, tm_T))))
        do_split = false;

      if(!first_split && new_position == previous_split_position)
        do_split = false;

#ifdef CGAL_PMP_SNAP_DEBUG_PP
      std::cout << " -.-.-. Splitting " << edge(h_to_split, tm_T) << " |||| "
                << " Vs " << source(h_to_split, tm_T) << " (" << tm_T.point(source(h_to_split, tm_T)) << ")"
                << " --- Vt " << target(h_to_split, tm_T) << " (" << tm_T.point(target(h_to_split, tm_T)) << ")" << std::endl;
      std::cout << "With point: " << vnp.second << std::endl;
      std::cout << "Actually split? " << do_split << std::endl;
#endif

      // Split and update positions
      vertex_descriptor new_v = boost::graph_traits<TriangleMesh>::null_vertex();
      if(do_split)
      {
        CGAL::Euler::split_edge(h_to_split, tm_T);
        new_v = source(h_to_split, tm_T);
        put(vpm_T, new_v, new_position); // position of the new point on the target mesh
      }

      if(!is_source_mesh_fixed)
        put(vpm_S, splitter_v, new_position);

      first_split = false;
      previous_split_position = new_position;
      ++snapped_n;

      // Everything below is choosing the diagonal to triangulate the quad formed by the edge split
      // So, it's only relevant if splitting has been performed
      if(!do_split)
        continue;

      /*          new_p
       *         /   \
       *    res /     \ h_to_split
       *       /       \
       *      /         \
       *    left       right
       *     |         /
       *     |        /
       *     |       /
       *     |      /
       *     |     /
       *     |    /
       *      opp
       */

      const halfedge_descriptor res = prev(h_to_split, tm_T);
      const Point_ref left_pt = get(vpm_T, source(res, tm_T));
      const Point_ref right_pt = get(vpm_T, target(h_to_split, tm_T));
      const Point_ref opp = get(vpm_T, target(next(opposite(res, tm_T), tm_T), tm_T));

      // Check if 'p' is "visible" from 'opp' (i.e. its projection on the plane 'Pl(left, opp, right)'
      // falls in the cone with apex 'opp' and sides given by 'left' and 'right')
      const Vector n = gt.construct_orthogonal_vector_3_object()(right_pt, left_pt, opp);
      const Point trans_left_pt = gt.construct_translated_point_3_object()(left_pt, n);
      const Point trans_right_pt = gt.construct_translated_point_3_object()(right_pt, n);

      const Point_ref new_p = get(vpm_T, new_v);
      const bool left_of_left = (gt.orientation_3_object()(trans_left_pt, left_pt, opp, new_p) == CGAL::POSITIVE);
      const bool right_of_right = (gt.orientation_3_object()(right_pt, trans_right_pt, opp, new_p) == CGAL::POSITIVE);

      const bool is_visible = (!left_of_left && !right_of_right);

#ifdef CGAL_PMP_SNAP_DEBUG_PP
      std::cout << "Left/Right: " << left_of_left << " " << right_of_right << std::endl;
      std::cout << "visible from " << opp << " ? " << is_visible << std::endl;
#endif

      // h_to_split is equal to 'next(res)' after splitting
      const halfedge_descriptor h_to_split_opp = opposite(h_to_split, tm_T);

      if(is_visible)
      {
        halfedge_descriptor new_hd = CGAL::Euler::split_face(h_to_split_opp,
                                                             prev(prev(h_to_split_opp, tm_T), tm_T), tm_T);
        h_to_split = opposite(prev(new_hd, tm_T), tm_T);
      }
      else
      {
        halfedge_descriptor new_hd = CGAL::Euler::split_face(opposite(res, tm_T),
                                                             prev(h_to_split_opp, tm_T), tm_T);
        h_to_split = opposite(next(new_hd, tm_T), tm_T);
      }
    }
  }

  return snapped_n;
}

// Collect border points that can be projected onto a border edge
//
// If multiple vertices on the source border have the same geometric position, they are moved
// all at once and we use the smallest tolerance from all source vertices with that position.
// This is done to solve situation such as:
//          ||
//          ||
//      ____||____
//      ----------
// (vertex-vertex snapping into non-conformal snapping) in one go.
template <typename ConcurrencyTag = CGAL::Sequential_tag,
          typename HalfedgeRange, typename TriangleMesh,
          typename LockedVertexMap, typename LockedHalfedgeMap, typename ToleranceMap,
          typename VertexPatchMap_S, typename FacePatchMap_T,
          typename SourceNamedParameters, typename TargetNamedParameters>
std::size_t snap_non_conformal_one_way(const HalfedgeRange& halfedge_range_S,
                                       TriangleMesh& tm_S,
                                       ToleranceMap tolerance_map_S,
                                       VertexPatchMap_S vertex_patch_map_S,
                                       LockedVertexMap locked_vertices_S,
                                       const HalfedgeRange& halfedge_range_T,
                                       TriangleMesh& tm_T,
                                       FacePatchMap_T face_patch_map_T,
                                       LockedHalfedgeMap locked_halfedges_T,
                                       const bool is_source_mesh_fixed,
                                       const SourceNamedParameters& snp,
                                       const TargetNamedParameters& tnp)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor           vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor         halfedge_descriptor;

  typedef typename GetGeomTraits<TriangleMesh, TargetNamedParameters>::type       GT;
  typedef typename GT::FT                                                         FT;

  typedef typename GetVertexPointMap<TriangleMesh, SourceNamedParameters>::type   VPMS;
  typedef typename GetVertexPointMap<TriangleMesh, TargetNamedParameters>::type   VPMT;

  typedef typename boost::property_traits<VPMT>::value_type                       Point;

  typedef CGAL::AABB_halfedge_graph_segment_primitive<TriangleMesh, VPMT>         Primitive;
  typedef CGAL::AABB_traits<GT, Primitive>                                        AABB_Traits;
  typedef CGAL::AABB_tree<AABB_Traits>                                            AABB_tree;

  typedef std::pair<halfedge_descriptor, Point>                                   Vertex_with_new_position;
  typedef std::vector<Vertex_with_new_position>                                   Vertices_with_new_position;

  using parameters::get_parameter;
  using parameters::choose_parameter;

  VPMS vpm_S = choose_parameter(get_parameter(snp, internal_np::vertex_point), get_property_map(vertex_point, tm_S));
  VPMT vpm_T = choose_parameter(get_parameter(tnp, internal_np::vertex_point), get_property_map(vertex_point, tm_T));
  const GT gt = choose_parameter<GT>(get_parameter(snp, internal_np::geom_traits));

#ifdef CGAL_PMP_SNAP_DEBUG
  std::cout << "Gather unique points in source range..." << std::endl;
#endif

  typedef std::pair<halfedge_descriptor, FT>                                      Vertex_with_tolerance;
  typedef std::vector<Vertex_with_tolerance>                                      Vertices_with_tolerance;

  Vertices_with_tolerance vertices_to_snap;
  vertices_to_snap.reserve(halfedge_range_S.size()); // ensures that iterators stay valid

  // Take the min tolerance for all points that have the same coordinates
  std::map<Point, FT> point_tolerance_map;

  for(halfedge_descriptor h : halfedge_range_S)
  {
    if(get(locked_vertices_S, target(h, tm_S)))
      continue;

    // Skip the source vertex if its two incident halfedges are geometrically identical (it means that
    // the two halfedges are already stitchable and we don't want this common vertex to be used
    // to split a halfedge somewhere else)
    if(get(vpm_S, source(h, tm_S)) == get(vpm_S, target(next(h, tm_S), tm_S)))
      continue;

    const vertex_descriptor v = target(h, tm_S);
    const FT tolerance = get(tolerance_map_S, v);

    vertices_to_snap.emplace_back(h, tolerance);

    std::pair<Point, FT> entry(get(vpm_S, v), tolerance);
    std::pair<typename std::map<Point, FT>::iterator, bool> is_insert_successful =
      point_tolerance_map.insert(entry);
    if(!is_insert_successful.second)
      is_insert_successful.first->second = (std::min)(is_insert_successful.first->second, tolerance);

    #ifdef CGAL_PMP_SNAP_DEBUG_PP
    std::cout << "Non-conformal query: " << v << " (" << get(vpm_S, v) << "), tolerance: " << tolerance << std::endl;
#endif
  }

  for(auto& p : vertices_to_snap)
    p.second = point_tolerance_map[get(vpm_S, target(p.first, tm_S))];

  // Since we're inserting primitives one by one, we can't pass this shared data in the constructor of the tree
  AABB_Traits aabb_traits;
  aabb_traits.set_shared_data(tm_T, vpm_T);
  AABB_tree aabb_tree(aabb_traits);

  for(halfedge_descriptor h : halfedge_range_T)
  {
    CGAL_precondition(is_border(edge(h, tm_T), tm_T));
    if(get(locked_halfedges_T, h))
    {
#ifdef CGAL_PMP_SNAP_DEBUG_PP
      std::cout << edge(h, tm_T) << " is locked and not a valid target" << std::endl;
#endif
      continue;
    }

    aabb_tree.insert(Primitive(edge(h, tm_T), tm_T, vpm_T));
  }

  // Now, check which edges of the target range ought to be split by source vertices
#ifdef CGAL_PMP_SNAP_DEBUG_PP
  std::cout << "Collect edges to split with " << vertices_to_snap.size() << " vertices" << std::endl;
#endif

#ifndef CGAL_LINKED_WITH_TBB
  CGAL_static_assertion_msg (!(std::is_convertible<ConcurrencyTag, Parallel_tag>::value),
                             "Parallel_tag is enabled but TBB is unavailable.");
#else
  if(std::is_convertible<ConcurrencyTag, Parallel_tag>::value)
  {
#ifdef CGAL_PMP_SNAP_DEBUG
    std::cout << "Parallel find splittable edges!" << std::endl;
#endif

    typedef tbb::concurrent_hash_map<halfedge_descriptor,
                                     Vertices_with_new_position>                 Concurrent_edge_to_split_container;
    typedef internal::Find_splittable_edge_for_parallel_for<
              Vertices_with_tolerance, TriangleMesh,
              Concurrent_edge_to_split_container, AABB_tree,
              VertexPatchMap_S, FacePatchMap_T, VPMS, VPMT, GT>                  Functor;

    CGAL::Real_timer timer;
    timer.start();

    Concurrent_edge_to_split_container edges_to_split;
    Functor f(vertices_to_snap, edges_to_split, &aabb_tree,
              tm_S, vertex_patch_map_S, vpm_S, tm_T, face_patch_map_T, vpm_T, gt);
    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, vertices_to_snap.size()), f);

    std::cout << "Time to gather edges: " << timer.time() << std::endl;

    return split_edges(edges_to_split, tm_S, vpm_S,  tm_T, vpm_T, gt, is_source_mesh_fixed);
  }
  else
#endif // CGAL_LINKED_WITH_TBB
  {
#ifdef CGAL_PMP_SNAP_DEBUG
    std::cout << "Sequential find splittable edges!" << std::endl;
#endif

    std::map<halfedge_descriptor, Vertices_with_new_position> edges_to_split;
    for(const Vertex_with_tolerance& vt : vertices_to_snap)
    {
      internal::find_splittable_edge(vt, edges_to_split, &aabb_tree,
                                     tm_S, vertex_patch_map_S, vpm_S,
                                     tm_T, face_patch_map_T, vpm_T, gt);
    }

    return split_edges(edges_to_split, tm_S, vpm_S, tm_T, vpm_T,  gt, is_source_mesh_fixed);
  }
}

// \ingroup PMP_repairing_grp
//
// Attempts to snap the vertices in `halfedge_range_A` onto edges of `halfedge_range_B`, and reciprocally.
// A vertex from the first range is only snapped to an edge of the second range if the distance to
// the edge is smaller than the tolerance prescribed at the vertex.
//
// \warning This function does not give any guarantee on the conformity between the two meshes after the snapping.
// \warning This function does not merge vertices or the meshes, it is purely geometric.
//
// \tparam TriangleMesh a model of `FaceListGraph` and `MutableFaceGraph`
// \tparam HalfedgeRange_A a model of `Range` with value type `boost::graph_traits<TriangleMesh>::%halfedge_descriptor`
// \tparam HalfedgeRange_B a model of `Range` with value type `boost::graph_traits<TriangleMesh>::%halfedge_descriptor`
// \tparam ToleranceMap_A a model of `ReadablePropertyMap` with key type `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
//                        and value type `GetGeomTraits<TriangleMesh, NamedParameters_A>::type::FT`
// \tparam ToleranceMap_B a model of `ReadablePropertyMap` with key type `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
//                        and value type `GetGeomTraits<TriangleMesh, NamedParameters_A>::type::FT`
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
//     an internal property map for `CGAL::vertex_point_t` must be available in `TriangleMesh`.
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
//     an internal property map for `CGAL::vertex_point_t` must be available in `TriangleMesh`.
//   \cgalParamEnd
// \cgalNamedParamsEnd
//
// \return the number of snapped vertex pairs
//
template <typename ConcurrencyTag = CGAL::Sequential_tag,
          typename HalfedgeRange, typename TriangleMesh,
          typename ToleranceMap_A, typename ToleranceMap_B,
          typename NamedParameters_A, typename NamedParameters_B>
std::size_t snap_non_conformal(HalfedgeRange& halfedge_range_A,
                               TriangleMesh& tm_A,
                               ToleranceMap_A tolerance_map_A,
                               HalfedgeRange& halfedge_range_B,
                               TriangleMesh& tm_B,
                               ToleranceMap_B tolerance_map_B,
                               const bool is_self_snapping, // == true if range and meshes are equal
                               const NamedParameters_A& np_A,
                               const NamedParameters_B& np_B)
{
#ifdef CGAL_PMP_SNAP_DEBUG
  std::cout.precision(17);
  std::cerr.precision(17);
#endif

  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor           vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor         halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor             face_descriptor;

  typedef CGAL::dynamic_vertex_property_t<bool>                                   Vertex_bool_tag;
  typedef CGAL::dynamic_vertex_property_t<std::size_t>                            Vertex_size_t_tag;
  typedef CGAL::dynamic_halfedge_property_t<bool>                                 Halfedge_bool_tag;

  typedef typename boost::property_map<TriangleMesh, Vertex_bool_tag>::type       Locked_vertices;
  typedef typename boost::property_map<TriangleMesh, Halfedge_bool_tag>::type     Locked_halfedges;
  typedef typename boost::property_map<TriangleMesh, Vertex_size_t_tag>::type     Vertex_patch;

  typedef typename internal_np::Lookup_named_param_def <
    internal_np::face_patch_t, NamedParameters_A,
      Constant_property_map<face_descriptor, std::size_t> /*default*/ >::type     Face_patch_map_A;
  typedef typename internal_np::Lookup_named_param_def <
    internal_np::face_patch_t, NamedParameters_B,
      Constant_property_map<face_descriptor, std::size_t> /*default*/ >::type     Face_patch_map_B;

  using CGAL::parameters::choose_parameter;
  using CGAL::parameters::is_default_parameter;
  using CGAL::parameters::get_parameter;

  const bool is_same_mesh = (&tm_A == &tm_B);
  const bool simplify_first_mesh = choose_parameter(get_parameter(np_A, internal_np::do_simplify_border), false);
  const bool simplify_second_mesh = choose_parameter(get_parameter(np_B, internal_np::do_simplify_border), false);
  const bool is_second_mesh_fixed = choose_parameter(get_parameter(np_B, internal_np::do_lock_mesh), false);

  // vertex-vertex and vertex-edge snapping is only considered within compatible patches
  Face_patch_map_A face_patch_map_A = choose_parameter(get_parameter(np_A, internal_np::face_patch),
                                                       Constant_property_map<face_descriptor, std::size_t>(-1));
  Face_patch_map_B face_patch_map_B = choose_parameter(get_parameter(np_B, internal_np::face_patch),
                                                       Constant_property_map<face_descriptor, std::size_t>(-1));

  Vertex_patch vertex_patch_map_A = get(Vertex_size_t_tag(), tm_A);
  Vertex_patch vertex_patch_map_B = get(Vertex_size_t_tag(), tm_B);

  for(const halfedge_descriptor h : halfedge_range_A)
  {
    CGAL_precondition(is_border(h, tm_A));
    halfedge_descriptor h_opp = opposite(h, tm_A);
    for(const vertex_descriptor v : vertices_around_face(h_opp, tm_A))
      put(vertex_patch_map_A, v, get(face_patch_map_A, face(h_opp, tm_A)));
  }

  // @todo avoid that when 'self_snapping' is true
  for(const halfedge_descriptor h : halfedge_range_B)
  {
    CGAL_precondition(is_border(h, tm_B));
    halfedge_descriptor h_opp = opposite(h, tm_B);
    for(const vertex_descriptor v : vertices_around_face(h_opp, tm_B))
      put(vertex_patch_map_B, v, get(face_patch_map_B, face(h_opp, tm_B)));
  }

#ifdef CGAL_PMP_SNAP_DEBUG
  std::cout << "Non-conformal snapping... Range sizes: "
            << std::distance(halfedge_range_A.begin(), halfedge_range_A.end()) << " and "
            << std::distance(halfedge_range_B.begin(), halfedge_range_B.end()) << std::endl;
#endif

  CGAL_expensive_precondition(is_valid_polygon_mesh(tm_S) && is_triangle_mesh(tm_S));
  CGAL_expensive_precondition(is_valid_polygon_mesh(tm_T) && is_triangle_mesh(tm_T));

  // Steps:
  // - #1 Simplify the source range
  // - #1bis Simplify the target range (if it's not locked)
  // - #2 two-way vertex-vertex snapping
  // - #3 two-way vertex-edge snapping

  //////////////////////////////////////////////////////////////////////////////////////////////////
  /// #1 and #1bis (Simplification of borders)
  //////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef CGAL_PMP_SNAP_DEBUG
  std::cout << "Simplify ranges (" << simplify_first_mesh << " " << simplify_second_mesh << ")..." << std::endl;
#endif

  if(simplify_first_mesh)
  {
    internal::simplify_range(halfedge_range_A, tm_A, tolerance_map_A, np_A);

#ifdef CGAL_PMP_SNAP_DEBUG_OUTPUT
    std::ofstream("results/simplified_A.off") << std::setprecision(17) << tm_A;
#endif
  }

  if(!is_self_snapping && !is_second_mesh_fixed && simplify_second_mesh)
  {
    internal::simplify_range(halfedge_range_B, tm_B, tolerance_map_B, np_B);

#ifdef CGAL_PMP_SNAP_DEBUG_OUTPUT
    std::ofstream("results/simplified_B.off") << std::setprecision(17) << tm_B;
#endif
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  /// #2 (Two-way vertex-vertex snapping)
  //////////////////////////////////////////////////////////////////////////////////////////////////

  // We keep in memory pairs of source/target edges that are stitchable after vertex-vertex snapping
  // --> these halfedges should not be considered as targets in non-conformal snapping
  // Similarly, matching vertices whose incident edges have matching directions are also locked
  Locked_vertices locked_vertices_A = get(Vertex_bool_tag(), tm_A);
  Locked_vertices locked_vertices_B = get(Vertex_bool_tag(), tm_B);
  Locked_halfedges locked_halfedges_A = get(Halfedge_bool_tag(), tm_A);
  Locked_halfedges locked_halfedges_B = get(Halfedge_bool_tag(), tm_B);

  std::size_t snapped_n = 0;

  std::vector<std::pair<vertex_descriptor, vertex_descriptor> > locked_vertices;
  std::vector<halfedge_descriptor> locked_halfedges_A_vector, locked_halfedges_B_vector;

  snapped_n += internal::snap_vertices_two_way<ConcurrencyTag>(
                 halfedge_range_A, tm_A, tolerance_map_A, vertex_patch_map_A,
                 halfedge_range_B, tm_B, tolerance_map_B, vertex_patch_map_B,
                 std::back_inserter(locked_vertices),
                 std::back_inserter(locked_halfedges_A_vector),
                 std::back_inserter(locked_halfedges_B_vector),
                 is_second_mesh_fixed, np_A, np_B);

  for(const auto& vpair : locked_vertices)
  {
    put(locked_vertices_A, vpair.first, true);
    put(locked_vertices_B, vpair.second, true);
    if(is_same_mesh)
    {
      put(locked_vertices_B, vpair.first, true);
      put(locked_vertices_A, vpair.second, true);
    }
  }

  for(const halfedge_descriptor h : locked_halfedges_A_vector)
    put(locked_halfedges_A, h, true);
  for(const halfedge_descriptor h : locked_halfedges_B_vector)
    put(locked_halfedges_B, h, true);

  if(is_same_mesh)
  {
    for(const halfedge_descriptor h : locked_halfedges_A_vector)
      put(locked_halfedges_B, h, true);
    for(const halfedge_descriptor h : locked_halfedges_B_vector)
      put(locked_halfedges_A, h, true);
  }

#ifdef CGAL_PMP_SNAP_DEBUG_OUTPUT
  std::ofstream("results/vertex_vertex_A.off") << std::setprecision(17) << tm_A;
  std::ofstream("results/vertex_vertex_B.off") << std::setprecision(17) << tm_B;
#endif

  //////////////////////////////////////////////////////////////////////////////////////////////////
  /// #3 (Two one-way vertex-edge snapping)
  //////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef CGAL_PMP_SNAP_DEBUG
  std::cout << " ///////////// Two one-way vertex-edge snapping (A --> B) " << std::endl;
#endif

  snapped_n += internal::snap_non_conformal_one_way<ConcurrencyTag>(
                 halfedge_range_A, tm_A, tolerance_map_A, vertex_patch_map_A, locked_vertices_A,
                 halfedge_range_B, tm_B, face_patch_map_B, locked_halfedges_B,
                 false /*source is never fixed*/, np_A, np_B);

#ifdef CGAL_PMP_SNAP_DEBUG_OUTPUT
  std::ofstream("results/vertex_edge_A.off") << std::setprecision(17) << tm_A;
  std::ofstream("results/vertex_edge_B.off") << std::setprecision(17) << tm_B;
#endif

  if(!is_self_snapping)
  {
#ifdef CGAL_PMP_SNAP_DEBUG
  std::cout << " ///////////// Two one-way vertex-edge snapping (B --> A) " << std::endl;
#endif

    snapped_n += internal::snap_non_conformal_one_way<ConcurrencyTag>(
                   halfedge_range_B, tm_B, tolerance_map_B, vertex_patch_map_B, locked_vertices_B,
                   halfedge_range_A, tm_A, face_patch_map_A, locked_halfedges_A,
                   is_second_mesh_fixed, np_B, np_A);
  }

  return snapped_n;
}

} // namespace internal

namespace experimental {

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Convenience overloads
////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename ConcurrencyTag = CGAL::Sequential_tag,
          typename TriangleMesh,
          typename ToleranceMap_A, typename ToleranceMap_B,
          typename NamedParameters_A, typename NamedParameters_B>
std::size_t snap_borders(TriangleMesh& tm_A,
                         ToleranceMap_A tolerance_map_A,
                         TriangleMesh& tm_B,
                         ToleranceMap_B tolerance_map_B,
                         const NamedParameters_A& np_A,
                         const NamedParameters_B& np_B)
{
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor        halfedge_descriptor;

  std::vector<halfedge_descriptor> border_vertices_A;
  border_halfedges(tm_A, std::back_inserter(border_vertices_A));
  std::vector<halfedge_descriptor> border_vertices_B;
  border_halfedges(tm_B, std::back_inserter(border_vertices_B));

  return internal::snap_non_conformal<ConcurrencyTag>(border_vertices_A, tm_A, tolerance_map_A,
                                                      border_vertices_B, tm_B, tolerance_map_B,
                                                      false /*not self snapping*/, np_A, np_B);
}

template <typename ConcurrencyTag = CGAL::Sequential_tag,
          typename TriangleMesh,
          typename NamedParameters_A, typename NamedParameters_B>
std::size_t snap_borders(TriangleMesh& tm_A,
                         TriangleMesh& tm_B,
                         const NamedParameters_A& np_A,
                         const NamedParameters_B& np_B)
{
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor        halfedge_descriptor;

  typedef typename GetGeomTraits<TriangleMesh>::type                             GT;
  typedef typename GT::FT                                                        FT;

  typedef CGAL::dynamic_vertex_property_t<FT>                                    Vertex_property_tag;
  typedef typename boost::property_map<TriangleMesh, Vertex_property_tag>::type  Tolerance_map;

  std::vector<halfedge_descriptor> border_vertices_A;
  std::vector<halfedge_descriptor> border_vertices_B;
  border_halfedges(tm_A, std::back_inserter(border_vertices_A));
  border_halfedges(tm_B, std::back_inserter(border_vertices_B));

  const FT tol_mx((std::numeric_limits<double>::max)());
  Tolerance_map tolerance_map_A = get(Vertex_property_tag(), tm_A);
  internal::assign_tolerance_with_local_edge_length_bound(border_vertices_A, tolerance_map_A, tol_mx, tm_A, np_A);
  Tolerance_map tolerance_map_B = get(Vertex_property_tag(), tm_B);
  internal::assign_tolerance_with_local_edge_length_bound(border_vertices_B, tolerance_map_B, tol_mx, tm_B, np_B);

  return internal::snap_non_conformal<ConcurrencyTag>(border_vertices_A, tm_A, tolerance_map_A,
                                                      border_vertices_B, tm_B, tolerance_map_B,
                                                      false /*no self snapping*/, np_A, np_B);
}

template <typename ConcurrencyTag = CGAL::Sequential_tag,
          typename TriangleMesh,
          typename ToleranceMap_A, typename ToleranceMap_B>
std::size_t snap_borders(TriangleMesh& tm_A,
                         ToleranceMap_A tolerance_map_A,
                         TriangleMesh& tm_B,
                         ToleranceMap_B tolerance_map_B)
{
  return snap_borders<ConcurrencyTag>(tm_A, tolerance_map_A, tm_B, tolerance_map_B,
                                      CGAL::parameters::all_default(), CGAL::parameters::all_default());
}

template <typename ConcurrencyTag = CGAL::Sequential_tag,
          typename TriangleMesh>
std::size_t snap_borders(TriangleMesh& tm_A,
                         TriangleMesh& tm_B)
{
  return snap_borders<ConcurrencyTag>(tm_A, tm_B, CGAL::parameters::all_default(), CGAL::parameters::all_default());
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// Same as above, but with a single mesh

template <typename ConcurrencyTag = CGAL::Sequential_tag,
          typename TriangleMesh,
          typename ToleranceMap,
          typename CGAL_PMP_NP_TEMPLATE_PARAMETERS>
std::size_t snap_borders(TriangleMesh& tm,
                         ToleranceMap tolerance_map,
                         const CGAL_PMP_NP_CLASS& np)
{
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor        halfedge_descriptor;

  std::vector<halfedge_descriptor> border_vertices;
  border_halfedges(tm, std::back_inserter(border_vertices));

  return internal::snap_non_conformal<ConcurrencyTag>(border_vertices, tm, tolerance_map,
                                                      border_vertices, tm, tolerance_map,
                                                      true /*self snapping*/, np, np);
}

template <typename ConcurrencyTag = CGAL::Sequential_tag,
          typename TriangleMesh,
          typename ToleranceMap>
std::size_t snap_borders(TriangleMesh& tm,
                         ToleranceMap tolerance_map)
{
  return snap_borders<ConcurrencyTag>(tm, tolerance_map, CGAL::parameters::all_default());
}

template <typename ConcurrencyTag = CGAL::Sequential_tag,
          typename TriangleMesh,
          typename CGAL_PMP_NP_TEMPLATE_PARAMETERS>
std::size_t snap_borders(TriangleMesh& tm,
                         const CGAL_PMP_NP_CLASS& np)
{
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor        halfedge_descriptor;

  typedef typename GetGeomTraits<TriangleMesh>::type                             GT;
  typedef typename GT::FT                                                        FT;

  typedef CGAL::dynamic_vertex_property_t<FT>                                    Vertex_property_tag;
  typedef typename boost::property_map<TriangleMesh, Vertex_property_tag>::type  Tolerance_map;

  std::vector<halfedge_descriptor> border_vertices;
  border_halfedges(tm, std::back_inserter(border_vertices));

  const FT tol_mx((std::numeric_limits<double>::max)());
  Tolerance_map tolerance_map = get(Vertex_property_tag(), tm);
  internal::assign_tolerance_with_local_edge_length_bound(border_vertices, tolerance_map, tol_mx, tm, np);

  return internal::snap_non_conformal<ConcurrencyTag>(border_vertices, tm, tolerance_map,
                                                      border_vertices, tm, tolerance_map,
                                                      true /*self snapping*/, np, np);
}

template <typename ConcurrencyTag = CGAL::Sequential_tag,
          typename TriangleMesh>
std::size_t snap_borders(TriangleMesh& tm)
{
  return snap_borders<ConcurrencyTag>(tm, CGAL::parameters::all_default());
}

} // end namespace experimental
} // end namespace Polygon_mesh_processing
} // end namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_SNAPPING_SNAP_H
