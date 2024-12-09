// Copyright (c) 2019 GeometryFactory SARL (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s) : Sebastien Loriot, Martin Skrodzki, Dmitry Anisimov

#ifndef CGAL_PMP_INTERNAL_AABB_TRAVERSAL_TRAITS_WITH_HAUSDORFF_DISTANCE
#define CGAL_PMP_INTERNAL_AABB_TRAVERSAL_TRAITS_WITH_HAUSDORFF_DISTANCE

#include <CGAL/license/Polygon_mesh_processing/distance.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Bbox_3.h>

#include <iostream>
#include <queue>
#include <utility>
#include <vector>

namespace CGAL {
namespace Polygon_mesh_processing {
namespace internal {

template <class Kernel,
          class Face_handle_1,
          class Face_handle_2>
struct Bounds
{
  using FT = typename Kernel::FT;

  Bounds(const FT infinity_value) : lower(infinity_value), upper(infinity_value) { }
  Bounds(const FT lower, const FT upper) : lower(lower), upper(upper) { }

  FT lower;
  FT upper;
};

template <class Kernel,
          class Face_handle_1,
          class Face_handle_2>
struct Local_bounds
  : public Bounds<Kernel, Face_handle_1, Face_handle_2>
{
  using Base = Bounds<Kernel, Face_handle_1, Face_handle_2>;
  using FT = typename Kernel::FT;

  Local_bounds(const FT infinity_value) : Base(infinity_value) { }
  Local_bounds(const FT lower, const FT upper) : Base(lower, upper) { }

  Face_handle_2 tm2_lface = Face_handle_2();
  Face_handle_2 tm2_uface = Face_handle_2();
};

template <class Kernel,
          class Face_handle_1,
          class Face_handle_2>
struct Global_bounds
  : public Bounds<Kernel, Face_handle_1, Face_handle_2>
{
  using Base = Bounds<Kernel, Face_handle_1, Face_handle_2>;
  using FT = typename Kernel::FT;

  Global_bounds(const FT infinity_value) : Base(infinity_value) { }
  Global_bounds(const FT lower, const FT upper) : Base(lower, upper) { }

  std::pair<Face_handle_1, Face_handle_2> lpair = default_face_pair();
  std::pair<Face_handle_1, Face_handle_2> upair = default_face_pair();

  constexpr std::pair<Face_handle_1, Face_handle_2> default_face_pair() const
  {
    return std::make_pair(Face_handle_1(), Face_handle_2());
  }
};

// Candidate triangle.
template <class Kernel,
          class Face_handle_1,
          class Face_handle_2>
struct Candidate_triangle
{
  using FT = typename Kernel::FT;
  using Triangle_3 = typename Kernel::Triangle_3;
  using Candidate_bounds = Local_bounds<Kernel, Face_handle_1, Face_handle_2>;

  Candidate_triangle(const Triangle_3& triangle,
                     const Candidate_bounds& bounds,
                     const Face_handle_1 fh)
    : triangle(triangle), bounds(bounds), tm1_face(fh)
  { }

  Triangle_3 triangle;
  Candidate_bounds bounds;
  Face_handle_1 tm1_face;

  // Comparator for the priority queue.
  // Provide std::less for Candidate_triangle to have the largest 'upper' value at the top of the PQ
  bool operator<(const Candidate_triangle& other) const
  {
    CGAL_precondition(bounds.upper >= FT(0));
    CGAL_precondition(other.bounds.upper >= FT(0));

    return bounds.upper < other.bounds.upper;
  }
};

// Hausdorff primitive traits on TM2.
template <class Query,
          class Kernel,
          class TriangleMesh1,
          class TriangleMesh2,
          class VPM2>
class Hausdorff_primitive_traits_tm2
{
  using FT         = typename Kernel::FT;
  using Point_3    = typename Kernel::Point_3;
  using Vector_3   = typename Kernel::Vector_3;
  using Triangle_3 = typename Kernel::Triangle_3;

  using Project_point_3 = typename Kernel::Construct_projected_point_3;
  using Face_handle_1 = typename boost::graph_traits<TriangleMesh1>::face_descriptor;
  using Face_handle_2 = typename boost::graph_traits<TriangleMesh2>::face_descriptor;

  using Local_bounds = internal::Local_bounds<Kernel, Face_handle_1, Face_handle_2>;
  using Global_bounds = internal::Global_bounds<Kernel, Face_handle_1, Face_handle_2>;

  using TM2_face_to_triangle_map = Triangle_from_face_descriptor_map<TriangleMesh2, VPM2>;

public:
  using Priority = FT;

private:
  const Bbox_3& m_t1_bbox;
  const TriangleMesh2& m_tm2;
  const VPM2 m_vpm2;
  const TM2_face_to_triangle_map m_face_to_triangle_map;

  Local_bounds m_local_bounds; // local Hausdorff bounds for the query triangle
  const Global_bounds& m_global_bounds;
  FT m_v0_lower, m_v1_lower, m_v2_lower;

  bool m_early_exit;

public:
  Hausdorff_primitive_traits_tm2(const Bbox_3& t1_bbox,
                                 const TriangleMesh2& tm2, const VPM2 vpm2,
                                 const Local_bounds& initial_bounds,
                                 const Global_bounds& global_bounds,
                                 const FT infinity_value)
    : m_t1_bbox(t1_bbox),
      m_tm2(tm2), m_vpm2(vpm2),
      m_face_to_triangle_map(&m_tm2, m_vpm2),
      m_local_bounds(initial_bounds),
      m_global_bounds(global_bounds),
      m_v0_lower(infinity_value),
      m_v1_lower(infinity_value),
      m_v2_lower(infinity_value),
      m_early_exit(false)
  { }

  // Return the local Hausdorff bounds computed for the passed query triangle.
  Local_bounds& get_local_bounds() { return m_local_bounds; }
  const Local_bounds& get_local_bounds() const { return m_local_bounds; }

  // Because
  //   h(TM1, TM2) := max_{query in TM1} h(query, TM2),
  // it is pointless to continue trying to find a smaller bound if the value is already known
  // to be below the current max computed through another TM1 face
  bool go_further() const { return !m_early_exit; }

  // Compute the explicit Hausdorff distance to the given primitive.
  template<class Primitive>
  void intersection(const Query& query, const Primitive& primitive)
  {
    if(m_early_exit)
      return;

#ifdef CGAL_HAUSDORFF_DEBUG_TM2_TRAVERSAL
    std::cout << "Intersection with TM2's " << primitive.id() << std::endl;
    std::cout << "Initial local bounds " << m_local_bounds.lower << " " << m_local_bounds.upper << std::endl;
#endif

    CGAL_assertion(m_local_bounds.lower >= FT(0));
    CGAL_assertion(m_local_bounds.upper >= FT(0));

    /* Have reached a single triangle, process it.
    / Determine the upper distance according to
    /   min_{b \in primitive} ( max_{vertex in query} ( d(vertex, b) ) )
    /
    / Here, we only have one triangle in B, i.e. tm2. Thus, it suffices to
    / compute the distance of the vertices of the query triangle to the
    / primitive triangle and use the maximum of the obtained distances.
    */

    // The query object is a triangle from TM1, get its vertices.
    const Point_3& v0 = query.vertex(0);
    const Point_3& v1 = query.vertex(1);
    const Point_3& v2 = query.vertex(2);

    CGAL_assertion(primitive.id() != Face_handle_2());
    const Triangle_3 triangle = get(m_face_to_triangle_map, primitive.id());

#ifdef CGAL_HAUSDORFF_DEBUG_TM2_TRAVERSAL
    std::cout << "Geometry: " << triangle << std::endl;
#endif

    // Compute distances of the vertices to the primitive triangle in TM2.
    const FT v0_dist = CGAL::squared_distance(v0, triangle);
    if(v0_dist < m_v0_lower)
      m_v0_lower = v0_dist;

    const FT v1_dist = CGAL::squared_distance(v1, triangle);
    if(v1_dist < m_v1_lower)
      m_v1_lower = v1_dist;

    const FT v2_dist = CGAL::squared_distance(v2, triangle);
    if(v2_dist < m_v2_lower)
      m_v2_lower = v2_dist;

    // Get the distance as maximizers over all vertices.
    //
    // Since we are at the level of a single triangle in TM2, distance_upper is actually
    // the Hausdorff distance from the query triangle in TM1 to the primitive triangle in TM2.

    // max_{v in query} (v, primitive), used in h_upper_i(query, TM2)
    const FT distance_upper = (CGAL::max)((CGAL::max)(v0_dist, v1_dist), v2_dist);

    // h_lower_i(query, TM2) := max_{v in query} min_{1<=j<=i} d(v, primitive_j)
    const FT distance_lower = (CGAL::max)((CGAL::max)(m_v0_lower, m_v1_lower), m_v2_lower);

#ifdef CGAL_HAUSDORFF_DEBUG_TM2_TRAVERSAL
    std::cout << "Distance from vertices of t1 to t2: " << v0_dist << " " << v1_dist << " " << v2_dist << std::endl;
#endif

    CGAL_assertion(distance_lower >= FT(0));
    CGAL_assertion(distance_upper >= distance_lower);

    // With each new TM2 face, the min value m_v{k}_lower can become smaller,
    // and thus also the value max_{v in query} min_{1<=j<=i} d(v, primitive_j)
    if(distance_lower < m_local_bounds.lower)
    {
#ifdef CGAL_HAUSDORFF_DEBUG_TM2_TRAVERSAL
      std::cout << "new best lower (" << distance_lower << ") with TM2 face: " << triangle << std::endl;
#endif
      m_local_bounds.lower = distance_lower;
      m_local_bounds.tm2_lface = primitive.id();
    }

    // This is the 'min_{1<=j<=i}' part in:
    //   h_upper_i(query, TM2) = min_{1<=j<=i} max_{v in query} (v, primitive_j), Equation (10)
    if(distance_upper < m_local_bounds.upper)
    {
#ifdef CGAL_HAUSDORFF_DEBUG_TM2_TRAVERSAL
      std::cout << "new best upper (" << distance_upper << ") with TM2 face: " << triangle << std::endl;
#endif
      m_local_bounds.upper = distance_upper;
      m_local_bounds.tm2_uface = primitive.id();
    }

#ifdef CGAL_HAUSDORFF_DEBUG_TM2_TRAVERSAL
    std::cout << "Distance from vertices of t1 to t2: " << v0_dist << " " << v1_dist << " " << v2_dist << std::endl;
    std::cout << "Current local bounds " << m_local_bounds.lower << " " << m_local_bounds.upper << std::endl;
#endif

    CGAL_assertion(m_local_bounds.lower >= FT(0));
    CGAL_assertion(m_local_bounds.lower <= m_local_bounds.upper);

// #define CGAL_PMP_HDIST_NO_CULLING_DURING_TRAVERSAL
#ifndef CGAL_PMP_HDIST_NO_CULLING_DURING_TRAVERSAL
    // the lhs can only go down with every additional TM2 face,
    // whereas the rhs can only go up with every additional TM1 face
    if(m_local_bounds.upper < m_global_bounds.lower) // Section 4.1, first §
    {
#ifdef CGAL_HAUSDORFF_DEBUG_TM2_TRAVERSAL
      std::cout << "Quitting early (TM2 traversal), global lower " << m_global_bounds.lower << " greater than local upper " << m_local_bounds.upper << std::endl;
#endif
      m_early_exit = true;
    }
#endif
  }

  // Determine whether child nodes will still contribute to a smaller
  // Hausdorff distance and thus have to be entered.
  template<class Node>
  std::pair<bool, Priority>
  do_intersect_with_priority(const Query&, const Node& node) const
  {
    if(m_early_exit)
      return std::make_pair(false, FT(0));

#ifdef CGAL_HAUSDORFF_DEBUG_TM2_TRAVERSAL
    std::cout << "Do_intersect TM2 node with bbox: " << node.bbox() << std::endl;
#endif

    // Compute a lower bound between the query (face of TM1) and a group of TM2 faces.
    const Bbox_3 node_bbox = node.bbox();

    // Distance along the x-axis.
    double dist_x = 0.;
    if(m_t1_bbox.xmax() < node_bbox.xmin())
      dist_x = node_bbox.xmin() - m_t1_bbox.xmax();
    else if(node_bbox.xmax() < m_t1_bbox.xmin())
      dist_x = m_t1_bbox.xmin() - node_bbox.xmax();

    // Distance along the y-axis.
    double dist_y = 0.;
    if(m_t1_bbox.ymax() < node_bbox.ymin())
      dist_y = node_bbox.ymin() - m_t1_bbox.ymax();
    else if(node_bbox.ymax() < m_t1_bbox.ymin())
      dist_y = m_t1_bbox.ymin() - node_bbox.ymax();

    // Distance along the z-axis.
    double dist_z = 0.;
    if(m_t1_bbox.zmax() < node_bbox.zmin())
      dist_z = node_bbox.zmin() - m_t1_bbox.zmax();
    else if(node_bbox.zmax() < m_t1_bbox.zmin())
      dist_z = m_t1_bbox.zmin() - node_bbox.zmax();

    const FT sq_dist = square(dist_x) + square(dist_y) + square(dist_z);

    // Culling on TM2:
    // The value 'dist' is the distance between bboxes and thus a lower bound on the distance
    // between the query and the TM2 primitives that are children of this node.
    // If this lower bound is greater than the current upper bound for this query,
    // then none of these primitives will reduce the Hausdorff distance between the query and TM2.
#ifdef CGAL_HAUSDORFF_DEBUG_TM2_TRAVERSAL
    std::cout << "Culling TM2? dist vs local bound upper " << sq_dist << " " << m_local_bounds.upper << std::endl;
#endif
    CGAL_assertion(m_local_bounds.upper >= FT(0));
    if(sq_dist > m_local_bounds.upper)
      return std::make_pair(false, FT(0));
    else
      return std::make_pair(true , -sq_dist);
  }

  template<class Node>
  bool do_intersect(const Query& query, const Node& node) const
  {
    if(m_early_exit)
      return false;

#ifdef CGAL_PMP_HDIST_NO_CULLING_DURING_TRAVERSAL
    CGAL_USE(query); CGAL_USE(node);
    return true;
#else
    return this->do_intersect_with_priority(query, node).first;
#endif
  }

  template<class PrimitiveConstIterator>
  void traverse_group(const Query& query,
                      PrimitiveConstIterator group_begin,
                      PrimitiveConstIterator group_end)
  {
    for(PrimitiveConstIterator it = group_begin; it != group_end; ++it)
      this->intersection(query, *it);
  }
};

// Hausdorff primitive traits on TM1.
template <class Query,
          class Kernel,
          class TriangleMesh1,
          class TriangleMesh2,
          class VPM1,
          class VPM2>
class Hausdorff_primitive_traits_tm1
{
  using FT         = typename Kernel::FT;
  using Point_3    = typename Kernel::Point_3;
  using Vector_3   = typename Kernel::Vector_3;
  using Triangle_3 = typename Kernel::Triangle_3;

  using TM2_primitive = AABB_face_graph_triangle_primitive<TriangleMesh2, VPM2>;
  using TM2_traits    = AABB_traits_3<Kernel, TM2_primitive>;
  using TM2_tree      = AABB_tree<TM2_traits>;
  using TM2_hd_traits = Hausdorff_primitive_traits_tm2<Triangle_3, Kernel, TriangleMesh1, TriangleMesh2, VPM2>;

  using TM1_face_to_triangle_map = Triangle_from_face_descriptor_map<TriangleMesh1, VPM1>;

  using Face_handle_1 = typename boost::graph_traits<TriangleMesh1>::face_descriptor;
  using Face_handle_2 = typename boost::graph_traits<TriangleMesh2>::face_descriptor;

  using Global_bounds = internal::Global_bounds<Kernel, Face_handle_1, Face_handle_2>;
  using Candidate = Candidate_triangle<Kernel, Face_handle_1, Face_handle_2>;
  using Heap_type = std::priority_queue<Candidate>;

public:
  using Priority = FT;

private:
  // Input data.
  const TriangleMesh1& m_tm1;
  const TriangleMesh2& m_tm2;
  const VPM1 m_vpm1;
  const VPM2 m_vpm2;
  const TM2_tree& m_tm2_tree;
  const TM1_face_to_triangle_map m_face_to_triangle_map;

  // Internal bounds and values.
  const FT m_sq_initial_bound;
  const FT m_sq_distance_bound;
  const FT m_infinity_value;
  Global_bounds m_global_bounds;
  bool m_early_exit;

  // All candidate triangles.
  Heap_type m_candidate_triangles;

public:
  Hausdorff_primitive_traits_tm1(const TM2_tree& tree,
                                 const TriangleMesh1& tm1,
                                 const TriangleMesh2& tm2,
                                 const VPM1 vpm1,
                                 const VPM2 vpm2,
                                 const FT infinity_value,
                                 const FT sq_initial_bound,
                                 const FT sq_distance_bound)
    : m_tm1(tm1), m_tm2(tm2),
      m_vpm1(vpm1), m_vpm2(vpm2),
      m_tm2_tree(tree),
      m_face_to_triangle_map(&m_tm1, m_vpm1),
      m_sq_initial_bound(sq_initial_bound),
      m_sq_distance_bound(sq_distance_bound),
      m_infinity_value(infinity_value),
      m_global_bounds(m_infinity_value),
      m_early_exit(false)
  {
    CGAL_precondition(m_infinity_value >= FT(0));

    // Bounds grow with every face of TM1 (Equation (6)).
    // If we initialize to zero here, then we are very slow even for big input error bounds!
    // Instead, we can use the error bound as our initial guess to filter out all pairs
    // which are already within this bound. It makes the code faster for close meshes.
    m_global_bounds.lower = m_sq_initial_bound;
    m_global_bounds.upper = m_sq_initial_bound;
  }

  // Return those triangles from TM1, which are candidates for including a
  // point realizing the Hausdorff distance.
  Heap_type& get_candidate_triangles() { return m_candidate_triangles; }

  // Return the global Hausdorff bounds computed for the passed query triangle.
  Global_bounds get_global_bounds()
  {
    CGAL_assertion(m_global_bounds.lower >= FT(0));
    CGAL_assertion(m_global_bounds.upper >= m_global_bounds.lower);

    update_global_bounds();
    return m_global_bounds;
  }

  // The maximum distance from one of the face corners to the second mesh, and the face realizing this max
  std::pair<FT, Face_handle_2> get_maximum_distance(const Face_handle_1 tm1_face) const
  {
    const Triangle_3 triangle = get(m_face_to_triangle_map, tm1_face);
    const Point_3& v0 = triangle.vertex(0);
    const Point_3& v1 = triangle.vertex(1);
    const Point_3& v2 = triangle.vertex(2);

    const auto pair0 = m_tm2_tree.closest_point_and_primitive(v0);
    const auto pair1 = m_tm2_tree.closest_point_and_primitive(v1);
    const auto pair2 = m_tm2_tree.closest_point_and_primitive(v2);

    const FT sq_dist0 = CGAL::squared_distance(v0, pair0.first);
    const FT sq_dist1 = CGAL::squared_distance(v1, pair1.first);
    const FT sq_dist2 = CGAL::squared_distance(v2, pair2.first);

    if(sq_dist0 > sq_dist1)
    {
      if(sq_dist0 > sq_dist2)
        return std::make_pair(sq_dist0, pair0.second);
      else
        return std::make_pair(sq_dist2, pair2.second);
    }
    else
    {
      if(sq_dist1 > sq_dist2)
        return std::make_pair(sq_dist1, pair1.second);
      else
        return std::make_pair(sq_dist2, pair2.second);
    }
  }

  // In case, we did not enter any loop, we set the realizing triangles here.
  void update_global_bounds()
  {
    if(m_candidate_triangles.size() > 0)
    {
      const Candidate& top = m_candidate_triangles.top();

      if(m_global_bounds.lpair.first == Face_handle_1())
        m_global_bounds.lpair.first = top.tm1_face;
      if(m_global_bounds.lpair.second == Face_handle_2())
        m_global_bounds.lpair.second = top.bounds.tm2_lface;

      if(m_global_bounds.upair.first == Face_handle_1())
        m_global_bounds.upair.first = top.tm1_face;
      if(m_global_bounds.upair.second == Face_handle_2())
        m_global_bounds.upair.second = top.bounds.tm2_uface;
    }
    else
    {
      Face_handle_1 tm1_f = *(faces(m_tm1).begin());
      const std::pair<FT, Face_handle_2> max_dist = get_maximum_distance(tm1_f);

      if(m_global_bounds.lpair.first == Face_handle_1())
        m_global_bounds.lpair.first = tm1_f;
      if(m_global_bounds.lpair.second == Face_handle_2())
        m_global_bounds.lpair.second = max_dist.second;

      if(m_global_bounds.upair.first == Face_handle_1())
        m_global_bounds.upair.first = tm1_f;
      if(m_global_bounds.upair.second == Face_handle_2())
        m_global_bounds.upair.second = max_dist.second;
    }
  }

  // Traversal-related
  bool early_exit() const { return m_early_exit; }

  // If the distance is already larger than the user-defined bound, traversal can stop
  bool go_further() const { return !m_early_exit; }

  // Compute Hausdorff distance bounds between a TM1 face and TM2
  template<class Primitive>
  void intersection(const Query&, const Primitive& primitive)
  {
    if(m_early_exit)
      return;

#ifdef CGAL_HAUSDORFF_DEBUG_TM1_TRAVERSAL
    std::cout << "Intersection with TM1's " << primitive.id() << std::endl;
    std::cout << "Initial global bounds " << m_global_bounds.lower << " " << m_global_bounds.upper << std::endl;
#endif

    // Set initial tight bounds.
    CGAL_assertion(primitive.id() != Face_handle_1());
    const Face_handle_1 tm1_face = primitive.id();
    const Triangle_3 triangle = get(m_face_to_triangle_map, tm1_face);

#ifdef CGAL_HAUSDORFF_DEBUG_TM1_TRAVERSAL
    std::cout << "Geometry: " << triangle << std::endl;
#endif

    // Call culling on TM2 with the TM1 triangle.
    const Bbox_3 t1_bbox = triangle.bbox();
    Local_bounds<Kernel, Face_handle_1, Face_handle_2> initial_bounds(m_infinity_value);
    TM2_hd_traits traversal_traits_tm2(t1_bbox, m_tm2, m_vpm2, initial_bounds, m_global_bounds, m_infinity_value);
    m_tm2_tree.traversal_with_priority(triangle, traversal_traits_tm2);

    // Post traversal, we have computed h_lower(query, TM2) and h_upper(query, TM2)
    const auto& local_bounds = traversal_traits_tm2.get_local_bounds();
#ifdef CGAL_HAUSDORFF_DEBUG_TM1_TRAVERSAL
    std::cout << "Bounds for TM1 primitive: " << local_bounds.lower << " " << local_bounds.upper << std::endl;
#endif

    CGAL_assertion(local_bounds.lower >= FT(0));
    CGAL_assertion(local_bounds.upper >= local_bounds.lower);
    CGAL_assertion(local_bounds.tm2_lface != boost::graph_traits<TriangleMesh2>::null_face());
    CGAL_assertion(local_bounds.tm2_uface != boost::graph_traits<TriangleMesh2>::null_face());

    // Update global Hausdorff bounds according to the obtained local bounds.
    // h_lower(TM1, TM2) = max_{query in TM1} h_lower(query, TM2)
    CGAL_assertion(m_global_bounds.lower >= FT(0));
    if(local_bounds.lower > m_global_bounds.lower) // Equation (6) in the paper, see also Algorithm 1, L.5
    {
      m_global_bounds.lower = local_bounds.lower;
      m_global_bounds.lpair.first = tm1_face;
      m_global_bounds.lpair.second = local_bounds.tm2_lface;
    }

    // h_upper(TM1, TM2) = max_{query in TM1} h_upper(query, TM2)
    CGAL_assertion(m_global_bounds.upper >= FT(0));
    if(local_bounds.upper > m_global_bounds.upper) // Equation (6) in the paper, see also Algorithm 1, L.8
    {
      m_global_bounds.upper = local_bounds.upper;
      m_global_bounds.upair.first = tm1_face;
      m_global_bounds.upair.second = local_bounds.tm2_uface;
    }

    CGAL_postcondition(m_global_bounds.upper >= m_global_bounds.lower);

    // Bounds only grow with each additional face of TM1 considered (Eq. (6)),
    // if the lower bound is already larger than the user-defined upper bound, we can stop
    if(is_positive(m_sq_distance_bound) &&
       m_sq_distance_bound <= m_global_bounds.lower)
    {
      m_early_exit = true;
      return;
    }

    // Store the TM1 triangle given as primitive in this function as a candidate triangle,
    // together with the local bounds it obtained to send it to subdivision later.
    m_candidate_triangles.emplace(triangle, local_bounds, tm1_face);
  }

  // Determine whether child nodes will still contribute to a larger
  // Hausdorff distance and thus have to be entered.
  template<class Node>
  std::pair<bool, Priority>
  do_intersect_with_priority(const Query&, const Node& node)
  {
    if(m_early_exit)
      return std::make_pair(false, FT(0));

#ifdef CGAL_HAUSDORFF_DEBUG_TM1_TRAVERSAL
    std::cout << "Do_intersect TM1 node with bbox " << node.bbox() << std::endl;
#endif

    // Compute an upper bound on the distance between the closest point in TM2 and
    // the corner of the bbox farthest from the closest point. This is an upper bound
    // on the Hausdorff distance between any children primitive of the node and TM2.
    //
    // @todo could find the bbox vertex that is closest to 'closest' to reduce a bit the bound
    const Bbox_3 bbox = node.bbox();
    const Point_3 bp(bbox.xmin(), bbox.ymin(), bbox.zmin());
    const Point_3 closest = m_tm2_tree.closest_point(bp);
    const Vector_3 difference(bp, closest);
    const Vector_3 diag(bbox.x_span(), bbox.y_span(), bbox.z_span());

    // @todo something better to avoid the sqrt
    const FT sq_dist = square(CGAL::approximate_sqrt(difference.squared_length())
                             + CGAL::approximate_sqrt(diag.squared_length()));

    // The Hausdorff distance grows with every TM1 face.
    // If the upper bound is smaller than the current global lower bound,
    // it is pointless to visit this node (and its children) because a larger distance
    // has been found somewhere else.
#ifdef CGAL_HAUSDORFF_DEBUG_TM1_TRAVERSAL
    std::cout << "Culling TM1? dist & global lower bound: " << sq_dist << " " << m_global_bounds.lower << std::endl;
#endif
    CGAL_assertion(m_global_bounds.lower >= FT(0));
    if(sq_dist < m_global_bounds.lower)
      return std::make_pair(false, FT(0));
    else
      return std::make_pair(true , +sq_dist);
  }

  template<class Node>
  bool do_intersect(const Query& query, const Node& node)
  {
    if(m_early_exit)
      return false;

#ifdef CGAL_PMP_HDIST_NO_CULLING_DURING_TRAVERSAL
    CGAL_USE(query); CGAL_USE(node);
    return true;
#else
    return this->do_intersect_with_priority(query, node).first;
#endif
  }

  template<class PrimitiveConstIterator>
  void traverse_group(const Query& query,
                      PrimitiveConstIterator group_begin,
                      PrimitiveConstIterator group_end)
  {
    CGAL_assertion_msg(false, "ERROR: we should not call the group traversal on TM1!");

    for(PrimitiveConstIterator it = group_begin; it != group_end; ++it)
      this->intersection(query, *it);
  }
};

} // namespace internal
} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif // CGAL_PMP_INTERNAL_AABB_TRAVERSAL_TRAITS_WITH_HAUSDORFF_DISTANCE
