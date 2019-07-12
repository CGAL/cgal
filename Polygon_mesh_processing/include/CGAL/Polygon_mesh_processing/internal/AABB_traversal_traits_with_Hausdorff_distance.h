// Copyright (c) 2019 GeometryFactory (France).
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
// Author(s) : Sebastien Loriot, Martin Skrodzki
//

#ifndef CGAL_PMP_INTERNAL_AABB_TRAVERSAL_TRAITS_WITH_HAUSDORFF_DISTANCE
#define CGAL_PMP_INTERNAL_AABB_TRAVERSAL_TRAITS_WITH_HAUSDORFF_DISTANCE

#include <iostream>

namespace CGAL {

  typedef std::pair<double, double> Hausdorff_bounds;

  /**
   * @class Hausdorff_primitive_traits_tm2
   */
  template<typename AABBTraits, typename Query, typename Kernel, typename TriangleMesh, typename VPM2>
  class Hausdorff_primitive_traits_tm2
  {
    typedef typename AABBTraits::Primitive Primitive;
    typedef ::CGAL::AABB_node<AABBTraits> Node;
    typedef typename AABBTraits::Bounding_box Bounding_box;
    typename Kernel::Construct_projected_point_3 project_point;
    typedef typename Kernel::Point_3 Point_3;

  public:
    Hausdorff_primitive_traits_tm2(
      const AABBTraits& traits,
      const TriangleMesh& tm2,
      const VPM2& vpm2,
      const double h_lower_init, const double h_upper_init,
      const double h_v0_lower_init, const double h_v1_lower_init, const double h_v2_lower_init
    )
      : m_traits(traits), m_tm2(tm2), m_vpm2(vpm2) {
        // Initialize the global and local bounds with the given values
        h_local_lower = h_lower_init;
        h_local_upper = h_upper_init;
        h_v0_lower = h_v0_lower_init;
        h_v1_lower = h_v1_lower_init;
        h_v2_lower = h_v2_lower_init;
      }

    // Explore the whole tree, i.e. always enter children if the methods
    // do_intersect() below determines that it is worthwhile.
    bool go_further() const { return true; }

    // Compute the explicit Hausdorff distance to the given primitive
    void intersection(const Query& query, const Primitive& primitive)
    {
      /* Have reached a single triangle, process it */

      /*
      / Determine the distance accroding to
      /       min_{b \in primitive} ( max_{vertex in query} ( d(vertex, b)))
      /
      / Here, we only have one triangle in B, i.e. tm2. Thus, it suffices to
      / compute the distance of the vertices of the query triangles to the
      / primitive triangle and use the maximum of the obtained distances.
      */

      // The query object is a triangle from TM1, get its vertices
      Point_3 v0 = query.vertex(0);
      Point_3 v1 = query.vertex(1);
      Point_3 v2 = query.vertex(2);

      // Compute distances of the vertices to the primitive triangle in TM2
      Triangle_from_face_descriptor_map<TriangleMesh, VPM2> face_to_triangle_map(&m_tm2, m_vpm2);
      double v0_dist = approximate_sqrt(squared_distance(
        project_point( get(face_to_triangle_map, primitive.id()), v0), v0 ) );
      if (v0_dist < h_v0_lower) h_v0_lower = v0_dist;

      double v1_dist = approximate_sqrt(squared_distance(
        project_point( get(face_to_triangle_map, primitive.id()), v1), v1 ) );
      if (v1_dist < h_v1_lower) h_v1_lower = v1_dist;

      double v2_dist = approximate_sqrt(squared_distance(
        project_point( get(face_to_triangle_map, primitive.id()), v2), v2 ) );
      if (v2_dist < h_v2_lower) h_v2_lower = v2_dist;

      // Get the distance as maximizers over all vertices
      double distance_upper = std::max(std::max(v0_dist,v1_dist),v2_dist);
      double distance_lower = std::max(std::max(h_v0_lower,h_v1_lower),h_v2_lower);

      // Since we are at the level of a single triangle in TM2, distance_upper is
      // actually the correct Hausdorff distance from the query triangle in
      // TM1 to the primitive triangle in TM2
      if ( distance_upper < h_local_upper ) h_local_upper = distance_upper;
      if ( distance_lower < h_local_lower ) h_local_lower = distance_lower;
    }

    // Determine whether child nodes will still contribute to a smaller
    // Hausdorff distance and thus have to be entered
    bool do_intersect(const Query& query, const Node& node) const
    {
      /* Have reached a node, determine whether or not to enter it */

      // Get the bounding box of the nodes
      Bounding_box bbox = node.bbox();
      // Compute its center
      Point_3 bb_center = Point_3(
        (bbox.min(0) + bbox.max(0)) / 2,
        (bbox.min(1) + bbox.max(1)) / 2,
        (bbox.min(2) + bbox.max(2)) / 2);

      // Get the center of the query triangle
      // The query object is a triangle from TM1, get its vertices
      Point_3 v0 = query.vertex(0);
      Point_3 v1 = query.vertex(1);
      Point_3 v2 = query.vertex(2);
      // Compute the centroid of the triangle
      Point_3 tri_center = centroid( query );

      // Compute the distance of the center to the closest point in tm2
      double dist = approximate_sqrt(squared_distance(bb_center, tri_center));

      // Compute the radius of the circumsphere of the bounding boxes
      double bb_radius = approximate_sqrt(squared_distance(
        Point_3(bbox.min(0),bbox.min(1),bbox.min(2)),
        Point_3(bbox.max(0),bbox.max(1),bbox.max(2)))
      )/2.;

      // Compute the radius of the circumcircle of the triangle
      double tri_radius = approximate_sqrt( std::max(
        std::max(
          squared_distance(tri_center, v0),
          squared_distance(tri_center, v1)
        ),
        squared_distance(tri_center, v2)
      ));

      // If triangles in the node can still lower the upper bound, enter it
      if (tri_radius + bb_radius >= dist) {
        return true;
      } else {
        return false;
      }
    }

    // Return the local Hausdorff bounds computed for the passed query triangle
    Hausdorff_bounds get_local_bounds()
    {
      return Hausdorff_bounds( h_local_lower, h_local_upper );
    }

  private:
    const AABBTraits& m_traits;
    // the mesh of this tree
    const TriangleMesh& m_tm2;
    // its vertex point map
    const VPM2& m_vpm2;
    // Local Hausdorff bounds for the query triangle
    double h_local_upper;
    double h_local_lower;
    double h_v0_lower;
    double h_v1_lower;
    double h_v2_lower;
  };


  /**
   * @class Hausdorff_primitive_traits_tm1
   */
  template<typename AABBTraits, typename Query, typename Kernel, typename TriangleMesh, typename VPM1, typename VPM2>
  class Hausdorff_primitive_traits_tm1
  {
    typedef AABB_face_graph_triangle_primitive<TriangleMesh, VPM2> TM2_primitive;
    typedef typename AABB_tree< AABB_traits<Kernel, TM2_primitive> >::AABB_traits Tree_traits;
    typedef typename AABBTraits::Primitive Primitive;
    typedef typename AABBTraits::Bounding_box Bounding_box;
    typedef ::CGAL::AABB_node<AABBTraits> Node;
    typedef typename Kernel::Point_3 Point_3;
    typedef typename Kernel::Triangle_3 Triangle_3;
    typedef std::pair<Triangle_3, Hausdorff_bounds> Candidate_triangle;
    typedef typename std::vector<Candidate_triangle> Candidate_set;
    typedef AABB_tree< AABB_traits<Kernel, TM2_primitive> > TM2_tree;
    typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;

  public:
    Hausdorff_primitive_traits_tm1(const AABBTraits& traits, const TM2_tree& tree, const TriangleMesh& tm1, const TriangleMesh& tm2 , const VPM1& vpm1, const VPM2& vpm2 )
      : m_traits(traits), m_tm2_tree(tree), m_tm1(tm1), m_tm2(tm2), m_vpm1(vpm1), m_vpm2(vpm2) {
        // Initialize the global bounds with 0., they will only grow.
        h_lower = 0.;
        h_upper = 0.;
      }

    // Explore the whole tree, i.e. always enter children if the methods
    // do_intersect() below determines that it is worthwhile.
    bool go_further() const { return true; }

    // Compute the explicit Hausdorff distance to the given primitive
    void intersection(const Query& query, const Primitive& primitive)
    {
      /* Have reached a single triangle, process it */

      // Map to transform faces of TM1 to actual triangles
      Triangle_from_face_descriptor_map<TriangleMesh, VPM1> m_face_to_triangle_map( &m_tm1, m_vpm1 );
      Triangle_3 candidate_triangle = get(m_face_to_triangle_map, primitive.id());

      // Call Culling on B with the single triangle found.
      Hausdorff_primitive_traits_tm2<
        Tree_traits, Triangle_3, Kernel, TriangleMesh, VPM2
      > traversal_traits_tm2(
        m_tm2_tree.traits(), m_tm2, m_vpm2,
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity()
      );
      m_tm2_tree.traversal(candidate_triangle, traversal_traits_tm2);

      // Update global Hausdorff bounds according to the obtained local bounds
      Hausdorff_bounds local_bounds = traversal_traits_tm2.get_local_bounds();
      if (local_bounds.first > h_lower) {
        h_lower = local_bounds.first;
      }
      if (local_bounds.second > h_upper) {
        h_upper = local_bounds.second;
      }

      // Store the triangle given as primitive here as candidate triangle
      // together with the local bounds it obtained to sind it to subdivision
      // later
      m_candidiate_triangles.push_back(Candidate_triangle(candidate_triangle, local_bounds));
    }

    // Determine whether child nodes will still contribute to a larger
    // Hausdorff distance and thus have to be entered
    // TODO Document Query object, explain why I don't need it here.
    bool do_intersect(const Query& /*query*/, const Node& node) const
    {
      /* Have reached a node, determine whether or not to enter it */

      // TODO What's the closest distance of TM2 to the box given by node?
      //      Can we have a sharper bound on this than the one implemented below?

      // Get the bounding box of the nodes
      Bounding_box bbox = node.bbox();
      // Compute its center
      Point_3 center = Point_3(
        (bbox.min(0) + bbox.max(0)) / 2,
        (bbox.min(1) + bbox.max(1)) / 2,
        (bbox.min(2) + bbox.max(2)) / 2);
      // Find the point from TM2 closest to the center
      // TODO Insert a hint here to accelerate the query
      Point_3 closest = m_tm2_tree.closest_point(center);
      // Compute the distance of the center to the closest point in tm2
      double dist = approximate_sqrt(squared_distance(center, closest));
      // Compute the radius of the circumsphere of the bounding boxes
      double radius = approximate_sqrt(squared_distance(
        Point_3(bbox.min(0),bbox.min(1),bbox.min(2)),
        Point_3(bbox.max(0),bbox.max(1),bbox.max(2)))
      )/2.;
      // If the distance is larger than the global lower bound, enter the node, i.e. return true.
      if (dist + radius > h_lower) {
          return true;
      } else {
        return false;
      }
    }

    // Return those triangles from TM1 which are candidates for including a
    // point realizing the Hausdorff distance
    Candidate_set get_candidate_triangles() {
      return m_candidiate_triangles;
    }

    // Return the local Hausdorff bounds computed for the passed query triangle
    Hausdorff_bounds get_global_bounds()
    {
      return Hausdorff_bounds( h_lower, h_upper );
    }

  private:
    const AABBTraits& m_traits;
    // The two triangle meshes
    const TriangleMesh& m_tm1;
    const TriangleMesh& m_tm2;
    // Their vertex-point-maps
    const VPM1 m_vpm1;
    const VPM2 m_vpm2;
    // AABB tree for the second triangle meshes
    const TM2_tree& m_tm2_tree;
    // Global Hausdorff bounds to be taken track of during the traversal
    double h_lower;
    double h_upper;
    // List of candidate triangles with their Hausdorff bounds attached
    Candidate_set m_candidiate_triangles;
  };
}

#endif //CGAL_PMP_INTERNAL_AABB_TRAVERSAL_TRAITS_WITH_HAUSDORFF_DISTANCE
