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
  template<typename AABBTraits, typename Query, typename Kernel, typename TriangleMesh, typename VPM1, typename VPM2>
  class Hausdorff_primitive_traits_tm2
  {
    typedef typename AABBTraits::Primitive Primitive;
    typedef ::CGAL::AABB_node<AABBTraits> Node;
    typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;
    typename Kernel::Construct_projected_point_3 project_point;

  public:
    Hausdorff_primitive_traits_tm2(const AABBTraits& traits, const TriangleMesh& tm1, const TriangleMesh& tm2, const VPM1& vpm1, const VPM2& vpm2)
      : m_traits(traits), m_tm1(tm1), m_tm2(tm2), m_vpm1(vpm1), m_vpm2(vpm2) {
        // Initialize the global bounds with 0., they will only grow.
        h_local_upper = std::numeric_limits<double>::infinity();
        h_local_lower_0 = std::numeric_limits<double>::infinity();
        h_local_lower_1 = std::numeric_limits<double>::infinity();
        h_local_lower_2 = std::numeric_limits<double>::infinity();
      }

    // Explore the whole tree, i.e. always enter children if the methods
    // do_intersect() below determines that it is worthwhile.
    bool go_further() const { return true; }

    // Compute the explicit Hausdorff distance to the given primitive
    void intersection(const Query& query, const Primitive& primitive)
    {
      // Have reached a single triangle
      std::cout << "Reached Triangle in TM2:" << primitive.id() << '\n';

      /*
      / Determine the distance accroding to
      /       min_{b \in primitive} ( max_{vertex in query} ( d(vertex, b)))
      /
      / Here, we only have one triangle in B, i.e. tm2. Thus, it suffices to
      / compute the distance of the vertices of the query triangles to the
      / primitive triangle and use the maximum of the obtained distances.
      */

      // The query object is a triangle from TM1, get its vertices
      halfedge_descriptor hd = halfedge(query.id(), m_tm1);
      vertex_descriptor v0 = source(hd,m_tm1);
      vertex_descriptor v1 = target(hd,m_tm1);
      vertex_descriptor v2 = target(next(hd, m_tm1), m_tm1);

      // Compute distances of the vertices to the primitive triangle in TM2
      Triangle_from_face_descriptor_map<TriangleMesh, VPM2> face_to_triangle_map(&m_tm2, m_vpm2);
      double v0_dist = squared_distance(
        project_point(get(face_to_triangle_map, primitive.id()), get(m_vpm1, v0)),
        get(m_vpm1, v0)
      );
      double v1_dist = squared_distance(
        project_point(get(face_to_triangle_map, primitive.id()), get(m_vpm1, v1)),
        get(m_vpm1, v1)
      );
      double v2_dist = squared_distance(
        project_point(get(face_to_triangle_map, primitive.id()), get(m_vpm1, v2)),
        get(m_vpm1, v2)
      );

      // Get the distance as maximizers over all vertices
      double distance = approximate_sqrt(std::max(std::max(v0_dist,v1_dist),v2_dist));

      // Update local upper bound
      if ( distance < h_local_upper ) h_local_upper = distance;

      double vertex_distance = 0.;
      /*
      /  TODO For all vertices v of the query triangle, determine the distance by
      /       min_{b \in primitive} ( d(v, b) )
      /       Update h_local_lower_i accordingly
      */

      h_local_lower = std::max( std::max (h_local_lower_0, h_local_lower_1), h_local_lower_2 );
    }

    // Determine whether child nodes will still contribute to a smaller
    // Hausdorff distance and thus have to be entered
    bool do_intersect(const Query& query, const Node& node) const
    {
      // Have reached a node, determine whether or not to enter it
      double distance = 0.;

      /* TODO Determine distance of the node's bounding box (node.bbox()) to
      /       the triangle given by the query object.
      */

      if (distance <= h_local_upper) return true;
      else return false;
    }

    Hausdorff_bounds get_local_bounds()
    {
      return Hausdorff_bounds( h_local_lower, h_local_upper );
    }

  private:
    const AABBTraits& m_traits;
    // The two triangle meshes
    const TriangleMesh& m_tm1;
    const TriangleMesh& m_tm2;
    // Their respective vertex-point Maps
    const VPM1& m_vpm1;
    const VPM2& m_vpm2;
    // Local Hausdorff bounds for the query triangle
    double h_local_upper;
    double h_local_lower;
    // Local Hausdorff bounds for the query triangle's vertices
    double h_local_lower_0;
    double h_local_lower_1;
    double h_local_lower_2;
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
    typedef AABB_tree< AABB_traits<Kernel, TM2_primitive> > TM2_tree;
    typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;

  public:
    Hausdorff_primitive_traits_tm1(const AABBTraits& traits, const TM2_tree& tree, const TriangleMesh& tm1, const TriangleMesh& tm2 , const VPM1& vpm1, const VPM2& vpm2 )
      : m_traits(traits), tm2_tree(tree), m_tm1(tm1), m_tm2(tm2), m_vpm1(vpm1), m_vpm2(vpm2) {
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
      // Have reached a single triangle
      std::cout << "Reached Triangle in TM1: " << primitive.id() << '\n';

      // Call Culling on B with the single triangle found.
      Hausdorff_primitive_traits_tm2<Tree_traits, Primitive, Kernel, TriangleMesh, VPM1, VPM2> traversal_traits_tm2( tm2_tree.traits(), m_tm1, m_tm2, m_vpm1, m_vpm2 );
      tm2_tree.traversal(primitive, traversal_traits_tm2);

      // Update global Hausdorff bounds according to the obtained local bounds
      Hausdorff_bounds local_bounds = traversal_traits_tm2.get_local_bounds();
      if (local_bounds.first > h_lower) {
        h_lower = local_bounds.first;
      }
      if (local_bounds.second > h_upper) {
        h_upper = local_bounds.second;
      }
    }

    // Determine whether child nodes will still contribute to a larger
    // Hausdorff distance and thus have to be entered
    bool do_intersect(const Query& query, const Node& node) const
    {
      /* Have reached a node, determine whether or not to enter it */

      // Get the bounding box of the nodes
      Bounding_box bbox = node.bbox();
      // Compute its center
      Point_3 center = Point_3(
        (bbox.min(0) + bbox.max(0)) / 2,
        (bbox.min(1) + bbox.max(1)) / 2,
        (bbox.min(2) + bbox.max(2)) / 2);
      // Find the point from TM2 closest to the center
      Point_3 closest = tm2_tree.closest_point(center);
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

  private:
    const AABBTraits& m_traits;
    // The two triangle meshes
    const TriangleMesh& m_tm1;
    const TriangleMesh& m_tm2;
    // Their vertex-point-maps
    const VPM1 m_vpm1;
    const VPM2 m_vpm2;
    // AABB tree for the second triangle meshes
    const TM2_tree& tm2_tree;
    // Global Hausdorff bounds to be taken track of during the traversal
    double h_lower;
    double h_upper;
  };
}

#endif //CGAL_PMP_INTERNAL_AABB_TRAVERSAL_TRAITS_WITH_HAUSDORFF_DISTANCE
