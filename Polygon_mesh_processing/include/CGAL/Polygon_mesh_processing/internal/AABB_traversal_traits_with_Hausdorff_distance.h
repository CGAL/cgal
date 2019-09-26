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
   * @struct Candidate_triangle
   */
  template<typename Kernel>
  struct Candidate_triangle {
    typedef typename Kernel::Triangle_3 Triangle_3;

    Candidate_triangle(const Triangle_3& triangle, const Hausdorff_bounds& bounds)
      : m_triangle(triangle), m_bounds(bounds) {}

    Triangle_3 m_triangle;
    Hausdorff_bounds m_bounds;

    #if BOOST_VERSION >= 105000
        bool operator<(const Candidate_triangle& other) const { return m_bounds.second < other.m_bounds.second; }
        bool operator>(const Candidate_triangle& other) const { return m_bounds.second > other.m_bounds.second; }
    #else
        bool operator>(const Candidate_triangle& other) const { return m_bounds.second < other.m_bounds.second; }
        bool operator<(const Candidate_triangle& other) const { return m_bounds.second > other.m_bounds.second; }
    #endif
  };

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
    typedef typename Kernel::Vector_3 Vector_3;

  public:
    typedef int Priority;
    
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
      // Get the bounding box of the nodes
      Bounding_box bbox = node.bbox();
      // Get the vertices of the query triangle
      Point_3 v0 = query.vertex(0);
      Point_3 v1 = query.vertex(1);
      Point_3 v2 = query.vertex(2);
      // Find the axis aligned bbox of the triangle
      Point_3 tri_min = Point_3 (
        std::min(std::min( v0.x(), v1.x()), v2.x() ),
        std::min(std::min( v0.y(), v1.y()), v2.y() ),
        std::min(std::min( v0.z(), v1.z()), v2.z() )
      );
      Point_3 tri_max = Point_3 (
        std::max(std::max( v0.x(), v1.x()), v2.x() ),
        std::max(std::max( v0.y(), v1.y()), v2.y() ),
        std::max(std::max( v0.z(), v1.z()), v2.z() )
      );

      // Compute distance of the bounding boxes
      // Distance along the x-axis
      double dist_x = 0.;
      if ( tri_max.x() < bbox.min(0) ) {
        dist_x = bbox.min(0) - tri_max.x();
      } else if ( bbox.max(0) < tri_min.x() ) {
        dist_x = tri_min.x() - bbox.max(0);
      }
      // Distance along the y-axis
      double dist_y = 0.;
      if ( tri_max.y() < bbox.min(1) ) {
        dist_y = bbox.min(1) - tri_max.y();
      } else if ( bbox.max(1) < tri_min.y() ) {
        dist_y = tri_min.y() - bbox.max(1);
      }
      // Distance along the y-axis
      double dist_z = 0.;
      if ( tri_max.z() < bbox.min(2) ) {
        dist_z = bbox.min(2) - tri_max.z();
      } else if ( bbox.max(2) < tri_min.z() ) {
        dist_z = tri_min.z() - bbox.max(2);
      }

      // Lower bound on the distance between the two bounding boxes is given
      // as the length of the diagonal of the bounding box between them
      double dist = approximate_sqrt( Vector_3(dist_x,dist_y,dist_z).squared_length() );

      // Check whether investigating the bbox can still lower the Hausdorff
      // distance. If so, enter the box.
      if ( dist <= h_local_lower) {
        return true;
      } else {
        return false;
      }
    }

    std::pair<bool,Priority> do_intersect_with_priority(const Query& query, const Node& node) const
    {
      bool b = do_intersect(query, node);
      return std::make_pair(b, Priority(0));
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
    typedef typename Kernel::Vector_3 Vector_3;
    typedef typename Kernel::Triangle_3 Triangle_3;
    typedef typename std::vector<Candidate_triangle<Kernel>> Candidate_set;
    typedef AABB_tree< AABB_traits<Kernel, TM2_primitive> > TM2_tree;
    typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;
    typedef
    #if BOOST_VERSION >= 105000
          boost::heap::priority_queue< Candidate_triangle<Kernel>, boost::heap::compare< std::greater<Candidate_triangle<Kernel> > > >
    #else
          std::priority_queue< Candidate_triangle<Kernel> >
    #endif
          Heap_type;

  public:
    typedef int Priority;
    Hausdorff_primitive_traits_tm1(
      const AABBTraits& traits, const TM2_tree& tree, const TriangleMesh& tm1,
      const TriangleMesh& tm2 , const VPM1& vpm1, const VPM2& vpm2,
      const Point_3 hint )
      : m_traits(traits), m_tm2_tree(tree), m_tm1(tm1), m_tm2(tm2),
      m_vpm1(vpm1), m_vpm2(vpm2) {
        // Initialize the global bounds with 0., they will only grow.
        h_lower = 0.;
        h_upper = 0.;
        // Initialize number of culled triangles
        m_investigated_on_tm1 = 0;
      }

    // Explore the whole tree, i.e. always enter children if the methods
    // do_intersect() below determines that it is worthwhile.
    bool go_further() const { return true; }

    // Compute the explicit Hausdorff distance to the given primitive
    void intersection(const Query& query, const Primitive& primitive)
    {
      /* Have reached a single triangle, process it */
      m_investigated_on_tm1++;

      // Map to transform faces of TM1 to actual triangles
      Triangle_from_face_descriptor_map<TriangleMesh, VPM1> m_face_to_triangle_map( &m_tm1, m_vpm1 );
      Triangle_3 candidate_triangle = get(m_face_to_triangle_map, primitive.id());

      // TODO Can we initialize the bounds here, s.t. we don't start with infty?
      // Can we initialize the bounds depending on the closest points in tm2
      // seen from the three vertices?

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
      m_candidiate_triangles.push( Candidate_triangle<Kernel>(candidate_triangle, local_bounds) );
      pq.push( Candidate_triangle<Kernel>(candidate_triangle, local_bounds) );
    }

    // Determine whether child nodes will still contribute to a larger
    // Hausdorff distance and thus have to be entered
    // TODO Document Query object, explain why I don't need it here.
    bool do_intersect(const Query& /*query*/, const Node& node) const
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
      Point_3 closest = m_tm2_tree.closest_point(center);
      // Compute the difference vector between the bbox center and the closest
      // point in tm2
      Vector_3 difference = Vector_3( closest, center );
      // Shift the vector to be the difference between the farthest corner
      // of the bounding box away from the closest point on TM2
      double diff_x = (bbox.max(0) - bbox.min(0)) / 2;
      if (difference.x() < 0) diff_x = diff_x * -1.;
      double diff_y = (bbox.max(1) - bbox.min(1)) / 2;
      if (difference.y() < 0) diff_y = diff_y * -1.;
      double diff_z = (bbox.max(2) - bbox.min(2)) / 2;
      if (difference.z() < 0) diff_z = diff_z * -1.;
      difference = difference + Vector_3( diff_x, diff_y, diff_z );
      // Compute distance from the farthest corner of the bbox to the closest
      // point in TM2
      double dist = approximate_sqrt( difference.squared_length() );

      // If the distance is larger than the global lower bound, enter the node, i.e. return true.
      if (dist > h_lower) {
          return true;
      } else {
        return false;
      }
    }

    std::pair<bool,Priority> do_intersect_with_priority(const Query& query, const Node& node) const
    {
      bool b = do_intersect(query, node);
      return std::make_pair(b, Priority(0));
    }
          
    // Return those triangles from TM1 which are candidates for including a
    // point realizing the Hausdorff distance
    Heap_type get_candidate_triangles() {
      return m_candidiate_triangles;
    }

    // Return the local Hausdorff bounds computed for the passed query triangle
    Hausdorff_bounds get_global_bounds()
    {
      return Hausdorff_bounds( h_lower, h_upper );
    }

    int get_num_culled_triangles()
    {
      return (m_tm1.num_faces() - m_investigated_on_tm1);
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
    Heap_type m_candidiate_triangles;
    // Heap of candidate triangles with their Hausdorff bounds attached
    Heap_type pq;
    // Number of triangles investigated in the procedure
    int m_investigated_on_tm1;
  };
}

#endif //CGAL_PMP_INTERNAL_AABB_TRAVERSAL_TRAITS_WITH_HAUSDORFF_DISTANCE
