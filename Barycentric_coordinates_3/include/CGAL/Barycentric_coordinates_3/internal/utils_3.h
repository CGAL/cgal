// Copyright (c) 2021 GeometryFactory SARL (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Antonio Gomes, Dmitry Anisimov
//

#ifndef CGAL_BARYCENTRIC_INTERNAL_UTILS_3_H
#define CGAL_BARYCENTRIC_INTERNAL_UTILS_3_H

#include <CGAL/license/Barycentric_coordinates_3.h>

// STL includes
#include <tuple>

// Internal includes
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/convexity_check_3.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/Barycentric_coordinates_2/triangle_coordinates_2.h>
#include <CGAL/Barycentric_coordinates_2/Wachspress_coordinates_2.h>

namespace CGAL {
namespace Barycentric_coordinates {
namespace internal {

enum class Edge_case {

  EXTERIOR = 0, // exterior part of the polyhedron
  INTERIOR = 1,  // interior part of the polyhedron
  BOUNDARY = 2, // boundary part of the polyhedron
  EXTERIOR_BOUNDARY = 3, // extension of the boundary
};

  template<typename FT>
  FT get_tolerance() {
    return 0.000000000001;
  }

// Compute barycentric coordinates in the space.
  template<
  typename OutputIterator,
  typename GeomTraits>
  OutputIterator tetrahedron_coordinates_impl(
    const typename GeomTraits::Point_3& p0,
    const typename GeomTraits::Point_3& p1,
    const typename GeomTraits::Point_3& p2,
    const typename GeomTraits::Point_3& p3,
    const typename GeomTraits::Point_3& query,
    OutputIterator coordinates,
    GeomTraits& traits) {

    // Number type.
    using FT = typename GeomTraits::FT;

    // Functions.
    auto volume_3 = traits.compute_volume_3_object();
    const FT total_volume = volume_3(p0, p1, p2, p3);

    CGAL_precondition(total_volume != FT(0));
    if (total_volume == FT(0)) {
      get_default(4, coordinates);
      return coordinates;
    }

    // Compute some related sub-volumes.
    const FT V1 = volume_3(p1, p3, p2, query);
    const FT V2 = volume_3(p2, p3, p0, query);
    const FT V3 = volume_3(p3, p1, p0, query);

    // Compute the inverted total volume of the tetrahedron.
    CGAL_assertion(total_volume != FT(0));
    const FT inverted_total_volume = FT(1) / total_volume;

    // Compute coordinates.
    const FT b0 = V1 * inverted_total_volume;
    const FT b1 = V2 * inverted_total_volume;
    const FT b2 = V3 * inverted_total_volume;
    const FT b3 = FT(1) - b0 - b1 - b2;

    // Return coordinates.
    *(coordinates++) = b0;
    *(coordinates++) = b1;
    *(coordinates++) = b2;
    *(coordinates++) = b3;

    return coordinates;
  }

  // Compute normal vector of the face (not normalized).
  template<
    typename Face,
    typename VertexPointMap,
    typename TriangleMesh,
    typename GeomTraits>
  typename GeomTraits::Vector_3 get_face_normal(
    const Face& face,
    const VertexPointMap& vpm,
    const TriangleMesh& tmesh,
    GeomTraits& traits){

    using Vector_3 = typename GeomTraits::Vector_3;
    auto cross_3 = traits.construct_cross_product_vector_3_object();
    auto vector_3 = traits.construct_vector_3_object();

    const auto hedge = halfedge(face, tmesh);
    const auto vertices = vertices_around_face(hedge, tmesh);
    CGAL_precondition(vertices.size() >= 3);

    auto vertex = vertices.begin();
    auto v0 = *vertex; ++vertex;
    auto v1 = *vertex; ++vertex;
    auto v2 = *vertex;

    const Vector_3 u = vector_3(get(vpm, v0), get(vpm, v1));
    const Vector_3 v = vector_3(get(vpm, v0), get(vpm, v2));
    Vector_3 face_normal = cross_3(u, v);

    return face_normal;
  }

  //Compute cotangent of dihedral angle between two faces
  template<typename GeomTraits>
  typename GeomTraits::FT cot_dihedral_angle(
    const typename GeomTraits::Vector_3& vec_1,
    const typename GeomTraits::Vector_3& vec_2,
    GeomTraits& traits){

    using FT = typename GeomTraits::FT;
    auto dot_3 = traits.compute_scalar_product_3_object();
    auto cross_3 = traits.construct_cross_product_vector_3_object();
    auto squared_length_3 = traits.compute_squared_length_3_object();
    auto sqrt(Get_sqrt<GeomTraits>::sqrt_object(traits));

    assert(squared_length_3(vec_1) != FT(0));
    assert(squared_length_3(vec_2) != FT(0));

    const FT approximate_dot_3 = dot_3(vec_1, vec_2);

    const FT approximate_cross_3_length = sqrt(squared_length_3(cross_3(vec_1, vec_2)));

    assert(approximate_cross_3_length != FT(0));

    return approximate_dot_3/approximate_cross_3_length;
  }

  template<typename VertexRange>
  inline bool are_vertices_distinct(VertexRange& vertices_face){

    // Check if the vertices are distinct
    for(auto itr1 = vertices_face.begin(); itr1 != vertices_face.end(); itr1++)
      for(auto itr2 = std::next(itr1, 1); itr2 != vertices_face.end(); itr2++)
        if(*itr1 == *itr2) return false;

    return true;
  }

  template<
    typename VertexRange,
    typename VertexPointMap,
    typename TriangleMesh,
    typename OutputIterator,
    typename GeomTraits>
  OutputIterator boundary_coordinates_3(
    VertexRange& vertices_face,
    const VertexPointMap& vertex_point_map,
    const TriangleMesh& tmesh,
    const typename GeomTraits::Point_3& query,
    OutputIterator coordinates,
    GeomTraits& traits,
    bool use_wp_flag){

    CGAL_USE(traits);
    using FT = typename GeomTraits::FT;
    using Plane_3 = typename GeomTraits::Plane_3;
    using Point_2 = typename GeomTraits::Point_2;

    const FT tol = get_tolerance<FT>();
    const std::size_t num_sides_face = vertices_face.size();

    // Check if the vertices are distinct
    CGAL_assertion(are_vertices_distinct(vertices_face));

    // Create plane
    auto vertex_itr = vertices_face.begin();
    const auto v0 = *vertex_itr; vertex_itr++;
    const auto v1 = *vertex_itr; vertex_itr++;
    const auto v2 = *vertex_itr;
    typename GeomTraits::Construct_plane_3 construct_plane_3 =
      traits.construct_plane_3_object();
    const Plane_3 face_plane = construct_plane_3(get(vertex_point_map, v0),
     get(vertex_point_map, v1), get(vertex_point_map, v2));

    typename GeomTraits::Construct_projected_xy_point_2 project_onto_plane =
      traits.construct_projected_xy_point_2_object();

    // Store 2d vertices
    std::vector<Point_2> polygon;
    polygon.reserve(num_sides_face);
    auto polygon_itr = std::back_inserter(polygon);
    Point_2 query_2 = project_onto_plane(face_plane, query);
    for(auto v : vertices_face){

      *polygon_itr = project_onto_plane(face_plane, get(vertex_point_map, v));
      polygon_itr++;
    }

    // Check convexity
    if constexpr (internal_kernel_traits::Has_nested_R<GeomTraits>::value)
      CGAL_assertion(is_convex_2(polygon.begin(), polygon.end(), traits));

    // Store 2d barycentric coordinates
    std::vector<FT> bar_coords_2;
    bar_coords_2.reserve(num_sides_face);

    // Use wp_2 or triangle coordinates
    if(use_wp_flag){
      wachspress_coordinates_2(polygon, query_2, std::back_inserter(bar_coords_2), traits);
    }
    else{
      CGAL_assertion(polygon.size() == 3);
      triangle_coordinates_2(polygon[0], polygon[1], polygon[2], query_2, std::back_inserter(bar_coords_2), traits);
    }

    // Fill coordinates
    CGAL_assertion(bar_coords_2.size() == num_sides_face);
    for(auto vertex_polyhedron : vertices(tmesh)){

      bool found_vertex = false;
      auto bar_coords_itr = bar_coords_2.begin();
      for(auto vertex_face : vertices_face){

        if(vertex_polyhedron == vertex_face){

          *coordinates = CGAL::abs(*bar_coords_itr) < tol? FT(0): *bar_coords_itr;
          found_vertex = true;
          break;
        }
        bar_coords_itr++;
      }

      if(!found_vertex)
        *coordinates = FT(0);

      coordinates++;
    }

    return coordinates;
  }

  // Determine if the query point is on the interior, exterior or boundary
  template<
    typename VertexPointMap,
    typename PolygonMesh,
    typename OutputIterator,
    typename GeomTraits>
  Edge_case locate_wrt_polyhedron(
    const VertexPointMap& vertex_point_map,
    const PolygonMesh& polygon_mesh,
    const typename GeomTraits::Point_3& query,
    OutputIterator coordinates,
    const GeomTraits& traits,
    const bool use_wp_flag = false) {

    using Vector_3 = typename GeomTraits::Vector_3;
    using FT = typename GeomTraits::FT;
    typename GeomTraits::Compute_scalar_product_3 dot_3 = traits.compute_scalar_product_3_object();
    typename GeomTraits::Construct_vector_3 construct_vector_3 = traits.construct_vector_3_object();
    typename GeomTraits::Compute_squared_length_3 squared_length_3 = traits.compute_squared_length_3_object();
    typename GeomTraits::Construct_divided_vector_3 construct_divided_vector_3 = traits.construct_divided_vector_3_object();
    const auto sqrt(Get_sqrt<GeomTraits>::sqrt_object(traits));

    // Flags that indicates position of the query point
    bool exterior_flag = false;
    bool boundary_flag = false;

    const FT tol = get_tolerance<FT>();
    auto face_range = faces(polygon_mesh);

    for(auto face : face_range){

      const auto hedge = halfedge(face, polygon_mesh);
      const auto vertices_face = vertices_around_face(hedge, polygon_mesh);
      CGAL_precondition(vertices_face.size() >= 3);

      auto vertex = vertices_face.begin();
      const auto vertex_val = get(vertex_point_map, *vertex);

      // Vector connecting query point to vertex;
      const Vector_3 query_vertex = construct_vector_3(query, vertex_val);

      // Calculate normals of faces
      Vector_3 face_normal_i = get_face_normal(
        face, vertex_point_map, polygon_mesh, traits);
      face_normal_i = construct_divided_vector_3(face_normal_i, sqrt(squared_length_3(face_normal_i)));

      // Distance of query to face
      const FT perp_dist_i = dot_3(query_vertex, face_normal_i);

      // Verify location of query point;
      if(CGAL::abs(perp_dist_i) < tol) {
        if(!boundary_flag)
          boundary_coordinates_3(vertices_face, vertex_point_map, polygon_mesh,
           query, coordinates, traits, use_wp_flag);
        boundary_flag = true;
      }
      else if(perp_dist_i < 0)
        exterior_flag = true;
    }

    // Choose location
    if(boundary_flag && exterior_flag)
      return Edge_case::EXTERIOR_BOUNDARY;
    else if(boundary_flag)
      return Edge_case::BOUNDARY;
    else if(exterior_flag)
      return Edge_case::EXTERIOR;
    else
      return Edge_case::INTERIOR;
  }

} // namespace internal
} // namespace Barycentric_coordinates
} // namespace CGAL

#endif
