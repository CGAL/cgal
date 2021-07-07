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

#ifndef CGAL_BARYCENTRIC_MEAN_VALUE_COORDINATES_3_H
#define CGAL_BARYCENTRIC_MEAN_VALUE_COORDINATES_3_H

// #include <CGAL/license/Barycentric_coordinates_3.h>

// Internal includes.
#include <CGAL/Barycentric_coordinates_3/internal/utils_3.h>
#include <CGAL/Barycentric_coordinates_3/barycentric_enum_3.h>

namespace CGAL {
namespace Barycentric_coordinates {

  template<
  typename PolygonMesh,
  typename GeomTraits,
  typename VertexToPointMap = typename property_map_selector<PolygonMesh,
    CGAL::vertex_point_t>::const_type>
  class Mean_value_coordinates_3 {

  public:
    using Polygon_mesh = PolygonMesh;
    using Geom_Traits = GeomTraits;
    using Vertex_to_point_map = VertexToPointMap;

    using Construct_vec_3 = typename GeomTraits::Construct_vector_3;
    using Cross_3 = typename GeomTraits::Construct_cross_product_vector_3;
    using Dot_3 = typename GeomTraits::Compute_scalar_product_3;
    using Sqrt = typename internal::Get_sqrt<GeomTraits>::Sqrt;
    using Approximate_angle_3 = typename GeomTraits::Compute_approximate_angle_3;

	  typedef typename GeomTraits::FT FT;
    typedef typename GeomTraits::Point_3 Point_3;
    typedef typename GeomTraits::Vector_3 Vector_3;

    Mean_value_coordinates_3(
      const PolygonMesh& polygon_mesh,
      const Computation_policy_3 policy,
      const VertexToPointMap vertex_to_point_map,
      const GeomTraits traits = GeomTraits()) :
    m_polygon_mesh(polygon_mesh),
    m_computation_policy(policy),
    m_vertex_to_point_map(vertex_to_point_map),
    m_traits(traits),
    m_construct_vector_3(m_traits.construct_vector_3_object()),
    m_cross_3(m_traits.construct_cross_product_vector_3_object()),
    m_dot_3(m_traits.compute_scalar_product_3_object()),
    m_approximate_angle_3(m_traits.compute_approximate_angle_3_object()),
    sqrt(internal::Get_sqrt<GeomTraits>::sqrt_object(m_traits)){

      // Check if polyhedron is strongly convex
      CGAL_assertion(is_strongly_convex_3(m_polygon_mesh, m_traits));
      m_weights.resize(vertices(m_polygon_mesh).size());
      query_vertex_vectors.resize(3);
      unit_vectors.resize(3);
      m_vectors.resize(3);
      angles.resize(3);
    }

    Mean_value_coordinates_3(
      const PolygonMesh& polygon_mesh,
      const Computation_policy_3 policy =
      Computation_policy_3::DEFAULT,
      const GeomTraits traits = GeomTraits()) :
    Mean_value_coordinates_3(
      polygon_mesh,
      policy,
      get_const_property_map(CGAL::vertex_point, polygon_mesh),
      traits) { }

    template<typename OutIterator>
    OutIterator operator()(const Point_3& query, OutIterator c_begin) {
      return compute(query, c_begin);
    }

    VertexToPointMap get_vertex_to_point_map(){
      return m_vertex_to_point_map;
    }

  private:
    const PolygonMesh& m_polygon_mesh;
  	const Computation_policy_3 m_computation_policy;
  	const VertexToPointMap m_vertex_to_point_map; // use it to map vertex to Point_3
  	const GeomTraits m_traits;

    const Construct_vec_3 m_construct_vector_3;
    const Cross_3 m_cross_3;
    const Dot_3 m_dot_3;
    const Sqrt sqrt;
    const Approximate_angle_3 m_approximate_angle_3;

    // Store useful information
    std::vector<FT> m_weights;
    std::vector<Vector_3> query_vertex_vectors;
    std::vector<Vector_3> unit_vectors;
    std::vector<Vector_3> m_vectors;
    std::vector<FT> angles;

  	template<typename OutputIterator>
    OutputIterator compute(
      const Point_3& query, OutputIterator coordinates) {

      // Compute weights.
      const FT sum = compute_weights(query);
      CGAL_assertion(sum != FT(0));

      // The coordinates must be saved in the same order as vertices in the vertex range.
      const auto vd = vertices(m_polygon_mesh);
      CGAL_assertion(m_weights.size() == vd.size());

      for (std::size_t vi = 0; vi < vd.size(); vi++) {

        CGAL_assertion(vi < m_weights.size());
        const FT coordinate = m_weights[vi]/sum;
        *(coordinates++) = coordinate;
      }

      return coordinates;
    }

    FT compute_weights(const Point_3& query) {

      // Sum of weights to normalize them later.
      FT sum = FT(0);

	    // Vertex index.
      std::size_t vi = 0;
      const auto vd = vertices(m_polygon_mesh);
      for (const auto& vertex : vd) {

        // Call function to calculate coordinates
        const FT weight = compute_mv_vertex_query(vertex, query);

    	  CGAL_assertion(vi < m_weights.size());
    	  m_weights[vi] = weight;
    	  sum += weight;
    	  ++vi; // update vi
      }

      CGAL_assertion(sum != FT(0));
      return sum;
    }

    // Compute wp coordinates for a given vertex v and a query q
    template<typename Vertex>
    FT compute_mv_vertex_query(const Vertex& vertex, const Point_3& query){

      // Map vertex descriptor to point_3
      const Point_3& vertex_val = get(m_vertex_to_point_map, vertex);

      // Circulator of faces around the vertex
      CGAL::Face_around_target_circulator<Polygon_mesh>
      face_circulator(m_polygon_mesh.halfedge(vertex), m_polygon_mesh);

      CGAL::Face_around_target_circulator<Polygon_mesh>
      face_done(face_circulator);

      // Compute weight w_v
      FT weight = FT(0);

      // Iterate using the circulator
      do{

        // Vertices around face iterator
        const auto hedge = halfedge(*face_circulator, m_polygon_mesh);
        const auto vertices = vertices_around_face(hedge, m_polygon_mesh);
        auto vertex_itr = vertices.begin();
        CGAL_precondition(vertices.size() == 3);

        // Weight of vertex for this particular face
        FT partial_weight = FT(0);
        int vertex_idx = -1;

        for(std::size_t i = 0; i < 3; i++){

          if(*vertex_itr == vertex)
            vertex_idx = i;

          const Vector_3 p = m_construct_vector_3(query, get(m_vertex_to_point_map, *vertex_itr));
          query_vertex_vectors[i] = p;
          vertex_itr++;
        }

        // Current vertex should be present in face
        assert(vertex_idx != -1);

        for(std::size_t i = 0; i < 3; i++){

          assert(query_vertex_vectors[i].squared_length() > 0);
          unit_vectors[i] = query_vertex_vectors[i]/sqrt(query_vertex_vectors[i].squared_length());

          m_vectors[i] = m_cross_3(query_vertex_vectors[i], query_vertex_vectors[(i+1)%3]);
          assert(m_vectors[i].squared_length() > 0);
          m_vectors[i] /= sqrt(m_vectors[i].squared_length());

          angles[i] = m_approximate_angle_3(query_vertex_vectors[i], query_vertex_vectors[(i+1)%3]);
        }

        partial_weight += angles[0] * m_dot_3(m_vectors[0], m_vectors[(vertex_idx+1)%3]);
        partial_weight += angles[1] * m_dot_3(m_vectors[1], m_vectors[(vertex_idx+1)%3]);
        partial_weight += angles[2] * m_dot_3(m_vectors[2], m_vectors[(vertex_idx+1)%3]);

        const FT dot_unit_m = m_dot_3(unit_vectors[vertex_idx], m_vectors[(vertex_idx+1)%3]);
        assert(dot_unit_m != 0);
        partial_weight /= dot_unit_m;

        weight += partial_weight;

        face_circulator++;

      }while(face_circulator!=face_done);

      const FT vertex_query_squared_len = m_construct_vector_3(vertex_val, query).squared_length();
      assert(vertex_query_squared_len != 0);

      return weight/sqrt(vertex_query_squared_len);
    }
  };

  template<
  typename Point_3,
  typename Mesh,
  typename OutIterator>
  OutIterator mean_value_coordinates_3(
    const Mesh& surface_mesh,
    const Point_3& query,
    OutIterator c_begin,
    const Computation_policy_3 policy =
    Computation_policy_3::DEFAULT) {

    using Geom_Traits = typename Kernel_traits<Point_3>::Kernel;

    Mean_value_coordinates_3<Mesh, Geom_Traits> mean_value(surface_mesh, policy);
    return mean_value(query, c_begin);
  }

} // namespace Barycentric_coordinates
} // namespace CGAL

#endif // CGAL_BARYCENTRIC_MEAN_VALUE_COORDINATES_3_H
