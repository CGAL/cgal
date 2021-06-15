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

#ifndef CGAL_BARYCENTRIC_WACHSPRESS_COORDINATES_3_H
#define CGAL_BARYCENTRIC_WACHSPRESS_COORDINATES_3_H

// Internal includes.
#include <CGAL/Barycentric_coordinates_3/internal/utils_3.h>

namespace CGAL {
namespace Barycentric_coordinates {

  template<
  typename PolygonMesh,
  typename GeomTraits,
  typename VertexToPointMap = typename property_map_selector<PolygonMesh, CGAL::vertex_point_t>::const_type>
  class Wachspress_coordinates_3 {

  public:
    using Polygon_mesh = PolygonMesh;
    using Geom_traits = GeomTraits;
    using Vertex_to_point_map = VertexToPointMap;

    using Dot_3 = typename GeomTraits::Compute_scalar_product_3;
    using Det_3 = typename GeomTraits::Compute_determinant_3;
    using Cross_3 = typename GeomTraits::Construct_cross_product_vector_3;

	typedef typename GeomTraits::FT FT;
    typedef typename GeomTraits::Point_3 Point_3;
    typedef typename GeomTraits::Vector_3 Vector_3;

  public:
    Wachspress_coordinates_3(
    	const PolygonMesh& polygon_mesh,
    	const Computation_policy_3 policy,
    	const VertexToPointMap vertex_to_point_map,
    	const GeomTraits traits = GeomTraits()) :
    m_polygon_mesh(polygon_mesh),
    m_computation_policy(policy),
    m_vertex_to_point_map(vertex_to_point_map),
    m_traits(traits),
    m_dot3(m_traits.compute_scalar_product_3_object()),
    m_det3(m_traits.compute_determinant_3_object()),
    m_cross3(m_traits.construct_cross_product_vector_3_object()){

        //put here conditions to ensure that the polyhedron is convex and combinatorial consistent
      // REVIEW: what about this: https://doc.cgal.org/latest/Convex_hull_3/group__PkgConvexityChecking.html#ga1e42ab960e30ee450cdebf0b5e58f9e3 ?
      // Also: faces(m_polygon_mesh).size() >= 4 ?
    	CGAL_precondition(vertices(m_polygon_mesh).size()>=4);
    	m_weights.resize(vertices(m_polygon_mesh).size());
    }

    Wachspress_coordinates_3(
    	const PolygonMesh& polygon_mesh,
    	const Computation_policy_3 policy =
    	Computation_policy_3::DEFAULT,
    	const GeomTraits traits = GeomTraits()) :
    Wachspress_coordinates_3(
      polygon_mesh,
      policy,
      get_const_property_map(CGAL::vertex_point, polygon_mesh),
      traits) { }

    template<typename OutIterator>
    OutIterator operator()(const Point_3& query, OutIterator c_begin) {
      return compute(query, c_begin);
    }

  private:
  	const PolygonMesh& m_polygon_mesh;
  	const Computation_policy_3 m_computation_policy;
  	const VertexToPointMap m_vertex_to_point_map; // use it to map vertex to Point_3
  	const GeomTraits m_traits;

    // REVIEW: change to m_dot_3 or m_dot_product_3, always use underscores, this is CGAL style.
    const Dot_3 m_dot3;
    const Det_3 m_det3;
    const Cross_3 m_cross3;

  	std::vector<FT> m_weights;

  	template<typename OutputIterator>
    OutputIterator compute(
      const Point_3& query, OutputIterator coordinates) {

        // Compute weights.
        const FT sum = compute_weights(query);
        CGAL_assertion(sum != FT(0));

        // The coordinates must be saved in the same order as vertices in the vertex range.
        std::size_t vi = 0;
        const auto vd = vertices(m_polygon_mesh);
        // REVIEW: CGAL_assertion(m_weights.size() == vd.size()) ?
        // Try to add as many as possible assertions.
        for (const auto vertex : vd) { // REVIEW: I get warning: vertex is unused variable

            CGAL_assertion(vi < m_weights.size());
            const FT coordinate = m_weights[vi]/sum;
            *(coordinates++) = coordinate;
            ++vi;
        }
        // REVIEW: you forgot return! Do you compile with -Wall and -Wextra ?
    }

    FT compute_weights(const Point_3& query) {

    	// Sum of weights to normalize them later.
    	FT sum = FT(0);

		// Vertex index.
    	std::size_t vi = 0;
    	const auto vd = vertices(m_polygon_mesh);
    	for (const auto vertex : vd) { // REVIEW: warning, here vertex should be passed by ref

            // Call function to calculate wp coordinates
            const FT weight = compute_wp_vertex_query(vertex, query);

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
    FT compute_wp_vertex_query(const Vertex& vertex, const Point_3& query){

        // Map vertex descriptor to point_3
        // REVIEW: use Point_3& that is with the reference.
        // Even if the property_map returns by value, the life time of this ref
        // will be extended so it is ok.
        const Point_3 vertex_val = get(m_vertex_to_point_map, vertex);

        // Circulator of faces around the vertex
        CGAL::Face_around_target_circulator<Polygon_mesh>
        face_circulator(m_polygon_mesh.halfedge(vertex), m_polygon_mesh);

        CGAL::Face_around_target_circulator<Polygon_mesh>
        done(face_circulator);

        // Vector connecting query point to vertex;
        // REVIEW: Use const Vector_3 query_vertex....
        // Or even better: https://doc.cgal.org/latest/Kernel_23/classKernel.html#abf20ed4375d501fdf4d0e139838ca0dc
        Vector_3 query_vertex(query, vertex_val);

        // First face. p_1 is negated because the order of the circulator is reversed
        // REVIEW: const const const
        Vector_3 face_normal_1 = get_face_normal(*face_circulator);
        FT dist_perp_1 = m_dot3(query_vertex, face_normal_1);
        CGAL_assertion(dist_perp_1 != 0); // REVIEW: != FT(0) because now you compare against integer but FT is double or whatsoever
        Vector_3 p_1 = -face_normal_1/dist_perp_1;
        face_circulator++;

        // Compute weight w_v
        FT weight = 0; // REVIEW: FT(0), better keep using FT(value), explicit conversion instead of implicit

        // Iterate using the circulator
        do{
            // Calculate normals of faces
            Vector_3 face_normal_i = get_face_normal(*face_circulator); face_circulator++;
            Vector_3 face_normal_i_1 = get_face_normal(*face_circulator); face_circulator++;

            // Distance of query to face
            FT perp_dist_i = m_dot3(query_vertex, face_normal_i);
            CGAL_assertion(perp_dist_i != 0);
            FT perp_dist_i_1 = m_dot3(query_vertex, face_normal_i_1);
            CGAL_assertion(perp_dist_i_1 != 0);

            // pf vector
            Vector_3 p_i = face_normal_i/perp_dist_i;
            Vector_3 p_i_1 = face_normal_i_1/perp_dist_i_1;

            // Sum partial result to weight
            weight += m_det3(p_1, p_i, p_i_1);

        }while(face_circulator!=done);

        return weight;
    }

    // Compute normal vector of the face (not normalized).
    template<typename Face>
    Vector_3 get_face_normal(const Face& face) {

      const auto hedge = halfedge(face, m_polygon_mesh);
      const auto vertices = vertices_around_face(hedge, m_polygon_mesh);
      CGAL_precondition(vertices.size() >= 3);

      auto vertex = vertices.begin();
      const Point_3& point1 = get(m_vertex_to_point_map, *vertex); ++vertex;
      const Point_3& point2 = get(m_vertex_to_point_map, *vertex); ++vertex;
      const Point_3& point3 = get(m_vertex_to_point_map, *vertex);

      const Vector_3 u = point2 - point1;
      const Vector_3 v = point3 - point1;
      const Vector_3 face_normal = m_cross3(u, v);

      return face_normal;
    }

  };

  template<
  typename Point_3,
  typename OutIterator,
  typename GeomTraits>
  OutIterator wachspress_coordinates_3(
    const CGAL::Surface_mesh<Point_3>& surface_mesh,
    const Point_3& query,
    OutIterator c_begin,
    const Computation_policy_3 policy =
    Computation_policy_3::DEFAULT) {

    using Geom_Traits = typename Kernel_traits<Point_3>::Kernel;
    using SM = CGAL::Surface_mesh<Point_3>;

    Wachspress_coordinates_3<SM, Geom_Traits> wachspress(surface_mesh, policy);
    return wachspress(query, c_begin);
  }

} // namespace Barycentric_coordinates
} // namespace CGAL

#endif // CGAL_BARYCENTRIC_WACHSPRESS_COORDINATES_3_H
