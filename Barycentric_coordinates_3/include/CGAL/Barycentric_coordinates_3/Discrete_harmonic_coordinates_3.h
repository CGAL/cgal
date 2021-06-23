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

#ifndef CGAL_BARYCENTRIC_DISCRETE_HARMONIC_COORDINATES_3_H
#define CGAL_BARYCENTRIC_DISCRETE_HARMONIC_COORDINATES_3_H

// #include <CGAL/license/Barycentric_coordinates_3.h>

// Internal includes.
#include <CGAL/Barycentric_coordinates_3/internal/utils_3.h>
#include <CGAL/Barycentric_coordinates_3/barycentric_enum_3.h>

namespace CGAL {
namespace Barycentric_coordinates {

  //Default sqrt
  template<class Traits>
  class Default_sqrt{
      typedef typename Traits::FT FT;

  public:
      FT operator()(const FT &value) const{
          return FT(CGAL::sqrt(CGAL::to_double(value)));
      }
  };

  BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Has_nested_type_Sqrt, Sqrt, false)

  // Case: do_not_use_default = false.
  template<class Traits, bool do_not_use_default = Has_nested_type_Sqrt<Traits>::value>
      class Get_sqrt
  {
  public:
      typedef Default_sqrt<Traits> Sqrt;

      static Sqrt sqrt_object(const Traits&)
      {
          return Sqrt();
      }
  };

  // Case: do_not_use_default = true.
  template<class Traits>
      class Get_sqrt<Traits, true>
  {
  public:
      typedef typename Traits::Sqrt Sqrt;

      static Sqrt sqrt_object(const Traits &traits)
      {
          return traits.sqrt_object();
      }
  };

  template<
  typename PolygonMesh,
  typename GeomTraits,
  typename VertexToPointMap = typename property_map_selector<PolygonMesh,
    CGAL::vertex_point_t>::const_type>
  class Discrete_harmonic_coordinates_3 {

  public:
    using Polygon_mesh = PolygonMesh;
    using Geom_traits = GeomTraits;
    using Vertex_to_point_map = VertexToPointMap;

    using Dihedral_angle_3 = typename GeomTraits::Compute_approximate_dihedral_angle_3;
    using Construct_vec_3 = typename GeomTraits::Construct_vector_3;
    using Cross_3 = typename GeomTraits::Construct_cross_product_vector_3;
    using Dot_3 = typename GeomTraits::Compute_scalar_product_3;
    using Sqrt = typename Get_sqrt<GeomTraits>::Sqrt;

	  typedef typename GeomTraits::FT FT;
    typedef typename GeomTraits::Point_3 Point_3;
    typedef typename GeomTraits::Vector_3 Vector_3;

    Discrete_harmonic_coordinates_3(
      const PolygonMesh& polygon_mesh,
      const Computation_policy_3 policy,
      const VertexToPointMap vertex_to_point_map,
      const GeomTraits traits = GeomTraits()) :
    m_polygon_mesh(polygon_mesh),
    m_computation_policy(policy),
    m_vertex_to_point_map(vertex_to_point_map),
    m_traits(traits),
    m_dihedral_angle_3(m_traits.compute_approximate_dihedral_angle_3_object()),
    m_construct_vector_3(m_traits.construct_vector_3_object()),
    m_cross_3(m_traits.construct_cross_product_vector_3_object()),
    m_dot_3(m_traits.compute_scalar_product_3_object()),
    sqrt(Get_sqrt<GeomTraits>::sqrt_object(m_traits)){

      // Check if polyhedron is strongly convex
      CGAL_assertion(is_strongly_convex_3(m_polygon_mesh, m_traits));
      m_weights.resize(vertices(m_polygon_mesh).size());
    }

    Discrete_harmonic_coordinates_3(
      const PolygonMesh& polygon_mesh,
      const Computation_policy_3 policy =
      Computation_policy_3::DEFAULT,
      const GeomTraits traits = GeomTraits()) :
    Discrete_harmonic_coordinates_3(
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

    const Dihedral_angle_3 m_dihedral_angle_3;
    const Construct_vec_3 m_construct_vector_3;
    const Cross_3 m_cross_3;
    const Dot_3 m_dot_3;
    const Sqrt sqrt;

  	std::vector<FT> m_weights;

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

        // Call function to calculate wp coordinates
        const FT weight = compute_dh_vertex_query(vertex, query);

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
    FT compute_dh_vertex_query(const Vertex& vertex, const Point_3& query){

      // Map vertex descriptor to point_3
      const Point_3& vertex_val = get(m_vertex_to_point_map, vertex);

      // Circulator of faces around the vertex
      CGAL::Face_around_target_circulator<Polygon_mesh>
      face_circulator(m_polygon_mesh.halfedge(vertex), m_polygon_mesh);

      CGAL::Face_around_target_circulator<Polygon_mesh>
      done(face_circulator);

      // Compute weight w_v
      FT weight = FT(0);

      std::cout << query << " \n";

      // Iterate using the circulator
      do{

        const auto hedge = halfedge(*face_circulator, m_polygon_mesh);
        const auto vertices = vertices_around_face(hedge, m_polygon_mesh);
        CGAL_precondition(vertices.size() >= 3);

        auto vertex_itr = vertices.begin();
        std::vector<Point_3> points;

        for(std::size_t i = 0; i < 3; i++){

          if(*vertex_itr != vertex)
            points.push_back(get(m_vertex_to_point_map, *vertex_itr));

          ++vertex_itr;
        }

        const Point_3& point1 = points[0];
        const Point_3& point2 = points[1];
        Vector_3 opposite_edge = m_construct_vector_3(point2, point1);
        FT edge_length = sqrt(opposite_edge.squared_length ());

        Vector_3 normal_query = m_cross_3(m_construct_vector_3(point2, query),
         m_construct_vector_3(query, point1));

        FT cot_dihedral = cot_dihedral_angle(get_face_normal(*face_circulator), normal_query);

        weight += (cot_dihedral * edge_length) / 2;
        face_circulator++;
        std::cout << (cot_dihedral) << " a\n";

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
      const Vector_3 face_normal = m_cross_3(u, v);

      return face_normal;
    }

    //Compute cotangent of dihedral angle between two faces
    FT cot_dihedral_angle(const Vector_3& normal_1, const Vector_3& normal_2){

      FT approximate_cos = m_dot_3(normal_1, normal_2) /
        sqrt(normal_1.squared_length()*normal_2.squared_length());

      FT approximate_sin = sqrt(m_cross_3(normal_1, normal_2).squared_length()) /
        sqrt(normal_1.squared_length()*normal_2.squared_length());

      assert(approximate_sin != FT(0));

      return approximate_cos/approximate_sin;
    }

  };

  template<
  typename Point_3,
  typename OutIterator,
  typename GeomTraits>
  OutIterator discrete_harmonic_coordinates_3(
    const CGAL::Surface_mesh<Point_3>& surface_mesh,
    const Point_3& query,
    OutIterator c_begin,
    const Computation_policy_3 policy =
    Computation_policy_3::DEFAULT) {

    using Geom_Traits = typename Kernel_traits<Point_3>::Kernel;
    using SM = CGAL::Surface_mesh<Point_3>;

    Discrete_harmonic_coordinates_3<SM, Geom_Traits> discrete_harmonic(surface_mesh, policy);
    return discrete_harmonic(query, c_begin);
  }

} // namespace Barycentric_coordinates
} // namespace CGAL

#endif // CGAL_BARYCENTRIC_DISCRETE_HARMONIC_COORDINATES_3_H
