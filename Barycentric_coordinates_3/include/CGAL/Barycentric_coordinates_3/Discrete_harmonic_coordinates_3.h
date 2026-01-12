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

#include <CGAL/license/Barycentric_coordinates_3.h>

// Internal includes.
#include <CGAL/Barycentric_coordinates_3/internal/utils_3.h>
#include <CGAL/Barycentric_coordinates_3/barycentric_enum_3.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/boost/graph/property_maps.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Named_function_parameters.h>

namespace CGAL {
namespace Barycentric_coordinates {

/*!
  \ingroup PkgBarycentricCoordinates3RefAnalytic

  \brief computes 3D discrete harmonic coordinates with respect to a closed convex triangle mesh.

  This class implements 3D discrete harmonic coordinates \cite cgal:bc:jlw-ggcccsp-07, which can be computed
  at any point inside a convex polyhedron with triangular faces.

  Discrete harmonic coordinates are well-defined in the closure of a convex polyhedron
  with triangular faces but they are not necessarily positive. The coordinates are
  computed analytically.

  \tparam TriangleMesh
  must be a model of the concept `FaceListGraph`.

  \tparam VertexPointMap
  a property map with boost::graph_traits<TriangleMesh>::vertex_descriptor as
  key type and `GeomTraits::Point_3` as value type.

  \tparam GeomTraits
  a model of `BarycentricTraits_3`.
*/
template<typename TriangleMesh,
         typename VertexPointMap = typename boost::property_map<TriangleMesh, CGAL::vertex_point_t>::const_type,
         typename GeomTraits = typename Kernel_traits<typename boost::property_traits<VertexPointMap>::value_type>::Kernel>
class Discrete_harmonic_coordinates_3 {
public:

  /// \name Types
  /// @{

  /// \cond SKIP_IN_MANUAL
  using Triangle_mesh = TriangleMesh;
  using Geom_traits = GeomTraits;
  using Vertex_point_map = VertexPointMap;

  using Compute_squared_length_3 = typename GeomTraits::Compute_squared_length_3;
  using Construct_vec_3 = typename GeomTraits::Construct_vector_3;
  using Construct_scaled_vector_3 = typename GeomTraits::Construct_scaled_vector_3;
  using Cross_3 = typename GeomTraits::Construct_cross_product_vector_3;
  using Dot_3 = typename GeomTraits::Compute_scalar_product_3;
  using Sqrt = typename internal::Get_sqrt<GeomTraits>::Sqrt;
  /// \endcond

  /// Number type.
  using FT = typename GeomTraits::FT;

  /// Point type.
  using Point_3 = typename GeomTraits::Point_3;

  /// %Vector type.
  using Vector_3 = typename GeomTraits::Vector_3;

  /// @}

  /// \name Initialization
  /// @{

  /*!
    \brief initializes all internal data structures.

    This class implements the behavior of discrete harmonic coordinates
    for 3D query points.

    \param tmesh
    an instance of `TriangleMesh`

    \param policy
    one of the `CGAL::Barycentric_coordinates::Computation_policy_3`;
    the default is `Computation_policy_3::FAST_WITH_EDGE_CASES`

    \param traits
    a traits class with geometric objects, predicates, and constructions;
    the default initialization is provided

    \param vertex_point_map
    an instance of `VertexPointMap` that maps a vertex from `tmesh` to `Point_3`;
    the default initialization is provided

    \pre boost::num_vertices(`tmesh`) >= 4.
    \pre CGAL::is_triangle_mesh(`tmesh`).
    \pre CGAL::is_closed(`tmesh`).
    \pre CGAL::is_strongly_convex_3(`tmesh`).
  */
  Discrete_harmonic_coordinates_3(const TriangleMesh& tmesh,
                                  const Computation_policy_3 policy,
                                  const VertexPointMap vertex_point_map,
                                  const GeomTraits traits = GeomTraits())
    : m_tmesh(tmesh)
    , m_computation_policy(policy)
    , m_vertex_point_map(vertex_point_map)
    , m_traits(traits)
    , m_compute_squared_length_3(m_traits.compute_squared_length_3_object())
    , m_construct_vector_3(m_traits.construct_vector_3_object())
    , m_construct_scaled_vector_3(m_traits.construct_scaled_vector_3_object())
    , m_cross_3(m_traits.construct_cross_product_vector_3_object())
    , m_dot_3(m_traits.compute_scalar_product_3_object())
    , sqrt(internal::Get_sqrt<GeomTraits>::sqrt_object(m_traits))
  {
    // Check if polyhedron is strongly convex
    if constexpr (internal_kernel_traits::Has_nested_R<GeomTraits>::value)
      CGAL_assertion(is_strongly_convex_3(m_tmesh, m_traits));

    m_weights.resize(vertices(m_tmesh).size());
  }

  /// @}

  Discrete_harmonic_coordinates_3(const TriangleMesh& tmesh,
                                  const Computation_policy_3 policy =
                                  Computation_policy_3::FAST_WITH_EDGE_CASES,
                                  const GeomTraits traits = GeomTraits())
    : Discrete_harmonic_coordinates_3(tmesh,
                                      policy,
                                      get_const_property_map(CGAL::vertex_point, tmesh),
                                      traits)
  {}

  /// \name Access
  /// @{


  /*!
    \brief computes 3D discrete harmonic coordinates with respect to a closed convex triangle mesh.

    This function fills `oi` with 3D discrete harmonic coordinates computed
    at the `query` point with respect to the vertices of the input polyhedron.

    The number of returned coordinates equals the number of vertices of `tmesh`.

    After the coordinates \f$b_i\f$ with \f$i = 0\dots n-1\f$ are computed, where
    \f$n\f$ is the number of vertices, the query point \f$q\f$ can be obtained
    as \f$q = \sum_{i = 0}^{n-1}b_ip_i\f$, where \f$p_i\f$ are the polyhedron vertices.

    \tparam OutputIterator
    a model of `OutputIterator` that accepts values of type `FT`

    \param query
    a query point

    \param oi
    the beginning of the destination range with the computed coordinates

    \return an output iterator to the element in the destination range,
    one past the last coordinate stored
  */
  template<typename OutputIterator>
  OutputIterator
  operator()(const Point_3& query, OutputIterator oi)
  {
    return compute(query, oi);
  }

  /// @}

private:
  const TriangleMesh& m_tmesh;
  const Computation_policy_3 m_computation_policy;
  const VertexPointMap m_vertex_point_map; // use it to map vertex to Point_3
  GeomTraits m_traits;

  Compute_squared_length_3 m_compute_squared_length_3;
  Construct_vec_3 m_construct_vector_3;
  Construct_scaled_vector_3 m_construct_scaled_vector_3;
  Cross_3 m_cross_3;
  Dot_3 m_dot_3;
  Sqrt sqrt;

  std::vector<FT> m_weights;

  template<typename OutputIterator>
  OutputIterator
  compute(const Point_3& query, OutputIterator coordinates)
  {
    switch(m_computation_policy){

      case Computation_policy_3::FAST:{
        return compute_coords(query, coordinates);
      }

      case Computation_policy_3::FAST_WITH_EDGE_CASES:{
        // Calculate query position relative to the polyhedron
        const auto edge_case = internal::locate_wrt_polyhedron(
          m_vertex_point_map, m_tmesh, query, coordinates, m_traits);

        if(edge_case == internal::Edge_case::BOUNDARY) {
          return coordinates;
        }
        if(edge_case == internal::Edge_case::EXTERIOR_BOUNDARY){
#ifdef CGAL_BARYCENTRIC_COORDINATES_3_VERBOSE
          std::cerr << std::endl <<
          "WARNING: query does not belong to the polyhedron!" << std::endl;
#endif
          return coordinates;
        }
        if(edge_case == internal::Edge_case::EXTERIOR) {
#ifdef CGAL_BARYCENTRIC_COORDINATES_3_VERBOSE
          std::cerr << std::endl <<
          "WARNING: query does not belong to the polyhedron!" << std::endl;
#endif
        }

        return compute_coords(query, coordinates);
      }

      default:{
        internal::get_default(vertices(m_tmesh).size(), coordinates);
        return coordinates;
      }
    }
    return coordinates;
  }

  template<typename OutputIterator>
  OutputIterator
  compute_coords(const Point_3& query, OutputIterator coordinates)
  {
    // Compute weights.
    const FT sum = compute_weights(query);
    CGAL_assertion(sum != FT(0));

    // The coordinates must be saved in the same order as vertices in the vertex range.
    const auto vd = vertices(m_tmesh);
    CGAL_assertion(m_weights.size() == vd.size());

    for (std::size_t vi = 0; vi < vd.size(); vi++) {

      CGAL_assertion(vi < m_weights.size());
      const FT coordinate = m_weights[vi]/sum;
      *(coordinates++) = coordinate;
    }

    return coordinates;
  }

  FT compute_weights(const Point_3& query)
  {
    // Sum of weights to normalize them later.
    FT sum = FT(0);

    // Vertex index.
    std::size_t vi = 0;
    const auto vd = vertices(m_tmesh);

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
  FT compute_dh_vertex_query(const Vertex& vertex, const Point_3& query)
  {
    const Point_3 vertex_val = get(m_vertex_point_map, vertex);

    // Circulator of faces around the vertex
    CGAL::Face_around_target_circulator<Triangle_mesh>
    face_circulator(halfedge(vertex, m_tmesh), m_tmesh);

    CGAL::Face_around_target_circulator<Triangle_mesh>
    face_done(face_circulator);

    // Compute weight w_v
    FT weight = FT(0);

    // Iterate using the circulator
    do{

      //Vertices around face iterator
      const auto hedge = halfedge(*face_circulator, m_tmesh);
      const auto vertices = vertices_around_face(hedge, m_tmesh);
      auto vertex_itr = vertices.begin();
      CGAL_precondition(vertices.size() == 3);

      int vertex_parity = 1;
      std::vector<Point_3> points;
      points.resize(2);
      int point_count = 0;

      for(std::size_t i = 0; i < 3; i++){

        if(*vertex_itr!=vertex){

          points[point_count] = get(m_vertex_point_map, *vertex_itr);
          point_count++;
        }
        else
          vertex_parity *= (i & 1)? -1 : 1;

        vertex_itr++;
      }

      const Point_3& point2 = points[0];
      const Point_3& point1 = points[1];

      const Vector_3 opposite_edge = m_construct_vector_3(point2, point1);
      const FT edge_length = sqrt(m_compute_squared_length_3(opposite_edge));

      const Vector_3 normal_query = m_construct_scaled_vector_3(m_cross_3(m_construct_vector_3(query, point2),
       m_construct_vector_3(query, point1)), vertex_parity);

      const Vector_3 face_normal = internal::get_face_normal(
        *face_circulator, m_vertex_point_map, m_tmesh, m_traits);

      FT cot_dihedral = internal::cot_dihedral_angle(
        face_normal, normal_query, m_traits);

      const Vector_3 vertex_query = m_construct_vector_3(vertex_val, query);

      // Treat case when the point is outside
      if(m_dot_3(face_normal, vertex_query) > 0)
        cot_dihedral *= -1;

      weight += (cot_dihedral * edge_length) / 2;
      face_circulator++;

    }while(face_circulator!=face_done);

    return weight;
  }

};

/*!
  \ingroup PkgBarycentricCoordinates3RefFunctions

  \brief computes 3D discrete harmonic coordinates with respect to a closed convex triangle mesh.

  This function computes 3D discrete harmonic coordinates at a given `query` point
  with respect to the vertices of a convex `polyhedron` with triangular faces, that is one
  coordinate per vertex. The coordinates are stored in a destination range
  beginning at `oi`.

  Internally, the class `Discrete_harmonic_coordinates_3` is used. If one wants to process
  multiple query points, it is better to use that class. When using the free function,
  internal memory is allocated for each query point, while when using the class,
  it is allocated only once, which is much more efficient. However, for a few query
  points, it is easier to use this function. It can also be used when the processing
  time is not a concern.

  \tparam TriangleMesh
  must be a model of the concept `FaceListGraph`

  \tparam Point
  A model of `GeomTraits::Point_3` with `GeomTraits` being the type of the named parameter `geom_traits`.

  \tparam OutputIterator
  must be an output iterator accepting `GeomTraits::FT` with `GeomTraits` being the type of the named parameter `geom_traits`.

  \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

  \param tmesh
  an instance of `TriangleMesh`

  \param query
  a query point

  \param oi
  the beginning of the destination range with the computed coordinates

  \param np
  an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

  \return an output iterator to the element in the destination range,
  one past the last coordinate stored

  \cgalNamedParamsBegin
    \cgalParamNBegin{computation_policy}
      \cgalParamDescription{The computation policy for handling points near the boundary of the `tmesh`.}
      \cgalParamType{`CGAL::Barycentric_coordinates::Computation_policy_3`}
      \cgalParamDefault{`Computation_policy_3::FAST_WITH_EDGE_CASES`}
    \cgalParamNEnd
    \cgalParamNBegin{vertex_point_map}
      \cgalParamDescription{a property map associating points to the vertices of `tmesh`}
      \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
                     as key type and `%Point` as value type}
      \cgalParamDefault{`boost::get(CGAL::vertex_point, tmesh)`}
      \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
                      must be available in `TriangleMesh`.}
    \cgalParamNEnd
    \cgalParamNBegin{geom_traits}
      \cgalParamDescription{an instance of a geometric traits class}
      \cgalParamType{a class model of `BarycentricTraits_3`}
      \cgalParamDefault{a \cgal Kernel deduced from the value type of the vertex-point map, using `CGAL::Kernel_traits`}
      \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
    \cgalParamNEnd
  \cgalNamedParamsEnd

  \pre boost::num_vertices(`tmesh`) >= 4.
  \pre CGAL::is_triangle_mesh(`tmesh`).
  \pre CGAL::is_closed(`tmesh`).
  \pre CGAL::is_strongly_convex_3(`tmesh`).
*/
template<typename TriangleMesh,
         typename Point,
         typename OutputIterator,
         typename NamedParameters = parameters::Default_named_parameters>
OutputIterator
discrete_harmonic_coordinates_3(const TriangleMesh& tmesh,
                                const Point& query,
                                OutputIterator oi,
                                const NamedParameters& np = parameters::default_values())
{
  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type Geom_traits;
  static_assert(std::is_same_v<Geom_traits, typename Kernel_traits<Point>::Kernel>);

  const Computation_policy_3 policy = parameters::choose_parameter(parameters::get_parameter(np, internal_np::computation_policy), Computation_policy_3::FAST_WITH_EDGE_CASES);

  using VPM = typename GetVertexPointMap<TriangleMesh, NamedParameters>::const_type;
  VPM vpm = parameters::choose_parameter(parameters::get_parameter(np, internal_np::vertex_point), get_const_property_map(CGAL::vertex_point, tmesh));

  if constexpr (std::is_same_v<typename GetGeomTraits<TriangleMesh, NamedParameters>::GT_from_NP, internal_np::Param_not_found>)
    return Discrete_harmonic_coordinates_3<TriangleMesh, VPM, Geom_traits>(tmesh, policy, vpm, Geom_traits())(query, oi);
  else
    return Discrete_harmonic_coordinates_3<TriangleMesh, VPM, Geom_traits>(tmesh, policy, vpm, parameters::get_parameter(np, internal_np::geom_traits))(query, oi);
}

} // namespace Barycentric_coordinates
} // namespace CGAL

#endif // CGAL_BARYCENTRIC_DISCRETE_HARMONIC_COORDINATES_3_H
