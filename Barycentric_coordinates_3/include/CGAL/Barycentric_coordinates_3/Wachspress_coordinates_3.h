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

#include <CGAL/license/Barycentric_coordinates_3.h>

// Internal includes.
#include <CGAL/Barycentric_coordinates_3/internal/utils_3.h>
#include <CGAL/Barycentric_coordinates_3/barycentric_enum_3.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/boost/graph/property_maps.h>
#include <CGAL/Named_function_parameters.h>

namespace CGAL {
namespace Barycentric_coordinates {


/*!
  \ingroup PkgBarycentricCoordinates3RefAnalytic

  \brief computes 3D Wachspress coordinates with respect to a closed convex triangle mesh.

  This class implements 3D Wachspress coordinates ( \cite cgal:bc:f-wmvc-14,
  \cite cgal:bc:jlw-ggcccsp-07 ), which can be computed at any point inside a convex
  polyhedron with triangular faces.

  Wachspress coordinates are well-defined and non-negative in the closure of a convex
  polyhedron with triangular faces. The coordinates are computed analytically.

  \tparam TriangleMesh
  must be a model of the concept `FaceListGraph`

  \tparam GeomTraits
  a model of `BarycentricTraits_3`

  \tparam VertexPointMap
  a property map with boost::graph_traits<TriangleMesh>::vertex_descriptor as
  key type and `GeomTraits::Point_3` as value type
*/
template<typename TriangleMesh,
         typename GeomTraits,
         typename VertexPointMap = typename boost::property_map<TriangleMesh, CGAL::vertex_point_t>::const_type>
class Wachspress_coordinates_3 {
public:

  /// \name Types
  /// @{

  /// \cond SKIP_IN_MANUAL
  using Triangle_mesh = TriangleMesh;
  using Geom_traits = GeomTraits;
  using Vertex_point_map = VertexPointMap;

  using Dot_3 = typename GeomTraits::Compute_scalar_product_3;
  using Det_3 = typename GeomTraits::Compute_determinant_3;
  using Cross_3 = typename GeomTraits::Construct_cross_product_vector_3;
  using Compute_squared_length_3 = typename GeomTraits::Compute_squared_length_3;
  using Construct_vec_3 = typename GeomTraits::Construct_vector_3;
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

    This class implements the behavior of Wachspress coordinates
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
  Wachspress_coordinates_3(const TriangleMesh& tmesh,
                           const Computation_policy_3 policy,
                           const VertexPointMap vertex_point_map,
                           const GeomTraits traits = GeomTraits())
    : m_tmesh(tmesh)
    , m_computation_policy(policy)
    , m_vertex_point_map(vertex_point_map)
    , m_traits(traits)
    , m_dot_3(m_traits.compute_scalar_product_3_object())
    , m_det_3(m_traits.compute_determinant_3_object())
    , m_cross_3(m_traits.construct_cross_product_vector_3_object())
    , m_compute_squared_length_3(m_traits.compute_squared_length_3_object())
    , m_construct_vector_3(m_traits.construct_vector_3_object())
  {
    // Check if polyhedron is strongly convex
    if constexpr (internal_kernel_traits::Has_nested_R<GeomTraits>::value)
      CGAL_assertion(is_strongly_convex_3(m_tmesh, m_traits));

    m_weights.resize(vertices(m_tmesh).size());
  }

  /// @}

  Wachspress_coordinates_3(const TriangleMesh& tmesh,
                           const Computation_policy_3 policy =
                           Computation_policy_3::FAST_WITH_EDGE_CASES,
                           const GeomTraits traits = GeomTraits())
    : Wachspress_coordinates_3(tmesh,
                               policy,
                               get_const_property_map(CGAL::vertex_point, tmesh),
                               traits)
  {}

  /// \name Access
  /// @{

  /*!
    \brief computes 3D Wachspress coordinates with respect to a closed convex triangle mesh.

    This function fills `oi` with 3D Wachspress coordinates computed
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
  OutputIterator operator()(const Point_3& query, OutputIterator oi)
  {
    return compute(query, oi);
  }

  /// @}

private:
  const TriangleMesh& m_tmesh;
  const Computation_policy_3 m_computation_policy;
  const VertexPointMap m_vertex_point_map; // use it to map vertex to Point_3
  GeomTraits m_traits;

  Dot_3 m_dot_3;
  Det_3 m_det_3;
  Cross_3 m_cross_3;
  Compute_squared_length_3 m_compute_squared_length_3;
  Construct_vec_3 m_construct_vector_3;

  std::vector<FT> m_weights;

  template<typename OutputIterator>
  OutputIterator compute(const Point_3& query, OutputIterator coordinates)
  {
    switch(m_computation_policy){

      case Computation_policy_3::FAST:{
        return compute_coords(query, coordinates);
      }

      case Computation_policy_3::FAST_WITH_EDGE_CASES:{
        // Calculate query position relative to the polyhedron
        const auto edge_case = internal::locate_wrt_polyhedron(
          m_vertex_point_map, m_tmesh, query, coordinates, m_traits, true);

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
  OutputIterator compute_coords(
    const Point_3& query, OutputIterator coordinates)
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
    for (auto vertex : vertices(m_tmesh)) {

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
  // cf. Wachspress and mean value coordinates by Michael S. Floater Pages 17~18
  template<typename Vertex>
  FT compute_wp_vertex_query(const Vertex& vertex, const Point_3& query)
  {
    // Map vertex descriptor to point_3
    const Point_3& vertex_val = get(m_vertex_point_map, vertex);
    // Vector connecting query point to vertex;
    const Vector_3 query_vertex = m_construct_vector_3(query, vertex_val);

    // Loop on the faces the vertex
    using halfedge_descriptor = typename boost::graph_traits<TriangleMesh>::halfedge_descriptor;
    halfedge_descriptor first_h = halfedge(vertex, m_tmesh);

    typename Geom_traits::Construct_divided_vector_3 construct_divided_vector_3 =
      m_traits.construct_divided_vector_3_object();

    auto compute_pf_i = [&](halfedge_descriptor h)
    {
      Vector_3 nf = internal::get_face_normal(
        face(h, m_tmesh), m_vertex_point_map, m_tmesh, m_traits);
      nf = construct_divided_vector_3(nf, approximate_sqrt(m_compute_squared_length_3(nf)));
      const FT hfx = m_dot_3(query_vertex, nf);
      CGAL_assertion(hfx != FT(0));
      return construct_divided_vector_3(nf, hfx);
    };

    // First face.
    const Vector_3 pf_1 = compute_pf_i(first_h);
    halfedge_descriptor h_i=prev(opposite(first_h, m_tmesh), m_tmesh);
    Vector_3 pf_i=compute_pf_i(h_i);
    // Compute weight w_v
    FT weight = FT(0);

    // Iterate using the circulator
    do{
      halfedge_descriptor h_i_p_1=prev(opposite(h_i, m_tmesh), m_tmesh);
      if (h_i_p_1==first_h)
        break;
      const Vector_3 pf_i_p_1=compute_pf_i(h_i_p_1);

      // Sum partial result to weight
      weight += m_det_3(pf_1, pf_i, pf_i_p_1);
      h_i=h_i_p_1;
      pf_i = pf_i_p_1;
    }while(true);

    return weight;
  }

};

/*!
  \ingroup PkgBarycentricCoordinates3RefFunctions

  \brief computes 3D Wachspress coordinates with respect to a closed convex triangle mesh.

  This function computes 3D Wachspress coordinates at a given `query` point
  with respect to the vertices of `tmesh`, that is one
  coordinate per vertex. The coordinates are stored in a destination range
  beginning at `oi`.

  Internally, the class `Wachspress_coordinates_3` is used. If one wants to process
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
wachspress_coordinates_3(const TriangleMesh& tmesh,
                         const Point query,
                         OutputIterator oi,
                         const NamedParameters& np = parameters::default_values())
{
  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type Geom_traits;
  static_assert(std::is_same_v<Geom_traits, typename Kernel_traits<Point>::Kernel>);

  const Computation_policy_3 policy = parameters::choose_parameter(parameters::get_parameter(np, internal_np::computation_policy), Computation_policy_3::FAST_WITH_EDGE_CASES);

  using VPM = typename GetVertexPointMap<TriangleMesh, NamedParameters>::const_type;
  VPM vpm = parameters::choose_parameter(parameters::get_parameter(np, internal_np::vertex_point), get_const_property_map(CGAL::vertex_point, tmesh));

  if constexpr (std::is_same_v<typename GetGeomTraits<TriangleMesh, NamedParameters>::GT_from_NP, internal_np::Param_not_found>)
    return Wachspress_coordinates_3<TriangleMesh, Geom_traits, VPM>(tmesh, policy, vpm, Geom_traits())(query, oi);
  else
    return Wachspress_coordinates_3<TriangleMesh, Geom_traits, VPM>(tmesh, policy, vpm, parameters::get_parameter(np, internal_np::geom_traits))(query, oi);
}

} // namespace Barycentric_coordinates
} // namespace CGAL

#endif // CGAL_BARYCENTRIC_WACHSPRESS_COORDINATES_3_H
