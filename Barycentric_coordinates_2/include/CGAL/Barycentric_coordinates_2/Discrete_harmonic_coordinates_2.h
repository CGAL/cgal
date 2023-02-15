// Copyright (c) 2014 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Dmitry Anisimov, David Bommes, Kai Hormann, Pierre Alliez
//

#ifndef CGAL_BARYCENTRIC_DISCRETE_HARMONIC_COORDINATES_2_H
#define CGAL_BARYCENTRIC_DISCRETE_HARMONIC_COORDINATES_2_H

#include <CGAL/license/Barycentric_coordinates_2.h>

// Internal includes.
#include <CGAL/Weights/discrete_harmonic_weights.h>
#include <CGAL/Barycentric_coordinates_2/internal/utils_2.h>

// [1] Reference: "M. S. Floater, K. Hormann, and G. Kos.
// A general construction of barycentric coordinates over convex polygons.
// Advances in Computational Mathematics, 24(1-4):311-331, 2006.".

namespace CGAL {
namespace Barycentric_coordinates {

  /*!
    \ingroup PkgBarycentricCoordinates2RefAnalytic

    \brief 2D discrete harmonic coordinates.

    This class implements 2D discrete harmonic coordinates ( \cite cgal:bc:fhk-gcbcocp-06,
    \cite cgal:pp-cdmsc-93, \cite cgal:bc:eddhls-maam-95 ), which can be computed
    at any point inside a strictly convex polygon.

    Discrete harmonic coordinates are well-defined in the closure of a strictly
    convex polygon but they are not necessarily positive. The coordinates are
    computed analytically. See more details in the user manual \ref compute_dh_coord "here".

    \tparam VertexRange
    a model of `ConstRange` whose iterator type is `RandomAccessIterator`

    \tparam GeomTraits
    a model of `BarycentricTraits_2`

    \tparam PointMap
    a model of `ReadablePropertyMap` whose key type is `VertexRange::value_type` and
    value type is `Point_2`. The default is `CGAL::Identity_property_map`.
  */
  template<
  typename VertexRange,
  typename GeomTraits,
  typename PointMap = CGAL::Identity_property_map<typename GeomTraits::Point_2> >
  class Discrete_harmonic_coordinates_2 {

  public:

    /// \name Types
    /// @{

    /// \cond SKIP_IN_MANUAL
    using Vertex_range = VertexRange;
    using Geom_traits = GeomTraits;
    using Point_map = PointMap;

    using Area_2 = typename GeomTraits::Compute_area_2;
    using Squared_distance_2 = typename GeomTraits::Compute_squared_distance_2;

    using Discrete_harmonic_weights_2 =
      Weights::Discrete_harmonic_weights_2<VertexRange, GeomTraits, PointMap>;
    /// \endcond

    /// Number type.
    typedef typename GeomTraits::FT FT;

    /// Point type.
    typedef typename GeomTraits::Point_2 Point_2;

    /// @}

    /// \name Initialization
    /// @{

    /*!
      \brief initializes all internal data structures.

      This class implements the behavior of discrete harmonic coordinates
      for 2D query points.

      \param polygon
      an instance of `VertexRange` with the vertices of a strictly convex polygon

      \param policy
      one of the `Computation_policy_2`;
      the default is `Computation_policy_2::PRECISE_WITH_EDGE_CASES`

      \param traits
      a traits class with geometric objects, predicates, and constructions;
      the default initialization is provided

      \param point_map
      an instance of `PointMap` that maps a vertex from `polygon` to `Point_2`;
      the default initialization is provided

      \pre polygon.size() >= 3
      \pre polygon is simple
      \pre polygon is strictly convex
    */
    Discrete_harmonic_coordinates_2(
      const VertexRange& polygon,
      const Computation_policy_2 policy
      = Computation_policy_2::PRECISE_WITH_EDGE_CASES,
      const GeomTraits traits = GeomTraits(),
      const PointMap point_map = PointMap()) :
    m_polygon(polygon),
    m_computation_policy(policy),
    m_traits(traits),
    m_point_map(point_map),
    m_area_2(m_traits.compute_area_2_object()),
    m_squared_distance_2(m_traits.compute_squared_distance_2_object()),
    m_discrete_harmonic_weights_2(
      polygon, traits, point_map) {

      CGAL_precondition(
        polygon.size() >= 3);
      CGAL_precondition(
        internal::is_simple_2(polygon, traits, point_map));
      CGAL_precondition(
        internal::polygon_type_2(polygon, traits, point_map) ==
        internal::Polygon_type::STRICTLY_CONVEX);
      resize();
    }

    /// @}

    /// \name Access
    /// @{

    /*!
      \brief computes 2D discrete harmonic weights.

      This function fills `weights` with 2D discrete harmonic weights computed at the `query`
      point with respect to the vertices of the input polygon. If `query` belongs to
      the polygon boundary, the returned weights are normalized.

      The number of returned weights equals to the number of polygon vertices.

      \tparam OutIterator
      a model of `OutputIterator` that accepts values of type `FT`

      \param query
      a query point

      \param w_begin
      the beginning of the destination range with the computed weights

      \return an output iterator to the element in the destination range,
      one past the last weight stored
    */
    template<typename OutIterator>
    OutIterator weights(const Point_2& query, OutIterator w_begin) {
      const bool normalize = false;
      return compute(query, w_begin, normalize);
    }

    /*!
      \brief computes 2D discrete harmonic coordinates.

      This function fills `c_begin` with 2D discrete harmonic coordinates computed
      at the `query` point with respect to the vertices of the input polygon.

      The number of returned coordinates equals to the number of polygon vertices.

      After the coordinates \f$b_i\f$ with \f$i = 1\dots n\f$ are computed, where
      \f$n\f$ is the number of polygon vertices, the query point \f$q\f$ can be obtained
      as \f$q = \sum_{i = 1}^{n}b_ip_i\f$, where \f$p_i\f$ are the polygon vertices.

      \tparam OutIterator
      a model of `OutputIterator` that accepts values of type `FT`

      \param query
      a query point

      \param c_begin
      the beginning of the destination range with the computed coordinates

      \return an output iterator to the element in the destination range,
      one past the last coordinate stored
    */
    template<typename OutIterator>
    OutIterator operator()(const Point_2& query, OutIterator c_begin) {
      const bool normalize = true;
      return compute(query, c_begin, normalize);
    }

    /// @}

  private:

    // Fields.
    const VertexRange& m_polygon;
    const Computation_policy_2 m_computation_policy;
    const GeomTraits m_traits;
    const PointMap m_point_map;

    const Area_2 m_area_2;
    const Squared_distance_2 m_squared_distance_2;

    Discrete_harmonic_weights_2 m_discrete_harmonic_weights_2;

    std::vector<FT> r;
    std::vector<FT> A;
    std::vector<FT> B;
    std::vector<FT> w;

    // Functions.
    void resize() {
      r.resize(m_polygon.size());
      A.resize(m_polygon.size());
      B.resize(m_polygon.size());
      w.resize(m_polygon.size());
    }

    template<typename OutputIterator>
    OutputIterator compute(
      const Point_2& query, OutputIterator output, const bool normalize) {

      switch (m_computation_policy) {

        case Computation_policy_2::PRECISE: {
          if (normalize) {
            return max_precision_coordinates(query, output);
          } else {
            std::cerr << "WARNING: you can't use the precise version of unnormalized weights! ";
            std::cerr << "They are not valid weights!" << std::endl;
            internal::get_default(m_polygon.size(), output);
            return output;
          }
        }

        case Computation_policy_2::PRECISE_WITH_EDGE_CASES: {
          const auto edge_case = verify(query, output);
          if (edge_case == internal::Edge_case::BOUNDARY) {
            return output;
          }
          if (edge_case == internal::Edge_case::EXTERIOR) {
            std::cerr << std::endl <<
            "WARNING: query does not belong to the polygon!" << std::endl;
          }
          if (normalize) {
            return max_precision_coordinates(query, output);
          } else {
            std::cerr << "WARNING: you can't use the precise version of unnormalized weights! ";
            std::cerr << "They are not valid weights!" << std::endl;
            internal::get_default(m_polygon.size(), output);
            return output;
          }
        }

        case Computation_policy_2::FAST: {
          return m_discrete_harmonic_weights_2(query, output, normalize);
        }

        case Computation_policy_2::FAST_WITH_EDGE_CASES: {
          const auto edge_case = verify(query, output);
          if (edge_case == internal::Edge_case::BOUNDARY) {
            return output;
          }
          if (edge_case == internal::Edge_case::EXTERIOR) {
            std::cerr << std::endl <<
            "WARNING: query does not belong to the polygon!" << std::endl;
          }
          return m_discrete_harmonic_weights_2(query, output, normalize);
        }

        default: {
          internal::get_default(m_polygon.size(), output);
          return output;
        }
      }
      return output;
    }

    template<typename OutputIterator>
    internal::Edge_case verify(
      const Point_2& query, OutputIterator output) const {

      const auto result = internal::locate_wrt_polygon_2(
        m_polygon, query, m_traits, m_point_map);
      if (!result) {
        return internal::Edge_case::EXTERIOR;
      }

      const auto location = (*result).first;
      const std::size_t index = (*result).second;
      if (location == internal::Query_point_location::ON_UNBOUNDED_SIDE) {
        return internal::Edge_case::EXTERIOR;
      }

      if (
        location == internal::Query_point_location::ON_VERTEX ||
        location == internal::Query_point_location::ON_EDGE ) {
        internal::boundary_coordinates_2(
          m_polygon, query, location, index, output, m_traits, m_point_map);
        return internal::Edge_case::BOUNDARY;
      }
      return internal::Edge_case::INTERIOR;
    }

    template<typename OutputIterator>
    OutputIterator max_precision_coordinates(
      const Point_2& query, OutputIterator coordinates) {

      // Get the number of vertices in the polygon.
      const std::size_t n = m_polygon.size();

      // Compute areas A, B, and distances r following the notation from [1].
      // Split the loop to make this computation faster.
      const auto& p1 = get(m_point_map, *(m_polygon.begin() + 0));
      const auto& p2 = get(m_point_map, *(m_polygon.begin() + 1));
      const auto& pn = get(m_point_map, *(m_polygon.begin() + (n - 1)));

      r[0] = m_squared_distance_2(p1, query);
      A[0] = m_area_2(p1, p2, query);
      B[0] = m_area_2(pn, p2, query);

      for (std::size_t i = 1; i < n - 1; ++i) {
        const auto& pi0 = get(m_point_map, *(m_polygon.begin() + (i - 1)));
        const auto& pi1 = get(m_point_map, *(m_polygon.begin() + (i + 0)));
        const auto& pi2 = get(m_point_map, *(m_polygon.begin() + (i + 1)));

        r[i] = m_squared_distance_2(pi1, query);
        A[i] = m_area_2(pi1, pi2, query);
        B[i] = m_area_2(pi0, pi2, query);
      }

      const auto& pm = get(m_point_map, *(m_polygon.begin() + (n - 2)));
      r[n - 1] = m_squared_distance_2(pn, query);
      A[n - 1] = m_area_2(pn, p1, query);
      B[n - 1] = m_area_2(pm, p1, query);

      // Initialize weights with the numerator of the formula (25) with p = 2 from [1].
      // Then we multiply them by areas A as in the formula (5) in [1]. We also split the loop.
      w[0] = r[1] * A[n - 1] - r[0] * B[0] + r[n - 1] * A[0];
      for (std::size_t j = 1; j < n - 1; ++j) {
        w[0] *= A[j];
      }

      for (std::size_t i = 1; i < n - 1; ++i) {
        w[i] = r[i + 1] * A[i - 1] - r[i] * B[i] + r[i - 1] * A[i];

        for (std::size_t j = 0; j < i - 1; ++j) {
          w[i] *= A[j];
        }
        for (std::size_t j = i + 1; j < n; ++j) {
          w[i] *= A[j];
        }
      }

      w[n - 1] = r[0] * A[n - 2] - r[n - 1] * B[n - 1] + r[n - 2] * A[n - 1];
      for (std::size_t j = 0; j < n - 2; ++j) {
        w[n - 1] *= A[j];
      }

      // Return coordinates.
      internal::normalize(w);
      for (std::size_t i = 0; i < n; ++i) {
        *(coordinates++) = w[i];
      }
      return coordinates;
    }
  };

  /*!
    \ingroup PkgBarycentricCoordinates2RefFunctions

    \brief computes 2D discrete harmonic weights.

    This function computes 2D discrete harmonic weights at a given `query` point
    with respect to the vertices of a strictly convex `polygon`, that is one
    weight per vertex. The weights are stored in a destination range
    beginning at `w_begin`.

    Internally, the class `Discrete_harmonic_coordinates_2` is used. If one wants to process
    multiple query points, it is better to use that class. When using the free function,
    internal memory is allocated for each query point, while when using the class,
    it is allocated only once, which is much more efficient. However, for a few query
    points, it is easier to use this function. It can also be used when the processing
    time is not a concern.

    \tparam PointRange
    a model of `ConstRange` whose iterator type is `RandomAccessIterator`
    and value type is `GeomTraits::Point_2`

    \tparam OutIterator
    a model of `OutputIterator` that accepts values of type `GeomTraits::FT`

    \tparam GeomTraits
    a model of `BarycentricTraits_2`

    \param polygon
    an instance of `PointRange` with 2D points, which form a strictly convex polygon

    \param query
    a query point

    \param w_begin
    the beginning of the destination range with the computed weights

    \param traits
    a traits class with geometric objects, predicates, and constructions;
    this parameter can be omitted if the traits class can be deduced from the point type

    \param policy
    one of the `Computation_policy_2`;
    the default is `Computation_policy_2::FAST_WITH_EDGE_CASES`

    \return an output iterator to the element in the destination range,
    one past the last weight stored

    \pre polygon.size() >= 3
    \pre polygon is simple
    \pre polygon is strictly convex
  */
  template<
  typename PointRange,
  typename OutIterator,
  typename GeomTraits>
  OutIterator discrete_harmonic_weights_2(
    const PointRange& polygon,
    const typename GeomTraits::Point_2& query,
    OutIterator w_begin,
    const GeomTraits& traits,
    const Computation_policy_2 policy =
    Computation_policy_2::FAST_WITH_EDGE_CASES) {

    Discrete_harmonic_coordinates_2<PointRange, GeomTraits>
      discrete_harmonic(polygon, policy, traits);
    return discrete_harmonic.weights(query, w_begin);
  }

  /// \cond SKIP_IN_MANUAL
  template<
  typename PointRange,
  typename OutIterator>
  OutIterator discrete_harmonic_weights_2(
    const PointRange& polygon,
    const typename PointRange::value_type& query,
    OutIterator w_begin,
    const Computation_policy_2 policy =
    Computation_policy_2::FAST_WITH_EDGE_CASES) {

    using Point_2 = typename PointRange::value_type;
    using GeomTraits = typename Kernel_traits<Point_2>::Kernel;
    const GeomTraits traits;
    return discrete_harmonic_weights_2(
      polygon, query, w_begin, traits, policy);
  }
  /// \endcond

  /*!
    \ingroup PkgBarycentricCoordinates2RefFunctions

    \brief computes 2D discrete harmonic coordinates.

    This function computes 2D discrete harmonic coordinates at a given `query` point
    with respect to the vertices of a strictly convex `polygon`, that is one
    coordinate per vertex. The coordinates are stored in a destination range
    beginning at `c_begin`.

    Internally, the class `Discrete_harmonic_coordinates_2` is used. If one wants to process
    multiple query points, it is better to use that class. When using the free function,
    internal memory is allocated for each query point, while when using the class,
    it is allocated only once, which is much more efficient. However, for a few query
    points, it is easier to use this function. It can also be used when the processing
    time is not a concern.

    \tparam PointRange
    a model of `ConstRange` whose iterator type is `RandomAccessIterator`
    and value type is `GeomTraits::Point_2`

    \tparam OutIterator
    a model of `OutputIterator` that accepts values of type `GeomTraits::FT`

    \tparam GeomTraits
    a model of `BarycentricTraits_2`

    \param polygon
    an instance of `PointRange` with 2D points, which form a strictly convex polygon

    \param query
    a query point

    \param c_begin
    the beginning of the destination range with the computed coordinates

    \param traits
    a traits class with geometric objects, predicates, and constructions;
    this parameter can be omitted if the traits class can be deduced from the point type

    \param policy
    one of the `Computation_policy_2`;
    the default is `Computation_policy_2::PRECISE_WITH_EDGE_CASES`

    \return an output iterator to the element in the destination range,
    one past the last coordinate stored

    \pre polygon.size() >= 3
    \pre polygon is simple
    \pre polygon is strictly convex
  */
  template<
  typename PointRange,
  typename OutIterator,
  typename GeomTraits>
  OutIterator discrete_harmonic_coordinates_2(
    const PointRange& polygon,
    const typename GeomTraits::Point_2& query,
    OutIterator c_begin,
    const GeomTraits& traits,
    const Computation_policy_2 policy =
    Computation_policy_2::PRECISE_WITH_EDGE_CASES) {

    Discrete_harmonic_coordinates_2<PointRange, GeomTraits>
      discrete_harmonic(polygon, policy, traits);
    return discrete_harmonic(query, c_begin);
  }

  /// \cond SKIP_IN_MANUAL
  template<
  typename PointRange,
  typename OutIterator>
  OutIterator discrete_harmonic_coordinates_2(
    const PointRange& polygon,
    const typename PointRange::value_type& query,
    OutIterator c_begin,
    const Computation_policy_2 policy =
    Computation_policy_2::PRECISE_WITH_EDGE_CASES) {

    using Point_2 = typename PointRange::value_type;
    using GeomTraits = typename Kernel_traits<Point_2>::Kernel;
    const GeomTraits traits;
    return discrete_harmonic_coordinates_2(
      polygon, query, c_begin, traits, policy);
  }
  /// \endcond

} // namespace Barycentric_coordinates
} // namespace CGAL

#endif // CGAL_BARYCENTRIC_DISCRETE_HARMONIC_COORDINATES_2_H
