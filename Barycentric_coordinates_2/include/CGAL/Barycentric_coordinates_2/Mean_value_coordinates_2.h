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

#ifndef CGAL_BARYCENTRIC_MEAN_VALUE_COORDINATES_2_H
#define CGAL_BARYCENTRIC_MEAN_VALUE_COORDINATES_2_H

#include <CGAL/license/Barycentric_coordinates_2.h>

// Internal includes.
#include <CGAL/Weights/mean_value_weights.h>
#include <CGAL/Barycentric_coordinates_2/internal/utils_2.h>

// [1] Reference: "K. Hormann and M. Floater.
// Mean value coordinates for arbitrary planar polygons.
// ACM Transactions on Graphics, 25(4):1424-1441, 2006.".

// [2] Reference: "M. S. Floater.
// Wachspress and mean value coordinates.
// Proceedings of the 14th International Conference on Approximation Theory.
// G. Fasshauer and L. L. Schumaker (eds.)."

namespace CGAL {
namespace Barycentric_coordinates {

  /*!
    \ingroup PkgBarycentricCoordinates2RefAnalytic

    \brief 2D mean value coordinates.

    This class implements 2D mean value coordinates ( \cite cgal:bc:hf-mvcapp-06,
    \cite cgal:bc:fhk-gcbcocp-06, \cite cgal:f-mvc-03 ), which can be computed
    at any point in the plane.

    Mean value coordinates are well-defined everywhere in the plane and are
    non-negative in the kernel of a star-shaped polygon. The coordinates are
    computed analytically. See more details in the user manual \ref compute_mv_coord "here".

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
  class Mean_value_coordinates_2 {

  public:

    /// \name Types
    /// @{

    /// \cond SKIP_IN_MANUAL
    using Vertex_range = VertexRange;
    using Geom_traits = GeomTraits;
    using Point_map = PointMap;

    using Vector_2 = typename GeomTraits::Vector_2;
    using Area_2 = typename GeomTraits::Compute_area_2;
    using Construct_vector_2 = typename GeomTraits::Construct_vector_2;
    using Squared_length_2 = typename GeomTraits::Compute_squared_length_2;
    using Scalar_product_2 = typename GeomTraits::Compute_scalar_product_2;
    using Get_sqrt = internal::Get_sqrt<GeomTraits>;
    using Sqrt = typename Get_sqrt::Sqrt;

    using Mean_value_weights_2 =
      Weights::Mean_value_weights_2<VertexRange, GeomTraits, PointMap>;
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

      This class implements the behavior of mean value coordinates
      for 2D query points.

      \param polygon
      an instance of `VertexRange` with the vertices of a simple polygon

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
    */
    Mean_value_coordinates_2(
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
    m_construct_vector_2(m_traits.construct_vector_2_object()),
    m_squared_length_2(m_traits.compute_squared_length_2_object()),
    m_scalar_product_2(m_traits.compute_scalar_product_2_object()),
    m_sqrt(Get_sqrt::sqrt_object(m_traits)),
    m_mean_value_weights_2(
      polygon, traits, point_map) {

      CGAL_precondition(
        polygon.size() >= 3);
      CGAL_precondition(
        internal::is_simple_2(polygon, traits, point_map));
      resize();
    }

    /// @}

    /// \name Access
    /// @{

    /*!
      \brief computes 2D mean value weights.

      This function fills `weights` with 2D mean value weights computed at the `query`
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
      \brief computes 2D mean value coordinates.

      This function fills `c_begin` with 2D mean value coordinates computed
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
    const Construct_vector_2 m_construct_vector_2;
    const Squared_length_2 m_squared_length_2;
    const Scalar_product_2 m_scalar_product_2;
    const Sqrt m_sqrt;

    Mean_value_weights_2 m_mean_value_weights_2;

    std::vector<Vector_2> s;
    std::vector<FT> r;
    std::vector<FT> A;
    std::vector<FT> B;
    std::vector<FT> P;
    std::vector<FT> w;

    // Functions.
    void resize() {
      s.resize(m_polygon.size());
      r.resize(m_polygon.size());
      A.resize(m_polygon.size());
      B.resize(m_polygon.size());
      P.resize(m_polygon.size());
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
          return m_mean_value_weights_2(query, output, normalize);
        }

        case Computation_policy_2::FAST_WITH_EDGE_CASES: {
          const auto edge_case = verify(query, output);
          if (edge_case == internal::Edge_case::BOUNDARY) {
            return output;
          }
          return m_mean_value_weights_2(query, output, normalize);
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

      // Compute vectors s and its lengths r following the pseudo-code
      // in the Figure 10 from [1].
      const auto& p1 = get(m_point_map, *(m_polygon.begin() + 0));
      const auto& p2 = get(m_point_map, *(m_polygon.begin() + 1));
      const auto& pn = get(m_point_map, *(m_polygon.begin() + (n - 1)));

      s[0] = m_construct_vector_2(query, p1);
      r[0] = m_sqrt(m_squared_length_2(s[0]));

      // Compute areas A and B following the notation from [1] (see Figure 2).
      // Split the loop to make this computation faster.
      A[0] = m_area_2(p1, p2, query);
      B[0] = m_area_2(pn, p2, query);

      for (std::size_t i = 1; i < n - 1; ++i) {
        const auto& pi0 = get(m_point_map, *(m_polygon.begin() + (i - 1)));
        const auto& pi1 = get(m_point_map, *(m_polygon.begin() + (i + 0)));
        const auto& pi2 = get(m_point_map, *(m_polygon.begin() + (i + 1)));

        s[i] = m_construct_vector_2(query, pi1);
        r[i] = m_sqrt(m_squared_length_2(s[i]));

        A[i] = m_area_2(pi1, pi2, query);
        B[i] = m_area_2(pi0, pi2, query);
      }

      const auto& pm = get(m_point_map, *(m_polygon.begin() + (n - 2)));
      s[n - 1] = m_construct_vector_2(query, pn);
      r[n - 1] = m_sqrt(m_squared_length_2(s[n - 1]));

      A[n - 1] = m_area_2(pn, p1, query);
      B[n - 1] = m_area_2(pm, p1, query);

      // Following section 4.2 from [2] we denote P_j = r_j*r_{j+1} + dot_product(d_j, d_{j+1}).
      // Vector s_i from [1] corresponds to that one with the name d_i in [2].
      for (std::size_t j = 0; j < n - 1; ++j) {
        P[j] = (CGAL::max)(r[j] * r[j + 1] + m_scalar_product_2(s[j], s[j + 1]), FT(0));
      }
      P[n - 1] = (CGAL::max)(r[n - 1] * r[0] + m_scalar_product_2(s[n - 1], s[0]), FT(0));

      // Compute mean value weights using the formula (16) from [2].
      // Since the formula (16) always gives positive values,
      // we have to add a proper sign to all the weight functions.
      w[0] = r[n - 1] * r[1] - m_scalar_product_2(s[n - 1], s[1]);
      for (std::size_t j = 1; j < n - 1; ++j) {
        w[0] *= P[j];
      }
      w[0] = sign_of_weight(A[n - 1], A[0], B[0]) * m_sqrt(w[0]);

      for (std::size_t i = 1; i < n - 1; ++i) {
        w[i] = r[i - 1] * r[i + 1] - m_scalar_product_2(s[i - 1], s[i + 1]);

        for (std::size_t j = 0; j < i - 1; ++j) {
          w[i] *= P[j];
        }
        for (std::size_t j = i + 1; j < n; ++j) {
          w[i] *= P[j];
        }

        w[i] = sign_of_weight(A[i - 1], A[i], B[i]) * m_sqrt(w[i]);
      }

      w[n - 1] = r[n - 2] * r[0] - m_scalar_product_2(s[n - 2], s[0]);
      for (std::size_t j = 0; j < n - 2; ++j) {
        w[n - 1] *= P[j];
      }
      w[n - 1] = sign_of_weight(A[n - 2], A[n - 1], B[n - 1]) * m_sqrt(w[n - 1]);

      // Return coordinates.
      internal::normalize(w);
      for (std::size_t i = 0; i < n; ++i) {
        *(coordinates++) = w[i];
      }
      return coordinates;
    }

    // Return the sign of a mean value weight function.
    // We can have 3 different values: 0 if the weight = 0,
    // -1 if the weight is negative, and +1 if the weight is positive.
    FT sign_of_weight(const FT& A_prev, const FT& A, const FT& B) const {

      if (A_prev > FT(0) && A > FT(0) && B <= FT(0)) return  FT(1);
      if (A_prev < FT(0) && A < FT(0) && B >= FT(0)) return -FT(1);
      if (B > FT(0)) return  FT(1);
      if (B < FT(0)) return -FT(1);

      return FT(0);
    }
  };

  /*!
    \ingroup PkgBarycentricCoordinates2RefFunctions

    \brief computes 2D mean value weights.

    This function computes 2D mean value weights at a given `query` point
    with respect to the vertices of a simple `polygon`, that is one
    weight per vertex. The weights are stored in a destination range
    beginning at `w_begin`.

    Internally, the class `Mean_value_coordinates_2` is used. If one wants to process
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
    an instance of `PointRange` with 2D points, which form a simple polygon

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
  */
  template<
  typename PointRange,
  typename OutIterator,
  typename GeomTraits>
  OutIterator mean_value_weights_2(
    const PointRange& polygon,
    const typename GeomTraits::Point_2& query,
    OutIterator w_begin,
    const GeomTraits& traits,
    const Computation_policy_2 policy =
    Computation_policy_2::FAST_WITH_EDGE_CASES) {

    Mean_value_coordinates_2<PointRange, GeomTraits>
      mean_value(polygon, policy, traits);
    return mean_value.weights(query, w_begin);
  }

  /// \cond SKIP_IN_MANUAL
  template<
  typename PointRange,
  typename OutIterator>
  OutIterator mean_value_weights_2(
    const PointRange& polygon,
    const typename PointRange::value_type& query,
    OutIterator w_begin,
    const Computation_policy_2 policy =
    Computation_policy_2::FAST_WITH_EDGE_CASES) {

    using Point_2 = typename PointRange::value_type;
    using GeomTraits = typename Kernel_traits<Point_2>::Kernel;
    const GeomTraits traits;
    return mean_value_weights_2(
      polygon, query, w_begin, traits, policy);
  }
  /// \endcond

  /*!
    \ingroup PkgBarycentricCoordinates2RefFunctions

    \brief computes 2D mean value coordinates.

    This function computes 2D mean value coordinates at a given `query` point
    with respect to the vertices of a simple `polygon`, that is one
    coordinate per vertex. The coordinates are stored in a destination range
    beginning at `c_begin`.

    Internally, the class `Mean_value_coordinates_2` is used. If one wants to process
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
    an instance of `PointRange` with 2D points, which form a simple polygon

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
  */
  template<
  typename PointRange,
  typename OutIterator,
  typename GeomTraits>
  OutIterator mean_value_coordinates_2(
    const PointRange& polygon,
    const typename GeomTraits::Point_2& query,
    OutIterator c_begin,
    const GeomTraits& traits,
    const Computation_policy_2 policy =
    Computation_policy_2::PRECISE_WITH_EDGE_CASES) {

    Mean_value_coordinates_2<PointRange, GeomTraits>
      mean_value(polygon, policy, traits);
    return mean_value(query, c_begin);
  }

  /// \cond SKIP_IN_MANUAL
  template<
  typename PointRange,
  typename OutIterator>
  OutIterator mean_value_coordinates_2(
    const PointRange& polygon,
    const typename PointRange::value_type& query,
    OutIterator c_begin,
    const Computation_policy_2 policy =
    Computation_policy_2::PRECISE_WITH_EDGE_CASES) {

    using Point_2 = typename PointRange::value_type;
    using GeomTraits = typename Kernel_traits<Point_2>::Kernel;
    const GeomTraits traits;
    return mean_value_coordinates_2(
      polygon, query, c_begin, traits, policy);
  }
  /// \endcond

} // namespace Barycentric_coordinates
} // namespace CGAL

#endif // CGAL_BARYCENTRIC_MEAN_VALUE_COORDINATES_2_H
