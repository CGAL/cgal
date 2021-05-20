// Copyright (c) 2020 GeometryFactory SARL (France).
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
// Author(s)     : Dmitry Anisimov, David Bommes, Kai Hormann, Pierre Alliez
//

#ifndef CGAL_BARYCENTRIC_INTERNAL_MEAN_VALUE_WEIGHTS_2_H
#define CGAL_BARYCENTRIC_INTERNAL_MEAN_VALUE_WEIGHTS_2_H

#include <CGAL/license/Barycentric_coordinates_2.h>

// Internal includes.
#include <CGAL/Barycentric_coordinates_2/internal/utils_2.h>

namespace CGAL {
namespace Barycentric_coordinates {
namespace internal {

  /*!
    \ingroup PkgWeightsRefBarycentricMeanValueWeights

    \brief 2D mean value weights for polygons.

    This class implements 2D mean value weights ( \cite cgal:bc:hf-mvcapp-06,
    \cite cgal:bc:fhk-gcbcocp-06, \cite cgal:f-mvc-03 ), which can be computed
    at any point inside and outside a simple polygon.

    Mean value weights are well-defined inside and outside a simple polygon and are
    non-negative in the kernel of a star-shaped polygon. These weights are computed
    analytically using the formulation from the `tangent_weight()`.

    \tparam VertexRange
    a model of `ConstRange` whose iterator type is `RandomAccessIterator`

    \tparam GeomTraits
    a model of `AnalyticWeightTraits_2`

    \tparam PointMap
    a model of `ReadablePropertyMap` whose key type is `VertexRange::value_type` and
    value type is `Point_2`. The default is `CGAL::Identity_property_map`.

    \cgalModels `BarycentricWeights_2`
  */
  template<
  typename VertexRange,
  typename GeomTraits,
  typename PointMap = CGAL::Identity_property_map<typename GeomTraits::Point_2> >
  class Mean_value_weights_2 {

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

      This class implements the behavior of mean value weights
      for 2D query points inside simple polygons.

      \param polygon
      an instance of `VertexRange` with the vertices of a simple polygon

      \param traits
      a traits class with geometric objects, predicates, and constructions;
      the default initialization is provided

      \param point_map
      an instance of `PointMap` that maps a vertex from `polygon` to `Point_2`;
      the default initialization is provided

      \pre polygon.size() >= 3
      \pre polygon is simple
    */
    Mean_value_weights_2(
      const VertexRange& polygon,
      const GeomTraits traits = GeomTraits(),
      const PointMap point_map = PointMap()) :
    m_polygon(polygon),
    m_traits(traits),
    m_point_map(point_map),
    m_area_2(m_traits.compute_area_2_object()),
    m_construct_vector_2(m_traits.construct_vector_2_object()),
    m_squared_length_2(m_traits.compute_squared_length_2_object()),
    m_scalar_product_2(m_traits.compute_scalar_product_2_object()),
    m_sqrt(Get_sqrt::sqrt_object(m_traits))  {

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

      This function fills a destination range with 2D mean value weights computed at
      the `query` point with respect to the vertices of the input polygon.

      The number of computed weights equals to the number of polygon vertices.

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
    OutIterator operator()(
      const Point_2& query,
      OutIterator w_begin) {

      const bool normalize = false;
      return operator()(query, w_begin, normalize);
    }

    /// @}

    /// \cond SKIP_IN_MANUAL
    template<typename OutIterator>
    OutIterator operator()(
      const Point_2& query,
      OutIterator weights,
      const bool normalize) {

      return optimal_weights(
        query, weights, normalize);
    }
    /// \endcond

  private:

    // Fields.
    const VertexRange& m_polygon;
    const GeomTraits m_traits;
    const PointMap m_point_map;

    const Area_2 m_area_2;
    const Construct_vector_2 m_construct_vector_2;
    const Squared_length_2 m_squared_length_2;
    const Scalar_product_2 m_scalar_product_2;
    const Sqrt m_sqrt;

    std::vector<Vector_2> s;
    std::vector<FT> r;
    std::vector<FT> A;
    std::vector<FT> D;
    std::vector<FT> t;
    std::vector<FT> w;

    // Functions.
    void resize() {
      s.resize(m_polygon.size());
      r.resize(m_polygon.size());
      A.resize(m_polygon.size());
      D.resize(m_polygon.size());
      t.resize(m_polygon.size());
      w.resize(m_polygon.size());
    }

    template<typename OutputIterator>
    OutputIterator optimal_weights(
      const Point_2& query,
      OutputIterator weights,
      const bool normalize) {

      // Get the number of vertices in the polygon.
      const std::size_t n = m_polygon.size();

      // Compute vectors s following the pseudo-code in the Figure 10 from [1].
      for (std::size_t i = 0; i < n; ++i) {
        const auto& pi = get(m_point_map, *(m_polygon.begin() + i));
        s[i] = m_construct_vector_2(query, pi);
      }

      // Compute lengths r, areas A, and dot products D following the pseudo-code
      // in the Figure 10 from [1]. Split the loop to make this computation faster.
      const auto& p1 = get(m_point_map, *(m_polygon.begin() + 0));
      const auto& p2 = get(m_point_map, *(m_polygon.begin() + 1));

      r[0] = m_sqrt(m_squared_length_2(s[0]));
      A[0] = m_area_2(p1, p2, query);
      D[0] = m_scalar_product_2(s[0], s[1]);

      for (std::size_t i = 1; i < n - 1; ++i) {
        const auto& pi1 = get(m_point_map, *(m_polygon.begin() + (i + 0)));
        const auto& pi2 = get(m_point_map, *(m_polygon.begin() + (i + 1)));

        r[i] = m_sqrt(m_squared_length_2(s[i]));
        A[i] = m_area_2(pi1, pi2, query);
        D[i] = m_scalar_product_2(s[i], s[i + 1]);
      }

      const auto& pn = get(m_point_map, *(m_polygon.begin() + (n - 1)));
      r[n - 1] = m_sqrt(m_squared_length_2(s[n - 1]));
      A[n - 1] = m_area_2(pn, p1, query);
      D[n - 1] = m_scalar_product_2(s[n - 1], s[0]);

      // Compute intermediate values t using the formulas from slide 19 here
      // - http://www.inf.usi.ch/hormann/nsfworkshop/presentations/Hormann.pdf
      for (std::size_t i = 0; i < n - 1; ++i) {
        CGAL_assertion((r[i] * r[i + 1] + D[i]) != FT(0));
        t[i] = FT(2) * A[i] / (r[i] * r[i + 1] + D[i]);
      }

      CGAL_assertion((r[n - 1] * r[0] + D[n - 1]) != FT(0));
      t[n - 1] = FT(2) * A[n - 1] / (r[n - 1] * r[0] + D[n - 1]);

      // Compute mean value weights using the same pseudo-code as before.
      CGAL_assertion(r[0] != FT(0));
      w[0] = FT(2) * (t[n - 1] + t[0]) / r[0];

      for (std::size_t i = 1; i < n - 1; ++i) {
        CGAL_assertion(r[i] != FT(0));
        w[i] = FT(2) * (t[i - 1] + t[i]) / r[i];
      }

      CGAL_assertion(r[n - 1] != FT(0));
      w[n - 1] = FT(2) * (t[n - 2] + t[n - 1]) / r[n - 1];

      // Normalize if necessary.
      if (normalize)
        internal::normalize(w);

      // Return weights.
      for (std::size_t i = 0; i < n; ++i)
        *(weights++) = w[i];
      return weights;
    }
  };

} // namespace internal
} // namespace Barycentric_coordinates
} // namespace CGAL

#endif // CGAL_BARYCENTRIC_INTERNAL_MEAN_VALUE_WEIGHTS_2_H
