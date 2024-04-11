// Copyright (c) 2020 GeometryFactory SARL (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Dmitry Anisimov, David Bommes, Kai Hormann, Pierre Alliez
//

#ifndef CGAL_MEAN_VALUE_WEIGHTS_H
#define CGAL_MEAN_VALUE_WEIGHTS_H

#include <CGAL/Weights/internal/utils.h>
#include <CGAL/Weights/internal/polygon_utils_2.h>

#include <CGAL/Kernel_traits.h>
#include <CGAL/Point_2.h>
#include <CGAL/Point_3.h>
#include <CGAL/property_map.h>

#include <vector>

namespace CGAL {
namespace Weights {

/// \cond SKIP_IN_MANUAL

namespace mean_value_ns {

template<typename FT>
FT sign_of_weight(const FT A0, const FT A2, const FT B)
{
  if (A0 > FT(0) && A2 > FT(0) && B <= FT(0))
    return +FT(1);

  if (A0 < FT(0) && A2 < FT(0) && B >= FT(0))
    return -FT(1);

  if (B > FT(0))
    return +FT(1);

  if (B < FT(0))
    return -FT(1);

  return FT(0);
}

template<typename GeomTraits>
typename GeomTraits::FT weight(const typename GeomTraits::FT d0,
                               const typename GeomTraits::FT d2,
                               const typename GeomTraits::FT d,
                               const typename GeomTraits::FT D0,
                               const typename GeomTraits::FT D2,
                               const typename GeomTraits::FT D,
                               const typename GeomTraits::FT sign,
                               const GeomTraits& traits)
{
  using FT = typename GeomTraits::FT;

  using Get_sqrt = internal::Get_sqrt<GeomTraits>;
  auto sqrt = Get_sqrt::sqrt_object(traits);

  const FT P1 = d * d0 + D0;
  const FT P2 = d * d2 + D2;

  FT w = FT(0);
  CGAL_precondition(!is_zero(P1) && !is_zero(P2));
  const FT prod = P1 * P2;
  if (!is_zero(prod))
  {
    w = FT(2) * (d0 * d2 - D) / prod;
    CGAL_assertion(w >= FT(0));
    w = sqrt(w);
  }

  w *= sign * FT(2);
  return w;
}

} // namespace mean_value_ns

/// \endcond

// 2D ==============================================================================================

/*!
  \ingroup PkgWeightsRefMeanValueWeights
  \brief computes the mean value weight in 2D at `q` using the points `p0`, `p1`, and `p2`.
  \tparam GeomTraits a model of `AnalyticWeightTraits_2`
*/
template<typename GeomTraits>
typename GeomTraits::FT mean_value_weight(const typename GeomTraits::Point_2& p0,
                                          const typename GeomTraits::Point_2& p1,
                                          const typename GeomTraits::Point_2& p2,
                                          const typename GeomTraits::Point_2& q,
                                          const GeomTraits& traits)
{
  using FT = typename GeomTraits::FT;
  using Vector_2 = typename GeomTraits::Vector_2;

  auto vector_2 = traits.construct_vector_2_object();
  auto dot_product_2 = traits.compute_scalar_product_2_object();
  auto area_2 = traits.compute_area_2_object();

  const Vector_2 v1 = vector_2(q, p0);
  const Vector_2 v = vector_2(q, p1);
  const Vector_2 v2 = vector_2(q, p2);

  const FT d0 = internal::length_2(v1, traits);
  const FT d = internal::length_2(v, traits);
  const FT d2 = internal::length_2(v2, traits);

  const FT D0 = dot_product_2(v1, v);
  const FT D2 = dot_product_2(v, v2);
  const FT D = dot_product_2(v1, v2);

  const FT A0 = area_2(p1, q, p0);
  const FT A2 = area_2(p2, q, p1);
  const FT B = area_2(p2, q, p0);

  const FT sign = mean_value_ns::sign_of_weight(A0, A2, B);
  return mean_value_ns::weight(d0, d2, d, D0, D2, D, sign, traits);
}

/*!
  \ingroup PkgWeightsRefMeanValueWeights
  \brief computes the mean value weight in 2D at `q` using the points `p0`, `p1`, and `p2`.
  \tparam Kernel a model of `Kernel`
*/
template<typename Kernel>
typename Kernel::FT mean_value_weight(const CGAL::Point_2<Kernel>& p0,
                                      const CGAL::Point_2<Kernel>& p1,
                                      const CGAL::Point_2<Kernel>& p2,
                                      const CGAL::Point_2<Kernel>& q)
{
  const Kernel traits;
  return mean_value_weight(p0, p1, p2, q, traits);
}

// 3D ==============================================================================================

/// \cond SKIP_IN_MANUAL

template<typename GeomTraits>
typename GeomTraits::FT mean_value_weight(const typename GeomTraits::Point_3& p0,
                                          const typename GeomTraits::Point_3& p1,
                                          const typename GeomTraits::Point_3& p2,
                                          const typename GeomTraits::Point_3& q,
                                          const GeomTraits& traits)
{
  using Point_2 = typename GeomTraits::Point_2;

  Point_2 p0f, p1f, p2f, qf;
  internal::flatten(p0, p1, p2, q,
                    p0f, p1f, p2f, qf,
                    traits);
  return CGAL::Weights::mean_value_weight(p0f, p1f, p2f, qf, traits);
}

template<typename Kernel>
typename Kernel::FT mean_value_weight(const CGAL::Point_3<Kernel>& p0,
                                      const CGAL::Point_3<Kernel>& p1,
                                      const CGAL::Point_3<Kernel>& p2,
                                      const CGAL::Point_3<Kernel>& q)
{
  const Kernel traits;
  return mean_value_weight(p0, p1, p2, q, traits);
}

/// \endcond

/*!
  \ingroup PkgWeightsRefBarycentricMeanValueWeights

  \brief 2D mean value weights for polygons.

  This class implements 2D mean value weights (\cite cgal:bc:fhk-gcbcocp-06, \cite cgal:f-mvc-03,
  \cite cgal:bc:hf-mvcapp-06) which can be computed at any point inside and outside a simple polygon.

  Mean value weights are well-defined inside and outside a simple polygon and are
  non-negative in the kernel of a star-shaped polygon. These weights are computed
  analytically using the formulation from `tangent_weight()`.

  \tparam VertexRange a model of `ConstRange` whose iterator type is `RandomAccessIterator`
  \tparam GeomTraits a model of `AnalyticWeightTraits_2`
  \tparam PointMap a model of `ReadablePropertyMap` whose key type is `VertexRange::value_type` and
                   value type is `Point_2`. The default is `CGAL::Identity_property_map`.

  \cgalModels{BarycentricWeights_2}
*/
template<typename VertexRange,
         typename GeomTraits,
         typename PointMap = CGAL::Identity_property_map<typename GeomTraits::Point_2> >
class Mean_value_weights_2
{
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

    \param polygon an instance of `VertexRange` with the vertices of a simple polygon
    \param traits a traits class with geometric objects, predicates, and constructions;
                  the default initialization is provided
    \param point_map an instance of `PointMap` that maps a vertex from `polygon` to `Point_2`;
                     the default initialization is provided

    \pre polygon.size() >= 3
    \pre polygon is simple
  */
  Mean_value_weights_2(const VertexRange& polygon,
                       const GeomTraits traits = GeomTraits(),
                       const PointMap point_map = PointMap())
    : m_polygon(polygon),
      m_traits(traits),
      m_point_map(point_map),
      m_area_2(m_traits.compute_area_2_object()),
      m_construct_vector_2(m_traits.construct_vector_2_object()),
      m_squared_length_2(m_traits.compute_squared_length_2_object()),
      m_scalar_product_2(m_traits.compute_scalar_product_2_object()),
      m_sqrt(Get_sqrt::sqrt_object(m_traits))
  {
    CGAL_precondition(polygon.size() >= 3);
    CGAL_precondition(internal::is_simple_2(polygon, traits, point_map));
    resize();
  }

  /// @}

  /// \name Access
  /// @{

  /*!
    \brief computes 2D mean value weights.

    This function fills a destination range with 2D mean value weights computed at
    the `query` point with respect to the vertices of the input polygon.

    The number of computed weights is equal to the number of polygon vertices.

    \tparam OutIterator a model of `OutputIterator` whose value type is `FT`

    \param query a query point
    \param w_begin the beginning of the destination range with the computed weights

    \return an output iterator to the element in the destination range, one past the last weight stored
  */
  template<typename OutIterator>
  OutIterator operator()(const Point_2& query, OutIterator w_begin)
  {
    const bool normalize = false;
    return operator()(query, w_begin, normalize);
  }

  /// @}

  /// \cond SKIP_IN_MANUAL

  template<typename OutIterator>
  OutIterator operator()(const Point_2& query,
                         OutIterator weights,
                         const bool normalize)
  {
    return optimal_weights(query, weights, normalize);
  }

  /// \endcond

private:
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

  void resize()
  {
    s.resize(m_polygon.size());
    r.resize(m_polygon.size());
    A.resize(m_polygon.size());
    D.resize(m_polygon.size());
    t.resize(m_polygon.size());
    w.resize(m_polygon.size());
  }

  template<typename OutputIterator>
  OutputIterator optimal_weights(const Point_2& query,
                                 OutputIterator weights,
                                 const bool normalize)
  {
    const std::size_t n = m_polygon.size();

    // Compute vectors s following the pseudo-code in the Figure 10 from [1].
    for (std::size_t i = 0; i < n; ++i)
    {
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

    for (std::size_t i = 1; i < n - 1; ++i)
    {
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
    for (std::size_t i = 0; i < n - 1; ++i)
    {
      CGAL_assertion((r[i] * r[i + 1] + D[i]) != FT(0));
      t[i] = FT(2) * A[i] / (r[i] * r[i + 1] + D[i]);
    }

    CGAL_assertion((r[n - 1] * r[0] + D[n - 1]) != FT(0));
    t[n - 1] = FT(2) * A[n - 1] / (r[n - 1] * r[0] + D[n - 1]);

    // Compute mean value weights using the same pseudo-code as before.
    CGAL_assertion(r[0] != FT(0));
    w[0] = FT(2) * (t[n - 1] + t[0]) / r[0];

    for (std::size_t i = 1; i < n - 1; ++i)
    {
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

/*!
  \ingroup PkgWeightsRefBarycentricMeanValueWeights

  \brief computes 2D mean value weights for polygons.

  This function computes 2D mean value weights at a given `query` point
  with respect to the vertices of a simple `polygon`, that is one
  weight per vertex. The weights are stored in a destination range
  beginning at `w_begin`.

  Internally, the class `Mean_value_weights_2` is used. If one wants to process
  multiple query points, it is better to use that class. When using the free function,
  internal memory is allocated for each query point, while when using the class,
  it is allocated only once which is much more efficient. However, for a few query
  points, it is easier to use this function. It can also be used when the processing
  time is not a concern.

  \tparam PointRange a model of `ConstRange` whose iterator type is `RandomAccessIterator`
                     and value type is `GeomTraits::Point_2`
  \tparam OutIterator a model of `OutputIterator` whose value type is `GeomTraits::FT`
  \tparam GeomTraits a model of `AnalyticWeightTraits_2`

  \param polygon an instance of `PointRange` with 2D points which form a simple polygon
  \param query a query point
  \param w_begin the beginning of the destination range with the computed weights
  \param traits a traits class with geometric objects, predicates, and constructions;
                this parameter can be omitted if the traits class can be deduced from the point type

  \return an output iterator to the element in the destination range, one past the last weight stored

  \pre `polygon.size() >= 3`
  \pre `polygon` is simple
*/
template<typename PointRange,
         typename OutIterator,
         typename GeomTraits>
OutIterator mean_value_weights_2(const PointRange& polygon,
                                 const typename GeomTraits::Point_2& query,
                                 OutIterator w_begin,
                                 const GeomTraits& traits)
{
  Mean_value_weights_2<PointRange, GeomTraits> mean_value(polygon, traits);
  return mean_value(query, w_begin);
}

/// \cond SKIP_IN_MANUAL

template<typename PointRange,
         typename OutIterator>
OutIterator mean_value_weights_2(const PointRange& polygon,
                                 const typename PointRange::value_type& query,
                                 OutIterator w_begin)
{
  using Point_2 = typename PointRange::value_type;
  using GeomTraits = typename Kernel_traits<Point_2>::Kernel;

  const GeomTraits traits;
  return mean_value_weights_2(polygon, query, w_begin, traits);
}

/// \endcond

} // namespace Weights
} // namespace CGAL

#endif // CGAL_MEAN_VALUE_WEIGHTS_H
