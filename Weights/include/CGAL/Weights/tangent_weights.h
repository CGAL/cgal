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
// Author(s)     : Dmitry Anisimov
//

#ifndef CGAL_TANGENT_WEIGHTS_H
#define CGAL_TANGENT_WEIGHTS_H

#include <CGAL/Weights/utils.h>

#include <CGAL/Point_2.h>
#include <CGAL/Point_3.h>

#include <cmath>

namespace CGAL {
namespace Weights {

/// \cond SKIP_IN_MANUAL

namespace tangent_ns {

template<typename FT>
FT half_weight(const FT t, const FT r)
{
  FT w = FT(0);
  CGAL_precondition(!is_zero(r));
  if (!is_zero(r))
    w = FT(2) * t / r;

  return w;
}

template<typename FT>
FT weight(const FT t0, const FT t2, const FT r)
{
  FT w = FT(0);
  CGAL_precondition(r != FT(0));
  if (r != FT(0))
    w = FT(2) * (t0 + t2) / r;

  return w;
}

template<typename FT>
FT weight(const FT d0, const FT d2, const FT d,
          const FT A0, const FT A2,
          const FT D0, const FT D2)
{
  const FT P0 = d * d0 + D0;
  const FT P2 = d * d2 + D2;

  FT w = FT(0);
  CGAL_precondition(!is_zero(P0) && !is_zero(P2));
  if (!is_zero(P0) && !is_zero(P2))
  {
    const FT t0 = FT(2) * A0 / P0;
    const FT t2 = FT(2) * A2 / P2;
    w = weight(t0, t2, d);
  }

  return w;
}

// This is positive case only.
// This version is based on the positive area.
// This version is more precise for all positive cases.
template<typename GeomTraits>
typename GeomTraits::FT tangent_weight_v1(const typename GeomTraits::Point_3& p0,
                                          const typename GeomTraits::Point_3& p1,
                                          const typename GeomTraits::Point_3& p2,
                                          const typename GeomTraits::Point_3& q,
                                          const GeomTraits& traits)
{
  using FT = typename GeomTraits::FT;
  using Vector_3 = typename GeomTraits::Vector_3;

  auto dot_product_3 = traits.compute_scalar_product_3_object();
  auto vector_3 = traits.construct_vector_3_object();

  const Vector_3 v0 = vector_3(q, p0);
  const Vector_3 v = vector_3(q, p1);
  const Vector_3 v2 = vector_3(q, p2);

  const FT d0 = internal::length_3(v0, traits);
  const FT d = internal::length_3(v, traits);
  const FT d2 = internal::length_3(v2, traits);

  const FT A0 = internal::positive_area_3(p1, q, p0, traits);
  const FT A2 = internal::positive_area_3(p2, q, p1, traits);

  const FT D0 = dot_product_3(v0, v);
  const FT D2 = dot_product_3(v, v2);

  return weight(d0, d2, d, A0, A2, D0, D2);
}

// This version handles both positive and negative cases.
// However, it is less precise.
template<typename GeomTraits>
typename GeomTraits::FT tangent_weight_v2(const typename GeomTraits::Point_3& p0,
                                          const typename GeomTraits::Point_3& p1,
                                          const typename GeomTraits::Point_3& p2,
                                          const typename GeomTraits::Point_3& q,
                                          const GeomTraits& traits)
{
  using FT = typename GeomTraits::FT;
  using Vector_3 = typename GeomTraits::Vector_3;

  auto vector_3 = traits.construct_vector_3_object();

  Vector_3 v0 = vector_3(q, p0);
  Vector_3 v = vector_3(q, p1);
  Vector_3 v2 = vector_3(q, p2);

  const FT l2 = internal::length_3(v, traits);

  const double ha_rad_1 = internal::angle_3(v0, v, traits) / 2.0;
  const double ha_rad_2 = internal::angle_3(v, v2, traits) / 2.0;
  const FT t0 = static_cast<FT>(std::tan(ha_rad_1));
  const FT t2 = static_cast<FT>(std::tan(ha_rad_2));

  return weight(t0, t2, l2);
}

} // namespace tangent_ns

/// \endcond

// 2D ==============================================================================================

/*!
  \ingroup PkgWeightsRefTangentWeights

  \brief computes the half value of the tangent weight.

  This function constructs the half of the tangent weight using the precomputed
  half angle tangent and distance values. The returned value is
  \f$\frac{2\textbf{tan05}}{\textbf{d}}\f$.

  \tparam FT a model of `FieldNumberType`

  \param tan05 the half angle tangent value
  \param d the distance value

  \pre d != 0

  \sa `tangent_half_angle()`
  \sa `tangent_weight()`
*/
template<typename FT>
FT half_tangent_weight(const FT tan05, const FT d)
{
  return tangent_ns::half_weight(tan05, d);
}

/*!
  \ingroup PkgWeightsRefTangentWeights

  \brief computes the tangent of the half angle.

  This function computes the tangent of the half angle using the precomputed
  distance, area, and dot product values. The returned value is
  \f$\frac{2\textbf{A}}{\textbf{d}\textbf{l} + \textbf{D}}\f$.

  \tparam FT a model of `FieldNumberType`

  \param d the distance value
  \param l the distance value
  \param A the area value
  \param D the dot product value

  \pre (d * l + D) != 0

  \sa `half_tangent_weight()`
*/
template<typename FT>
FT tangent_half_angle(const FT d, const FT l, const FT A, const FT D)
{
  // tan(theta/2) = sin(theta) / ( 1 + cos(theta) ), also = (1 - cos(theta)) / sin(theta).
  //              = ( 2*A / |v1|*|v2| ) / ( 1 + v1.v2 / |v1|*|v2| )
  //              = 2*A / ( |v1|*|v2| + v1.v2 )

  FT t = FT(0);
  const FT P = d * l + D;
  CGAL_precondition(!is_zero(P));
  if (!is_zero(P))
    t = FT(2) * A / P;

  return t;
}

/*!
  \ingroup PkgWeightsRefTangentWeights

  \brief computes the half value of the tangent weight.

  This function constructs the half of the tangent weight using the precomputed
  distance, area, and dot product values. The returned value is
  \f$\frac{2\textbf{t}}{\textbf{d}}\f$ where
  \f$\textbf{t} = \frac{2\textbf{A}}{\textbf{d}\textbf{l} + \textbf{D}}\f$.

  \tparam FT a model of `FieldNumberType`

  \param d the distance value
  \param l the distance value
  \param A the area value
  \param D the dot product value

  \pre (d * l + D) != 0 && d != 0

  \sa `tangent_weight()`
*/
template<typename FT>
FT half_tangent_weight(const FT d, const FT l, const FT A, const FT D)
{
  const FT tan05 = tangent_half_angle(d, l, A, D);
  return tangent_ns::half_weight(tan05, d);
}

/// \cond SKIP_IN_MANUAL

template<typename GeomTraits>
typename GeomTraits::FT half_tangent_weight(const typename GeomTraits::Point_2& p0,
                                            const typename GeomTraits::Point_2& q,
                                            const typename GeomTraits::Point_2& p2,
                                            const GeomTraits& traits)
{
  using FT = typename GeomTraits::FT;
  using Vector_2 = typename GeomTraits::Vector_2;

  auto vector_2 = traits.construct_vector_2_object();
  auto dot_product_2 = traits.compute_scalar_product_2_object();
  auto area_2 = traits.compute_area_2_object();

  const Vector_2 v0 = vector_2(q, p0);
  const Vector_2 v2 = vector_2(q, p2);

  const FT l0 = internal::length_2(v0, traits);
  const FT l2 = internal::length_2(v2, traits);
  const FT A = area_2(p2, q, p0);
  const FT D = dot_product_2(v0, v2);

  return half_tangent_weight(l0, l2, A, D);
}

template<typename GeomTraits>
typename GeomTraits::FT half_tangent_weight(const typename GeomTraits::Point_3& p0,
                                            const typename GeomTraits::Point_3& q,
                                            const typename GeomTraits::Point_3& p2,
                                            const GeomTraits& traits)
{
  using FT = typename GeomTraits::FT;
  using Vector_3 = typename GeomTraits::Vector_3;

  auto vector_3 = traits.construct_vector_3_object();
  auto dot_product_3 = traits.compute_scalar_product_3_object();

  const Vector_3 v0 = vector_3(q, p0);
  const Vector_3 v2 = vector_3(q, p2);

  const FT l0 = internal::length_3(v0, traits);
  const FT l2 = internal::length_3(v2, traits);
  const FT A = internal::area_3(p2, q, p0, traits);
  const FT D = dot_product_3(v0, v2);

  return half_tangent_weight(l0, l2, A, D);
}

/// \endcond

// 2D ==============================================================================================

/*!
  \ingroup PkgWeightsRefTangentWeights
  \brief computes the tangent weight in 2D at `q` using the points `p0`, `p1`, and `p2`
  \tparam GeomTraits a model of `AnalyticWeightTraits_2`
*/
template<typename GeomTraits>
typename GeomTraits::FT tangent_weight(const typename GeomTraits::Point_2& p0,
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

  const Vector_2 v0 = vector_2(q, p0);
  const Vector_2 v = vector_2(q, p1);
  const Vector_2 v2 = vector_2(q, p2);

  const FT l0 = internal::length_2(v0, traits);
  const FT l = internal::length_2(v, traits);
  const FT l2 = internal::length_2(v2, traits);

  const FT A0 = area_2(p1, q, p0);
  const FT A2 = area_2(p2, q, p1);

  const FT D0 = dot_product_2(v0, v);
  const FT D2 = dot_product_2(v, v2);

  return tangent_ns::weight(l0, l2, l, A0, A2, D0, D2);
}

/*!
  \ingroup PkgWeightsRefTangentWeights
  \brief computes the tangent weight in 2D at `q` using the points `p0`, `p1`, and `p2`
  \tparam Kernel a model of `Kernel`
*/
template<typename Kernel>
typename Kernel::FT tangent_weight(const CGAL::Point_2<Kernel>& p0,
                                   const CGAL::Point_2<Kernel>& p1,
                                   const CGAL::Point_2<Kernel>& p2,
                                   const CGAL::Point_2<Kernel>& q)
{
  const Kernel traits;
  return tangent_weight(p0, p1, p2, q, traits);
}

// 3D ==============================================================================================

/*!
  \ingroup PkgWeightsRefTangentWeights
  \brief computes the tangent weight in 3D at `q` using the points `p0`, `p1`, and `p2`
  \tparam GeomTraits a model of `AnalyticWeightTraits_3`
*/
template<typename GeomTraits>
typename GeomTraits::FT tangent_weight(const typename GeomTraits::Point_3& p0,
                                       const typename GeomTraits::Point_3& p1,
                                       const typename GeomTraits::Point_3& p2,
                                       const typename GeomTraits::Point_3& q,
                                       const GeomTraits& traits)
{
//  return tangent_ns::tangent_weight_v1(p0, p1, p2, q, traits);
  return tangent_ns::tangent_weight_v2(p0, p1, p2, q, traits);
}

/*!
  \ingroup PkgWeightsRefTangentWeights
  \brief computes the tangent weight in 3D at `q` using the points `p0`, `p1`, and `p2`
  \tparam Kernel a model of `Kernel`
*/
template<typename Kernel>
typename Kernel::FT tangent_weight(const CGAL::Point_3<Kernel>& p0,
                                   const CGAL::Point_3<Kernel>& p1,
                                   const CGAL::Point_3<Kernel>& p2,
                                   const CGAL::Point_3<Kernel>& q)
{
  const Kernel traits;
  return tangent_weight(p0, p1, p2, q, traits);
}

/// \cond SKIP_IN_MANUAL

// Undocumented tangent weight class.
//
// Its constructor takes a polygon mesh and a vertex to point map
// and its operator() is defined based on the halfedge_descriptor only.
// This version is currently used in:
// Surface_mesh_parameterizer -> Iterative_authalic_parameterizer_3.h
template<
    typename PolygonMesh,
    typename VertexPointMap,
    typename GeomTraits>
class Edge_tangent_weight
{
  using vertex_descriptor = typename boost::graph_traits<PolygonMesh>::vertex_descriptor;
  using halfedge_descriptor = typename boost::graph_traits<PolygonMesh>::halfedge_descriptor;

  using Point_ref = typename boost::property_traits<VertexPointMap>::reference;
  using FT = typename GeomTraits::FT;

  const PolygonMesh& m_pmesh;
  const VertexPointMap m_pmap;
  const GeomTraits m_traits;

public:

  Edge_tangent_weight(const PolygonMesh& pmesh,
                      const VertexPointMap pmap,
                      const GeomTraits& traits)
    : m_pmesh(pmesh), m_pmap(pmap), m_traits(traits)
  { }

  FT operator()(halfedge_descriptor he) const
  {
    if(is_border(he, m_pmesh))
      return FT(0);

    FT weight = FT(0);
    if (is_border_edge(he, m_pmesh)) // ie, opp(he, pmesh) is a border halfedge
    {
      const halfedge_descriptor h1 = next(he, m_pmesh);

      const vertex_descriptor v0 = target(he, m_pmesh);
      const vertex_descriptor v1 = source(he, m_pmesh);
      const vertex_descriptor v2 = target(h1, m_pmesh);

      const Point_ref p0 = get(m_pmap, v0);
      const Point_ref p1 = get(m_pmap, v1);
      const Point_ref p2 = get(m_pmap, v2);

      weight = half_tangent_weight(p1, p0, p2, m_traits) / FT(2);
    }
    else
    {
      const halfedge_descriptor h1 = next(he, m_pmesh);
      const halfedge_descriptor h2 = prev(opposite(he, m_pmesh), m_pmesh);

      const vertex_descriptor v0 = target(he, m_pmesh);
      const vertex_descriptor v1 = source(he, m_pmesh);
      const vertex_descriptor v2 = target(h1, m_pmesh);
      const vertex_descriptor v3 = source(h2, m_pmesh);

      const Point_ref p0 = get(m_pmap, v0);
      const Point_ref p1 = get(m_pmap, v1);
      const Point_ref p2 = get(m_pmap, v2);
      const Point_ref p3 = get(m_pmap, v3);

      weight = tangent_weight(p2, p1, p3, p0, m_traits) / FT(2);
    }
    return weight;
  }
};

// Undocumented tangent weight class.
// Returns - std::tan(theta/2); uses positive areas.
//
// Its constructor takes three points either in 2D or 3D.
// This version is currently used in:
// Surface_mesh_parameterizer -> MVC_post_processor_3.h
// Surface_mesh_parameterizer -> Orbifold_Tutte_parameterizer_3.h
template<typename FT>
class Tangent_weight
{
  FT m_d_r, m_d_p, m_w_base;

public:
  template<typename Kernel>
  Tangent_weight(const CGAL::Point_2<Kernel>& p,
                 const CGAL::Point_2<Kernel>& q,
                 const CGAL::Point_2<Kernel>& r)
  {
    const Kernel traits;

    using Vector_2 = typename Kernel::Vector_2;

    auto vector_2 = traits.construct_vector_2_object();
    auto scalar_product_2 = traits.compute_scalar_product_2_object();

    m_d_r = internal::distance_2(q, r, traits);
    CGAL_assertion(is_positive(m_d_r)); // two points are identical!
    m_d_p = internal::distance_2(q, p, traits);
    CGAL_assertion(is_positive(m_d_p)); // two points are identical!

    const Vector_2 v1 = vector_2(q, r);
    const Vector_2 v2 = vector_2(q, p);

    const FT A = internal::positive_area_2(p, q, r, traits);
    CGAL_assertion(!is_zero(A));

    const FT S = scalar_product_2(v1, v2);
    m_w_base = -tangent_half_angle(m_d_r, m_d_p, A, S);
  }

  template<typename Kernel>
  Tangent_weight(const CGAL::Point_3<Kernel>& p,
                 const CGAL::Point_3<Kernel>& q,
                 const CGAL::Point_3<Kernel>& r)
  {
    const Kernel traits;

    using Vector_3 = typename Kernel::Vector_3;

    auto vector_3 = traits.construct_vector_3_object();
    auto scalar_product_3 = traits.compute_scalar_product_3_object();

    m_d_r = internal::distance_3(q, r, traits);
    CGAL_assertion(is_positive(m_d_r)); // two points are identical!
    m_d_p = internal::distance_3(q, p, traits);
    CGAL_assertion(is_positive(m_d_p)); // two points are identical!

    const Vector_3 v1 = vector_3(q, r);
    const Vector_3 v2 = vector_3(q, p);

    const FT A = internal::positive_area_3(p, q, r, traits);
    CGAL_assertion(is_positive(A));

    const FT S = scalar_product_3(v1, v2);
    m_w_base = -tangent_half_angle(m_d_r, m_d_p, A, S);
  }

  FT get_w_r() const
  {
    return half_tangent_weight(m_w_base, m_d_r) / FT(2);
  }

  FT get_w_p() const
  {
    return half_tangent_weight(m_w_base, m_d_p) / FT(2);
  }
};

/// \endcond

} // namespace Weights
} // namespace CGAL

#endif // CGAL_TANGENT_WEIGHTS_H
