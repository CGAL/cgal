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

// Internal includes.
#include <CGAL/Weights/internal/utils.h>

namespace CGAL {
namespace Weights {

  /// \cond SKIP_IN_MANUAL
  namespace tangent_ns {

    template<typename FT>
    FT half_angle_tangent(const FT r, const FT d, const FT A, const FT D) {

      FT t = FT(0);
      const FT P = r * d + D;
      CGAL_precondition(P != FT(0));
      if (P != FT(0)) {
        const FT inv = FT(2) / P;
        t = A * inv;
      }
      return t;
    }

    template<typename FT>
    FT half_weight(const FT t, const FT r) {

      FT w = FT(0);
      CGAL_precondition(r != FT(0));
      if (r != FT(0)) {
        const FT inv = FT(2) / r;
        w = t * inv;
      }
      return w;
    }

    template<typename FT>
    FT weight(const FT t1, const FT t2, const FT r) {

      FT w = FT(0);
      CGAL_precondition(r != FT(0));
      if (r != FT(0)) {
        const FT inv = FT(2) / r;
        w = (t1 + t2) * inv;
      }
      return w;
    }

    template<typename FT>
    FT weight(
      const FT d1, const FT r, const FT d2,
      const FT A1, const FT A2,
      const FT D1, const FT D2) {

      const FT P1 = d1 * r + D1;
      const FT P2 = d2 * r + D2;

      FT w = FT(0);
      CGAL_precondition(P1 != FT(0) && P2 != FT(0));
      if (P1 != FT(0) && P2 != FT(0)) {
        const FT inv1 = FT(2) / P1;
        const FT inv2 = FT(2) / P2;
        const FT t1 = A1 * inv1;
        const FT t2 = A2 * inv2;
        w = weight(t1, t2, r);
      }
      return w;
    }

    // This is positive case only.
    // This version is based on the positive area.
    // This version is more precise for all positive cases.
    template<typename GeomTraits>
    typename GeomTraits::FT tangent_weight_v1(
      const typename GeomTraits::Point_3& t,
      const typename GeomTraits::Point_3& r,
      const typename GeomTraits::Point_3& p,
      const typename GeomTraits::Point_3& q,
      const GeomTraits& traits) {

      using FT = typename GeomTraits::FT;
      const auto dot_product_3 =
        traits.compute_scalar_product_3_object();
      const auto construct_vector_3 =
        traits.construct_vector_3_object();

      const auto v1 = construct_vector_3(q, t);
      const auto v2 = construct_vector_3(q, r);
      const auto v3 = construct_vector_3(q, p);

      const FT l1 = internal::length_3(traits, v1);
      const FT l2 = internal::length_3(traits, v2);
      const FT l3 = internal::length_3(traits, v3);

      const FT A1 = internal::positive_area_3(traits, r, q, t);
      const FT A2 = internal::positive_area_3(traits, p, q, r);

      const FT D1 = dot_product_3(v1, v2);
      const FT D2 = dot_product_3(v2, v3);

      return weight(l1, l2, l3, A1, A2, D1, D2);
    }

    // This version handles both positive and negative cases.
    // However, it is less precise.
    template<typename GeomTraits>
    typename GeomTraits::FT tangent_weight_v2(
      const typename GeomTraits::Point_3& t,
      const typename GeomTraits::Point_3& r,
      const typename GeomTraits::Point_3& p,
      const typename GeomTraits::Point_3& q,
      const GeomTraits& traits) {

      using FT = typename GeomTraits::FT;
      const auto construct_vector_3 =
        traits.construct_vector_3_object();

      auto v1 = construct_vector_3(q, t);
      auto v2 = construct_vector_3(q, r);
      auto v3 = construct_vector_3(q, p);

      const FT l2 = internal::length_3(traits, v2);

      internal::normalize_3(traits, v1);
      internal::normalize_3(traits, v2);
      internal::normalize_3(traits, v3);

      const double ha_rad_1 = internal::angle_3(traits, v1, v2) / 2.0;
      const double ha_rad_2 = internal::angle_3(traits, v2, v3) / 2.0;
      const FT t1 = static_cast<FT>(std::tan(ha_rad_1));
      const FT t2 = static_cast<FT>(std::tan(ha_rad_2));

      return weight(t1, t2, l2);
    }
  }
  /// \endcond

  /*!
    \ingroup PkgWeightsRefTangentWeights

    \brief computes the tangent of the half angle.

    This function computes the tangent of the half angle using the precomputed
    distance, area, and dot product values. The returned value is
    \f$\frac{2\textbf{A}}{\textbf{d}\textbf{l} + \textbf{D}}\f$.

    \tparam FT
    a model of `FieldNumberType`

    \param d
    the distance value

    \param l
    the distance value

    \param A
    the area value

    \param D
    the dot product value

    \pre (d * l + D) != 0

    \sa `half_tangent_weight()`
  */
  template<typename FT>
  FT tangent_half_angle(const FT d, const FT l, const FT A, const FT D) {
    return tangent_ns::half_angle_tangent(d, l, A, D);
  }

  /*!
    \ingroup PkgWeightsRefTangentWeights

    \brief computes the half value of the tangent weight.

    This function constructs the half of the tangent weight using the precomputed
    half angle tangent and distance values. The returned value is
    \f$\frac{2\textbf{tan05}}{\textbf{d}}\f$.

    \tparam FT
    a model of `FieldNumberType`

    \param tan05
    the half angle tangent value

    \param d
    the distance value

    \pre d != 0

    \sa `tangent_half_angle()`
    \sa `tangent_weight()`
  */
  template<typename FT>
  FT half_tangent_weight(const FT tan05, const FT d) {
    return tangent_ns::half_weight(tan05, d);
  }

  /*!
    \ingroup PkgWeightsRefTangentWeights

    \brief computes the half value of the tangent weight.

    This function constructs the half of the tangent weight using the precomputed
    distance, area, and dot product values. The returned value is
    \f$\frac{2\textbf{t}}{\textbf{d}}\f$ where
    \f$\textbf{t} = \frac{2\textbf{A}}{\textbf{d}\textbf{l} + \textbf{D}}\f$.

    \tparam FT
    a model of `FieldNumberType`

    \param d
    the distance value

    \param l
    the distance value

    \param A
    the area value

    \param D
    the dot product value

    \pre (d * l + D) != 0 && d != 0

    \sa `tangent_weight()`
  */
  template<typename FT>
  FT half_tangent_weight(const FT d, const FT l, const FT A, const FT D) {
    const FT tan05 = tangent_half_angle(d, l, A, D);
    return half_tangent_weight(tan05, d);
  }

  #if defined(DOXYGEN_RUNNING)

  /*!
    \ingroup PkgWeightsRefTangentWeights

    \brief computes the tangent weight in 2D at `q` using the points `p0`, `p1`,
    and `p2`, given a traits class `traits` with geometric objects, predicates, and constructions.
  */
  template<typename GeomTraits>
  typename GeomTraits::FT tangent_weight(
    const typename GeomTraits::Point_2& p0,
    const typename GeomTraits::Point_2& p1,
    const typename GeomTraits::Point_2& p2,
    const typename GeomTraits::Point_2& q,
    const GeomTraits& traits) { }

  /*!
    \ingroup PkgWeightsRefTangentWeights

    \brief computes the tangent weight in 3D at `q` using the points `p0`, `p1`,
    and `p2`, given a traits class `traits` with geometric objects, predicates, and constructions.
  */
  template<typename GeomTraits>
  typename GeomTraits::FT tangent_weight(
    const typename GeomTraits::Point_3& p0,
    const typename GeomTraits::Point_3& p1,
    const typename GeomTraits::Point_3& p2,
    const typename GeomTraits::Point_3& q,
    const GeomTraits& traits) { }

  /*!
    \ingroup PkgWeightsRefTangentWeights

    \brief computes the tangent weight in 2D at `q` using the points `p0`, `p1`,
    and `p2` which are parameterized by a `Kernel` K.
  */
  template<typename K>
  typename K::FT tangent_weight(
    const CGAL::Point_2<K>& p0,
    const CGAL::Point_2<K>& p1,
    const CGAL::Point_2<K>& p2,
    const CGAL::Point_2<K>& q) { }

  /*!
    \ingroup PkgWeightsRefTangentWeights

    \brief computes the tangent weight in 3D at `q` using the points `p0`, `p1`,
    and `p2` which are parameterized by a `Kernel` K.
  */
  template<typename K>
  typename K::FT tangent_weight(
    const CGAL::Point_3<K>& p0,
    const CGAL::Point_3<K>& p1,
    const CGAL::Point_3<K>& p2,
    const CGAL::Point_3<K>& q) { }

  #endif // DOXYGEN_RUNNING

  /// \cond SKIP_IN_MANUAL
  template<typename GeomTraits>
  typename GeomTraits::FT tangent_weight(
    const typename GeomTraits::Point_2& t,
    const typename GeomTraits::Point_2& r,
    const typename GeomTraits::Point_2& p,
    const typename GeomTraits::Point_2& q,
    const GeomTraits& traits) {

    using FT = typename GeomTraits::FT;
    const auto dot_product_2 =
      traits.compute_scalar_product_2_object();
    const auto construct_vector_2 =
      traits.construct_vector_2_object();

    const auto v1 = construct_vector_2(q, t);
    const auto v2 = construct_vector_2(q, r);
    const auto v3 = construct_vector_2(q, p);

    const FT l1 = internal::length_2(traits, v1);
    const FT l2 = internal::length_2(traits, v2);
    const FT l3 = internal::length_2(traits, v3);

    const FT A1 = internal::area_2(traits, r, q, t);
    const FT A2 = internal::area_2(traits, p, q, r);

    const FT D1 = dot_product_2(v1, v2);
    const FT D2 = dot_product_2(v2, v3);

    return tangent_ns::weight(
      l1, l2, l3, A1, A2, D1, D2);
  }

  template<typename GeomTraits>
  typename GeomTraits::FT tangent_weight(
    const CGAL::Point_2<GeomTraits>& t,
    const CGAL::Point_2<GeomTraits>& r,
    const CGAL::Point_2<GeomTraits>& p,
    const CGAL::Point_2<GeomTraits>& q) {

    const GeomTraits traits;
    return tangent_weight(t, r, p, q, traits);
  }

  template<typename GeomTraits>
  typename GeomTraits::FT tangent_weight(
    const typename GeomTraits::Point_3& t,
    const typename GeomTraits::Point_3& r,
    const typename GeomTraits::Point_3& p,
    const typename GeomTraits::Point_3& q,
    const GeomTraits& traits) {

    // return tangent_ns::tangent_weight_v1(t, r, p, q, traits);
       return tangent_ns::tangent_weight_v2(t, r, p, q, traits);
  }

  template<typename GeomTraits>
  typename GeomTraits::FT tangent_weight(
    const CGAL::Point_3<GeomTraits>& t,
    const CGAL::Point_3<GeomTraits>& r,
    const CGAL::Point_3<GeomTraits>& p,
    const CGAL::Point_3<GeomTraits>& q) {

    const GeomTraits traits;
    return tangent_weight(t, r, p, q, traits);
  }

  // Undocumented tangent weight class.
  // Its constructor takes a polygon mesh and a vertex to point map
  // and its operator() is defined based on the halfedge_descriptor only.
  // This version is currently used in:
  // Surface_mesh_parameterizer -> Iterative_authalic_parameterizer_3.h
  template<
  typename PolygonMesh,
  typename VertexPointMap = typename boost::property_map<PolygonMesh, CGAL::vertex_point_t>::type>
  class Edge_tangent_weight {

    using GeomTraits = typename CGAL::Kernel_traits<
        typename boost::property_traits<VertexPointMap>::value_type>::type;
    using FT = typename GeomTraits::FT;

    const PolygonMesh& m_pmesh;
    const VertexPointMap m_pmap;
    const GeomTraits m_traits;

  public:
    using vertex_descriptor = typename boost::graph_traits<PolygonMesh>::vertex_descriptor;
    using halfedge_descriptor = typename boost::graph_traits<PolygonMesh>::halfedge_descriptor;

    Edge_tangent_weight(const PolygonMesh& pmesh, const VertexPointMap pmap) :
    m_pmesh(pmesh), m_pmap(pmap), m_traits() { }

    FT operator()(const halfedge_descriptor he) const {

      FT weight = FT(0);
      if (is_border_edge(he, m_pmesh)) {
        const auto h1 = next(he, m_pmesh);

        const auto v0 = target(he, m_pmesh);
        const auto v1 = source(he, m_pmesh);
        const auto v2 = target(h1, m_pmesh);

        const auto& p0 = get(m_pmap, v0);
        const auto& p1 = get(m_pmap, v1);
        const auto& p2 = get(m_pmap, v2);

        weight = internal::tangent_3(m_traits, p0, p2, p1);

      } else {
        const auto h1 = next(he, m_pmesh);
        const auto h2 = prev(opposite(he, m_pmesh), m_pmesh);

        const auto v0 = target(he, m_pmesh);
        const auto v1 = source(he, m_pmesh);
        const auto v2 = target(h1, m_pmesh);
        const auto v3 = source(h2, m_pmesh);

        const auto& p0 = get(m_pmap, v0);
        const auto& p1 = get(m_pmap, v1);
        const auto& p2 = get(m_pmap, v2);
        const auto& p3 = get(m_pmap, v3);

        weight = tangent_weight(p2, p1, p3, p0) / FT(2);
      }
      return weight;
    }
  };

  // Undocumented tangent weight class.
  // Its constructor takes three points either in 2D or 3D.
  // This version is currently used in:
  // Surface_mesh_parameterizer -> MVC_post_processor_3.h
  // Surface_mesh_parameterizer -> Orbifold_Tutte_parameterizer_3.h
  template<typename FT>
  class Tangent_weight {
    FT m_d_r, m_d_p, m_w_base;

  public:
    template<typename GeomTraits>
    Tangent_weight(
      const CGAL::Point_2<GeomTraits>& p,
      const CGAL::Point_2<GeomTraits>& q,
      const CGAL::Point_2<GeomTraits>& r) {

      const GeomTraits traits;
      const auto scalar_product_2 =
        traits.compute_scalar_product_2_object();
      const auto construct_vector_2 =
        traits.construct_vector_2_object();

      m_d_r = internal::distance_2(traits, q, r);
      CGAL_assertion(m_d_r != FT(0)); // two points are identical!
      m_d_p = internal::distance_2(traits, q, p);
      CGAL_assertion(m_d_p != FT(0)); // two points are identical!

      const auto v1 = construct_vector_2(q, r);
      const auto v2 = construct_vector_2(q, p);

      const auto A = internal::area_2(traits, p, q, r);
      CGAL_assertion(A != FT(0)); // three points are identical!
      const auto S = scalar_product_2(v1, v2);
      m_w_base = -tangent_half_angle(m_d_r, m_d_p, A, S);
    }

    template<typename GeomTraits>
    Tangent_weight(
      const CGAL::Point_3<GeomTraits>& p,
      const CGAL::Point_3<GeomTraits>& q,
      const CGAL::Point_3<GeomTraits>& r) {

      const GeomTraits traits;
      const auto scalar_product_3 =
        traits.compute_scalar_product_3_object();
      const auto construct_vector_3 =
        traits.construct_vector_3_object();

      m_d_r = internal::distance_3(traits, q, r);
      CGAL_assertion(m_d_r != FT(0)); // two points are identical!
      m_d_p = internal::distance_3(traits, q, p);
      CGAL_assertion(m_d_p != FT(0)); // two points are identical!

      const auto v1 = construct_vector_3(q, r);
      const auto v2 = construct_vector_3(q, p);

      const auto A = internal::positive_area_3(traits, p, q, r);
      CGAL_assertion(A != FT(0)); // three points are identical!
      const auto S = scalar_product_3(v1, v2);
      m_w_base = -tangent_half_angle(m_d_r, m_d_p, A, S);
    }

    FT get_w_r() const {
      return half_tangent_weight(m_w_base, m_d_r) / FT(2);
    }

    FT get_w_p() const {
      return half_tangent_weight(m_w_base, m_d_p) / FT(2);
    }
  };

  /// \endcond

} // namespace Weights
} // namespace CGAL

#endif // CGAL_TANGENT_WEIGHTS_H
