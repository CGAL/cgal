// Copyright (c) 2020 GeometryFactory SARL (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Dmitry Anisimov
//

#ifndef CGAL_WEIGHTS_INTERNAL_TOOLS_H
#define CGAL_WEIGHTS_INTERNAL_TOOLS_H

// #include <CGAL/license/Weights.h>

// CGAL and Boost includes.
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Kernel/global_functions.h>

// Internal includes.
#include <CGAL/Weights.h>

namespace CGAL {
namespace Weights {
namespace internal {

// Computes a secure cotanget between two 3D vectors.
template<typename Point_3>
decltype(auto) cotangent_3_secure(
  const Point_3& p, const Point_3& q, const Point_3& r) {

  using GeomTraits = typename Kernel_traits<Point_3>::Kernel;
  using FT = typename GeomTraits::FT;
  using Get_sqrt = Get_sqrt<GeomTraits>;
  const GeomTraits traits;
  const auto sqrt = Get_sqrt::sqrt_object(traits);

  const auto dot_product_3 =
    traits.compute_scalar_product_3_object();
  const auto construct_vector_3 =
    traits.construct_vector_3_object();

  const auto v1 = construct_vector_3(q, r);
  const auto v2 = construct_vector_3(q, p);

  const FT dot = dot_product_3(v1, v2);
  const FT length_v1 = length_3(traits, v1);
  const FT length_v2 = length_3(traits, v2);

  const FT lb = -FT(999) / FT(1000), ub = FT(999) / FT(1000);
  FT cosine = dot / length_v1 / length_v2;
  cosine = (cosine < lb) ? lb : cosine;
  cosine = (cosine > ub) ? ub : cosine;
  const FT sine = sqrt(FT(1) - cosine * cosine);

  CGAL_assertion(sine != FT(0));
  if (sine != FT(0)) {
    return cosine / sine;
  }
  return FT(0); // undefined
}

template<typename FT>
class Tangent_weight_wrapper {
  FT m_d_r, m_d_p, m_w_base;

public:

  template<typename Point>
  Tangent_weight_wrapper(
    const Point& p, const Point& q, const Point& r) {

    m_d_r = CGAL::Weights::distance(q, r);
    CGAL_assertion(m_d_r != FT(0)); // two points are identical!
    m_d_p = CGAL::Weights::distance(q, p);
    CGAL_assertion(m_d_p != FT(0)); // two points are identical!

    const auto A = CGAL::Weights::area(p, q, r);
    CGAL_assertion(A != FT(0)); // three points are identical!
    const auto S = CGAL::Weights::scalar_product(p, q, r);
    m_w_base = -CGAL::Weights::tangent_half_angle(m_d_r, m_d_p, A, S);
  }

  FT get_w_r() const {
    return CGAL::Weights::half_tangent_weight(m_w_base, m_d_r) / FT(2);
  }

  FT get_w_p() const {
    return CGAL::Weights::half_tangent_weight(m_w_base, m_d_p) / FT(2);
  }
};

template<typename PolygonMesh>
class Cotangent_weight_wrapper {
  bool m_use_secure_version;

public:
  using halfedge_descriptor = typename boost::graph_traits<PolygonMesh>::halfedge_descriptor;
  using vertex_descriptor   = typename boost::graph_traits<PolygonMesh>::vertex_descriptor;

  Cotangent_weight_wrapper(
    const bool use_secure_version = false) :
  m_use_secure_version(use_secure_version)
  { }

  template<class VertexPointMap>
  decltype(auto) operator()(
    const halfedge_descriptor he,
    const PolygonMesh& pmesh,
    const VertexPointMap& pmap) const {

    const auto v0 = target(he, pmesh);
    const auto v1 = source(he, pmesh);

    const auto& p0 = get(pmap, v0);
    const auto& p1 = get(pmap, v1);

    using Kernel = typename CGAL::Kernel_traits<
      typename boost::property_traits<VertexPointMap>::value_type>::type;
    using FT = typename Kernel::FT;

    FT weight = FT(0);
    if (is_border_edge(he, pmesh)) {
      const auto he_cw = opposite(next(he, pmesh), pmesh);
      auto v2 = source(he_cw, pmesh);

      if (is_border_edge(he_cw, pmesh)) {
        const auto he_ccw = prev(opposite(he, pmesh), pmesh);
        v2 = source(he_ccw, pmesh);

        const auto& p2 = get(pmap, v2);
        if (m_use_secure_version) {
          weight = cotangent_3_secure(p1, p2, p0);
        } else {
          weight = CGAL::Weights::cotangent(p1, p2, p0);
        }
        weight = (CGAL::max)(FT(0), weight);
        weight /= FT(2);
      } else {
        const auto& p2 = get(pmap, v2);
        if (m_use_secure_version) {
          weight = cotangent_3_secure(p0, p2, p1);
        } else {
          weight = CGAL::Weights::cotangent(p0, p2, p1);
        }
        weight = (CGAL::max)(FT(0), weight);
        weight /= FT(2);
      }

    } else {
      const auto he_cw = opposite(next(he, pmesh), pmesh);
      const auto v2 = source(he_cw, pmesh);
      const auto he_ccw = prev(opposite(he, pmesh), pmesh);
      const auto v3 = source(he_ccw, pmesh);

      const auto& p2 = get(pmap, v2);
      const auto& p3 = get(pmap, v3);
      FT cot_beta = FT(0), cot_gamma = FT(0);

      if (m_use_secure_version) {
        cot_beta = cotangent_3_secure(p0, p2, p1);
      } else {
        cot_beta = CGAL::Weights::cotangent(p0, p2, p1);
      }

      if (m_use_secure_version) {
        cot_gamma = cotangent_3_secure(p1, p3, p0);
      } else {
        cot_gamma = CGAL::Weights::cotangent(p1, p3, p0);
      }

      cot_beta  = (CGAL::max)(FT(0), cot_beta);  cot_beta  /= FT(2);
      cot_gamma = (CGAL::max)(FT(0), cot_gamma); cot_gamma /= FT(2);
      weight = cot_beta + cot_gamma;
    }
    return weight;
  }
};

template<
typename PolygonMesh,
typename VertexPointMap = typename boost::property_map<PolygonMesh, vertex_point_t>::type>
class Cotangent_weight_secure_with_voronoi_wrapper {

  using GeomTraits = typename CGAL::Kernel_traits<
      typename boost::property_traits<VertexPointMap>::value_type>::type;
  using FT = typename GeomTraits::FT;

  const PolygonMesh& m_pmesh;
  const VertexPointMap m_pmap;

public:
  using halfedge_descriptor = typename boost::graph_traits<PolygonMesh>::halfedge_descriptor;
  using vertex_descriptor   = typename boost::graph_traits<PolygonMesh>::vertex_descriptor;

  Cotangent_weight_secure_with_voronoi_wrapper(
    const PolygonMesh& pmesh, const VertexPointMap pmap) :
  m_pmesh(pmesh), m_pmap(pmap) { }

  FT w_i(const vertex_descriptor v_i) const {
    return FT(1) / (FT(2) * voronoi(v_i));
  }

  FT w_ij(const halfedge_descriptor he) const {
    return cotangent_secure(he);
  }

private:
  FT cotangent_secure(const halfedge_descriptor he) const {

    const auto v0 = target(he, m_pmesh);
    const auto v1 = source(he, m_pmesh);

    const auto& p0 = get(m_pmap, v0);
    const auto& p1 = get(m_pmap, v1);

    FT weight = FT(0);
    if (is_border_edge(he, m_pmesh)) {
      const auto he_cw = opposite(next(he, m_pmesh), m_pmesh);
      auto v2 = source(he_cw, m_pmesh);

      if (is_border_edge(he_cw, m_pmesh)) {
        const auto he_ccw = prev(opposite(he, m_pmesh), m_pmesh);
        v2 = source(he_ccw, m_pmesh);

        const auto& p2 = get(m_pmap, v2);
        weight = cotangent_3_secure(p1, p2, p0);
      } else {
        const auto& p2 = get(m_pmap, v2);
        weight = cotangent_3_secure(p0, p2, p1);
      }

    } else {
      const auto he_cw = opposite(next(he, m_pmesh), m_pmesh);
      const auto v2 = source(he_cw, m_pmesh);
      const auto he_ccw = prev(opposite(he, m_pmesh), m_pmesh);
      const auto v3 = source(he_ccw, m_pmesh);

      const auto& p2 = get(m_pmap, v2);
      const auto& p3 = get(m_pmap, v3);

      const FT cot_beta  = cotangent_3_secure(p0, p2, p1);
      const FT cot_gamma = cotangent_3_secure(p1, p3, p0);
      weight = cot_beta + cot_gamma;
    }
    return weight;
  }

  FT voronoi(const vertex_descriptor v0) const {

    FT voronoi_area = FT(0);
    CGAL_assertion(CGAL::is_triangle_mesh(m_pmesh));
    for (const auto& he : halfedges_around_target(halfedge(v0, m_pmesh), m_pmesh)) {
      CGAL_assertion(v0 == target(he, m_pmesh));
      if (is_border(he, m_pmesh)) {
        continue;
      }

      const auto v1 = source(he, m_pmesh);
      const auto v2 = target(next(he, m_pmesh), m_pmesh);

      const auto& p0 = get(m_pmap, v0);
      const auto& p1 = get(m_pmap, v1);
      const auto& p2 = get(m_pmap, v2);

      const auto angle0 = CGAL::angle(p1, p0, p2);
      const auto angle1 = CGAL::angle(p2, p1, p0);
      const auto angle2 = CGAL::angle(p0, p2, p1);

      const bool obtuse =
        (angle0 == CGAL::OBTUSE) ||
        (angle1 == CGAL::OBTUSE) ||
        (angle2 == CGAL::OBTUSE);

      if (!obtuse) {
        const FT cot_p1 = CGAL::Weights::cotangent(p2, p1, p0);
        const FT cot_p2 = CGAL::Weights::cotangent(p0, p2, p1);

        const FT t1 = cot_p1 * (p2 - p0).squared_length();
        const FT t2 = cot_p2 * (p1 - p0).squared_length();
        voronoi_area += (t1 + t2) / FT(8);

      } else {

        const FT A = static_cast<FT>(CGAL::sqrt(CGAL::to_double(
          CGAL::squared_area(p0, p1, p2))));
        if (angle0 == CGAL::OBTUSE) {
          voronoi_area += A / FT(2);
        } else {
          voronoi_area += A / FT(4);
        }
      }
    }
    CGAL_assertion(voronoi_area != FT(0));
    return voronoi_area;
  }
};

template<
typename PolygonMesh,
typename VertexPointMap = typename boost::property_map<PolygonMesh, vertex_point_t>::type>
class Edge_cotangent_weight_wrapper {

  using GeomTraits = typename CGAL::Kernel_traits<
      typename boost::property_traits<VertexPointMap>::value_type>::type;
  using FT = typename GeomTraits::FT;

  const PolygonMesh& m_pmesh;
  const VertexPointMap m_pmap;

public:
  using halfedge_descriptor = typename boost::graph_traits<PolygonMesh>::halfedge_descriptor;
  using vertex_descriptor   = typename boost::graph_traits<PolygonMesh>::vertex_descriptor;

  Edge_cotangent_weight_wrapper(const PolygonMesh& pmesh, const VertexPointMap pmap) :
  m_pmesh(pmesh), m_pmap(pmap) { }

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

      weight = CGAL::Weights::cotangent(p0, p2, p1);

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

      weight = CGAL::Weights::cotangent_weight(p2, p1, p3, p0) / FT(2);
    }
    return weight;
  }
};

template<
typename PolygonMesh,
typename VertexPointMap = typename boost::property_map<PolygonMesh, vertex_point_t>::type>
class Mean_value_weight_wrapper {

  using GeomTraits = typename CGAL::Kernel_traits<
      typename boost::property_traits<VertexPointMap>::value_type>::type;
  using FT = typename GeomTraits::FT;
  using Vector = typename GeomTraits::Vector_3;

  const PolygonMesh& m_pmesh;
  const VertexPointMap m_pmap;

public:
  using halfedge_descriptor = typename boost::graph_traits<PolygonMesh>::halfedge_descriptor;
  using vertex_descriptor   = typename boost::graph_traits<PolygonMesh>::vertex_descriptor;

  Mean_value_weight_wrapper(const PolygonMesh& pmesh, const VertexPointMap pmap) :
  m_pmesh(pmesh), m_pmap(pmap) { }

  // Returns the mean-value coordinate of the specified halfedge_descriptor.
  // Returns different values for different edge orientations (which is normal
  // behaviour according to the formula).
  FT operator()(const halfedge_descriptor he) const {

    const vertex_descriptor v0 = target(he, m_pmesh);
    const vertex_descriptor v1 = source(he, m_pmesh);
    const Vector vec = get(m_pmap, v0) - get(m_pmap, v1);
    const FT norm = static_cast<FT>(
      CGAL::sqrt(CGAL::to_double(vec.squared_length())));

    // Only one triangle for border edges.
    if (is_border_edge(he, m_pmesh)) {
      const halfedge_descriptor he_cw = opposite(next(he, m_pmesh), m_pmesh);
      vertex_descriptor v2 = source(he_cw, m_pmesh);
      if (is_border_edge(he_cw, m_pmesh)) {
        const halfedge_descriptor he_ccw = prev(opposite(he, m_pmesh), m_pmesh);
        v2 = source(he_ccw, m_pmesh);
      }
      return half_tan_value_2(v1, v0, v2) / norm;
    } else {
      const halfedge_descriptor he_cw = opposite(next(he, m_pmesh), m_pmesh);
      const vertex_descriptor v2 = source(he_cw, m_pmesh);
      const halfedge_descriptor he_ccw = prev(opposite(he, m_pmesh), m_pmesh);
      const vertex_descriptor v3 = source(he_ccw, m_pmesh);
      return (
        half_tan_value_2(v1, v0, v2) / norm +
        half_tan_value_2(v1, v0, v3) / norm );
    }
  }

private:
  // The authors deviation built on Meyer_02.
  // See Iterative_authalic_parameterizer_3.h.
  FT half_tan_value_2(
    const vertex_descriptor v0,
    const vertex_descriptor v1,
    const vertex_descriptor v2) const {

    const Vector a = get(m_pmap, v0) - get(m_pmap, v1);
    const Vector b = get(m_pmap, v2) - get(m_pmap, v1);
    const FT dot_ab = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    const FT dot_aa = a.squared_length();
    const FT dot_bb = b.squared_length();
    const FT dot_aa_bb = dot_aa * dot_bb;

    const FT cos_rep = dot_ab;
    const FT sin_rep = static_cast<FT>(
      CGAL::sqrt(CGAL::to_double(dot_aa_bb  - dot_ab * dot_ab)));
    const FT normalizer = static_cast<FT>(
      CGAL::sqrt(CGAL::to_double(dot_aa_bb))); // |a| * |b|

    // The formula from [Floater04] page 4:
    // tan(Q / 2) = (1 - cos(Q)) / sin(Q).
    return (normalizer - cos_rep) / sin_rep;
  }
};

template<typename PolygonMesh>
class Single_cotangent_weight_wrapper {

public:
  using halfedge_descriptor = typename boost::graph_traits<PolygonMesh>::halfedge_descriptor;
  using vertex_descriptor   = typename boost::graph_traits<PolygonMesh>::vertex_descriptor;

  template<class VertexPointMap>
  decltype(auto) operator()(
    const halfedge_descriptor he,
    const PolygonMesh& pmesh,
    const VertexPointMap& pmap) const {

    using Kernel = typename CGAL::Kernel_traits<
      typename boost::property_traits<VertexPointMap>::value_type>::type;
    using FT = typename Kernel::FT;

    if (is_border(he, pmesh)) {
      return FT(0);
    }

    const vertex_descriptor v0 = target(he, pmesh);
    const vertex_descriptor v1 = source(he, pmesh);
    const vertex_descriptor v2 = target(next(he, pmesh), pmesh);

    const auto& p0 = get(pmap, v0);
    const auto& p1 = get(pmap, v1);
    const auto& p2 = get(pmap, v2);

    return CGAL::Weights::cotangent(p0, p2, p1);
  }
};

template<class PolygonMesh>
class Uniform_weight_wrapper {

public:
  using halfedge_descriptor = typename boost::graph_traits<PolygonMesh>::halfedge_descriptor;
  using vertex_descriptor = typename boost::graph_traits<PolygonMesh>::vertex_descriptor;
  double w_ij(halfedge_descriptor) { return 1.0; }
  double w_i(vertex_descriptor) { return 1.0; }
};

} // namespace internal
} // namespace Weights
} // namespace CGAL

#endif // CGAL_WEIGHTS_INTERNAL_TOOLS_H
