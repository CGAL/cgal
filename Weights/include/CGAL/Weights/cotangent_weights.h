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

#ifndef CGAL_COTANGENT_WEIGHTS_H
#define CGAL_COTANGENT_WEIGHTS_H

// Internal includes.
#include <CGAL/Weights/internal/utils.h>

namespace CGAL {
namespace Weights {

  /// \cond SKIP_IN_MANUAL
  namespace cotangent_ns {

    template<typename FT>
    FT half_weight(const FT cot) {
      return FT(2) * cot;
    }

    template<typename FT>
    FT weight(const FT cot_beta, const FT cot_gamma) {
      return FT(2) * (cot_beta + cot_gamma);
    }
  }
  /// \endcond

  /*!
    \ingroup PkgWeightsRefCotangentWeights

    \brief computes the half value of the cotangent weight.

    This function constructs the half of the cotangent weight using the precomputed
    cotangent value. The returned value is
    \f$2\textbf{cot}\f$.

    \tparam FT
    a model of `FieldNumberType`

    \param cot
    the cotangent value

    \sa `cotangent_weight()`
  */
  template<typename FT>
  FT half_cotangent_weight(const FT cot) {
    return cotangent_ns::half_weight(cot);
  }

  #if defined(DOXYGEN_RUNNING)

  /*!
    \ingroup PkgWeightsRefCotangentWeights

    \brief computes the cotangent weight in 2D at `q` using the points `p0`, `p1`,
    and `p2`, given a traits class `traits` with geometric objects, predicates, and constructions.
  */
  template<typename GeomTraits>
  typename GeomTraits::FT cotangent_weight(
    const typename GeomTraits::Point_2& p0,
    const typename GeomTraits::Point_2& p1,
    const typename GeomTraits::Point_2& p2,
    const typename GeomTraits::Point_2& q,
    const GeomTraits& traits) { }

  /*!
    \ingroup PkgWeightsRefCotangentWeights

    \brief computes the cotangent weight in 3D at `q` using the points `p0`, `p1`,
    and `p2`, given a traits class `traits` with geometric objects, predicates, and constructions.
  */
  template<typename GeomTraits>
  typename GeomTraits::FT cotangent_weight(
    const typename GeomTraits::Point_3& p0,
    const typename GeomTraits::Point_3& p1,
    const typename GeomTraits::Point_3& p2,
    const typename GeomTraits::Point_3& q,
    const GeomTraits& traits) { }

  /*!
    \ingroup PkgWeightsRefCotangentWeights

    \brief computes the cotangent weight in 2D at `q` using the points `p0`, `p1`,
    and `p2` which are parameterized by a `Kernel` K.
  */
  template<typename K>
  typename K::FT cotangent_weight(
    const CGAL::Point_2<K>& p0,
    const CGAL::Point_2<K>& p1,
    const CGAL::Point_2<K>& p2,
    const CGAL::Point_2<K>& q) { }

  /*!
    \ingroup PkgWeightsRefCotangentWeights

    \brief computes the cotangent weight in 3D at `q` using the points `p0`, `p1`,
    and `p2` which are parameterized by a `Kernel` K.
  */
  template<typename K>
  typename K::FT cotangent_weight(
    const CGAL::Point_3<K>& p0,
    const CGAL::Point_3<K>& p1,
    const CGAL::Point_3<K>& p2,
    const CGAL::Point_3<K>& q) { }

  #endif // DOXYGEN_RUNNING

  /// \cond SKIP_IN_MANUAL
  template<typename GeomTraits>
  typename GeomTraits::FT cotangent_weight(
    const typename GeomTraits::Point_2& t,
    const typename GeomTraits::Point_2& r,
    const typename GeomTraits::Point_2& p,
    const typename GeomTraits::Point_2& q,
    const GeomTraits& traits) {

    using FT = typename GeomTraits::FT;
    const FT cot_beta  = internal::cotangent_2(traits, q, t, r);
    const FT cot_gamma = internal::cotangent_2(traits, r, p, q);
    return cotangent_ns::weight(cot_beta, cot_gamma);
  }

  template<typename GeomTraits>
  typename GeomTraits::FT cotangent_weight(
    const CGAL::Point_2<GeomTraits>& t,
    const CGAL::Point_2<GeomTraits>& r,
    const CGAL::Point_2<GeomTraits>& p,
    const CGAL::Point_2<GeomTraits>& q) {

    const GeomTraits traits;
    return cotangent_weight(t, r, p, q, traits);
  }

  template<typename GeomTraits>
  typename GeomTraits::FT cotangent_weight(
    const typename GeomTraits::Point_3& t,
    const typename GeomTraits::Point_3& r,
    const typename GeomTraits::Point_3& p,
    const typename GeomTraits::Point_3& q,
    const GeomTraits& traits) {

    using FT = typename GeomTraits::FT;
    const FT cot_beta  = internal::cotangent_3(traits, q, t, r);
    const FT cot_gamma = internal::cotangent_3(traits, r, p, q);
    return cotangent_ns::weight(cot_beta, cot_gamma);
  }

  template<typename GeomTraits>
  typename GeomTraits::FT cotangent_weight(
    const CGAL::Point_3<GeomTraits>& t,
    const CGAL::Point_3<GeomTraits>& r,
    const CGAL::Point_3<GeomTraits>& p,
    const CGAL::Point_3<GeomTraits>& q) {

    const GeomTraits traits;
    return cotangent_weight(t, r, p, q, traits);
  }

  // Undocumented cotangent weight class.
  // Its constructor takes a polygon mesh and a vertex to point map
  // and its operator() is defined based on the halfedge_descriptor only.
  // This version is currently used in:
  // Polygon_mesh_processing -> curvature_flow_impl.h
  template<
  typename PolygonMesh,
  typename VertexPointMap = typename boost::property_map<PolygonMesh, CGAL::vertex_point_t>::type>
  class Edge_cotangent_weight {

    using GeomTraits = typename CGAL::Kernel_traits<
        typename boost::property_traits<VertexPointMap>::value_type>::type;
    using FT = typename GeomTraits::FT;

    const PolygonMesh& m_pmesh;
    const VertexPointMap m_pmap;
    const GeomTraits m_traits;

  public:
    using vertex_descriptor = typename boost::graph_traits<PolygonMesh>::vertex_descriptor;
    using halfedge_descriptor = typename boost::graph_traits<PolygonMesh>::halfedge_descriptor;

    Edge_cotangent_weight(const PolygonMesh& pmesh, const VertexPointMap pmap) :
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

        weight = internal::cotangent_3(m_traits, p0, p2, p1);

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

        weight = cotangent_weight(p2, p1, p3, p0) / FT(2);
      }
      return weight;
    }
  };

  // Undocumented cotangent weight class.
  // Returns a single cotangent weight, its operator() is defined based on the
  // halfedge_descriptor, polygon mesh, and vertex to point map.
  // For border edges it returns zero.
  // This version is currently used in:
  // Surface_mesh_deformation -> Surface_mesh_deformation.h
  template<typename PolygonMesh>
  class Single_cotangent_weight {

  public:
    using vertex_descriptor = typename boost::graph_traits<PolygonMesh>::vertex_descriptor;
    using halfedge_descriptor = typename boost::graph_traits<PolygonMesh>::halfedge_descriptor;

    template<class VertexPointMap>
    decltype(auto) operator()(const halfedge_descriptor he,
      const PolygonMesh& pmesh, const VertexPointMap pmap) const {

      using GeomTraits = typename CGAL::Kernel_traits<
        typename boost::property_traits<VertexPointMap>::value_type>::type;
      using FT = typename GeomTraits::FT;
      const GeomTraits traits;

      if (is_border(he, pmesh)) {
        return FT(0);
      }

      const vertex_descriptor v0 = target(he, pmesh);
      const vertex_descriptor v1 = source(he, pmesh);
      const vertex_descriptor v2 = target(next(he, pmesh), pmesh);

      const auto& p0 = get(pmap, v0);
      const auto& p1 = get(pmap, v1);
      const auto& p2 = get(pmap, v2);

      return internal::cotangent_3(traits, p0, p2, p1);
    }
  };

  // Undocumented cotangent weight class.
  // Its constructor takes a boolean flag to choose between default and clamped
  // versions of the cotangent weights and its operator() is defined based on the
  // halfedge_descriptor, polygon mesh, and vertex to point map.
  // This version is currently used in:
  // Surface_mesh_deformation -> Surface_mesh_deformation.h (default version)
  // Surface_mesh_parameterizer -> Orbifold_Tutte_parameterizer_3.h (default version)
  // Surface_mesh_skeletonization -> Mean_curvature_flow_skeletonization.h (clamped version)
  template<typename PolygonMesh>
  class Cotangent_weight {
    bool m_use_clamped_version;

  public:
    using vertex_descriptor = typename boost::graph_traits<PolygonMesh>::vertex_descriptor;
    using halfedge_descriptor = typename boost::graph_traits<PolygonMesh>::halfedge_descriptor;

    Cotangent_weight(const bool use_clamped_version = false) :
    m_use_clamped_version(use_clamped_version) { }

    template<class VertexPointMap>
    decltype(auto) operator()(const halfedge_descriptor he,
      const PolygonMesh& pmesh, const VertexPointMap pmap) const {

      using GeomTraits = typename CGAL::Kernel_traits<
        typename boost::property_traits<VertexPointMap>::value_type>::type;
      using FT = typename GeomTraits::FT;
      const GeomTraits traits;

      const auto v0 = target(he, pmesh);
      const auto v1 = source(he, pmesh);

      const auto& p0 = get(pmap, v0);
      const auto& p1 = get(pmap, v1);

      FT weight = FT(0);
      if (is_border_edge(he, pmesh)) {
        const auto he_cw = opposite(next(he, pmesh), pmesh);
        auto v2 = source(he_cw, pmesh);

        if (is_border_edge(he_cw, pmesh)) {
          const auto he_ccw = prev(opposite(he, pmesh), pmesh);
          v2 = source(he_ccw, pmesh);

          const auto& p2 = get(pmap, v2);
          if (m_use_clamped_version) {
            weight = internal::cotangent_3_clamped(traits, p1, p2, p0);
          } else {
            weight = internal::cotangent_3(traits, p1, p2, p0);
          }
          weight = (CGAL::max)(FT(0), weight);
          weight /= FT(2);
        } else {
          const auto& p2 = get(pmap, v2);
          if (m_use_clamped_version) {
            weight = internal::cotangent_3_clamped(traits, p0, p2, p1);
          } else {
            weight = internal::cotangent_3(traits, p0, p2, p1);
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

        if (m_use_clamped_version) {
          cot_beta = internal::cotangent_3_clamped(traits, p0, p2, p1);
        } else {
          cot_beta = internal::cotangent_3(traits, p0, p2, p1);
        }

        if (m_use_clamped_version) {
          cot_gamma = internal::cotangent_3_clamped(traits, p1, p3, p0);
        } else {
          cot_gamma = internal::cotangent_3(traits, p1, p3, p0);
        }

        cot_beta  = (CGAL::max)(FT(0), cot_beta);  cot_beta  /= FT(2);
        cot_gamma = (CGAL::max)(FT(0), cot_gamma); cot_gamma /= FT(2);
        weight = cot_beta + cot_gamma;
      }
      return weight;
    }
  };

  // Undocumented cotangent weight class.
  // Its constructor takes a polygon mesh and a vertex to point map
  // and its operator() is defined based on the halfedge_descriptor only.
  // This class is using a special clamped version of the cotangent weights.
  // This version is currently used in:
  // Polygon_mesh_processing -> fair.h
  // Polyhedron demo -> Hole_filling_plugin.cpp
  template<
  typename PolygonMesh,
  typename VertexPointMap = typename boost::property_map<PolygonMesh, CGAL::vertex_point_t>::type>
  class Secure_cotangent_weight_with_voronoi_area {

    using GeomTraits = typename CGAL::Kernel_traits<
        typename boost::property_traits<VertexPointMap>::value_type>::type;
    using FT = typename GeomTraits::FT;

    const PolygonMesh& m_pmesh;
    const VertexPointMap m_pmap;
    const GeomTraits m_traits;

  public:
    using vertex_descriptor = typename boost::graph_traits<PolygonMesh>::vertex_descriptor;
    using halfedge_descriptor = typename boost::graph_traits<PolygonMesh>::halfedge_descriptor;

    Secure_cotangent_weight_with_voronoi_area(const PolygonMesh& pmesh, const VertexPointMap pmap) :
    m_pmesh(pmesh), m_pmap(pmap), m_traits() { }

    FT w_i(const vertex_descriptor v_i) const {
      return FT(1) / (FT(2) * voronoi(v_i));
    }

    FT w_ij(const halfedge_descriptor he) const {
      return cotangent_clamped(he);
    }

  private:
    FT cotangent_clamped(const halfedge_descriptor he) const {

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
          weight = internal::cotangent_3_clamped(m_traits, p1, p2, p0);
        } else {
          const auto& p2 = get(m_pmap, v2);
          weight = internal::cotangent_3_clamped(m_traits, p0, p2, p1);
        }

      } else {
        const auto he_cw = opposite(next(he, m_pmesh), m_pmesh);
        const auto v2 = source(he_cw, m_pmesh);
        const auto he_ccw = prev(opposite(he, m_pmesh), m_pmesh);
        const auto v3 = source(he_ccw, m_pmesh);

        const auto& p2 = get(m_pmap, v2);
        const auto& p3 = get(m_pmap, v3);

        const FT cot_beta  = internal::cotangent_3_clamped(m_traits, p0, p2, p1);
        const FT cot_gamma = internal::cotangent_3_clamped(m_traits, p1, p3, p0);
        weight = cot_beta + cot_gamma;
      }
      return weight;
    }

    FT voronoi(const vertex_descriptor v0) const {

      const auto squared_length_3 =
        m_traits.compute_squared_length_3_object();
      const auto construct_vector_3 =
        m_traits.construct_vector_3_object();

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
          const FT cot_p1 = internal::cotangent_3(m_traits, p2, p1, p0);
          const FT cot_p2 = internal::cotangent_3(m_traits, p0, p2, p1);

          const auto v1 = construct_vector_3(p0, p1);
          const auto v2 = construct_vector_3(p0, p2);

          const FT t1 = cot_p1 * squared_length_3(v2);
          const FT t2 = cot_p2 * squared_length_3(v1);
          voronoi_area += (t1 + t2) / FT(8);

        } else {

          const FT A = internal::positive_area_3(m_traits, p0, p1, p2);
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

  /// \endcond

} // namespace Weights
} // namespace CGAL

#endif // CGAL_COTANGENT_WEIGHTS_H
