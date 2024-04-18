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

#include <CGAL/Weights/utils.h>

#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Point_2.h>
#include <CGAL/Point_3.h>
#include <CGAL/utility.h>

namespace CGAL {
namespace Weights {

/// \cond SKIP_IN_MANUAL

namespace cotangent_ns {

template<typename FT>
FT half_weight(const FT cot)
{
  return FT(2) * cot;
}

template<typename FT>
FT weight(const FT cot_beta, const FT cot_gamma)
{
  return FT(2) * (cot_beta + cot_gamma);
}

} // namespace cotangent_ns

/// \endcond

/*!
  \ingroup PkgWeightsRefCotangentWeights

  \brief computes the half value of the cotangent weight.

  This function constructs the half of the cotangent weight using the precomputed
  cotangent value. The returned value is \f$2\textbf{cot}\f$.

  \tparam FT a model of `FieldNumberType`

  \param cot the cotangent value

  \sa `cotangent_weight()`
*/
template<typename FT>
FT half_cotangent_weight(const FT cot)
{
  return cotangent_ns::half_weight(cot);
}

/*!
  \ingroup PkgWeightsRefCotangentWeights
  \brief computes the cotangent weight in 2D at `q` using the points `p0`, `p1`, and `p2`.
  \tparam GeomTraits a model of `AnalyticWeightTraits_2`
*/
template<typename GeomTraits>
typename GeomTraits::FT cotangent_weight(const typename GeomTraits::Point_2& p0,
                                         const typename GeomTraits::Point_2& p1,
                                         const typename GeomTraits::Point_2& p2,
                                         const typename GeomTraits::Point_2& q,
                                         const GeomTraits& traits)
{
  using FT = typename GeomTraits::FT;

  const FT cot_beta = cotangent_2(q, p0, p1, traits);
  const FT cot_gamma = cotangent_2(p1, p2, q, traits);

  return cotangent_ns::weight(cot_beta, cot_gamma);
}

/*!
  \ingroup PkgWeightsRefCotangentWeights
  \brief computes the cotangent weight in 2D at `q` using the points `p0`, `p1`, and `p2`.
  \tparam Kernel a model of `Kernel`
*/
template<typename Kernel>
typename Kernel::FT cotangent_weight(const CGAL::Point_2<Kernel>& p0,
                                     const CGAL::Point_2<Kernel>& p1,
                                     const CGAL::Point_2<Kernel>& p2,
                                     const CGAL::Point_2<Kernel>& q)
{
  Kernel traits;
  return cotangent_weight(p0, p1, p2, q, traits);
}

// 3D ==============================================================================================

/*!
  \ingroup PkgWeightsRefCotangentWeights
  \brief computes the cotangent weight in 3D at `q` using the points `p0`, `p1`, and `p2`.
  \tparam GeomTraits a model of `AnalyticWeightTraits_3`
*/
template<typename GeomTraits>
typename GeomTraits::FT cotangent_weight(const typename GeomTraits::Point_3& p0,
                                         const typename GeomTraits::Point_3& p1,
                                         const typename GeomTraits::Point_3& p2,
                                         const typename GeomTraits::Point_3& q,
                                         const GeomTraits& traits)
{
  using FT = typename GeomTraits::FT;

  const FT cot_beta = cotangent_3(q, p0, p1, traits);
  const FT cot_gamma = cotangent_3(p1, p2, q, traits);

  return cotangent_ns::weight(cot_beta, cot_gamma);
}

/*!
  \ingroup PkgWeightsRefCotangentWeights
  \brief computes the cotangent weight in 3D at `q` using the points `p0`, `p1`, and `p2`.
  \tparam Kernel a model of `Kernel`
*/
template<typename Kernel>
typename Kernel::FT cotangent_weight(const CGAL::Point_3<Kernel>& p0,
                                     const CGAL::Point_3<Kernel>& p1,
                                     const CGAL::Point_3<Kernel>& p2,
                                     const CGAL::Point_3<Kernel>& q)
{
  Kernel traits;
  return cotangent_weight(p0, p1, p2, q, traits);
}

/// \cond SKIP_IN_MANUAL

// Undocumented cotangent weight class.
// Returns: cot(beta)
//
// Returns a single cotangent weight, its operator() is defined based on the
// halfedge_descriptor, polygon mesh, and vertex to point map.
// For border edges it returns zero.
// This version is currently used in:
// Surface_mesh_deformation -> Surface_mesh_deformation.h
template<typename PolygonMesh>
class Single_cotangent_weight
{
  using vertex_descriptor = typename boost::graph_traits<PolygonMesh>::vertex_descriptor;
  using halfedge_descriptor = typename boost::graph_traits<PolygonMesh>::halfedge_descriptor;

public:
  // Returns the cotangent of the opposite angle of the edge
  // 0 for border edges (which does not have an opposite angle).
  template <typename VPM>
  auto operator()(halfedge_descriptor he,
                  PolygonMesh& pmesh,
                  VPM vpm)
  {
    using Point = typename boost::property_traits<VPM>::value_type;
    using Point_ref = typename boost::property_traits<VPM>::reference;

    using GeomTraits = typename Kernel_traits<Point>::type;
    using FT = typename GeomTraits::FT;

    if(is_border(he, pmesh))
      return FT{0};

    const vertex_descriptor v0 = target(he, pmesh);
    const vertex_descriptor v1 = source(he, pmesh);
    const vertex_descriptor v2 = target(next(he, pmesh), pmesh);

    const Point_ref p0 = get(vpm, v0);
    const Point_ref p1 = get(vpm, v1);
    const Point_ref p2 = get(vpm, v2);

    return cotangent(p0, p2, p1);
  }
};

// Undocumented cotangent weight class.
// Returns: 0.5 * (cot(beta) + cot(gamma))
//
// Its constructor takes a boolean flag to choose between default and clamped
// versions of the cotangent weights and its operator() is defined based on the
// halfedge_descriptor, polygon mesh, and vertex to point map.
// This version is currently used in:
// Polygon_mesh_processing -> curvature_flow_impl.h (no clamping, no bounding)
// Surface_mesh_deformation -> Surface_mesh_deformation.h (default version)
// Surface_mesh_parameterizer -> Orbifold_Tutte_parameterizer_3.h (default version)
// Surface_mesh_skeletonization -> Mean_curvature_flow_skeletonization.h (clamped version)
//
// The API is a bit awkward: the template parameters VertexPointMap and GeomTraits
// are only meaningful in the API that calls the operator with a single parameter.
template<typename PolygonMesh,
         typename VertexPointMap = typename GetVertexPointMap<PolygonMesh>::type,
         typename GeomTraits = typename Kernel_traits<
                                 typename boost::property_traits<VertexPointMap>::value_type>::type>
class Cotangent_weight
{
  using vertex_descriptor = typename boost::graph_traits<PolygonMesh>::vertex_descriptor;
  using halfedge_descriptor = typename boost::graph_traits<PolygonMesh>::halfedge_descriptor;

private:
  // These class members are used only when the constructor initializing them
  // is used, but Surface_mesh_deformation has its own weight API locked
  // by the concept SurfaceMeshDeformationWeights.
  // A bit awkward, but better than duplicating code...
  PolygonMesh const * const m_pmesh_ptr;
  const VertexPointMap m_vpm;
  const GeomTraits m_traits;

  bool m_use_clamped_version;
  bool m_bound_from_below;

public:
  Cotangent_weight()
    : m_pmesh_ptr(nullptr), m_vpm(), m_traits(), m_use_clamped_version(false), m_bound_from_below(true)
  { }

  // Common API whether mesh/vpm/traits are initialized in the constructor,
  // or passed in the operator()
  template <typename VPM, typename GT>
  typename GT::FT
  operator()(const halfedge_descriptor he,
             const PolygonMesh& pmesh,
             const VPM vpm,
             const GT& traits) const
  {
    using Point_ref = typename boost::property_traits<VPM>::reference;
    using FT = typename GT::FT;

    if(is_border(he, pmesh))
      return FT{0};

    auto half_weight = [&] (const halfedge_descriptor he) -> FT
    {
      if(is_border(he, pmesh))
        return FT{0};

      const vertex_descriptor v0 = target(he, pmesh);
      const vertex_descriptor v1 = source(he, pmesh);
      const vertex_descriptor v2 = target(next(he, pmesh), pmesh);

      const Point_ref p0 = get(vpm, v0);
      const Point_ref p1 = get(vpm, v1);
      const Point_ref p2 = get(vpm, v2);

      FT weight = 0;
      if (m_use_clamped_version)
        weight = cotangent_3_clamped(p1, p2, p0, traits);
      else
        weight = cotangent_3(p1, p2, p0, traits);

      if(m_bound_from_below)
        weight = (CGAL::max)(FT(0), weight);

      return weight / FT(2);
    };

    FT weight = half_weight(he) + half_weight(opposite(he, pmesh));
    return weight;
  }

  // That is the API called by Surface_mesh_deformation
  template <typename VPM>
  auto // kernel_traits<VPM::value_type>::type::FT
  operator()(const halfedge_descriptor he,
             const PolygonMesh& pmesh,
             const VPM vpm) const
  {
    using Point = typename boost::property_traits<VPM>::value_type;
    using GT = typename Kernel_traits<Point>::type;
    return this->operator()(he, pmesh, vpm, GT());
  }

public:
  // This is the "normal" API: give all info to the constructor, and operator()(halfedge)
  Cotangent_weight(const PolygonMesh& pmesh,
                   const VertexPointMap vpm,
                   const GeomTraits& traits = GeomTraits(),
                   const bool use_clamped_version = false,
                   const bool bound_from_below = true)
    : m_pmesh_ptr(&pmesh), m_vpm(vpm), m_traits(traits),
      m_use_clamped_version(use_clamped_version),
      m_bound_from_below(bound_from_below)
  { }

  typename GeomTraits::FT operator()(const halfedge_descriptor he) const
  {
    CGAL_precondition(m_pmesh_ptr != nullptr);
    return this->operator()(he, *m_pmesh_ptr, m_vpm, m_traits);
  }
};

// Undocumented cotangent weight class.
//
// Its constructor takes a polygon mesh and a vertex to point map
// and its operator() is defined based on the halfedge_descriptor only.
// This class is using a special clamped version of the cotangent weights.
// This version is currently used in:
// Polygon_mesh_processing -> fair.h
// CGAL Lab -> Hole_filling_plugin.cpp
template<typename PolygonMesh,
         typename VertexPointMap,
         typename GeomTraits>
class Secure_cotangent_weight_with_voronoi_area
{
  using vertex_descriptor = typename boost::graph_traits<PolygonMesh>::vertex_descriptor;
  using halfedge_descriptor = typename boost::graph_traits<PolygonMesh>::halfedge_descriptor;

  using Point_ref = typename boost::property_traits<VertexPointMap>::reference;
  using FT = typename GeomTraits::FT;
  using Vector_3 = typename GeomTraits::Vector_3;

private:
  const PolygonMesh& m_pmesh;
  const VertexPointMap m_vpm;
  GeomTraits m_traits;

  Cotangent_weight<PolygonMesh, VertexPointMap, GeomTraits> cotangent_weight_calculator;

public:
  Secure_cotangent_weight_with_voronoi_area(const PolygonMesh& pmesh,
                                            const VertexPointMap vpm,
                                            const GeomTraits& traits = GeomTraits())
    : m_pmesh(pmesh), m_vpm(vpm), m_traits(traits),
      cotangent_weight_calculator(m_pmesh, m_vpm, m_traits,
                                  true /*clamp*/, true /*bound from below*/)
  { }

  FT w_i(const vertex_descriptor v_i) const
  {
    return FT(1) / (FT(2) * voronoi(v_i));
  }

  FT w_ij(const halfedge_descriptor he) const
  {
    return cotangent_weight_calculator(he);
  }

private:
  FT voronoi(const vertex_descriptor v0) const
  {
    auto squared_length_3 = m_traits.compute_squared_length_3_object();
    auto vector_3 = m_traits.construct_vector_3_object();

    FT voronoi_area = FT(0);
    for (const halfedge_descriptor he : halfedges_around_target(halfedge(v0, m_pmesh), m_pmesh))
    {
      CGAL_assertion(v0 == target(he, m_pmesh));
      CGAL_assertion(CGAL::is_triangle(he, m_pmesh));

      if (is_border(he, m_pmesh))
        continue;

      const vertex_descriptor v1 = source(he, m_pmesh);
      const vertex_descriptor v2 = target(next(he, m_pmesh), m_pmesh);

      const Point_ref p0 = get(m_vpm, v0);
      const Point_ref p1 = get(m_vpm, v1);
      const Point_ref p2 = get(m_vpm, v2);

      const CGAL::Angle angle0 = CGAL::angle(p1, p0, p2);
      if((angle0 == CGAL::OBTUSE) ||
         (CGAL::angle(p2, p1, p0) == CGAL::OBTUSE) ||
         (CGAL::angle(p0, p2, p1) == CGAL::OBTUSE))
      {
        const FT A = internal::positive_area_3(p0, p1, p2, m_traits);
        if (angle0 == CGAL::OBTUSE)
          voronoi_area += A / FT(2);
         else
          voronoi_area += A / FT(4);
      }
      else
      {
        const FT cot_p1 = cotangent_3_clamped(p2, p1, p0, m_traits);
        const FT cot_p2 = cotangent_3_clamped(p0, p2, p1, m_traits);

        const Vector_3 v1 = vector_3(p0, p1);
        const Vector_3 v2 = vector_3(p0, p2);

        const FT t1 = cot_p1 * squared_length_3(v2);
        const FT t2 = cot_p2 * squared_length_3(v1);
        voronoi_area += (t1 + t2) / FT(8);
      }
    }

    CGAL_assertion(!is_zero(voronoi_area));
    return voronoi_area;
  }
};

/// \endcond

} // namespace Weights
} // namespace CGAL

#endif // CGAL_COTANGENT_WEIGHTS_H
