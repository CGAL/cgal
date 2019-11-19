// Copyright (c) 2017-2018 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Pierre Alliez and Lingjie Zhu


#ifndef CGAL_L2_METRIC_PLANE_PROXY_H
#define CGAL_L2_METRIC_PLANE_PROXY_H

#include <CGAL/license/Surface_mesh_approximation.h>

#include <CGAL/Kernel/global_functions.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Dynamic_property_map.h>

#include <boost/graph/graph_traits.hpp>
#include <boost/unordered_map.hpp>

#include <list>

namespace CGAL {
namespace Surface_mesh_approximation {

/// \ingroup PkgTSMARef
/// @brief Approximation L2 metric of plane proxy.
///
/// \cgalModels `ErrorMetricProxy`
///
/// @tparam TriangleMesh a triangle `FaceGraph`
/// @tparam VertexPointMap a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
///    as key type, GeomTraits::Point_3 as value type
/// @tparam GeomTraits a model of Kernel
template <typename TriangleMesh,
  typename VertexPointMap
    = typename boost::property_map<TriangleMesh, boost::vertex_point_t>::type,
  typename GeomTraits
    = typename Kernel_traits<typename boost::property_traits<VertexPointMap>::value_type>::Kernel>
class L2_metric_plane_proxy {
  typedef typename GeomTraits::FT FT;
  typedef typename GeomTraits::Point_3 Point_3;
  typedef typename GeomTraits::Triangle_3 Triangle_3;

  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;

  typedef CGAL::dynamic_face_property_t<FT> Face_area_tag;
  typedef typename boost::property_map<TriangleMesh, Face_area_tag>::type Face_area_map;

  typedef TriangleMesh Triangle_mesh;

public:
  /// \name Types
  /// @{

  /// Proxy type
  typedef typename GeomTraits::Plane_3 Proxy;
  /// @}

  /// \name Constructor
  /// @{
  /*!
   * @brief Constructor
   * @param tm triangle mesh
   * @param vpmap vertex point map
   */
  L2_metric_plane_proxy(const TriangleMesh &tm, const VertexPointMap &vpmap)
    : m_tm(&tm), m_vpmap(vpmap), m_famap( get(Face_area_tag(), const_cast<TriangleMesh &>(*m_tm)) )
  {
    for(face_descriptor f : faces(*m_tm)) {
      const halfedge_descriptor he = halfedge(f, *m_tm);
      const Point_3 &p0 = m_vpmap[source(he, *m_tm)];
      const Point_3 &p1 = m_vpmap[target(he, *m_tm)];
      const Point_3 &p2 = m_vpmap[target(next(he, *m_tm), *m_tm)];
      put(m_famap, f, CGAL::approximate_sqrt(CGAL::squared_area(p0, p1, p2)));
    }
  }
  /// @}

  /*!
   * @brief computes the L21 error from a face to a proxy,
   * using integral (closed-form) computation.
   * @param tm input triangle mesh
   * @param f face_descriptor of a face
   * @param px proxy
   * @return computed error
   */
  FT compute_error(const face_descriptor f, const TriangleMesh &tm, const Proxy &px) const {
    (void)(tm);
    halfedge_descriptor he = halfedge(f, *m_tm);
    const Point_3 &p0 = m_vpmap[source(he, *m_tm)];
    const Point_3 &p1 = m_vpmap[target(he, *m_tm)];
    const Point_3 &p2 = m_vpmap[target(next(he, *m_tm), *m_tm)];
    const FT sq_d0 = CGAL::squared_distance(p0, px);
    const FT sq_d1 = CGAL::squared_distance(p1, px);
    const FT sq_d2 = CGAL::squared_distance(p2, px);
    const FT d0 = CGAL::approximate_sqrt(sq_d0);
    const FT d1 = CGAL::approximate_sqrt(sq_d1);
    const FT d2 = CGAL::approximate_sqrt(sq_d2);

    return (sq_d0 + sq_d1 + sq_d2 +
            d0 * d1 + d1 * d2 + d2 * d0) * get(m_famap, f) / FT(6.0);
  }

  /*!
   * @brief fits a proxy from a range of faces, in the L2 sense, with an
   * integral (closed-form) formulation. The best-fit plane passes
   * through the center of mass and is defined by the two principal
   * components of the integral covariance matrix.
   * @tparam FaceRange range of face descriptors, model of Range.
   * @param faces the range of faces to be fitted
   * @param tm input triangle mesh
   * @return fitted proxy
   */
  template <typename FaceRange>
  Proxy fit_proxy(const FaceRange &faces, const TriangleMesh &tm) const {
    (void)(tm);
    CGAL_assertion(!faces.empty());

    std::list<Triangle_3> tris;
    for(const face_descriptor f : faces) {
      const halfedge_descriptor he = halfedge(f, *m_tm);
      const Point_3 &p0 = m_vpmap[source(he, *m_tm)];
      const Point_3 &p1 = m_vpmap[target(he, *m_tm)];
      const Point_3 &p2 = m_vpmap[target(next(he, *m_tm), *m_tm)];
      tris.push_back(Triangle_3(p0, p1, p2));
    }

    // construct and fit proxy plane
    Proxy proxy;
    CGAL::linear_least_squares_fitting_3(
      tris.begin(),
      tris.end(),
      proxy,
      CGAL::Dimension_tag<2>());

    // TODO: check and flip plane normal?

    return proxy;
  }

private:
  const TriangleMesh *m_tm;
  const VertexPointMap m_vpmap;
  Face_area_map m_famap;
};

} // namespace Surface_mesh_approximation
} // namespace CGAL

#endif // CGAL_L2_METRIC_PLANE_PROXY_H
