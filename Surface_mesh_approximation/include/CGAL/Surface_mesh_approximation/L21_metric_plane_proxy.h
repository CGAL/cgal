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


#ifndef CGAL_L21_METRIC_PLANE_PROXY_H
#define CGAL_L21_METRIC_PLANE_PROXY_H

#include <CGAL/license/Surface_mesh_approximation.h>

#include <CGAL/Kernel/global_functions.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/Dynamic_property_map.h>

#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/unordered_map.hpp>

namespace CGAL {
namespace Surface_mesh_approximation {

/// \ingroup PkgTSMARef
/// @brief Approximation L21 metric of vector proxy.
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
class L21_metric_plane_proxy {
  typedef typename GeomTraits::FT FT;
  typedef typename GeomTraits::Vector_3 Vector_3;
  typedef typename GeomTraits::Point_3 Point_3;
  typedef typename GeomTraits::Construct_scaled_vector_3 Construct_scaled_vector_3;
  typedef typename GeomTraits::Construct_sum_of_vectors_3 Construct_sum_of_vectors_3;
  typedef typename GeomTraits::Compute_scalar_product_3 Compute_scalar_product_3;
  typedef typename GeomTraits::Collinear_3 Collinear_3;

  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;

  typedef CGAL::dynamic_face_property_t<Vector_3> Face_normal_tag;
  typedef CGAL::dynamic_face_property_t<FT> Face_area_tag;
  typedef typename boost::property_map<TriangleMesh, Face_normal_tag>::type Face_normal_map;
  typedef typename boost::property_map<TriangleMesh, Face_area_tag>::type Face_area_map;

  typedef TriangleMesh Triangle_mesh;

public:
  /// \name Types
  /// @{

  /// Proxy type
  typedef typename GeomTraits::Vector_3 Proxy;
  /// @}

  /// \name Constructor
  /// @{
  /*!
   * @brief Constructor
   * @param tm triangle mesh
   * @param vpmap vertex point map
   */
  L21_metric_plane_proxy(const TriangleMesh &tm, const VertexPointMap &vpmap)
    : m_fnmap( get(Face_normal_tag(), const_cast<TriangleMesh &>(tm)) )
    , m_famap( get(Face_area_tag(), const_cast<TriangleMesh &>(tm)) )
  {
    GeomTraits traits;
    m_scalar_product_functor = traits.compute_scalar_product_3_object();
    m_sum_functor = traits.construct_sum_of_vectors_3_object();
    m_scale_functor = traits.construct_scaled_vector_3_object();
    m_collinear_functor = traits.collinear_3_object();

    // construct internal face normal & area map
    for(face_descriptor f : faces(tm)) {
      const halfedge_descriptor he = halfedge(f, tm);
      const Point_3 &p0 = vpmap[source(he, tm)];
      const Point_3 &p1 = vpmap[target(he, tm)];
      const Point_3 &p2 = vpmap[target(next(he, tm), tm)];
      if (CGAL::collinear(p0, p1, p2))
        put(m_fnmap, f, CGAL::NULL_VECTOR);
      else
        put(m_fnmap, f, CGAL::unit_normal(p0, p1, p2));
      put(m_famap, f, CGAL::approximate_sqrt(CGAL::squared_area(p0, p1, p2)));
    }
  }
  /// @}

  /*!
   * @brief computes the L2,1 error from a face to a proxy.
   * @param tm input triangle mesh
   * @param f face_descriptor of a face
   * @param px proxy
   * @return computed error
   */
  FT compute_error(const face_descriptor f, const TriangleMesh &tm, const Proxy &px) const {
    (void)(tm);
    Vector_3 v = m_sum_functor(get(m_fnmap, f), m_scale_functor(px, FT(-1.0)));
    return get(m_famap, f) * m_scalar_product_functor(v, v);
  }

  /*!
   * @brief fits a proxy to a range of faces.
   * @tparam FaceRange range of face descriptors, model of Range.
   * @param faces the range of faces to be fitted
   * @param tm input triangle mesh
   * @return fitted proxy
   */
  template <typename FaceRange>
  Proxy fit_proxy(const FaceRange &faces, const TriangleMesh &tm) const {
    (void)(tm);
    CGAL_assertion(!faces.empty());

    // fitting normal
    Vector_3 norm = CGAL::NULL_VECTOR;
    for(const face_descriptor f : faces) {
      norm = m_sum_functor(norm,
        m_scale_functor(get(m_fnmap, f), get(m_famap, f)));
    }
    if (norm.squared_length() > FT(0.0))
      norm = m_scale_functor(norm,
        FT(1.0) / CGAL::approximate_sqrt(norm.squared_length()));

    return norm;
  }

private:
  Face_normal_map m_fnmap;
  Face_area_map m_famap;
  Construct_scaled_vector_3 m_scale_functor;
  Compute_scalar_product_3 m_scalar_product_functor;
  Construct_sum_of_vectors_3 m_sum_functor;
  Collinear_3 m_collinear_functor;
};

} // namespace Surface_mesh_approximation
} // namespace CGAL

#endif // CGAL_L21_METRIC_PLANE_PROXY_H
