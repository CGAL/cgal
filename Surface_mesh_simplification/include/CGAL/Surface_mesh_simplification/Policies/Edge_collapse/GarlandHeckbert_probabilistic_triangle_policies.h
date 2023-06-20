// Copyright (c) 2019  GeometryFactory (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Baskin Burak Senbaslar,
//                 Mael Rouxel-Labbé,
//                 Julian Komaromy

#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_GARLANDHECKBERT_PROBABILISTIC_TRIANGLE_POLICIES_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_GARLANDHECKBERT_PROBABILISTIC_TRIANGLE_POLICIES_H

#include <CGAL/license/Surface_mesh_simplification.h>

#include <CGAL/Surface_mesh_simplification/internal/Common.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/internal/GarlandHeckbert_policy_base.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/internal/GarlandHeckbert_functions.h>

#include <CGAL/Default.h>
#include <CGAL/property_map.h>

#include <utility>

namespace CGAL {
namespace Surface_mesh_simplification {
namespace internal {

template <typename TriangleMesh, typename GeomTraits, typename FaceVarianceMap>
class Probabilistic_triangle_quadric_calculator
{
  typedef typename GeomTraits::FT                                              FT;
  typedef typename GeomTraits::Point_3                                         Point_3;
  typedef typename GeomTraits::Vector_3                                        Vector_3;

  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor          face_descriptor;

  typedef Constant_property_map<face_descriptor, FT>                           Default_FVM;
  typedef typename Default::Get<FaceVarianceMap, Default_FVM>::type            Face_variance_map;

  typedef typename GarlandHeckbert_matrix_types<GeomTraits>::Mat_4             Mat_4;
  typedef typename GarlandHeckbert_matrix_types<GeomTraits>::Col_4             Col_4;

private:
  // same meaning as for probabilistic plane quadrics
  static constexpr FT default_variance_unit = 0.05;

  // this is only used when we fall back to probabilistic planes for the discontinuous edges,
  // the actual triangle quadric calculation only uses the normal variance
  static constexpr FT position_variance_factor = 1;

  Face_variance_map m_face_variance_map;

public:
  Probabilistic_triangle_quadric_calculator() = delete;

  template <typename FVM>
  Probabilistic_triangle_quadric_calculator(const FVM fvm)
      : m_face_variance_map(fvm)
  { }

  Probabilistic_triangle_quadric_calculator(TriangleMesh& tmesh,
                                            typename boost::enable_if<std::is_same<Face_variance_map, Default_FVM> >::type* = nullptr)
  {
    // try to initialize the face variance map using the estimated variance
    // parameters are constants defined for this class
    FT variance, discard_position;
    std::tie(variance, discard_position) =
      estimate_variances(tmesh, GeomTraits(), default_variance_unit, position_variance_factor);

    // see probabilistic plane quadrics
    m_face_variance_map = Default_FVM { variance };
  }

public:
  // we don't have a sensible way to construct a triangle quadric
  // from an edge, so we fall back to probabilistic plane quadrics here
  template<typename VertexPointMap>
  Mat_4 construct_quadric_from_edge(typename boost::graph_traits<TriangleMesh>::halfedge_descriptor he,
                                    const TriangleMesh& tmesh,
                                    const VertexPointMap point_map,
                                    const GeomTraits& gt) const
  {
    // same as in probabilistic plane policy
    const Vector_3 normal = construct_edge_normal(he, tmesh, point_map, gt);
    const Point_3 p = get(point_map, source(he, tmesh));

    const FT variance = get(m_face_variance_map, face(he, tmesh));

    // @fixme plane?
    return construct_prob_plane_quadric_from_normal(normal, p, gt, variance,
                                                    position_variance_factor * variance);
  }

  template<typename VertexPointMap>
  Mat_4 construct_quadric_from_face(typename boost::graph_traits<TriangleMesh>::face_descriptor f,
                                    const TriangleMesh& tmesh,
                                    const VertexPointMap point_map,
                                    const GeomTraits& gt) const
  {
    const FT variance = get(m_face_variance_map, f);

    return construct_prob_triangle_quadric_from_face(f, variance, tmesh, point_map, gt);
  }

  Col_4 construct_optimal_point(const Mat_4& quadric, const Col_4& /*p0*/, const Col_4& /*p1*/) const
  {
    // @fixme check this
    return construct_optimal_point_invertible<GeomTraits>(quadric);
  }
};

} // namespace internal

// implements probabilistic triangle faces and optionally takes a face variance map
// analogously to probabilistic plane quadrics
template<typename TriangleMesh,
         typename GeomTraits,
         typename FaceVarianceMap = CGAL::Default>
class GarlandHeckbert_probabilistic_triangle_policies
  : public internal::GarlandHeckbert_cost_and_placement<
             internal::Probabilistic_triangle_quadric_calculator<TriangleMesh, GeomTraits, FaceVarianceMap>,
             TriangleMesh, GeomTraits>
{
public:
  typedef internal::Probabilistic_triangle_quadric_calculator<
            TriangleMesh, GeomTraits, FaceVarianceMap>                         Quadric_calculator;

private:
  typedef internal::GarlandHeckbert_cost_and_placement<
            Quadric_calculator, TriangleMesh, GeomTraits>                      Base;
  typedef GarlandHeckbert_probabilistic_triangle_policies<
            TriangleMesh, GeomTraits>                                          Self;

public:
  typedef Self                                                                 Get_cost;
  typedef Self                                                                 Get_placement;

  typedef typename GeomTraits::FT                                              FT;

public:
  // Only available if the quadric calculator is using the default (constant) variance property map
  GarlandHeckbert_probabilistic_triangle_policies(TriangleMesh& tmesh,
                                                  const FT dm = FT(100))
    : Base(tmesh, Quadric_calculator(tmesh), dm)
  { }

  template <typename FVM>
  GarlandHeckbert_probabilistic_triangle_policies(TriangleMesh& tmesh,
                                                  const FT dm,
                                                  const FVM fvm)
    : Base(tmesh, Quadric_calculator(fvm), dm)
  { }

public:
  const Get_cost& get_cost() const { return *this; }
  const Get_placement& get_placement() const { return *this; }

  using Base::operator();
};

} // namespace Surface_mesh_simplification
} // namespace CGAL

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_GARLANDHECKBERT_PROBABILISTIC_TRIANGLE_POLICIES_H
