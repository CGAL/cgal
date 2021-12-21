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

#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_GARLANDHECKBERT_PROBABILISTIC_TRI_POLICIES_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_GARLANDHECKBERT_PROBABILISTIC_TRI_POLICIES_H

#include <CGAL/license/Surface_mesh_simplification.h>

#include <CGAL/Surface_mesh_simplification/internal/Common.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/internal/GarlandHeckbert_policy_base.h>

#include <Eigen/Dense>

#include <CGAL/boost/graph/Named_function_parameters.h>

namespace CGAL {
namespace Surface_mesh_simplification {

// derived class implements functions used in the base class that
// takes the derived class as template argument - see "CRTP"
//
// derives from cost_base and placement_base
// implements probabilistic triangle faces and optionally takes a face variance map
// analogously to probabilistic plane quadrics
template<typename TriangleMesh,
         typename GeomTraits,
         typename FaceVarianceMap =
           Constant_property_map<typename boost::graph_traits<TriangleMesh>::face_descriptor,
                                 typename GeomTraits::FT> >
class GarlandHeckbert_probabilistic_tri_policies
  : public internal::GarlandHeckbert_placement_base<
             typename boost::property_map<
               TriangleMesh,
               CGAL::dynamic_vertex_property_t<Eigen::Matrix<typename GeomTraits::FT, 4, 4, Eigen::DontAlign>>
             >::type,
             GeomTraits,
             GarlandHeckbert_probabilistic_tri_policies<TriangleMesh, GeomTraits, FaceVarianceMap>
           >,
  public internal::GarlandHeckbert_cost_base<
           typename boost::property_map<
             TriangleMesh,
             CGAL::dynamic_vertex_property_t<Eigen::Matrix<typename GeomTraits::FT, 4, 4, Eigen::DontAlign>>
           >::type,
           GeomTraits,
           GarlandHeckbert_probabilistic_tri_policies<TriangleMesh, GeomTraits, FaceVarianceMap> >
{
  typedef typename boost::property_traits<FaceVarianceMap>::value_type         Face_variance;

public:
  typedef typename GeomTraits::FT                                              FT;

  typedef typename Eigen::Matrix<FT, 3, 3, Eigen::DontAlign>                   Mat_3;
  typedef typename Eigen::Matrix<FT, 4, 4, Eigen::DontAlign>                   GH_matrix;
  typedef CGAL::dynamic_vertex_property_t<GH_matrix>                           Cost_property;
  typedef typename boost::property_map<TriangleMesh, Cost_property>::type      Vertex_cost_map;

  typedef internal::GarlandHeckbert_placement_base<
    Vertex_cost_map, GeomTraits, GarlandHeckbert_probabilistic_tri_policies<TriangleMesh, GeomTraits>
    >                                                                          Placement_base;

  typedef internal::GarlandHeckbert_cost_base<
    Vertex_cost_map, GeomTraits, GarlandHeckbert_probabilistic_tri_policies<TriangleMesh, GeomTraits>
    >                                                                          Cost_base;

  // both types are the same, this is so we avoid casting back to the base class in
  // get_cost() or get_placement()
  typedef GarlandHeckbert_probabilistic_tri_policies                           Get_cost;
  typedef GarlandHeckbert_probabilistic_tri_policies                           Get_placement;

  // these using directives are needed to choose between the definitions of these types
  // in Cost_base and Placement_base (even though they are the same)
  using typename Cost_base::Mat_4;
  using typename Cost_base::Col_4;
  using typename Cost_base::Point_3;
  using typename Cost_base::Vector_3;

private:
  Vertex_cost_map vcm_;
  FaceVarianceMap face_variance_map;

  // same meaning as for probabilistic plane quadrics
  static constexpr FT good_default_variance_unit = 0.05;

  // this is only used when we fall back to probabilistic planes for the discontinuous edges,
  // the actual triangle quadric calculation only uses the normal variance
  static constexpr FT position_variance_factor = 1;

public:
  // default discontinuity multiplier is 100
  GarlandHeckbert_probabilistic_tri_policies(TriangleMesh& tmesh)
    : GarlandHeckbert_probabilistic_tri_policies(tmesh, 100)
  { }

  GarlandHeckbert_probabilistic_tri_policies(TriangleMesh& tmesh, FT dm)
    : Cost_base(dm)
  {
    // initialize the private variable vcm so it's lifetime is bound to that of the policy's
    vcm_ = get(Cost_property(), tmesh);

    // initialize both vcms
    Cost_base::init_vcm(vcm_);
    Placement_base::init_vcm(vcm_);

    FT variance;
    FT discard_position;

    std::tie(variance, discard_position) = internal::estimate_variances(tmesh, GeomTraits(),
        good_default_variance_unit, position_variance_factor);

    // see probabilistic plane quadrics
    face_variance_map = FaceVarianceMap { variance };
  }

  GarlandHeckbert_probabilistic_tri_policies(TriangleMesh& tmesh,
                                             FT dm,
                                             const FaceVarianceMap* fvm)
    : Cost_base(dm), face_variance_map(fvm)
  {
    // we need positive variances so that we always get an invertible matrix
//    CGAL_precondition(sdn > 0. && sdp > 0.); // @fixme what was that, check history

    // initialize the private variable vcm so it's lifetime is bound to that of the policy's
    vcm_ = get(Cost_property(), tmesh);

    // initialize both vcms
    Cost_base::init_vcm(vcm_);
    Placement_base::init_vcm(vcm_);
  }

  const Get_cost& get_cost() const { return *this; }
  const Get_placement& get_placement() const { return *this; }

  // so that operator() gets overloaded, this is needed because now Get_cost and Get_placement are the same
  using Cost_base::operator();
  using Placement_base::operator();

public:
  template<typename VPM, typename TM>
  Mat_4 construct_quadric_from_face(typename boost::graph_traits<TM>::face_descriptor f,
                                    const TM& tmesh,
                                    const VPM point_map,
                                    const GeomTraits& gt) const
  {
    FT variance = get(face_variance_map, f);

    return internal::construct_prob_triangle_quadric_from_face(f, variance, tmesh, point_map, gt);
  }

  // we don't have a sensible way to construct a triangle quadric
  // from an edge, so we fall back to probabilistic plane quadrics here
  template<typename VPM, typename TM>
  Mat_4 construct_quadric_from_edge(typename boost::graph_traits<TM>::halfedge_descriptor he,
                                    const TM& tmesh,
                                    const VPM point_map,
                                    const GeomTraits& gt) const
  {
    // same as in probabilistic plane policy
    const Vector_3 normal = internal::construct_edge_normal(he, tmesh, point_map, gt);
    const Point_3 p = get(point_map, source(he, tmesh));

    FT variance = get(face_variance_map, face(he, tmesh));

    return internal::construct_prob_plane_quadric_from_normal(normal, p, gt, variance,
                                                              position_variance_factor * variance);
  }

  Col_4 construct_optimal_point(const Mat_4 quadric, const Col_4& p0, const Col_4& p1) const
  {
    return internal::construct_optimal_point_invertible<GeomTraits>(quadric);
  }
};

} // namespace Surface_mesh_simplification
} // namespace CGAL

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_GARLANDHECKBERT_PROBABILISTIC_TRI_POLICIES_H
