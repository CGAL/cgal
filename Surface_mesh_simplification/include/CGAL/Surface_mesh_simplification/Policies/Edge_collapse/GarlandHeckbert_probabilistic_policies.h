// Copyright (c) 2019  GeometryFactory (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Baskin Burak Senbaslar,
//                 Mael Rouxel-Labbé,
//                 Julian Komaromy

#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_GARLANDHECKBERT_PROBABILISTIC_POLICIES_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_GARLANDHECKBERT_PROBABILISTIC_POLICIES_H

#include <CGAL/license/Surface_mesh_simplification.h>
#include <CGAL/Surface_mesh_simplification/internal/Common.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/internal/GarlandHeckbert_policy_base.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/internal/GarlandHeckbert_functions.h>
#include <Eigen/Dense>
  
namespace CGAL {
namespace Surface_mesh_simplification {

// derived class implements functions used in the base class that
// takes the derived class as template argument - see "CRTP"
//
// implements probabilistic plane quadrics, optionally takes a face variance map
// giving a per-face variance 
template<typename TriangleMesh, typename GeomTraits, typename 
  FaceVarianceMap = Constant_property_map<
    typename boost::graph_traits<TriangleMesh>::face_descriptor, 
             std::pair<typename GeomTraits::FT, typename GeomTraits::FT>
  >> 
class GarlandHeckbert_probabilistic_policies :
  public internal::GarlandHeckbert_placement_base<
    typename boost::property_map<
      TriangleMesh,
      CGAL::dynamic_vertex_property_t<Eigen::Matrix<typename GeomTraits::FT, 4, 4, Eigen::DontAlign>>
    >::type,
    GeomTraits, 
    GarlandHeckbert_probabilistic_policies<TriangleMesh, GeomTraits, FaceVarianceMap>
  >,
  public internal::GarlandHeckbert_cost_base<
    typename boost::property_map<
      TriangleMesh,
      CGAL::dynamic_vertex_property_t<Eigen::Matrix<typename GeomTraits::FT, 4, 4, Eigen::DontAlign>>
    >::type,
    GeomTraits,
    GarlandHeckbert_probabilistic_policies<TriangleMesh, GeomTraits, FaceVarianceMap>
  >
{
  
  typedef typename GeomTraits::FT FT;
  
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor; 
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor; 
  
  typedef typename boost::property_traits<FaceVarianceMap>::value_type Face_variance;

  public:

  typedef typename Eigen::Matrix<FT, 4, 4, Eigen::DontAlign> GH_matrix;
  typedef CGAL::dynamic_vertex_property_t<GH_matrix> Cost_property;

  typedef typename boost::property_map<TriangleMesh, Cost_property>::type Vertex_cost_map;

  typedef internal::GarlandHeckbert_placement_base<
    Vertex_cost_map, GeomTraits, GarlandHeckbert_probabilistic_policies<TriangleMesh, GeomTraits>
    > Placement_base;

  typedef internal::GarlandHeckbert_cost_base<
    Vertex_cost_map, GeomTraits, GarlandHeckbert_probabilistic_policies<TriangleMesh, GeomTraits>
    > Cost_base;

  // both types are the same, this is so we avoid casting back to the base class in
  // get_cost() or get_placement()
  typedef GarlandHeckbert_probabilistic_policies Get_cost;
  typedef GarlandHeckbert_probabilistic_policies Get_placement;

  // so that operator() gets overloaded, this is needed because now Get_cost and Get_placement
  // are the same
  using Cost_base::operator();
  using Placement_base::operator();

  // these using directives are needed to choose between the definitions of these types
  // in Cost_base and Placement_base (even though they are the same)
  using typename Cost_base::Mat_4;
  using typename Cost_base::Col_4;
  using typename Cost_base::Point_3;
  using typename Cost_base::Vector_3;

  // default discontinuity multiplier is 100
  GarlandHeckbert_probabilistic_policies(TriangleMesh& tmesh)
    : GarlandHeckbert_probabilistic_policies(tmesh, 100) 
  { }
  
  GarlandHeckbert_probabilistic_policies(TriangleMesh& tmesh, FT dm)
    : Cost_base(dm)
  { 
    // initialize the private variable vcm so it's lifetime is bound to that of the policy's
    vcm_ = get(Cost_property(), tmesh);

    // initialize both vcms
    Cost_base::init_vcm(vcm_);
    Placement_base::init_vcm(vcm_);
    
    // try to initialize the face variance map using the estimated variance
    // parameters are constants defined for this class
    face_variance_map = FaceVarianceMap { 
      internal::estimate_variances(tmesh, GeomTraits(), good_default_variance_unit, 
          position_variance_factor) };  
  }
  
  GarlandHeckbert_probabilistic_policies(TriangleMesh& tmesh, FT dm, const FaceVarianceMap& fvm)
    : Cost_base(dm), face_variance_map(fvm)
  { 
     // initialize the private variable vcm so it's lifetime is bound to that of the policy's
    vcm_ = get(Cost_property(), tmesh);

    // initialize both vcms
    Cost_base::init_vcm(vcm_);
    Placement_base::init_vcm(vcm_);
  }
  
  template<typename VPM>
  Mat_4 construct_quadric_from_face(
      const VPM& point_map,
      const TriangleMesh& tmesh, 
      face_descriptor f, 
      const GeomTraits& gt) const
  {
    const Vector_3 normal = internal::construct_unit_normal_from_face<
      VPM, TriangleMesh, GeomTraits>(point_map, tmesh, f, gt);
    
    const Point_3 p = get(point_map, source(halfedge(f, tmesh), tmesh));

    FT n_variance;
    FT p_variance;

    std::tie(n_variance, p_variance) = get(face_variance_map, f);
    
    return internal::construct_prob_plane_quadric_from_normal(normal, p, gt, n_variance, p_variance);
  }
  
  template<typename VPM>
  Mat_4 construct_quadric_from_edge(
      const VPM& point_map,
      const TriangleMesh& tmesh, 
      halfedge_descriptor he, 
      const GeomTraits& gt) const
  {
    const Vector_3 normal = internal::construct_edge_normal(point_map, tmesh, he, gt);

    const Point_3 p = get(point_map, source(he, tmesh));

    FT n_variance;
    FT p_variance;

    std::tie(n_variance, p_variance) = get(face_variance_map, face(he, tmesh));
    
    return internal::construct_prob_plane_quadric_from_normal(normal, p, gt, n_variance, p_variance);
  }
  
  Col_4 construct_optimal_point(const Mat_4& quadric, const Col_4& p0, 
      const Col_4& p1) const
  {
    return internal::construct_optimal_point_invertible<GeomTraits>(quadric);
  }
  
  const Get_cost& get_cost() const { return *this; }
  
  const Get_placement& get_placement() const { return *this; }
  
  private:
  Vertex_cost_map vcm_;
  
  // magic number determined by some testing 
  static constexpr FT good_default_variance_unit = 0.05;

  // magic number - for most use cases, there is no input variance, so it makes sense to
  // set the positional variance to a smaller value than the normal variance
  static constexpr FT position_variance_factor = 0.1;
    
  FaceVarianceMap face_variance_map;
};

} // namespace Surface_mesh_simplification
} // namespace CGAL

#endif //CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_GARLANDHECKBERT_PROBABILISTIC_POLICIES_H
