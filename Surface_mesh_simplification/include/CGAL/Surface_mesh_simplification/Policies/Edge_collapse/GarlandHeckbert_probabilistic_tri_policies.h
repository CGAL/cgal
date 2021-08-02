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
template<typename TriangleMesh, typename GeomTraits, typename
  FaceVarianceMap = Constant_property_map<
  typename boost::graph_traits<TriangleMesh>::face_descriptor, typename GeomTraits::FT>>
class GarlandHeckbert_probabilistic_tri_policies :
  public internal::GarlandHeckbert_placement_base<
    typename boost::property_map<
      TriangleMesh,
      CGAL::dynamic_vertex_property_t<Eigen::Matrix<typename GeomTraits::FT, 4, 4, Eigen::DontAlign>>
    >::type,
    GeomTraits,
    GarlandHeckbert_probabilistic_tri_policies<TriangleMesh, GeomTraits>
  >,
  public internal::GarlandHeckbert_cost_base<
    typename boost::property_map<
      TriangleMesh,
      CGAL::dynamic_vertex_property_t<Eigen::Matrix<typename GeomTraits::FT, 4, 4, Eigen::DontAlign>>
    >::type,
    GeomTraits,
    GarlandHeckbert_probabilistic_tri_policies<TriangleMesh, GeomTraits>
  >
{

  typedef typename boost::property_traits<FaceVarianceMap>::value_type Face_variance;

  public:
    typedef typename GeomTraits::FT FT;

    typedef typename Eigen::Matrix<FT, 3, 3, Eigen::DontAlign> Mat_3;
    typedef typename Eigen::Matrix<FT, 4, 4, Eigen::DontAlign> GH_matrix;
    typedef CGAL::dynamic_vertex_property_t<GH_matrix> Cost_property;

    typedef typename boost::property_map<TriangleMesh, Cost_property>::type Vertex_cost_map;

    typedef internal::GarlandHeckbert_placement_base<
      Vertex_cost_map, GeomTraits, GarlandHeckbert_probabilistic_tri_policies<TriangleMesh, GeomTraits>
      > Placement_base;

    typedef internal::GarlandHeckbert_cost_base<
      Vertex_cost_map, GeomTraits, GarlandHeckbert_probabilistic_tri_policies<TriangleMesh, GeomTraits>
      > Cost_base;

    // both types are the same, this is so we avoid casting back to the base class in
    // get_cost() or get_placement()
    typedef GarlandHeckbert_probabilistic_tri_policies Get_cost;
    typedef GarlandHeckbert_probabilistic_tri_policies Get_placement;

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

    // TODO use property maps here as well
    // default discontinuity multiplier is 100
    GarlandHeckbert_probabilistic_tri_policies(TriangleMesh& tmesh)
      : GarlandHeckbert_probabilistic_tri_policies(tmesh, 100) 
    { }
    
    GarlandHeckbert_probabilistic_tri_policies(TriangleMesh& tmesh, FT dm)
      : Cost_base(dm)
    { 
      // initialize the private variable vcm so it's lifetime is bound to that of the policy's
      vcm_ = get(Cost_property(), tmesh);

      FT variance;
      std::tie(variance, variance) = estimate_variances(tmesh, GeomTraits());

      // initialize both vcms
      Cost_base::init_vcm(vcm_);
      Placement_base::init_vcm(vcm_);

      face_variance_map = FaceVarianceMap { variance };
    }

    GarlandHeckbert_probabilistic_tri_policies(TriangleMesh& tmesh,
                                           FT dm,
                                           const FaceVarianceMap* fvm) 
      : Cost_base(dm), face_variance_map(fvm)
    {

      // we need positive variances so that we always get an invertible matrix
      CGAL_precondition(sdn > 0.0 && sdp > 0.0);
      
      // initialize the private variable vcm so it's lifetime is bound to that of the policy's
      vcm_ = get(Cost_property(), tmesh);

      // initialize both vcms
      Cost_base::init_vcm(vcm_);
      Placement_base::init_vcm(vcm_);
    }
        
    template<typename VPM, typename TM>
    Mat_4 construct_quadric_from_face(
        const VPM& point_map,
        const TM& tmesh, 
        typename boost::graph_traits<TM>::face_descriptor f, 
        const GeomTraits& gt) const
    {
        return internal::construct_prob_triangle_quadric_from_face(
            point_map, tmesh, f, get(face_variance_map, f), gt);
    }

    template<typename VPM, typename TM>
    Mat_4 construct_quadric_from_edge(
        const VPM& point_map,
        const TM& tmesh, 
        typename boost::graph_traits<TM>::halfedge_descriptor he, 
        const GeomTraits& gt) const
    {
      // same as in probabilistic plane policy
      const Vector_3 normal = internal::construct_edge_normal(point_map, tmesh, he, gt);

      const Point_3 p = get(point_map, source(he, tmesh));

      FT variance = get(face_variance_map, face(he, tmesh));
      
      return internal::construct_prob_plane_quadric_from_normal(normal, p, gt, variance, variance);
    }
    
    Col_4 construct_optimal_point(const Mat_4 quadric, const Col_4& p0, const Col_4& p1) const
    {
      return internal::construct_optimal_point_invertible<GeomTraits>(quadric);
    }
    const Get_cost& get_cost() const { return *this; }
    const Get_placement& get_placement() const { return *this; }
    
  private:
    Vertex_cost_map vcm_;
    
    FaceVarianceMap face_variance_map;
    
   // give a very rough estimate of a decent variance for both parameters
    static std::pair<FT, FT> estimate_variances(const TriangleMesh& mesh, 
        const GeomTraits& gt)
    {
      typedef typename TriangleMesh::Vertex_index vertex_descriptor;
      typedef typename TriangleMesh::Edge_index edge_descriptor;
      
      FT average_edge_length = 0;
      
      auto construct_vector = gt.construct_vector_3_object();
      
      for (edge_descriptor e : edges(mesh))
      {
        vertex_descriptor v1 = mesh.vertex(e, 0);
        vertex_descriptor v2 = mesh.vertex(e, 1);

        const Point_3& p1 = mesh.point(v1); 
        const Point_3& p2 = mesh.point(v2); 

        const Vector_3 vec = construct_vector(p1, p2);
        average_edge_length += sqrt(vec.squared_length());
      }
      
      average_edge_length = average_edge_length / mesh.number_of_edges();
      const FT n2 = 0.05 * average_edge_length;
      const FT p2 = 0.05 * average_edge_length;
      
      return std::make_pair(n2, p2);
    }

    // magic number determined by some testing, this is a good variance for models that
    // fit inside a [-1, 1]^3 unit cube
    static constexpr FT good_default_variance_unit = 0.05;
};

} // namespace Surface_mesh_simplification
} // namespace CGAL

#endif //CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_GARLANDHECKBERT_PROBABILISTIC_TRI_POLICIES_H
