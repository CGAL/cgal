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
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/internal/GarlandHeckbert_optimizers.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/internal/GarlandHeckbert_common.h>
#include <Eigen/Dense>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/internal/GarlandHeckbert_plane_quadrics.h>

namespace CGAL {
namespace Surface_mesh_simplification {

// derived class implements functions used in the base class that
// takes the derived class as template argument - see "CRTP"
//
// derives from cost_base and placement_base
template<typename TriangleMesh, typename GeomTraits>
class GarlandHeckbert_probabilistic_policies :
  public internal::GarlandHeckbert_plane_edges<TriangleMesh, GeomTraits>,
  public internal::GarlandHeckbert_invertible_optimizer<GeomTraits>,
  public internal::GarlandHeckbert_placement_base<
    typename boost::property_map<
      TriangleMesh,
      CGAL::dynamic_vertex_property_t<Eigen::Matrix<typename GeomTraits::FT, 4, 4, Eigen::DontAlign>>
    >::type,
    GeomTraits,
    GarlandHeckbert_probabilistic_policies<TriangleMesh, GeomTraits>
  >,
  public internal::GarlandHeckbert_cost_base<
    typename boost::property_map<
      TriangleMesh,
      CGAL::dynamic_vertex_property_t<Eigen::Matrix<typename GeomTraits::FT, 4, 4, Eigen::DontAlign>>
    >::type,
    GeomTraits,
    GarlandHeckbert_probabilistic_policies<TriangleMesh, GeomTraits>
  >
{

  public:
    typedef typename GeomTraits::FT FT;

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

      std::tie(normal_variance, mean_variance) = estimate_variances(tmesh);
      
      // initialize both vcms
      Cost_base::init_vcm(vcm_);
      Placement_base::init_vcm(vcm_);
    }

    GarlandHeckbert_probabilistic_policies(TriangleMesh& tmesh,
                                           FT dm,
                                           FT sdn,
                                           FT sdp) : Cost_base(dm)
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
      const Vector_3 normal = internal::common::construct_unit_normal_from_face<
        GeomTraits, VPM, TM>(point_map, tmesh, f, gt);
      const Point_3 p = get(point_map, source(halfedge(f, tmesh), tmesh));

      return construct_quadric_from_normal(normal, p, gt);
    }
    
    const Get_cost& get_cost() const { return *this; }
    const Get_placement& get_placement() const { return *this; }
    
  private:
    Vertex_cost_map vcm_;
   
    Mat_4 construct_quadric_from_normal(
        const Vector_3& mean_normal,
        const Point_3& point, 
        const GeomTraits& gt) const
    {
      auto squared_length = gt.compute_squared_length_3_object();
      auto dot_product = gt.compute_scalar_product_3_object();
      auto construct_vec_3 = gt.construct_vector_3_object();

      const Vector_3 mean_vec = construct_vec_3(ORIGIN, point);
      const FT dot_mnmv = dot_product(mean_normal, mean_vec);

      // Eigen column vector of length 3
      const Eigen::Matrix<FT, 3, 1> mean_n_col{mean_normal.x(), mean_normal.y(), mean_normal.z()};

      // start by setting values along the diagonal
      Mat_4 mat = normal_variance * Mat_4::Identity();

      // add outer product of the mean normal with itself
      // to the upper left 3x3 block
      mat.block(0, 0, 3, 3) += mean_n_col * mean_n_col.transpose();

      // set the first 3 values of the last row and the first
      // 3 values of the last column
      // the negative sign comes from the fact that in the paper,
      // the b column and row appear with a negative sign
      const auto b1 = -(dot_mnmv * mean_normal + normal_variance * mean_vec);

      const Eigen::Matrix<FT, 3, 1> b {b1.x(), b1.y(), b1.z()};

      mat.col(3).head(3) = b;
      mat.row(3).head(3) = b.transpose();

      // set the value in the bottom right corner, we get this by considering
      // that we only have single variances given instead of covariance matrices
      mat(3, 3) = CGAL::square(dot_mnmv)
        + normal_variance * squared_length(mean_vec)
        + mean_variance * squared_length(mean_normal)
        + 3 * normal_variance * mean_variance;

      return mat;
    }
    // give a very rough estimate of a decent variance for both parameters
    static std::pair<FT, FT> estimate_variances(const TriangleMesh& mesh)
    {
      typedef typename TriangleMesh::Vertex_index vertex_descriptor;
      // get the bounding box of the mesh
      CGAL::Bbox_3 bbox { };

      for (vertex_descriptor v : vertices(mesh))
      {
        bbox += mesh.point(v).bbox();
      }

      // calculate geometric mean of edge lengths
      FT geometric_mean = 1;

      for (int i = 0; i < 3; ++i)
      {
        geometric_mean *= abs(bbox.max(i) - bbox.min(i));
      }

      geometric_mean = pow(geometric_mean, 1/3);

      // set the variances based on the default value for [-1, 1]^3 cubes
      // here we have the maximum edge length, so if we scale
      // our mesh down by max_edge_length it fits inside a [0, 1]^3 cube, hence the "/ 2"
      const FT n2 = good_default_variance_unit * geometric_mean / 2.0;
      const FT p2 = good_default_variance_unit * geometric_mean / 2.0;
      
      return std::make_pair(n2, p2);
    }

    // magic number determined by some testing, this is a good variance for models that
    // fit inside a [-1, 1]^3 unit cube
    static constexpr FT good_default_variance_unit = 1e-4;
    
  private:
    FT mean_variance;
    FT normal_variance;
};

} // namespace Surface_mesh_simplification
} // namespace CGAL

#endif //CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_GARLANDHECKBERT_PROBABILISTIC_POLICIES_H
