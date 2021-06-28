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
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/internal/GarlandHeckbert_core.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/internal/GarlandHeckbert_policy_base.h>

#include <Eigen/Dense>

namespace CGAL {
namespace Surface_mesh_simplification {

// derived class implements functions used in the base class that
// takes the derived class as template argument - see "CRTP"
//
// derives from cost_base and placement_base
template<typename TriangleMesh, typename GeomTraits>
class GarlandHeckbert_probabilistic_policies : 
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
    // TODO alternatives - e.g. rename base class types so they don't clash
    using typename Cost_base::Mat_4;
    using typename Cost_base::Col_4;
    using typename Cost_base::Point_3;
    using typename Cost_base::Vector_3;


    // TODO good default values
    GarlandHeckbert_probabilistic_policies(
        TriangleMesh& tmesh, 
        FT dm,
        FT sdn,
        FT sdp) 
      : sdev_n_2(square(sdn)), sdev_p_2(square(sdp))
    {
      // initialize the private variable vcm so it's lifetime is bound to that of the policy's
      vcm_ = get(Cost_property(), tmesh);

      // initialize both vcms
      Cost_base::init_vcm(vcm_);
      Placement_base::init_vcm(vcm_);
    }

    Col_4 construct_optimal_point(const Mat_4& aQuadric, const Col_4& p0, const Col_4& p1) const 
    {
      Mat_4 X;
      X << aQuadric.block(0, 0, 3, 4), 0, 0, 0, 1;

      Col_4 opt_pt;

      opt_pt = X.inverse().col(3); // == X.inverse() * (0 0 0 1)
      
      return opt_pt;

    }

    Mat_4 construct_quadric_from_normal(const Vector_3& mean_normal, const Point_3& point,
        const GeomTraits& gt) const {
      auto squared_length = gt.compute_squared_length_3_object();
      auto dot_product = gt.compute_scalar_product_3_object();
      auto construct_vec_3 = gt.construct_vector_3_object();

      const Vector_3 mean_vec = construct_vec_3(ORIGIN, point);
      const FT dot_mnmv = dot_product(mean_normal, mean_vec);

      // Eigen column vector of length 3
      const Eigen::Matrix<FT, 3, 1> mean_n_col{mean_normal.x(), mean_normal.y(), mean_normal.z()};

      // start by setting values along the diagonal
      Mat_4 mat = sdev_n_2 * Mat_4::Identity();

      // add outer product of the mean normal with itself
      // to the upper left 3x3 block
      mat.block(0, 0, 3, 3) += mean_n_col * mean_n_col.transpose();

      // set the first 3 values of the last row and the first 
      // 3 values of the last column
      //TODO why do we have to flip this sign as well? Probably linked to
      // our normal orientation
      const auto b1 = -(dot_mnmv * mean_normal + sdev_n_2 * mean_vec);
      
      const Eigen::Matrix<FT, 3, 1> b {b1.x(), b1.y(), b1.z()};
      
      mat.col(3).head(3) = b;
      mat.row(3).head(3) = b.transpose();

      // set the value in the bottom right corner
      mat(3, 3) = CGAL::square(dot_mnmv) 
        + sdev_n_2 * squared_length(mean_vec) 
        + sdev_p_2 * squared_length(mean_normal)
        + 3 * sdev_n_2 + sdev_p_2;

      return mat;
    }

    const Get_cost& get_cost() {
      return *this;
    }

    const Get_placement& get_placement() {
      return *this;
    }
    
  private:
    // the only parameters we need are the normal and position variances 
    // - ie the squared standard deviations
    FT sdev_n_2;
    FT sdev_p_2;

    Vertex_cost_map vcm_;
};

} //namespace Surface_mesh_simplification
} //namespace CGAL

#endif //CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_GARLANDHECKBERT_PROBABILISTIC_POLICIES_H
