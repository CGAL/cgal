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

#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_GARLANDHECKBERT_POLICIES_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_GARLANDHECKBERT_POLICIES_H

#include <CGAL/license/Surface_mesh_simplification.h>

#include <CGAL/Surface_mesh_simplification/internal/Common.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/internal/GarlandHeckbert_policy_base.h>

#include <Eigen/Dense>

namespace CGAL {
namespace Surface_mesh_simplification {

// derived class implements functions used in the base class that
// takes the derived class as template argument - see "CRTP"
template<typename TriangleMesh, typename GeomTraits>
class GarlandHeckbert_policies : 
  public internal::GarlandHeckbert_placement_base<
    typename boost::property_map<
      TriangleMesh, 
      CGAL::dynamic_vertex_property_t<Eigen::Matrix<typename GeomTraits::FT, 4, 4, Eigen::DontAlign>>
    >::type,
    GeomTraits,
    GarlandHeckbert_policies<TriangleMesh, GeomTraits>
  >,
  public internal::GarlandHeckbert_cost_base<
    typename boost::property_map<
      TriangleMesh, 
      CGAL::dynamic_vertex_property_t<Eigen::Matrix<typename GeomTraits::FT, 4, 4, Eigen::DontAlign>>
    >::type,
    GeomTraits,
    GarlandHeckbert_policies<TriangleMesh, GeomTraits>
  >
{

  public:
    typedef typename GeomTraits::FT FT;
    
    typedef typename Eigen::Matrix<FT, 4, 4, Eigen::DontAlign> GH_matrix;
    typedef CGAL::dynamic_vertex_property_t<GH_matrix> Cost_property;
    
    typedef typename boost::property_map<TriangleMesh, Cost_property>::type Vertex_cost_map;

    typedef internal::GarlandHeckbert_placement_base<
      Vertex_cost_map, GeomTraits, GarlandHeckbert_policies<TriangleMesh, GeomTraits>
      > Placement_base;

    typedef internal::GarlandHeckbert_cost_base<
      Vertex_cost_map, GeomTraits, GarlandHeckbert_policies<TriangleMesh, GeomTraits>
      > Cost_base;
    
    // both types are the same, this is so we avoid casting back to the base class in
    // get_cost() or get_placement()
    typedef GarlandHeckbert_policies Get_cost;
    typedef GarlandHeckbert_policies Get_placement;
    
    // introduce both operators into scope so we get an overload operator()
    // this is needed since Get_cost and Get_placement are the same type
    using Cost_base::operator();
    using Placement_base::operator();
      
    // these using directives are needed to choose between the definitions of these types
    // in Cost_base and Placement_base (even though they are the same)
    using typename Cost_base::Mat_4;
    using typename Cost_base::Col_4;
    using typename Cost_base::Point_3;
    using typename Cost_base::Vector_3;

    GarlandHeckbert_policies(TriangleMesh& tmesh, FT dm = FT(100)) 
    {
      vcm_ = get(Cost_property(), tmesh);

      /**
       * initialize the two base class cost matrices (protected members)
       */
      Cost_base::init_vcm(vcm_);
      Placement_base::init_vcm(vcm_);
    }

    Col_4 construct_optimal_point(const Mat_4& aQuadric, const Col_4& p0, const Col_4& p1) const 
    {
      Mat_4 X;
      X << aQuadric.block(0, 0, 3, 4), 0, 0, 0, 1;

      Col_4 opt_pt;

      if(X.determinant() == 0)
      {
        // not invertible
        const Col_4 p1mp0 = std::move(p1 - p0);
        const FT a = (p1mp0.transpose() * aQuadric * p1mp0)(0, 0);
        const FT b = 2 * (p0.transpose() * aQuadric * p1mp0)(0, 0);

        if(a == 0)
        {
          if(b < 0)
            opt_pt = p1;
          else if(b == 0)
            opt_pt = 0.5 * (p0 + p1);
          else
            opt_pt = p0;
        }
        else
        {
          FT ext_t = -b/(2*a);
          if(ext_t < 0 || ext_t > 1 || a < 0)
          {
            // one of endpoints
            FT p0_cost = (p0.transpose() * aQuadric * p0)(0, 0);
            FT p1_cost = (p1.transpose() * aQuadric * p1)(0, 0);
            if(p0_cost > p1_cost)
              opt_pt = p1;
            else
              opt_pt = p0;
          }
          else
          {
            // extremum of the parabola
            opt_pt = p0 + ext_t * (p1 - p0);
          }
        }
      }
      else // invertible
      {
        opt_pt = X.inverse().col(3); // == X.inverse() * (0 0 0 1)
      }
      return opt_pt;
    }

    template<typename VPM, typename TM> 
    Mat_4 construct_quadric_from_face(
        const VPM& point_map,
        const TM& tmesh, 
        typename boost::graph_traits<TM>::face_descriptor f, 
        const GeomTraits& gt) const
    {
      // this-> is used because we are calling a base class function
      const Vector_3 normal = this->construct_unit_normal_from_face(point_map, tmesh, f, gt);

      // get any point of the face 
      const auto p = get(point_map, source(halfedge(f, tmesh), tmesh));
      
      return construct_quadric_from_normal(normal, p, gt);
    }
    
    template<typename VPM, typename TM> 
    Mat_4 construct_quadric_from_edge(
        const VPM& point_map,
        const TM& tmesh, 
        typename boost::graph_traits<TM>::halfedge_descriptor he, 
        const GeomTraits& gt) const
    {
      const Vector_3 normal = this->construct_edge_normal(point_map, tmesh, he, gt);
      
      return construct_quadric_from_normal(normal, get(point_map, source(he, tmesh)), gt);
    }
    
    const Get_cost& get_cost() const {
      return *this;
    }
    
    const Get_placement& get_placement() const {
      return *this;
    }

  private:
  Vertex_cost_map vcm_;
  
  // helper function to construct quadrics

  Mat_4 construct_quadric_from_normal(const Vector_3& normal, const Point_3& point, 
      const GeomTraits& gt) const 
  {

    auto dot_product = gt.compute_scalar_product_3_object();
    auto construct_vector = gt.construct_vector_3_object();

    // negative dot product between the normal and the position vector
    const FT d = - dot_product(normal, construct_vector(ORIGIN, point));

    // row vector given by d appended to the normal
    const Eigen::Matrix<FT, 1, 4> row(normal.x(), normal.y(), normal.z(), d);

    // outer product
    return row.transpose() * row;
  }
};

} //namespace Surface_mesh_simplification
} //namespace CGAL

#endif //CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_GARLANDHECKBERT_POLICIES_H
