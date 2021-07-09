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
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/internal/GarlandHeckbert_optimizers.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/internal/GarlandHeckbert_common.h>
#include <Eigen/Dense>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/internal/GarlandHeckbert_plane_quadrics.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>

#include <CGAL/boost/graph/Named_function_parameters.h>

namespace CGAL {
namespace Surface_mesh_simplification {

// derived class implements functions used in the base class that
// takes the derived class as template argument - see "CRTP"
//
// derives from cost_base and placement_base
template<typename TriangleMesh, typename GeomTraits>
class GarlandHeckbert_probabilistic_tri_policies :
  public internal::GarlandHeckbert_plane_edges<TriangleMesh, GeomTraits>,
  public internal::GarlandHeckbert_invertible_optimizer<GeomTraits>,
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

    // default discontinuity multiplier is 100
    GarlandHeckbert_probabilistic_tri_policies(TriangleMesh& tmesh)
      : GarlandHeckbert_probabilistic_tri_policies(tmesh, 100) 
    { }
    
    GarlandHeckbert_probabilistic_tri_policies(TriangleMesh& tmesh, FT dm)
      : Cost_base(dm)
    { 
      // initialize the private variable vcm so it's lifetime is bound to that of the policy's
      vcm_ = get(Cost_property(), tmesh);

      std::tie(normal_variance, mean_variance) = estimate_variances(tmesh, GeomTraits());
      
      // initialize both vcms
      Cost_base::init_vcm(vcm_);
      Placement_base::init_vcm(vcm_);
    }

    GarlandHeckbert_probabilistic_tri_policies(TriangleMesh& tmesh,
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
    
    // overload for default parameters
    template<typename VPM, typename TM>
    Mat_4 construct_quadric_from_face(
        const VPM& point_map,
        const TM& tmesh, 
        typename boost::graph_traits<TM>::face_descriptor f, 
        const GeomTraits& gt) const
    {
      return construct_quadric_from_face(point_map, tmesh, f, gt, 
          CGAL::parameters::all_default());
    }
    
    template<typename VPM, typename TM, typename NamedParameters>
    Mat_4 construct_quadric_from_face(
        const VPM& point_map,
        const TM& tmesh, 
        typename boost::graph_traits<TM>::face_descriptor f, 
        const GeomTraits& gt,
        const NamedParameters& np) const
    {
      // TODO much of this is the same as for classical triangle quadrics
      auto construct_vector = gt.construct_vector_3_object();
      auto cross_product = gt.construct_cross_product_vector_3_object();
      auto sum_vectors = gt.construct_sum_of_vectors_3_object();
      auto dot_product = gt.compute_scalar_product_3_object();

      typedef typename boost::property_traits<VPM>::reference Point_reference;
      typedef typename boost::graph_traits<TM>::halfedge_descriptor halfedge_descriptor;

      const FT var = mean_variance;
      const halfedge_descriptor h = halfedge(f, tmesh);

      // get all points and turn them into location vectors so we can use cross product on them
      const Point_reference p = get(point_map, source(h, tmesh));
      const Point_reference q = get(point_map, target(h, tmesh));
      const Point_reference r = get(point_map, target(next(h, tmesh), tmesh));

      Vector_3 a = construct_vector(ORIGIN, p);
      Vector_3 b = construct_vector(ORIGIN, q);
      Vector_3 c = construct_vector(ORIGIN, r);
       
      // calculate certain vectors used later
      const Vector_3 ab = cross_product(a, b);
      const Vector_3 bc = cross_product(b, c);
      const Vector_3 ca = cross_product(c, a);

      const Vector_3 a_minus_b = sum_vectors(a, -b);
      const Vector_3 b_minus_c = sum_vectors(b, -c);
      const Vector_3 c_minus_a = sum_vectors(c, -a);
      
      const Mat_3 cp_ab = skew_sym_mat_cross_product(a_minus_b);
      const Mat_3 cp_bc = skew_sym_mat_cross_product(b_minus_c);
      const Mat_3 cp_ca = skew_sym_mat_cross_product(c_minus_a);
      
      const Vector_3 sum_of_cross_product = sum_vectors(sum_vectors(ab, bc), ca);

      const Eigen::Matrix<FT, 3, 1, Eigen::DontAlign> 
        sum_cp_col{ sum_of_cross_product.x(), sum_of_cross_product.y(), sum_of_cross_product.z() };
      
      Mat_3 A = sum_cp_col * sum_cp_col.transpose();
      A += var * (cp_ab * cp_ab.transpose() + cp_bc * cp_bc.transpose() + cp_ca * cp_ca.transpose());

      // add the 3 simple cross inference matrix - components (we only have one 
      // variance here)
      A += 6 * var * var * Mat_3::Identity();

      // we need the determinant of matrix with columns a, b, c - we use the scalar triple product
      const FT det = dot_product(ab, c);

      // compute the b vector, this follows the formula directly - but we can factor
      // out the diagonal covariance matrices 
      const Eigen::Matrix<FT, 3, 1> res_b = det * sum_cp_col
        - var * (
          vector_to_col_vec(cross_product(a_minus_b, ab))
        + vector_to_col_vec(cross_product(b_minus_c, bc))
        + vector_to_col_vec(cross_product(c_minus_a, ca)))
        + 2 * vector_to_col_vec(sum_vectors(sum_vectors(a, b), c)) * var * var;

      const FT res_c = det * det 
        + var * (
          dot_product(ab, ab) 
        + dot_product(bc, bc) 
        + dot_product(ca, ca)) 
        + var * var * ( 
            2 * (
            dot_product(a, a)
          + dot_product(b, b)
          + dot_product(c, c))
        + 6 * var);

      Mat_4 ret = Mat_4::Zero();
      ret.block(0, 0, 3, 3) = A;
      ret.block(3, 0, 1, 3) = -res_b.transpose();
      ret.block(0, 3, 3, 1) = -res_b;
      ret(4, 4) = res_c;

      return ret;
    }
    
    const Get_cost& get_cost() const { return *this; }
    const Get_placement& get_placement() const { return *this; }
    
  private:
    Vertex_cost_map vcm_;
   
    // get the matrix representation of the cross product, i.e.
    // the matrix representing the linear form
    // x -> v cross x
    static Mat_3 skew_sym_mat_cross_product(const Vector_3& v) 
    {
      Mat_3 mat;

      mat << 0, -v.z(), v.y(),
          v.z(), 0, -v.x(),
          -v.y(), v.x(), 0;

      return mat;
    }

    static Eigen::Matrix<FT, 3, 1> vector_to_col_vec(const Vector_3& v) 
    {
      Eigen::Matrix<FT, 3, 1> col {v.x(), v.y(), v.z()};
      return col;
    }
    
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
    
  private:
    FT mean_variance;
    FT normal_variance;
};

} // namespace Surface_mesh_simplification
} // namespace CGAL

#endif //CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_GARLANDHECKBERT_PROBABILISTIC_TRI_POLICIES_H
