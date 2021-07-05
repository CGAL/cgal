// Copyright (c) 2019  GeometryFactory (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé,
//                 Julian Komaromy

#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_GARLANDHECKBERT_PROBABILISTIC_FACES_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_GARLANDHECKBERT_PROBABILISTIC_FACES_H

#include <CGAL/license/Surface_mesh_simplification.h>

#include <CGAL/Surface_mesh_simplification/internal/Common.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/internal/GarlandHeckbert_common.h>
#include <Eigen/Dense>

namespace CGAL {
namespace Surface_mesh_simplification {
namespace internal {
  
template<typename TriangleMesh, typename GeomTraits>
class GarlandHeckbert_probabilistic_faces
{
  typedef typename GeomTraits::FT FT;
  typedef typename GeomTraits::Vector_3 Vector_3;
  typedef typename GeomTraits::Point_3 Point_3;
  typedef typename Eigen::Matrix<FT, 4, 4, Eigen::DontAlign> Mat_4;

  public:
  GarlandHeckbert_probabilistic_faces(FT var): GarlandHeckbert_probabilistic_faces(var, var) { }
  
  GarlandHeckbert_probabilistic_faces(FT mean_var, FT normal_var): 
    mean_variance(mean_var), normal_variance(normal_var) { } 
  
  template<typename VPM, typename TM>
  Mat_4 construct_quadric_from_face(
      const VPM& point_map,
      const TM& tmesh, 
      typename boost::graph_traits<TM>::face_descriptor f, 
      const GeomTraits& gt) const
  {
    const Vector_3 normal = common::construct_unit_normal_from_face<
      GeomTraits, VPM, TM>(point_map, tmesh, f, gt);
    const Point_3 p = get(point_map, source(halfedge(f, tmesh), tmesh));

    return construct_quadric_from_normal(point_map, tmesh, f, gt);
  }

  private:
  template<typename VPM, typename TM> 
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

  private:
  FT normal_variance;
  FT mean_variance;
};

} //namespace internal
} //namespace Surface_mesh_simplification
} //namespace CGAL

#endif //CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_GARLANDHECKBERT_PROBABILISTIC_FACES_H
