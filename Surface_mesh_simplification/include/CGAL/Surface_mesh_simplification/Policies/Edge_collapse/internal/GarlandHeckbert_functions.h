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

#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_GARLANDHECKBERT_FUNCTIONS_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_GARLANDHECKBERT_FUNCTIONS_H

#include <CGAL/license/Surface_mesh_simplification.h>

#include <CGAL/Surface_mesh_simplification/internal/Common.h>

#include <Eigen/Dense>
#include <iostream>

namespace CGAL {
namespace Surface_mesh_simplification {
namespace internal {

// taken from https://stackoverflow.com/questions/1148309/inverting-a-4x4-matrix/
template<typename Matrix>
bool invert_matrix_4(const Matrix& m, Matrix& im)
{
  double det;
  
  double A2323 = m(2, 2) * m(3, 3) - m(2, 3) * m(3, 2);
  double A1323 = m(2, 1) * m(3, 3) - m(2, 3) * m(3, 1);
  double A1223 = m(2, 1) * m(3, 2) - m(2, 2) * m(3, 1);
  double A0323 = m(2, 0) * m(3, 3) - m(2, 3) * m(3, 0);
  double A0223 = m(2, 0) * m(3, 2) - m(2, 2) * m(3, 0);
  double A0123 = m(2, 0) * m(3, 1) - m(2, 1) * m(3, 0);
  double A2313 = m(1, 2) * m(3, 3) - m(1, 3) * m(3, 2);
  double A1313 = m(1, 1) * m(3, 3) - m(1, 3) * m(3, 1);
  double A1213 = m(1, 1) * m(3, 2) - m(1, 2) * m(3, 1);
  double A2312 = m(1, 2) * m(2, 3) - m(1, 3) * m(2, 2);
  double A1312 = m(1, 1) * m(2, 3) - m(1, 3) * m(2, 1);
  double A1212 = m(1, 1) * m(2, 2) - m(1, 2) * m(2, 1);
  double A0313 = m(1, 0) * m(3, 3) - m(1, 3) * m(3, 0);
  double A0213 = m(1, 0) * m(3, 2) - m(1, 2) * m(3, 0);
  double A0312 = m(1, 0) * m(2, 3) - m(1, 3) * m(2, 0);
  double A0212 = m(1, 0) * m(2, 2) - m(1, 2) * m(2, 0);
  double A0113 = m(1, 0) * m(3, 1) - m(1, 1) * m(3, 0);
  double A0112 = m(1, 0) * m(2, 1) - m(1, 1) * m(2, 0);

  det = m(0, 0) * ( m(1, 1) * A2323 - m(1, 2) * A1323 + m(1, 3) * A1223 )
      - m(0, 1) * ( m(1, 0) * A2323 - m(1, 2) * A0323 + m(1, 3) * A0223 )
      + m(0, 2) * ( m(1, 0) * A1323 - m(1, 1) * A0323 + m(1, 3) * A0123 )
      - m(0, 3) * ( m(1, 0) * A1223 - m(1, 1) * A0223 + m(1, 2) * A0123 );
  det = 1 / det;

  if (det == 0.0)
  {
    return false;
  }

  // we never actually use values other than those in the third column,
  // so might as well not calculate them
  //im(0, 0) = det *   ( m(1, 1) * A2323 - m(1, 2) * A1323 + m(1, 3) * A1223 );
  //im(0, 1) = det * - ( m(0, 1) * A2323 - m(0, 2) * A1323 + m(0, 3) * A1223 );
  //im(0, 2) = det *   ( m(0, 1) * A2313 - m(0, 2) * A1313 + m(0, 3) * A1213 );
  im(0, 3) = det * - ( m(0, 1) * A2312 - m(0, 2) * A1312 + m(0, 3) * A1212 );
  //im(1, 0) = det * - ( m(1, 0) * A2323 - m(1, 2) * A0323 + m(1, 3) * A0223 );
  //im(1, 1) = det *   ( m(0, 0) * A2323 - m(0, 2) * A0323 + m(0, 3) * A0223 );
  //im(1, 2) = det * - ( m(0, 0) * A2313 - m(0, 2) * A0313 + m(0, 3) * A0213 );
  im(1, 3) = det *   ( m(0, 0) * A2312 - m(0, 2) * A0312 + m(0, 3) * A0212 );
  //im(2, 0) = det *   ( m(1, 0) * A1323 - m(1, 1) * A0323 + m(1, 3) * A0123 );
  //im(2, 1) = det * - ( m(0, 0) * A1323 - m(0, 1) * A0323 + m(0, 3) * A0123 );
  //im(2, 2) = det *   ( m(0, 0) * A1313 - m(0, 1) * A0313 + m(0, 3) * A0113 );
  im(2, 3) = det * - ( m(0, 0) * A1312 - m(0, 1) * A0312 + m(0, 3) * A0112 );
  //im(3, 0) = det * - ( m(1, 0) * A1223 - m(1, 1) * A0223 + m(1, 2) * A0123 );
  //im(3, 1) = det *   ( m(0, 0) * A1223 - m(0, 1) * A0223 + m(0, 2) * A0123 );
  //im(3, 2) = det * - ( m(0, 0) * A1213 - m(0, 1) * A0213 + m(0, 2) * A0113 );
  im(3, 3) = det *   ( m(0, 0) * A1212 - m(0, 1) * A0212 + m(0, 2) * A0112 );

  return true;
}

template<typename GeomTraits>
Eigen::Matrix<typename GeomTraits::FT, 3, 1> vector_to_col_vec(
    const typename GeomTraits::Vector_3& v) 
{
  Eigen::Matrix<typename GeomTraits::FT, 3, 1> col {v.x(), v.y(), v.z()};
  return col;
}

// convenience alias declarations to make function return and argument types more readable
template<typename GeomTraits>
using Mat_4 = Eigen::Matrix<typename GeomTraits::FT, 4, 4, Eigen::DontAlign>;

template<typename GeomTraits>
using Col_4 = Eigen::Matrix<typename GeomTraits::FT, 4, 1>;

template<typename VertexPointMap, typename TriangleMesh, typename GeomTraits>
typename GeomTraits::Vector_3 construct_unit_normal_from_face(
    const VertexPointMap& point_map,
    const TriangleMesh& tmesh, 
    typename boost::graph_traits<TriangleMesh>::face_descriptor f, 
    const GeomTraits& gt) 
{      
  // initialize all necessary kernel functions
  auto unit_normal = gt.construct_unit_normal_3_object();

  // reference and descriptor types
  typedef typename boost::property_traits<VertexPointMap>::reference point_reference;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;

  const halfedge_descriptor h = halfedge(f, tmesh);

  // get the three points of the face and calculate their unit normal
  const point_reference p = get(point_map, source(h, tmesh));
  const point_reference q = get(point_map, target(h, tmesh));
  const point_reference r = get(point_map, target(next(h, tmesh), tmesh));

  return unit_normal(p, q, r);
}

template<typename VertexPointMap, typename TriangleMesh, typename GeomTraits> 
typename GeomTraits::Vector_3 construct_edge_normal(
    const VertexPointMap& point_map,
    const TriangleMesh& tmesh, 
    typename boost::graph_traits<TriangleMesh>::halfedge_descriptor he, 
    const GeomTraits& gt)
{
  typedef typename GeomTraits::Vector_3 Vector_3;
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;

  const Vector_3 face_normal = construct_unit_normal_from_face(point_map, tmesh, face(he, tmesh), gt);

  const vertex_descriptor vs = source(he, tmesh);
  const vertex_descriptor vt = target(he, tmesh);

  const Vector_3 edge_vector = Vector_3(get(point_map, vs), get(point_map, vt));
  const Vector_3 discontinuity_normal = cross_product(edge_vector, face_normal);

  // normalize
  const Vector_3 normal = discontinuity_normal 
    / sqrt(discontinuity_normal.squared_length());

  return normal;
}

template<typename GeomTraits>
Mat_4<GeomTraits> construct_classic_plane_quadric_from_normal(
    const typename GeomTraits::Vector_3& normal,
    const typename GeomTraits::Point_3& point,
    const GeomTraits& gt)
{
  typedef typename GeomTraits::FT FT;

  auto dot_product = gt.compute_scalar_product_3_object();
  auto construct_vector = gt.construct_vector_3_object();

  // negative dot product between the normal and the position vector
  const FT d = -dot_product(normal, construct_vector(ORIGIN, point));

  // row vector given by d appended to the normal
  const Eigen::Matrix<FT, 1, 4> row (normal.x(), normal.y(), normal.z(), d);

  // outer product
  return row.transpose() * row;
}

template<typename VertexPointMap, typename TriangleMesh, typename GeomTraits>
Mat_4<GeomTraits> construct_classic_plane_quadric_from_face(
    const VertexPointMap& point_map,
    const TriangleMesh& mesh,
    typename boost::graph_traits<TriangleMesh>::face_descriptor f,
    const GeomTraits& gt) 
{
  const typename GeomTraits::Vector_3 normal 
    = construct_unit_normal_from_face(point_map, mesh, f, gt);

  // get any point of the face 
  const auto p = get(point_map, source(halfedge(f, mesh), mesh));

  return construct_classic_plane_quadric_from_normal(normal, p, gt);
}

template<typename VertexPointMap, typename TriangleMesh, typename GeomTraits>
Mat_4<GeomTraits> construct_classic_plane_quadric_from_edge(
    const VertexPointMap& point_map,
    const TriangleMesh& mesh,
    typename boost::graph_traits<TriangleMesh>::halfedge_descriptor he,
    const GeomTraits& gt) 
{

  typedef typename GeomTraits::Vector_3 Vector_3;
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;

  Vector_3 normal = construct_edge_normal(point_map, mesh, he, gt);
  
  // use this normal to construct the quadric analogously to constructing quadric
  // from the normal of the face
  return construct_classic_plane_quadric_from_normal(normal, get(point_map, source(he, mesh)), gt);
}

template <typename GeomTraits>
Mat_4<GeomTraits> construct_prob_plane_quadric_from_normal(
    const typename GeomTraits::Vector_3& mean_normal,
    const typename GeomTraits::Point_3& point, 
    const GeomTraits& gt,
    typename GeomTraits::FT face_nv,
    typename GeomTraits::FT face_mv)
{
  auto squared_length = gt.compute_squared_length_3_object();
  auto dot_product = gt.compute_scalar_product_3_object();
  auto construct_vec_3 = gt.construct_vector_3_object();

  typedef typename GeomTraits::FT FT;
  typedef typename GeomTraits::Vector_3 Vector_3;
  typedef Eigen::Matrix<FT, 4, 4> Mat_4;
  
  const Vector_3 mean_vec = construct_vec_3(ORIGIN, point);
  const FT dot_mnmv = dot_product(mean_normal, mean_vec);

  // Eigen column vector of length 3
  const Eigen::Matrix<FT, 3, 1> mean_n_col {mean_normal.x(), mean_normal.y(), mean_normal.z()};

  // start by setting values along the diagonal
  Mat_4 mat = face_nv * Mat_4::Identity();

  // add outer product of the mean normal with itself
  // to the upper left 3x3 block
  mat.block(0, 0, 3, 3) += mean_n_col * mean_n_col.transpose();

  // set the first 3 values of the last row and the first
  // 3 values of the last column
  // the negative sign comes from the fact that in the paper,
  // the b column and row appear with a negative sign
  const auto b1 = -(dot_mnmv * mean_normal + face_nv * mean_vec);

  const Eigen::Matrix<FT, 3, 1> b {b1.x(), b1.y(), b1.z()};

  mat.col(3).head(3) = b;
  mat.row(3).head(3) = b.transpose();

  // set the value in the bottom right corner, we get this by considering
  // that we only have single variances given instead of covariance matrices
  mat(3, 3) = CGAL::square(dot_mnmv)
    + face_nv * squared_length(mean_vec)
    + face_mv * squared_length(mean_normal)
    + 3 * face_nv * face_mv;

  return mat;
}

template<typename VertexPointMap, typename TriangleMesh, typename GeomTraits>
std::array<typename GeomTraits::Vector_3, 3> vectors_from_face(
    const VertexPointMap& point_map,
    const TriangleMesh& tmesh, 
    typename boost::graph_traits<TriangleMesh>::face_descriptor f, 
    const GeomTraits& gt) 
{
  auto construct_vector = gt.construct_vector_3_object();
  
  typedef typename boost::property_traits<VertexPointMap>::reference Point_reference;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;
  
  std::array<typename GeomTraits::Vector_3, 3> arr { };
  
  const halfedge_descriptor h = halfedge(f, tmesh);

  // get all points and turn them into location vectors so we can use cross product on them
  const Point_reference p = get(point_map, source(h, tmesh));
  const Point_reference q = get(point_map, target(h, tmesh));
  const Point_reference r = get(point_map, target(next(h, tmesh), tmesh));

  arr[0] = construct_vector(ORIGIN, p);
  arr[1] = construct_vector(ORIGIN, q);
  arr[2] = construct_vector(ORIGIN, r);

  return arr; 
}

template<typename VertexPointMap, typename TriangleMesh, typename GeomTraits> 
Mat_4<GeomTraits> construct_classic_triangle_quadric_from_face(
    const VertexPointMap& point_map,
    const TriangleMesh& tmesh, 
    typename boost::graph_traits<TriangleMesh>::face_descriptor f, 
    const GeomTraits& gt)
{
  auto cross_product = gt.construct_cross_product_vector_3_object();
  auto sum_vectors = gt.construct_sum_of_vectors_3_object();
  auto dot_product = gt.compute_scalar_product_3_object();

  typedef typename GeomTraits::FT FT;

  auto vectors = vectors_from_face(point_map, tmesh, f, gt);
  
  Vector_3 a = vectors[0];
  Vector_3 b = vectors[1];
  Vector_3 c = vectors[2];

  const Vector_3 ab = cross_product(a, b);
  const Vector_3 bc = cross_product(b, c);
  const Vector_3 ca = cross_product(c, a);

  const Vector_3 sum_of_cross_product = sum_vectors(sum_vectors(ab, bc), ca);
  const FT scalar_triple_product = dot_product(ab, c);

  Eigen::Matrix<FT, 1, 4> row;

  row << sum_of_cross_product.x(), sum_of_cross_product.y(), 
      sum_of_cross_product.z(), -scalar_triple_product;

  // calculate the outer product of row^t*row
  return row.transpose() * row;
}

template<typename GeomTraits>
Eigen::Matrix<typename GeomTraits::FT, 3, 3> skew_sym_mat_cross_product(
    const typename GeomTraits::Vector_3& v) 
{
  Eigen::Matrix<typename GeomTraits::FT, 3, 3> mat;
  
  mat << 0, -v.z(), v.y(),
      v.z(), 0, -v.x(),
      -v.y(), v.x(), 0;

  return mat;
}

template<typename VertexPointMap, typename TriangleMesh, typename GeomTraits>
Mat_4<GeomTraits> construct_prob_triangle_quadric_from_face(
      const VertexPointMap& point_map,
      const TriangleMesh& tmesh, 
      typename boost::graph_traits<TriangleMesh>::face_descriptor f, 
      typename GeomTraits::FT var,
      const GeomTraits& gt)
{
  auto construct_vector = gt.construct_vector_3_object();
  auto cross_product = gt.construct_cross_product_vector_3_object();
  auto sum_vectors = gt.construct_sum_of_vectors_3_object();
  auto dot_product = gt.compute_scalar_product_3_object();
  
  // array containing the position vectors corresponding to
  // the vertices of the given face
  auto vectors = vectors_from_face(point_map, tmesh, f, gt);
  
  typedef typename GeomTraits::FT FT;
  typedef typename GeomTraits::Vector_3 Vector_3;
  typedef Eigen::Matrix<FT, 3, 3> Mat_3;
  typedef Mat_4<GeomTraits> Mat_4;
  
  Vector_3 a = vectors[0];
  Vector_3 b = vectors[1];
  Vector_3 c = vectors[2];

  // calculate certain vectors used later
  const Vector_3 ab = cross_product(a, b);
  const Vector_3 bc = cross_product(b, c);
  const Vector_3 ca = cross_product(c, a);

  // subtracting vectors using GeomTraits
  const Vector_3 a_minus_b = sum_vectors(a, -b);
  const Vector_3 b_minus_c = sum_vectors(b, -c);
  const Vector_3 c_minus_a = sum_vectors(c, -a);

  const Mat_3 cp_ab = skew_sym_mat_cross_product<GeomTraits>(a_minus_b);
  const Mat_3 cp_bc = skew_sym_mat_cross_product<GeomTraits>(b_minus_c);
  const Mat_3 cp_ca = skew_sym_mat_cross_product<GeomTraits>(c_minus_a);

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
        vector_to_col_vec<GeomTraits>(cross_product(a_minus_b, ab))
        + vector_to_col_vec<GeomTraits>(cross_product(b_minus_c, bc))
        + vector_to_col_vec<GeomTraits>(cross_product(c_minus_a, ca)))
    + 2 * vector_to_col_vec<GeomTraits>(sum_vectors(sum_vectors(a, b), c)) * var * var;

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
  ret(3, 3) = res_c;

  return ret;
}

template<typename GeomTraits>
Col_4<GeomTraits> construct_optimal_point_invertible(const Mat_4<GeomTraits>& quadric)
{
  Mat_4<GeomTraits> x;
  x << quadric.block(0, 0, 3, 4), 0, 0, 0, 1;

  Col_4<GeomTraits> opt_pt;

  opt_pt = x.inverse().col(3); // == X.inverse() * (0 0 0 1)
  return opt_pt;
}

template <typename GeomTraits>
Col_4<GeomTraits> construct_optimal_point_singular(
    const Mat_4<GeomTraits>& quadric,
    const Col_4<GeomTraits>& p0, 
    const Col_4<GeomTraits>& p1) 
{
  typedef typename GeomTraits::FT FT;
  
  // in this case, the matrix mat may no be invertible,
  // so we save the result to check
  Mat_4<GeomTraits> mat;
  mat << quadric.block(0, 0, 3, 4), 0, 0, 0, 1;

  Mat_4<GeomTraits> inverse;
  bool invertible = invert_matrix_4(mat, inverse);
  
  if (invertible)
  {
    return inverse.col(3);
  }
  else 
  {
    Col_4<GeomTraits> opt_pt;

    const Col_4<GeomTraits> p1mp0 = p1 - p0;
    const FT a = (p1mp0.transpose() * quadric * p1mp0)(0, 0);
    const FT b = 2 * (p0.transpose() * quadric * p1mp0)(0, 0);

    if (a == 0)
    {
      if (b < 0)
        opt_pt = p1;
      else if (b == 0)
        opt_pt = 0.5 * (p0 + p1);
      else
        opt_pt = p0;
    }
    else
    {
      FT ext_t = -b / (2 * a);
      if(ext_t < 0 || ext_t > 1 || a < 0)
      {
        // one of endpoints
        FT p0_cost = (p0.transpose() * quadric * p0)(0, 0);
        FT p1_cost = (p1.transpose() * quadric * p1)(0, 0);

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

    return opt_pt;
  }
}

template<typename TriangleMesh, typename GeomTraits>
std::pair<typename GeomTraits::FT, typename GeomTraits::FT> estimate_variances(
    const TriangleMesh& mesh, const GeomTraits& gt, 
    typename GeomTraits::FT variance, 
    typename GeomTraits::FT p_factor)
{
  typedef typename TriangleMesh::Vertex_index vertex_descriptor;
  typedef typename TriangleMesh::Edge_index edge_descriptor;
  typedef typename GeomTraits::FT FT;
  typedef typename GeomTraits::Point_3 Point_3;

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
  const FT n2 = variance * average_edge_length;
  const FT p2 = p_factor * variance * average_edge_length;

  return std::make_pair(n2, p2);
}

} // namespace internal
} // namespace Surface_mesh_simplification
} // namespace CGAL

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_GARLANDHECKBERT_FUNCTIONS_H
