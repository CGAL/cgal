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
#include <CGAL/boost/graph/helpers.h>

#include <CGAL/Origin.h>

#include <iostream>

namespace CGAL {
namespace Surface_mesh_simplification {
namespace internal {

// taken from https://stackoverflow.com/questions/1148309/inverting-a-4x4-matrix/
template <typename Matrix>
bool invert_matrix_4(const Matrix& m, Matrix& im)
{
  typedef typename Matrix::value_type                                          FT;

  FT det;

  FT A2323 = m(2, 2) * m(3, 3) - m(2, 3) * m(3, 2);
  FT A1323 = m(2, 1) * m(3, 3) - m(2, 3) * m(3, 1);
  FT A1223 = m(2, 1) * m(3, 2) - m(2, 2) * m(3, 1);
  FT A0323 = m(2, 0) * m(3, 3) - m(2, 3) * m(3, 0);
  FT A0223 = m(2, 0) * m(3, 2) - m(2, 2) * m(3, 0);
  FT A0123 = m(2, 0) * m(3, 1) - m(2, 1) * m(3, 0);
  // FT A2313 = m(1, 2) * m(3, 3) - m(1, 3) * m(3, 2);
  // FT A1313 = m(1, 1) * m(3, 3) - m(1, 3) * m(3, 1);
  // FT A1213 = m(1, 1) * m(3, 2) - m(1, 2) * m(3, 1);
  FT A2312 = m(1, 2) * m(2, 3) - m(1, 3) * m(2, 2);
  FT A1312 = m(1, 1) * m(2, 3) - m(1, 3) * m(2, 1);
  FT A1212 = m(1, 1) * m(2, 2) - m(1, 2) * m(2, 1);
  // FT A0313 = m(1, 0) * m(3, 3) - m(1, 3) * m(3, 0);
  // FT A0213 = m(1, 0) * m(3, 2) - m(1, 2) * m(3, 0);
  FT A0312 = m(1, 0) * m(2, 3) - m(1, 3) * m(2, 0);
  FT A0212 = m(1, 0) * m(2, 2) - m(1, 2) * m(2, 0);
  // FT A0113 = m(1, 0) * m(3, 1) - m(1, 1) * m(3, 0);
  FT A0112 = m(1, 0) * m(2, 1) - m(1, 1) * m(2, 0);

  det = m(0, 0) * ( m(1, 1) * A2323 - m(1, 2) * A1323 + m(1, 3) * A1223 )
      - m(0, 1) * ( m(1, 0) * A2323 - m(1, 2) * A0323 + m(1, 3) * A0223 )
      + m(0, 2) * ( m(1, 0) * A1323 - m(1, 1) * A0323 + m(1, 3) * A0123 )
      - m(0, 3) * ( m(1, 0) * A1223 - m(1, 1) * A0223 + m(1, 2) * A0123 );

  if(is_zero(det))
    return false;

  det = 1 / det;

  // we never actually use values other than those in the third column,
  // so might as well not calculate them
  // im(0, 0) = det *   ( m(1, 1) * A2323 - m(1, 2) * A1323 + m(1, 3) * A1223 );
  // im(0, 1) = det * - ( m(0, 1) * A2323 - m(0, 2) * A1323 + m(0, 3) * A1223 );
  // im(0, 2) = det *   ( m(0, 1) * A2313 - m(0, 2) * A1313 + m(0, 3) * A1213 );
  im(0, 3) = det * - ( m(0, 1) * A2312 - m(0, 2) * A1312 + m(0, 3) * A1212 );
  // im(1, 0) = det * - ( m(1, 0) * A2323 - m(1, 2) * A0323 + m(1, 3) * A0223 );
  // im(1, 1) = det *   ( m(0, 0) * A2323 - m(0, 2) * A0323 + m(0, 3) * A0223 );
  // im(1, 2) = det * - ( m(0, 0) * A2313 - m(0, 2) * A0313 + m(0, 3) * A0213 );
  im(1, 3) = det *   ( m(0, 0) * A2312 - m(0, 2) * A0312 + m(0, 3) * A0212 );
  // im(2, 0) = det *   ( m(1, 0) * A1323 - m(1, 1) * A0323 + m(1, 3) * A0123 );
  // im(2, 1) = det * - ( m(0, 0) * A1323 - m(0, 1) * A0323 + m(0, 3) * A0123 );
  // im(2, 2) = det *   ( m(0, 0) * A1313 - m(0, 1) * A0313 + m(0, 3) * A0113 );
  im(2, 3) = det * - ( m(0, 0) * A1312 - m(0, 1) * A0312 + m(0, 3) * A0112 );
  // im(3, 0) = det * - ( m(1, 0) * A1223 - m(1, 1) * A0223 + m(1, 2) * A0123 );
  // im(3, 1) = det *   ( m(0, 0) * A1223 - m(0, 1) * A0223 + m(0, 2) * A0123 );
  // im(3, 2) = det * - ( m(0, 0) * A1213 - m(0, 1) * A0213 + m(0, 2) * A0113 );
  im(3, 3) = det *   ( m(0, 0) * A1212 - m(0, 1) * A0212 + m(0, 2) * A0112 );

  return true;
}

template <typename GeomTraits>
typename GarlandHeckbert_matrix_types<GeomTraits>::Col_3
vector_to_col_3(const typename GeomTraits::Vector_3& v)
{
  typename GarlandHeckbert_matrix_types<GeomTraits>::Col_3 col { v.x(), v.y(), v.z() };
  return col;
}

template <typename TriangleMesh, typename VertexPointMap, typename GeomTraits>
typename GeomTraits::Vector_3
construct_unit_normal_from_face(typename boost::graph_traits<TriangleMesh>::face_descriptor f,
                                const TriangleMesh& tmesh,
                                const VertexPointMap vpm,
                                const GeomTraits& gt)
{
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor      halfedge_descriptor;
  typedef typename boost::property_traits<VertexPointMap>::reference           Point_ref;

  typedef typename GeomTraits::FT                                              FT;
  typedef typename GeomTraits::Vector_3                                        Vector_3;

  auto cross_product = gt.construct_cross_product_vector_3_object();
  auto squared_length = gt.compute_squared_length_3_object();

  const halfedge_descriptor h = halfedge(f, tmesh);

  const Point_ref p = get(vpm, target(h, tmesh));
  const Point_ref q = get(vpm, target(next(h, tmesh), tmesh));
  const Point_ref r = get(vpm, source(h, tmesh));

  Vector_3 normal = cross_product(q - p, r - p);

  const FT norm = sqrt(squared_length(normal));
  if(!is_zero(norm))
    normal = normal / norm;

  return normal;
}

template <typename TriangleMesh, typename VertexPointMap, typename GeomTraits>
typename GeomTraits::Vector_3
construct_edge_normal(typename boost::graph_traits<TriangleMesh>::halfedge_descriptor h,
                      const TriangleMesh& tmesh,
                      const VertexPointMap vpm,
                      const GeomTraits& gt)
{
  typedef typename boost::property_traits<VertexPointMap>::reference           Point_ref;

  typedef typename GeomTraits::FT                                              FT;
  typedef typename GeomTraits::Vector_3                                        Vector_3;

  auto vector = gt.construct_vector_3_object();
  auto cross_product = gt.construct_cross_product_vector_3_object();
  auto squared_length = gt.compute_squared_length_3_object();

  const Point_ref p = get(vpm, target(h, tmesh));
  const Point_ref q = get(vpm, target(next(h, tmesh), tmesh));
  const Point_ref r = get(vpm, source(h, tmesh));
  const Vector_3 face_normal = cross_product(q - p, r - p);

  const Vector_3 edge_vector = vector(r, p);
  Vector_3 normal = cross_product(edge_vector, face_normal);

  const FT norm = sqrt(squared_length(normal));
  if(!is_zero(norm))
    normal = normal / norm;

  return normal;
}

template <typename VertexPointMap, typename TriangleMesh, typename GeomTraits>
std::array<typename GeomTraits::Vector_3, 3>
vectors_from_face(typename boost::graph_traits<TriangleMesh>::face_descriptor f,
                  const TriangleMesh& tmesh,
                  const VertexPointMap vpm,
                  const GeomTraits& gt)
{
  typedef typename boost::property_traits<VertexPointMap>::reference           Point_reference;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor      halfedge_descriptor;

  auto vector = gt.construct_vector_3_object();

  const halfedge_descriptor h = halfedge(f, tmesh);

  // get all points and turn them into location vectors so we can use cross product on them
  const Point_reference p = get(vpm, target(h, tmesh));
  const Point_reference q = get(vpm, target(next(h, tmesh), tmesh));
  const Point_reference r = get(vpm, source(h, tmesh));

  std::array<typename GeomTraits::Vector_3, 3> arr { vector(ORIGIN, p),
                                                     vector(ORIGIN, q),
                                                     vector(ORIGIN, r) };

  return arr;
}

template <typename GeomTraits>
typename GarlandHeckbert_matrix_types<GeomTraits>::Mat_3
skew_sym_mat_cross_product(const typename GeomTraits::Vector_3& v)
{
  typename GarlandHeckbert_matrix_types<GeomTraits>::Mat_3 mat;

  mat <<       0, - v.z(),   v.y(),
           v.z(),       0, - v.x(),
         - v.y(),   v.x(),       0;

  return mat;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// OPTIMAL POINT
////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename GeomTraits>
typename GarlandHeckbert_matrix_types<GeomTraits>::Col_4
construct_optimal_point_invertible(const typename GarlandHeckbert_matrix_types<GeomTraits>::Mat_4& quadric)
{
  typedef typename GarlandHeckbert_matrix_types<GeomTraits>::Mat_4             Mat_4;
  typedef typename GarlandHeckbert_matrix_types<GeomTraits>::Col_4             Col_4;

  Mat_4 x;
  x << quadric.block(0, 0, 3, 4), 0, 0, 0, 1;

  Col_4 opt_pt;
  opt_pt = x.inverse().col(3); // == X.inverse() * (0 0 0 1)

  return opt_pt;
}

template <typename GeomTraits>
typename GarlandHeckbert_matrix_types<GeomTraits>::Col_4
construct_optimal_point_singular(const typename GarlandHeckbert_matrix_types<GeomTraits>::Mat_4& quadric,
                                 const typename GarlandHeckbert_matrix_types<GeomTraits>::Col_4& p0,
                                 const typename GarlandHeckbert_matrix_types<GeomTraits>::Col_4& p1)
{
  typedef typename GeomTraits::FT                                              FT;
  typedef typename GarlandHeckbert_matrix_types<GeomTraits>::Mat_4             Mat_4;
  typedef typename GarlandHeckbert_matrix_types<GeomTraits>::Col_4             Col_4;

  // In this case, the matrix mat may no be invertible, so we save the result to check
  Mat_4 mat;
  mat << quadric.block(0, 0, 3, 4), 0, 0, 0, 1;

  Mat_4 inverse;
  bool invertible = invert_matrix_4(mat, inverse);

  if(invertible)
  {
    return inverse.col(3);
  }
  else
  {
    Col_4 opt_pt;

    const Col_4 p1mp0 = p1 - p0;
    const FT a = (p1mp0.transpose() * quadric * p1mp0)(0, 0);
    const FT b = 2 * (p0.transpose() * quadric * p1mp0)(0, 0);

    if(is_zero(a))
    {
      if(b < 0)
        opt_pt = p1;
      else if(is_zero(b))
        opt_pt = 0.5 * (p0 + p1);
      else
        opt_pt = p0;
    }
    else
    {
      FT ext_t = - b / (FT(2) * a);
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


////////////////////////////////////////////////////////////////////////////////////////////////////
/// CLASSIC PLANE
////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename GeomTraits>
typename GarlandHeckbert_matrix_types<GeomTraits>::Mat_4
construct_classic_plane_quadric_from_normal(const typename GeomTraits::Vector_3& normal,
                                            const typename GeomTraits::Point_3& point,
                                            const GeomTraits& gt)
{
  typedef typename GeomTraits::FT                                              FT;

  typedef typename GarlandHeckbert_matrix_types<GeomTraits>::Row_4             Row_4;

  auto dot_product = gt.compute_scalar_product_3_object();
  auto vector = gt.construct_vector_3_object();

  // negative dot product between the normal and the position vector
  const FT d = - dot_product(normal, vector(ORIGIN, point));

  // row vector given by d appended to the normal
  Row_4 row { normal.x(), normal.y(), normal.z(), d };

  // outer product
  return row.transpose() * row;
}

template <typename TriangleMesh, typename VertexPointMap, typename GeomTraits>
typename GarlandHeckbert_matrix_types<GeomTraits>::Mat_4
construct_classic_plane_quadric_from_edge(typename boost::graph_traits<TriangleMesh>::halfedge_descriptor he,
                                          const TriangleMesh& mesh,
                                          const VertexPointMap vpm,
                                          const GeomTraits& gt)
{
  typedef typename GeomTraits::Vector_3                                        Vector_3;

  const Vector_3 normal = construct_edge_normal(he, mesh, vpm, gt);

  // use this normal to construct the quadric analogously to constructing quadric
  // from the normal of the face
  return construct_classic_plane_quadric_from_normal(normal, get(vpm, target(he, mesh)), gt);
}

template <typename TriangleMesh, typename VertexPointMap, typename GeomTraits>
typename GarlandHeckbert_matrix_types<GeomTraits>::Mat_4
construct_classic_plane_quadric_from_face(typename boost::graph_traits<TriangleMesh>::face_descriptor f,
                                          const TriangleMesh& mesh,
                                          const VertexPointMap vpm,
                                          const GeomTraits& gt)
{
  auto normal = construct_unit_normal_from_face(f, mesh, vpm, gt);

  // get any point of the face
  const auto p = get(vpm, target(halfedge(f, mesh), mesh));

  return construct_classic_plane_quadric_from_normal(normal, p, gt);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// CLASSIC TRIANGLE
////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename TriangleMesh, typename VertexPointMap, typename GeomTraits>
typename GarlandHeckbert_matrix_types<GeomTraits>::Mat_4
construct_classic_triangle_quadric_from_face(typename boost::graph_traits<TriangleMesh>::face_descriptor f,
                                             const TriangleMesh& tmesh,
                                             const VertexPointMap vpm,
                                             const GeomTraits& gt)
{
  typedef typename GeomTraits::FT                                              FT;
  typedef typename GeomTraits::Vector_3                                        Vector_3;

  typedef typename GarlandHeckbert_matrix_types<GeomTraits>::Row_4             Row_4;

  auto cross_product = gt.construct_cross_product_vector_3_object();
  auto dot_product = gt.compute_scalar_product_3_object();
  auto sum_vectors = gt.construct_sum_of_vectors_3_object();

  std::array<Vector_3, 3> vectors = vectors_from_face(f, tmesh, vpm, gt);

  const Vector_3& a = vectors[0];
  const Vector_3& b = vectors[1];
  const Vector_3& c = vectors[2];

  const Vector_3 ab = cross_product(a, b);
  const Vector_3 bc = cross_product(b, c);
  const Vector_3 ca = cross_product(c, a);

  const Vector_3 sum_of_cross_products = sum_vectors(sum_vectors(ab, bc), ca);
  const FT scalar_triple_product = dot_product(ab, c);

  Row_4 row;
  row << sum_of_cross_products.x(),   sum_of_cross_products.y(),
         sum_of_cross_products.z(), - scalar_triple_product;

  return row.transpose() * row;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// PROB PLANE
////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename GeomTraits>
typename GarlandHeckbert_matrix_types<GeomTraits>::Mat_4
construct_prob_plane_quadric_from_normal(const typename GeomTraits::Vector_3& mean_normal,
                                         const typename GeomTraits::Point_3& point,
                                         const GeomTraits& gt,
                                         typename GeomTraits::FT face_nv,
                                         typename GeomTraits::FT face_mv)
{
  typedef typename GeomTraits::FT                                              FT;
  typedef typename GeomTraits::Vector_3                                        Vector_3;

  typedef typename GarlandHeckbert_matrix_types<GeomTraits>::Col_3             Col_3;
  typedef typename GarlandHeckbert_matrix_types<GeomTraits>::Mat_4             Mat_4;

  auto vector = gt.construct_vector_3_object();
  auto dot_product = gt.compute_scalar_product_3_object();
  auto squared_length = gt.compute_squared_length_3_object();

  const Vector_3 mean_vec = vector(ORIGIN, point);
  const Col_3 mean_n_col { mean_normal.x(), mean_normal.y(), mean_normal.z() };

  // start by setting values along the diagonal
  Mat_4 mat = face_nv * Mat_4::Identity();

  // add outer product of the mean normal with itself to the upper left 3x3 block
  mat.block(0, 0, 3, 3) += mean_n_col * mean_n_col.transpose();

  // set the first 3 values of the last row and the first 3 values of the last column
  const FT dot_mnmv = dot_product(mean_normal, mean_vec);
  const Vector_3 b1 = dot_mnmv * mean_normal + face_nv * mean_vec;

  const Col_3 b { b1.x(), b1.y(), b1.z() };

  mat.col(3).head(3) = - b;
  mat.row(3).head(3) = - b.transpose();

  // set the value in the bottom right corner, we get this by considering
  // that we only have single variances given instead of covariance matrices
  mat(3, 3) = square(dot_mnmv)
              + face_nv * squared_length(mean_vec)
              + face_mv * squared_length(mean_normal)
              + 3 * face_nv * face_mv; // tr(Sigma_n * Sigma_m)

  return mat;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// PROB TRIANGLE
////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename TriangleMesh, typename VertexPointMap, typename GeomTraits>
typename GarlandHeckbert_matrix_types<GeomTraits>::Mat_4
construct_prob_triangle_quadric_from_face(typename boost::graph_traits<TriangleMesh>::face_descriptor f,
                                          typename GeomTraits::FT var,
                                          const TriangleMesh& tmesh,
                                          const VertexPointMap vpm,
                                          const GeomTraits& gt)
{
  typedef typename GeomTraits::FT                                              FT;
  typedef typename GeomTraits::Vector_3                                        Vector_3;

  typedef typename GarlandHeckbert_matrix_types<GeomTraits>::Mat_3             Mat_3;
  typedef typename GarlandHeckbert_matrix_types<GeomTraits>::Col_3             Col_3;
  typedef typename GarlandHeckbert_matrix_types<GeomTraits>::Mat_4             Mat_4;

  auto cross_product = gt.construct_cross_product_vector_3_object();
  auto sum_vectors = gt.construct_sum_of_vectors_3_object();
  auto dot_product = gt.compute_scalar_product_3_object();
  auto squared_length = gt.compute_squared_length_3_object();

  // array containing the position vectors corresponding to
  // the vertices of the given face
  auto vectors = vectors_from_face(f, tmesh, vpm, gt);

  const Vector_3& a = vectors[0];
  const Vector_3& b = vectors[1];
  const Vector_3& c = vectors[2];

  const Vector_3 ab = cross_product(a, b);
  const Vector_3 bc = cross_product(b, c);
  const Vector_3 ca = cross_product(c, a);

  const Vector_3 a_minus_b = sum_vectors(a, -b);
  const Vector_3 b_minus_c = sum_vectors(b, -c);
  const Vector_3 c_minus_a = sum_vectors(c, -a);

  const Mat_3 cp_ab = skew_sym_mat_cross_product<GeomTraits>(a_minus_b);
  const Mat_3 cp_bc = skew_sym_mat_cross_product<GeomTraits>(b_minus_c);
  const Mat_3 cp_ca = skew_sym_mat_cross_product<GeomTraits>(c_minus_a);

  const Vector_3 sum_of_cross_product = sum_vectors(sum_vectors(ab, bc), ca);
  const Col_3 sum_cp_col { sum_of_cross_product.x(), sum_of_cross_product.y(), sum_of_cross_product.z() };

  Mat_3 A = sum_cp_col * sum_cp_col.transpose();
  A += var * (cp_ab * cp_ab.transpose() + cp_bc * cp_bc.transpose() + cp_ca * cp_ca.transpose());

  // Add the 3 simple cross inference matrix - components (we only have one variance here)
  A += 6 * square(var) * Mat_3::Identity();

  // we need the determinant of matrix with columns a, b, c - we use the scalar triple product
  const FT det = dot_product(ab, c);

  // Compute the b vector, this follows the formula directly - but we can factor
  // out the diagonal covariance matrices
  const Col_3 res_b = det * sum_cp_col - var * (vector_to_col_3<GeomTraits>(cross_product(a_minus_b, ab))
                                                + vector_to_col_3<GeomTraits>(cross_product(b_minus_c, bc))
                                                + vector_to_col_3<GeomTraits>(cross_product(c_minus_a, ca)))
                      + 2 * square(var) * vector_to_col_3<GeomTraits>(sum_vectors(sum_vectors(a, b), c));

  const FT ab2 = squared_length(ab);
  const FT bc2 = squared_length(bc);
  const FT ca2 = squared_length(ca);
  const FT a2 = squared_length(a);
  const FT b2 = squared_length(b);
  const FT c2 = squared_length(c);

  const FT res_c = square(det)
                   + var * (ab2 + bc2 + ca2)
                   + square(var) * (2 * (a2 + b2 + c2) + 6 * var);

  Mat_4 ret = Mat_4::Zero();
  ret.block(0, 0, 3, 3) = A;
  ret.block(3, 0, 1, 3) = - res_b.transpose();
  ret.block(0, 3, 3, 1) = - res_b;
  ret(3, 3) = res_c;

  return ret;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// PROB VARIANCE
////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename TriangleMesh, typename GeomTraits>
std::pair<typename GeomTraits::FT, typename GeomTraits::FT>
estimate_variances(const TriangleMesh& mesh,
                   const GeomTraits& gt,
                   const typename GeomTraits::FT variance,
                   const typename GeomTraits::FT p_factor)
{
  typedef typename TriangleMesh::Vertex_index                                  vertex_descriptor;
  typedef typename TriangleMesh::Edge_index                                    edge_descriptor;

  typedef typename GeomTraits::FT                                              FT;
  typedef typename GeomTraits::Point_3                                         Point_3;
  typedef typename GeomTraits::Vector_3                                        Vector_3;

  CGAL_precondition(!CGAL::is_empty(mesh));

  auto construct_vector = gt.construct_vector_3_object();
  auto squared_length = gt.compute_squared_length_3_object();

  FT average_edge_length = 0;

  std::size_t ne = 0;
  for(edge_descriptor e : edges(mesh))
  {
    vertex_descriptor v1 = mesh.vertex(e, 0);
    vertex_descriptor v2 = mesh.vertex(e, 1);

    const Point_3& p1 = mesh.point(v1); // @fixme Surface_mesh API
    const Point_3& p2 = mesh.point(v2);

    const Vector_3 vec = construct_vector(p1, p2);
    average_edge_length += sqrt(squared_length(vec));

    ++ne; // edges(mesh).size() can be costly, might as well increment now
  }

  average_edge_length = average_edge_length / ne;

  const FT n2 = variance * average_edge_length;
  const FT p2 = p_factor * variance * average_edge_length;

  return std::make_pair(n2, p2);
}

} // namespace internal
} // namespace Surface_mesh_simplification
} // namespace CGAL

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_GARLANDHECKBERT_FUNCTIONS_H
