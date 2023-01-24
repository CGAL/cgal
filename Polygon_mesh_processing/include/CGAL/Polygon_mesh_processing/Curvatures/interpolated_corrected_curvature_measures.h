// Copyright (c) 2022 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Hossam Saeed
//

#ifndef CGAL_POLYGON_MESH_PROCESSING_INTERPOLATED_CORRECTED_CURVATURE_MEASURES_H
#define CGAL_POLYGON_MESH_PROCESSING_INTERPOLATED_CORRECTED_CURVATURE_MEASURES_H

#include <CGAL/license/Polygon_mesh_processing/interpolated_corrected_curvature_measures.h>

#include <CGAL/assertions.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/property_map.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <Eigen/Eigenvalues>

#include <numeric>
#include <queue>
#include <unordered_set>
#include <functional>

#define EXPANDING_RADIUS_EPSILON 1e-6

namespace CGAL {

namespace Polygon_mesh_processing {

/**
  * \ingroup PMP_corrected_curvatures_grp
  *
  * \brief a struct for storing principal curvatures and directions.
  *
  * @tparam GT is the geometric traits class.
*/
template<typename GT>
struct Principal_curvatures_and_directions {

  /// min curvature magnitude
  typename GT::FT min_curvature;

  /// max curvature magnitude
  typename GT::FT max_curvature;

  /// min curvature direction vector
  typename GT::Vector_3 min_direction;

  /// max curvature direction vector
  typename GT::Vector_3 max_direction;

  Principal_curvatures_and_directions() {
    min_curvature = 0;
    max_curvature = 0;
    min_direction = typename GT::Vector_3(0, 0, 0);
    max_direction = typename GT::Vector_3(0, 0, 0);
  }

  Principal_curvatures_and_directions(
    typename GT::FT min_curvature,
    typename GT::FT max_curvature,
    typename GT::Vector_3 min_direction,
    typename GT::Vector_3 max_direction
  ) {
    this->min_curvature = min_curvature;
    this->max_curvature = max_curvature;
    this->min_direction = min_direction;
    this->max_direction = max_direction;
  }
};

namespace internal {

template<typename PolygonMesh, typename GT>
typename GT::FT average_edge_length(const PolygonMesh& pmesh) {
  const std::size_t n = edges(pmesh).size();
  if (n == 0)
    return 0;

  typename GT::FT avg_edge_length = 0;
  for (auto e : edges(pmesh))
    avg_edge_length += edge_length(e, pmesh);

  avg_edge_length /= n;
  return avg_edge_length;
}

enum Curvature_measure_index {
  MU0_AREA_MEASURE, ///< corrected area density
  MU1_MEAN_CURVATURE_MEASURE, ///< corrected mean curvature density
  MU2_GAUSSIAN_CURVATURE_MEASURE ///< corrected gaussian curvature density
};

template<typename GT>
struct Vertex_measures {
  typename GT::FT area_measure = 0;
  typename GT::FT mean_curvature_measure = 0;
  typename GT::FT gaussian_curvature_measure = 0;
  std::array<typename GT::FT, 3 * 3> anisotropic_measure = { 0, 0, 0,
                                                             0, 0, 0,
                                                             0, 0, 0 };
};

template<typename GT>
typename GT::FT interpolated_corrected_area_measure_face(const std::vector<typename GT::Vector_3>& u,
  const std::vector<typename GT::Vector_3>& x)
{
  const std::size_t n = x.size();
  CGAL_precondition(u.size() == n);
  CGAL_precondition(n >= 3);

  typename GT::Construct_cross_product_vector_3 cross_product;

  // Triangle: use triangle formula
  if (n == 3)
  {
    const typename GT::Vector_3 um = (u[0] + u[1] + u[2]) / 3.0;
    return 0.5 * um * cross_product(x[1] - x[0], x[2] - x[0]);
  }
  // Quad: use bilinear interpolation formula
  else if (n == 4)
  {
    // for the formulas below, values of verices 2 & 3 are swapped (compared to paper) to correct order.
    // the indices in paper vs in here are: 00 = 0, 10 = 1, 11 = 2, 01 = 3
    return (1.0 / 36.0) * (
      (4 * u[0] + 2 * u[1] + 2 * u[3] + u[2]) * cross_product(x[1] - x[0], x[3] - x[0])
      + (2 * u[0] + 4 * u[1] + u[3] + 2 * u[2]) * cross_product(x[1] - x[0], x[2] - x[1])
      + (2 * u[0] + u[1] + 4 * u[3] + 2 * u[2]) * cross_product(x[2] - x[3], x[3] - x[0])
      + (u[0] + 2 * u[1] + 2 * u[3] + 4 * u[2]) * cross_product(x[2] - x[3], x[2] - x[1])
      );
  }
  // N-gon: split into n triangles by polygon center and use triangle formula for each
  else
  {
    typename GT::FT mu0 = 0;

    // getting center of points
    typename GT::Vector_3 xc =
      std::accumulate(x.begin(), x.end(), typename GT::Vector_3(0, 0, 0));
    xc /= n;

    // getting unit average normal of points
    typename GT::Vector_3 uc =
      std::accumulate(u.begin(), u.end(), typename GT::Vector_3(0, 0, 0));
    uc /= sqrt(uc * uc);

    // summing each triangle's measure after triangulation by barycenter split.
    for (std::size_t i = 0; i < n; i++)
    {
      mu0 += interpolated_corrected_area_measure_face<GT>(
        std::vector<typename GT::Vector_3> {u[i], u[(i + 1) % n], uc},
        std::vector<typename GT::Vector_3> {x[i], x[(i + 1) % n], xc}
      );
    }
    return mu0;
  }
}

template<typename GT>
typename GT::FT interpolated_corrected_mean_curvature_measure_face(const std::vector<typename GT::Vector_3>& u,
  const std::vector<typename GT::Vector_3>& x)
{
  const std::size_t n = x.size();
  CGAL_precondition(u.size() == n);
  CGAL_precondition(n >= 3);

  typename GT::Construct_cross_product_vector_3 cross_product;

  // Triangle: use triangle formula
  if (n == 3)
  {
    const typename GT::Vector_3 um = (u[0] + u[1] + u[2]) / 3.0;

    return 0.5 * um * (cross_product(u[2] - u[1], x[0])
      + cross_product(u[0] - u[2], x[1])
      + cross_product(u[1] - u[0], x[2]));
  }
  // Quad: use bilinear interpolation formula
  else if (n == 4)
  {
    // for the formulas below, values of verices 2 & 3 are swapped (compared to paper) to correct order.
    // the indices in paper vs in here are: 00 = 0, 10 = 1, 11 = 2, 01 = 3

    const typename GT::Vector_3 u02 = u[2] - u[0];
    const typename GT::Vector_3 u13 = u[3] - u[1];
    const typename GT::Vector_3 x0_cross = cross_product(u13, x[0]);
    const typename GT::Vector_3 x1_cross = -cross_product(u02, x[1]);
    const typename GT::Vector_3 x3_cross = cross_product(u02, x[3]);
    const typename GT::Vector_3 x2_cross = -cross_product(u13, x[2]);

    return (1.0 / 12.0) * (
      u[0] * (2 * x0_cross - cross_product((u[3] + u[2]), x[1]) + cross_product((u[1] + u[2]), x[3]) + x2_cross)
      + u[1] * (cross_product((u[3] + u[2]), x[0]) + 2 * x1_cross + x3_cross - cross_product((u[0] + u[3]), x[2]))
      + u[3] * (-cross_product((u[1] + u[2]), x[0]) + x1_cross + 2 * x3_cross + cross_product((u[0] + u[1]), x[2]))
      + u[2] * (x0_cross + cross_product((u[0] + u[3]), x[1]) - cross_product((u[0] + u[1]), x[3]) + 2 * x2_cross)
      );
  }
  // N-gon: split into n triangles by polygon center and use triangle formula for each
  else
  {
    typename GT::FT mu1 = 0;

    // getting center of points
    typename GT::Vector_3 xc =
      std::accumulate(x.begin(), x.end(), typename GT::Vector_3(0, 0, 0));
    xc /= n;

    // getting unit average normal of points
    typename GT::Vector_3 uc =
      std::accumulate(u.begin(), u.end(), typename GT::Vector_3(0, 0, 0));
    uc /= sqrt(uc * uc);

    // summing each triangle's measure after triangulation by barycenter split.
    for (std::size_t i = 0; i < n; i++)
    {
      mu1 += interpolated_corrected_mean_curvature_measure_face<GT>(
        std::vector<typename GT::Vector_3> {u[i], u[(i + 1) % n], uc},
        std::vector<typename GT::Vector_3> {x[i], x[(i + 1) % n], xc}
      );
    }
    return mu1;
  }
}

template<typename GT>
typename GT::FT interpolated_corrected_gaussian_curvature_measure_face(const std::vector<typename GT::Vector_3>& u,
  const std::vector<typename GT::Vector_3>& x = {})
{
  const std::size_t n = u.size();
  CGAL_precondition(n >= 3);

  typename GT::Construct_cross_product_vector_3 cross_product;

  // Triangle: use triangle formula
  if (n == 3)
  {
    return 0.5 * u[0] * cross_product(u[1], u[2]);
  }
  // Quad: use bilinear interpolation formula
  else if (n == 4)
  {
    // for the formulas below, values of verices 2 & 3 are swapped (compared to paper) to correct order.
    // the indices in paper vs in here are: 00 = 0, 10 = 1, 11 = 2, 01 = 3
    return (1.0 / 36.0) * (
      (4 * u[0] + 2 * u[1] + 2 * u[3] + u[2]) * cross_product(u[1] - u[0], u[3] - u[0])
      + (2 * u[0] + 4 * u[1] + u[3] + 2 * u[2]) * cross_product(u[1] - u[0], u[2] - u[1])
      + (2 * u[0] + u[1] + 4 * u[3] + 2 * u[2]) * cross_product(u[2] - u[3], u[3] - u[0])
      + (u[0] + 2 * u[1] + 2 * u[3] + 4 * u[2]) * cross_product(u[2] - u[3], u[2] - u[1])
      );
  }
  // N-gon: split into n triangles by polygon center and use triangle formula for each
  else
  {
    typename GT::FT mu2 = 0;

    // getting unit average normal of points
    typename GT::Vector_3 uc =
      std::accumulate(u.begin(), u.end(), typename GT::Vector_3(0, 0, 0));
    uc /= sqrt(uc * uc);

    // summing each triangle's measure after triangulation by barycenter split.
    for (std::size_t i = 0; i < n; i++)
    {
      mu2 += interpolated_corrected_gaussian_curvature_measure_face<GT>(
        std::vector<typename GT::Vector_3> {u[i], u[(i + 1) % n], uc}
      );
    }
    return mu2;
  }
}

template<typename GT>
std::array<typename GT::FT, 3 * 3> interpolated_corrected_anisotropic_measure_face(const std::vector<typename GT::Vector_3>& u,
  const std::vector<typename GT::Vector_3>& x)
{
  const std::size_t n = x.size();
  CGAL_precondition(u.size() == n);
  CGAL_precondition(n >= 3);

  typename GT::Construct_cross_product_vector_3 cross_product;
  std::array<typename GT::FT, 3 * 3> muXY{ 0 };

  // Triangle: use triangle formula
  if (n == 3)
  {
    const typename GT::Vector_3 u01 = u[1] - u[0];
    const typename GT::Vector_3 u02 = u[2] - u[0];
    const typename GT::Vector_3 x01 = x[1] - x[0];
    const typename GT::Vector_3 x02 = x[2] - x[0];
    const typename GT::Vector_3 um = (u[0] + u[1] + u[2]) / 3.0;

    for (std::size_t ix = 0; ix < 3; ix++)
    {
      typename GT::Vector_3 X;
      if (ix == 0)
        X = typename GT::Vector_3(1, 0, 0);
      if (ix == 1)
        X = typename GT::Vector_3(0, 1, 0);
      if (ix == 2)
        X = typename GT::Vector_3(0, 0, 1);

      for (std::size_t iy = 0; iy < 3; iy++)
        muXY[ix * 3 + iy] = 0.5 * um * (cross_product(u02[iy] * X, x01) - cross_product(u01[iy] * X, x02));
    }
  }
  // Quad: use bilinear interpolation formula
  else if (n == 4)
  {
    // for the formulas below, values of verices 2 & 3 are swapped (compared to paper) to correct order.
    // the indices in paper vs in here are: 00 = 0, 10 = 1, 11 = 2, 01 = 3
    for (std::size_t ix = 0; ix < 3; ix++)
    {
      typename GT::Vector_3 X;
      if (ix == 0)
        X = typename GT::Vector_3(1, 0, 0);
      if (ix == 1)
        X = typename GT::Vector_3(0, 1, 0);
      if (ix == 2)
        X = typename GT::Vector_3(0, 0, 1);

      const typename GT::Vector_3 u0xX = cross_product(u[0], X);
      const typename GT::Vector_3 u1xX = cross_product(u[1], X);
      const typename GT::Vector_3 u2xX = cross_product(u[2], X);
      const typename GT::Vector_3 u3xX = cross_product(u[3], X);

      for (std::size_t iy = 0; iy < 3; iy++)
        muXY[ix * 3 + iy] = (1.0 / 72.0) * (

          u[0][iy] * (u0xX * (-x[0] - 11 * x[1] + 13 * x[3] - x[2])
            + u1xX * (-5 * x[0] - 7 * x[1] + 11 * x[3] + x[2])
            + u3xX * (x[0] - 7 * x[1] + 11 * x[3] - 5 * x[2])
            + u2xX * (-x[0] - 5 * x[1] + 7 * x[3] - x[2])
            )
          + u[1][iy] * (u0xX * (13 * x[0] - x[1] - 7 * x[3] - 5 * x[2])
            + u1xX * (17 * x[0] - 5 * x[1] - 5 * x[3] - 7 * x[2])
            + u3xX * (5 * x[0] + x[1] + x[3] - 7 * x[2])
            + u2xX * (7 * x[0] - x[1] + 5 * x[3] - 11 * x[2])
            )
          + u[2][iy] * (u0xX * (-11 * x[0] + 5 * x[1] - x[3] + 7 * x[2])
            + u1xX * (-7 * x[0] + x[1] + x[3] + 5 * x[2])
            + u3xX * (-7 * x[0] - 5 * x[1] - 5 * x[3] + 17 * x[2])
            + u2xX * (-5 * x[0] - 7 * x[1] - x[3] + 13 * x[2])
            )
          + u[3][iy] * (u0xX * (-x[0] + 7 * x[1] - 5 * x[3] - x[2])
            + u1xX * (-5 * x[0] + 11 * x[1] - 7 * x[3] + x[2])
            + u3xX * (x[0] + 11 * x[1] - 7 * x[3] - 5 * x[2])
            + u2xX * (-x[0] + 13 * x[1] - 11 * x[3] - x[2])
            )

          );
    }
  }
  // N-gon: split into n triangles by polygon center and use triangle formula for each
  else
  {
    // getting center of points
    typename GT::Vector_3 xc =
      std::accumulate(x.begin(), x.end(), typename GT::Vector_3(0, 0, 0));
    xc /= n;

    // getting unit average normal of points
    typename GT::Vector_3 uc =
      std::accumulate(u.begin(), u.end(), typename GT::Vector_3(0, 0, 0));
    uc /= sqrt(uc * uc);

    // summing each triangle's measure after triangulation by barycenter split.
    for (std::size_t i = 0; i < n; i++)
    {
      std::array<typename GT::FT, 3 * 3> muXY_curr_triangle =
        interpolated_corrected_anisotropic_measure_face<GT>(
          std::vector<typename GT::Vector_3> {u[i], u[(i + 1) % n], uc},
          std::vector<typename GT::Vector_3> {x[i], x[(i + 1) % n], xc}
      );

      for (std::size_t ix = 0; ix < 3; ix++)
        for (std::size_t iy = 0; iy < 3; iy++)
          muXY[ix * 3 + iy] += muXY_curr_triangle[ix * 3 + iy];
    }
  }
  return muXY;
}

//template<typename GT>
//typename GT::FT triangle_in_ball_ratio(const typename GT::Vector_3 x1,
//                                       const typename GT::Vector_3 x2,
//                                       const typename GT::Vector_3 x3,
//                                       const typename GT::FT r,
//                                       const typename GT::Vector_3 c,
//                                       const std::size_t res = 3)
//{
//    const typename GT::FT R = r * r;
//    const typename GT::FT acc = 1.0 / res;
//    std::size_t samples_in = 0;
//    for (typename GT::FT alpha = acc / 3; alpha < 1; alpha += acc)
//        for (typename GT::FT beta = acc / 3; beta < 1 - alpha; beta += acc)
//        {
//            if ((alpha * x1 + beta * x2 + (1 - alpha - beta) * x3 - c).squared_length() < R)
//                samples_in++;
//        }
//    return samples_in / (typename GT::FT)(res * (res + 1) / 2);
//}

template<typename GT>
typename GT::FT face_in_ball_ratio(const std::vector<typename GT::Vector_3>& x,
  const typename GT::FT r,
  const typename GT::Vector_3 c)
{
  const std::size_t n = x.size();

  // getting center of points
  typename GT::Vector_3 xm =
    std::accumulate(x.begin(), x.end(), typename GT::Vector_3(0, 0, 0));
  xm /= n;

  typename GT::FT d_min = (xm - c).squared_length();
  typename GT::FT d_max = d_min;

  for (const typename GT::Vector_3 xi : x)
  {
    const typename GT::FT d_sq = (xi - c).squared_length();
    d_max = (std::max)(d_sq, d_max);
    d_min = (std::min)(d_sq, d_min);
  }

  if (d_max <= r * r) return 1.0;
  else if (r * r <= d_min) return 0.0;

  d_max = sqrt(d_max);
  d_min = sqrt(d_min);

  return (r - d_min) / (d_max - d_min);
}

template<typename GT>
Principal_curvatures_and_directions<GT> principal_curvatures_and_directions_from_anisotropic_measures(
  const std::array<typename GT::FT, 3 * 3> anisotropic_measure,
  const typename GT::FT v_mu0,
  const typename GT::Vector_3 u_GT
)
{
  Eigen::Matrix<typename GT::FT, 3, 3> v_muXY = Eigen::Matrix<typename GT::FT, 3, 3>::Zero();

  for (std::size_t ix = 0; ix < 3; ix++)
    for (std::size_t iy = 0; iy < 3; iy++)
      v_muXY(ix, iy) = anisotropic_measure[ix * 3 + iy];

  Eigen::Matrix<typename GT::FT, 3, 1> u(u_GT.x(), u_GT.y(), u_GT.z());
  const typename GT::FT K = 1000 * v_mu0;

  v_muXY = 0.5 * (v_muXY + v_muXY.transpose()) + K * u * u.transpose();

  Eigen::SelfAdjointEigenSolver<Eigen::Matrix <typename GT::FT, 3, 3>> eigensolver;

  eigensolver.computeDirect(v_muXY);

  if (eigensolver.info() != Eigen::Success)
    return Principal_curvatures_and_directions<GT>();

  const Eigen::Matrix<typename GT::FT, 3, 1> eig_vals = eigensolver.eigenvalues();
  const Eigen::Matrix<typename GT::FT, 3, 3> eig_vecs = eigensolver.eigenvectors();

  const typename GT::Vector_3 min_eig_vec(eig_vecs(0, 1), eig_vecs(1, 1), eig_vecs(2, 1));
  const typename GT::Vector_3 max_eig_vec(eig_vecs(0, 0), eig_vecs(1, 0), eig_vecs(2, 0));

  return Principal_curvatures_and_directions<GT>(
    (v_mu0 != 0.0) ? -eig_vals[1] / v_mu0 : 0.0,
    (v_mu0 != 0.0) ? -eig_vals[0] / v_mu0 : 0.0,
    min_eig_vec,
    max_eig_vec
    );
}

template<typename GT, typename PolygonMesh, typename VPM, typename VNM>
Vertex_measures<GT> interpolated_corrected_measures_one_vertex_no_radius(
  const PolygonMesh pmesh,
  const typename boost::graph_traits<PolygonMesh>::vertex_descriptor v,
  const bool is_mean_curvature_selected,
  const bool is_gaussian_curvature_selected,
  const bool is_principal_curvatures_and_directions_selected,
  const VPM vpm,
  const VNM vnm
)
{
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor Face_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor Vertex_descriptor;
  typedef typename GT::Point_3 Point_3;
  typedef typename GT::Vector_3 Vector_3;
  typedef typename GT::FT FT;

  std::queue<Face_descriptor> bfs_queue;
  std::unordered_set<Face_descriptor> bfs_visited;

  Vertex_measures<GT> vertex_measures;

  std::vector<Vector_3> x;
  std::vector<Vector_3> u;

  for (Face_descriptor f : faces_around_target(halfedge(v, pmesh), pmesh)) {
    if (f != boost::graph_traits<PolygonMesh>::null_face())
    {
      for (Vertex_descriptor vi : vertices_around_face(halfedge(f, pmesh), pmesh))
      {
        Point_3 pi = get(vpm, vi);
        Vector_3 ui = get(vnm, vi);
        x.push_back(Vector_3(pi.x(), pi.y(), pi.z()));
        u.push_back(ui);
      }

      vertex_measures.area_measure += interpolated_corrected_area_measure_face<GT>(u, x);

      if (is_mean_curvature_selected)
        vertex_measures.mean_curvature_measure += interpolated_corrected_mean_curvature_measure_face<GT>(u, x);

      if (is_gaussian_curvature_selected)
        vertex_measures.gaussian_curvature_measure += interpolated_corrected_gaussian_curvature_measure_face<GT>(u, x);

      if (is_principal_curvatures_and_directions_selected)
      {
        const std::array<FT, 3 * 3> face_anisotropic_measure = interpolated_corrected_anisotropic_measure_face<GT>(u, x);
        for (std::size_t i = 0; i < 3 * 3; i++)
          vertex_measures.anisotropic_measure[i] += face_anisotropic_measure[i];
      }
    }

    x.clear();
    u.clear();
  }
  return vertex_measures;
}


template<typename GT, typename PolygonMesh, typename VPM, typename VNM>
Vertex_measures<GT> interpolated_corrected_measures_one_vertex(
  const PolygonMesh pmesh,
  const typename boost::graph_traits<PolygonMesh>::vertex_descriptor v,
  const typename GT::FT radius,
  const bool is_mean_curvature_selected,
  const bool is_gaussian_curvature_selected,
  const bool is_principal_curvatures_and_directions_selected,
  const VPM vpm,
  const VNM vnm
)
{
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor Face_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor Vertex_descriptor;
  typedef typename GT::Point_3 Point_3;
  typedef typename GT::Vector_3 Vector_3;
  typedef typename GT::FT FT;

  std::queue<Face_descriptor> bfs_queue;
  std::unordered_set<Face_descriptor> bfs_visited;

  Vertex_measures<GT> vertex_measures;

  Point_3 vp = get(vpm, v);
  Vector_3 c = Vector_3(vp.x(), vp.y(), vp.z());

  for (Face_descriptor f : faces_around_target(halfedge(v, pmesh), pmesh)) {
    if (f != boost::graph_traits<PolygonMesh>::null_face())
    {
      bfs_queue.push(f);
      bfs_visited.insert(f);
    }
  }
  std::vector<Vector_3> x;
  std::vector<Vector_3> u;
  while (!bfs_queue.empty()) {
    Face_descriptor fi = bfs_queue.front();
    bfs_queue.pop();

    // looping over vertices in face to get point coordinates

    for (Vertex_descriptor vi : vertices_around_face(halfedge(fi, pmesh), pmesh))
    {
      Point_3 pi = get(vpm, vi);
      Vector_3 ui = get(vnm, vi);
      x.push_back(Vector_3(pi.x(), pi.y(), pi.z()));
      u.push_back(ui);
    }

    const FT f_ratio = face_in_ball_ratio<GT>(x, radius, c);

    if (f_ratio != 0.0)
    {
      vertex_measures.area_measure += f_ratio * interpolated_corrected_area_measure_face<GT>(u, x);

      if (is_mean_curvature_selected)
        vertex_measures.mean_curvature_measure += f_ratio * interpolated_corrected_mean_curvature_measure_face<GT>(u, x);

      if (is_gaussian_curvature_selected)
        vertex_measures.gaussian_curvature_measure += f_ratio * interpolated_corrected_gaussian_curvature_measure_face<GT>(u, x);

      if (is_principal_curvatures_and_directions_selected)
      {
        const std::array<FT, 3 * 3> face_anisotropic_measure = interpolated_corrected_anisotropic_measure_face<GT>(u, x);
        for (std::size_t i = 0; i < 3 * 3; i++)
          vertex_measures.anisotropic_measure[i] += f_ratio * face_anisotropic_measure[i];
      }

      for (Face_descriptor fj : faces_around_face(halfedge(fi, pmesh), pmesh))
      {
        if (bfs_visited.find(fj) == bfs_visited.end() && fj != boost::graph_traits<PolygonMesh>::null_face())
        {
          bfs_queue.push(fj);
          bfs_visited.insert(fj);
        }
      }
    }

    x.clear();
    u.clear();
  }
  return vertex_measures;
}

template<typename PolygonMesh,
  typename NamedParameters = parameters::Default_named_parameters>
  void interpolated_corrected_curvatures_one_vertex(
    const PolygonMesh pmesh,
    const typename boost::graph_traits<PolygonMesh>::vertex_descriptor v,
    const NamedParameters& np = parameters::default_values()
  )
{
  typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type GT;
  typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::const_type Vertex_position_map;

  typedef dynamic_vertex_property_t<typename GT::Vector_3> Vector_map_tag;
  typedef typename boost::property_map<PolygonMesh, Vector_map_tag>::const_type Default_vector_map;
  typedef typename internal_np::Lookup_named_param_def<internal_np::vertex_normal_map_t,
    NamedParameters,
    Default_vector_map>::type Vertex_normal_map;

  using parameters::choose_parameter;
  using parameters::get_parameter;
  using parameters::is_default_parameter;

  Vertex_position_map vpm = choose_parameter(get_parameter(np, CGAL::vertex_point),
    get_const_property_map(CGAL::vertex_point, pmesh));

  Vertex_normal_map vnm = choose_parameter(get_parameter(np, internal_np::vertex_normal_map),
    get(Vector_map_tag(), pmesh));

  if (is_default_parameter<NamedParameters, internal_np::vertex_normal_map_t>::value)
    compute_vertex_normals(pmesh, vnm, np);

  typename GT::FT radius = choose_parameter(get_parameter(np, internal_np::ball_radius), -1);
  if (radius == 0)
    radius = average_edge_length<PolygonMesh, GT>(pmesh) * EXPANDING_RADIUS_EPSILON;

  typename GT::FT* vertex_mean_curvature = choose_parameter(get_parameter(np, internal_np::vertex_mean_curvature), nullptr);
  typename GT::FT* vertex_gaussian_curvature = choose_parameter(get_parameter(np, internal_np::vertex_gaussian_curvature), nullptr);
  Principal_curvatures_and_directions<GT>* vertex_principal_curvatures_and_directions = choose_parameter(get_parameter(np, internal_np::vertex_principal_curvatures_and_directions), nullptr);

  const bool is_mean_curvature_selected = (vertex_mean_curvature != nullptr);
  const bool is_gaussian_curvature_selected = (vertex_gaussian_curvature != nullptr);
  const bool is_principal_curvatures_and_directions_selected = (vertex_principal_curvatures_and_directions != nullptr);

  Vertex_measures<GT> vertex_measures;

  if (radius < 0)
  {
    vertex_measures = interpolated_corrected_measures_one_vertex_no_radius<GT>(
      pmesh,
      v,
      is_mean_curvature_selected,
      is_gaussian_curvature_selected,
      is_principal_curvatures_and_directions_selected,
      vpm,
      vnm
      );
  }
  else
  {
    vertex_measures = interpolated_corrected_measures_one_vertex<GT>(
      pmesh,
      v,
      radius,
      is_mean_curvature_selected,
      is_gaussian_curvature_selected,
      is_principal_curvatures_and_directions_selected,
      vpm,
      vnm
      );
  }

  if (is_mean_curvature_selected) {
    *vertex_mean_curvature = vertex_measures.area_measure != 0 ?
      0.5 * vertex_measures.mean_curvature_measure / vertex_measures.area_measure : 0;
  }

  if (is_gaussian_curvature_selected) {
    *vertex_gaussian_curvature = vertex_measures.area_measure != 0 ?
      vertex_measures.gaussian_curvature_measure / vertex_measures.area_measure : 0;
  }

  if (is_principal_curvatures_and_directions_selected) {
    const typename GT::Vector_3  v_normal = get(vnm, v);
    const Principal_curvatures_and_directions<GT> principal_curvatures_and_directions = principal_curvatures_and_directions_from_anisotropic_measures<GT>(
      vertex_measures.anisotropic_measure,
      vertex_measures.area_measure,
      v_normal
      );
    *vertex_principal_curvatures_and_directions = principal_curvatures_and_directions;
  }
}


template<typename PolygonMesh, typename NamedParameters = parameters::Default_named_parameters>
class Interpolated_corrected_curvatures_computer
{
  typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type GT;

  typedef typename GT::FT FT;
  typedef typename GT::Point_3 Point_3;
  typedef typename GT::Vector_3 Vector_3;

  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor Halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::edge_descriptor Edge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor Face_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor Vertex_descriptor;

  typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::const_type Vertex_position_map;

  typedef dynamic_vertex_property_t<Vector_3> Vector_map_tag;
  typedef typename boost::property_map<PolygonMesh, Vector_map_tag>::const_type Default_vector_map;
  typedef typename internal_np::Lookup_named_param_def<internal_np::vertex_normal_map_t,
    NamedParameters,
    Default_vector_map>::type Vertex_normal_map;

  typedef Constant_property_map<Vertex_descriptor, FT> Default_scalar_map;
  typedef typename internal_np::Lookup_named_param_def<internal_np::vertex_mean_curvature_map_t,
    NamedParameters,
    Default_scalar_map>::type Vertex_mean_curvature_map;

  typedef typename internal_np::Lookup_named_param_def<internal_np::vertex_gaussian_curvature_map_t,
    NamedParameters,
    Default_scalar_map>::type Vertex_gaussian_curvature_map;

  typedef Constant_property_map<Vertex_descriptor, Principal_curvatures_and_directions<GT>> Default_principal_map;
  typedef typename internal_np::Lookup_named_param_def<internal_np::vertex_principal_curvatures_and_directions_map_t,
    NamedParameters,
    Default_principal_map>::type Vertex_principal_curvatures_and_directions_map;

  typedef typename boost::property_map<PolygonMesh,
    CGAL::dynamic_face_property_t<FT>>::const_type Face_scalar_measure_map;
  typedef typename boost::property_map<PolygonMesh,
    CGAL::dynamic_face_property_t<std::array<FT, 3 * 3>>>::const_type Face_anisotropic_measure_map;

private:
  const PolygonMesh& pmesh;
  Vertex_position_map vpm;
  Vertex_normal_map vnm;
  FT ball_radius;

  bool is_mean_curvature_selected;
  bool is_gaussian_curvature_selected;
  bool is_principal_curvatures_and_directions_selected;

  Vertex_mean_curvature_map mean_curvature_map;
  Vertex_gaussian_curvature_map gaussian_curvature_map;
  Vertex_principal_curvatures_and_directions_map principal_curvatures_and_directions_map;

  Face_scalar_measure_map mu0_map, mu1_map, mu2_map;
  Face_anisotropic_measure_map muXY_map;

  void set_property_maps() {
    mu0_map = get(CGAL::dynamic_face_property_t<FT>(), pmesh);
    mu1_map = get(CGAL::dynamic_face_property_t<FT>(), pmesh);
    mu2_map = get(CGAL::dynamic_face_property_t<FT>(), pmesh);
    muXY_map = get(CGAL::dynamic_face_property_t<std::array<FT, 3 * 3>>(), pmesh);

  }

  void set_named_params(const NamedParameters& np)
  {
    using parameters::choose_parameter;
    using parameters::get_parameter;
    using parameters::is_default_parameter;

    vpm = choose_parameter(get_parameter(np, CGAL::vertex_point),
      get_const_property_map(CGAL::vertex_point, pmesh));

    vnm = choose_parameter(get_parameter(np, internal_np::vertex_normal_map),
      get(Vector_map_tag(), pmesh));

    if (is_default_parameter<NamedParameters, internal_np::vertex_normal_map_t>::value)
      compute_vertex_normals(pmesh, vnm, np);

    const FT radius = choose_parameter(get_parameter(np, internal_np::ball_radius), -1);

    is_mean_curvature_selected = !is_default_parameter<NamedParameters, internal_np::vertex_mean_curvature_map_t>::value;
    is_gaussian_curvature_selected = !is_default_parameter<NamedParameters, internal_np::vertex_gaussian_curvature_map_t>::value;
    is_principal_curvatures_and_directions_selected = !is_default_parameter<NamedParameters, internal_np::vertex_principal_curvatures_and_directions_map_t>::value;

    mean_curvature_map = choose_parameter(get_parameter(np, internal_np::vertex_mean_curvature_map), Default_scalar_map());
    gaussian_curvature_map = choose_parameter(get_parameter(np, internal_np::vertex_gaussian_curvature_map), Default_scalar_map());
    principal_curvatures_and_directions_map = choose_parameter(get_parameter(np, internal_np::vertex_principal_curvatures_and_directions_map), Default_principal_map());

    set_ball_radius(radius);
  }

  void set_ball_radius(const FT radius) {
    if (radius == 0)
      ball_radius = average_edge_length<PolygonMesh, GT>(pmesh) * EXPANDING_RADIUS_EPSILON;
    else
      ball_radius = radius;
  }

public:

  Interpolated_corrected_curvatures_computer(const PolygonMesh& pmesh,
    const NamedParameters& np = parameters::default_values()
  ) :
    pmesh(pmesh)
  {
    set_named_params(np);

    if (is_mean_curvature_selected || is_gaussian_curvature_selected || is_principal_curvatures_and_directions_selected)
    {
      set_property_maps();

      compute_selected_curvatures();
    }
  }

private:

  void interpolated_corrected_selected_measures_all_faces()
  {
    std::vector<Vector_3> x;
    std::vector<Vector_3> u;

    for (Face_descriptor f : faces(pmesh))
    {
      for (Vertex_descriptor v : vertices_around_face(halfedge(f, pmesh), pmesh))
      {
        Point_3 p = get(vpm, v);
        x.push_back(Vector_3(p.x(), p.y(), p.z()));
        u.push_back(get(vnm, v));
      }
      put(mu0_map, f, interpolated_corrected_area_measure_face<GT>(u, x));

      if (is_mean_curvature_selected)
        put(mu1_map, f, interpolated_corrected_mean_curvature_measure_face<GT>(u, x));

      if (is_gaussian_curvature_selected)
        put(mu2_map, f, interpolated_corrected_gaussian_curvature_measure_face<GT>(u));

      if (is_principal_curvatures_and_directions_selected)
        put(muXY_map, f, interpolated_corrected_anisotropic_measure_face<GT>(u, x));

      x.clear();
      u.clear();
    }
  }

  Vertex_measures<GT> expand_interpolated_corrected_measure_vertex_no_radius(Vertex_descriptor v)
  {
    Vertex_measures<GT> vertex_measures;

    for (Face_descriptor f : faces_around_target(halfedge(v, pmesh), pmesh)) {
      if (f == boost::graph_traits<PolygonMesh>::null_face())
        continue;

      vertex_measures.area_measure += get(mu0_map, f);

      if (is_mean_curvature_selected)
        vertex_measures.mean_curvature_measure += get(mu1_map, f);

      if (is_gaussian_curvature_selected)
        vertex_measures.gaussian_curvature_measure += get(mu2_map, f);

      if (is_principal_curvatures_and_directions_selected)
      {
        const std::array<FT, 3 * 3> face_anisotropic_measure = get(muXY_map, f);
        for (std::size_t i = 0; i < 3 * 3; i++)
          vertex_measures.anisotropic_measure[i] += face_anisotropic_measure[i];
      }
    }

    return vertex_measures;
  }

  Vertex_measures<GT> expand_interpolated_corrected_measure_vertex(Vertex_descriptor v)
  {
    std::queue<Face_descriptor> bfs_queue;
    std::unordered_set<Face_descriptor> bfs_visited;

    Point_3 vp = get(vpm, v);
    Vector_3 c = Vector_3(vp.x(), vp.y(), vp.z());

    Vertex_measures<GT> vertex_measures;

    for (Face_descriptor f : faces_around_target(halfedge(v, pmesh), pmesh)) {
      if (f != boost::graph_traits<PolygonMesh>::null_face())
      {
        bfs_queue.push(f);
        bfs_visited.insert(f);
      }
    }
    while (!bfs_queue.empty()) {
      Face_descriptor fi = bfs_queue.front();
      bfs_queue.pop();

      // looping over vertices in face to get point coordinates
      std::vector<Vector_3> x;
      for (Vertex_descriptor vi : vertices_around_face(halfedge(fi, pmesh), pmesh))
      {
        Point_3 pi = get(vpm, vi);
        x.push_back(Vector_3(pi.x(), pi.y(), pi.z()));
      }

      const FT f_ratio = face_in_ball_ratio<GT>(x, ball_radius, c);

      if (f_ratio != 0.0)
      {
        vertex_measures.area_measure += f_ratio * get(mu0_map, fi);

        if (is_mean_curvature_selected)
          vertex_measures.mean_curvature_measure += f_ratio * get(mu1_map, fi);

        if (is_gaussian_curvature_selected)
          vertex_measures.gaussian_curvature_measure += f_ratio * get(mu2_map, fi);

        if (is_principal_curvatures_and_directions_selected)
        {
          const std::array<FT, 3 * 3> face_anisotropic_measure = get(muXY_map, fi);
          for (std::size_t i = 0; i < 3 * 3; i++)
            vertex_measures.anisotropic_measure[i] += f_ratio * face_anisotropic_measure[i];
        }

        for (Face_descriptor fj : faces_around_face(halfedge(fi, pmesh), pmesh))
        {
          if (bfs_visited.find(fj) == bfs_visited.end() && fj != boost::graph_traits<PolygonMesh>::null_face())
          {
            bfs_queue.push(fj);
            bfs_visited.insert(fj);
          }
        }
      }
    }
    return vertex_measures;
  }

  void compute_selected_curvatures() {
    interpolated_corrected_selected_measures_all_faces();

    for (Vertex_descriptor v : vertices(pmesh))
    {
      Vertex_measures<GT> vertex_measures = (ball_radius < 0) ?
        expand_interpolated_corrected_measure_vertex_no_radius(v) :
        expand_interpolated_corrected_measure_vertex(v);

      if (is_mean_curvature_selected) {
        vertex_measures.area_measure != 0 ?
          put(mean_curvature_map, v, 0.5 * vertex_measures.mean_curvature_measure / vertex_measures.area_measure) :
          put(mean_curvature_map, v, 0);
      }

      if (is_gaussian_curvature_selected) {
        vertex_measures.area_measure != 0 ?
          put(gaussian_curvature_map, v, vertex_measures.gaussian_curvature_measure / vertex_measures.area_measure) :
          put(gaussian_curvature_map, v, 0);
      }

      if (is_principal_curvatures_and_directions_selected) {
        const Vector_3  v_normal = get(vnm, v);
        const Principal_curvatures_and_directions<GT> principal_curvatures_and_directions = principal_curvatures_and_directions_from_anisotropic_measures<GT>(
          vertex_measures.anisotropic_measure,
          vertex_measures.area_measure,
          v_normal
          );
        put(principal_curvatures_and_directions_map, v, principal_curvatures_and_directions);
      }
    }
  }
};

} // namespace internal

/**
* \ingroup PMP_corrected_curvatures_grp
*
* Computes the interpolated corrected curvatures across the mesh, based on the provided property maps.
* By providing mean, gaussian and/or principal curvature property maps as named parameters, the user
* can choose which curvatures to compute.
*
* @tparam PolygonMesh a model of `FaceListGraph`.
* @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters".
*
* @param pmesh the polygon mesh.
* @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*
* \cgalNamedParamsBegin
*
*   \cgalParamNBegin{vertex_point_map}
*     \cgalParamDescription{a property map associating points to the vertices of `pmesh`}
*     \cgalParamType{a class model of `ReadablePropertyMap` with
*                    `boost::graph_traits<PolygonMesh>::%Vertex_descriptor`
*                    as key type and `%Point_3` as value type}
*     \cgalParamDefault{`boost::get(CGAL::vertex_point, pmesh)`}
*     \cgalParamExtra{If this parameter is omitted, an internal property map for
*                     `CGAL::vertex_point_t` must be available in `PolygonMesh`.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{vertex_normal_map}
*     \cgalParamDescription{a property map associating normal vectors to the vertices of `pmesh`}
*     \cgalParamType{a class model of `ReadablePropertyMap` with
*                    `boost::graph_traits<PolygonMesh>::%Vertex_descriptor`
*                    as key type and `%Vector_3` as value type}
*     \cgalParamDefault{`get(dynamic_vertex_property_t<GT::Vector_3>(), pmesh)`}
*     \cgalParamExtra{If this parameter is omitted, vertex normals will be
*                     computed using compute_vertex_normals()}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{vertex_mean_curvature_map}
*     \cgalParamDescription{a property map associating mean curvatures to the vertices of `pmesh`}
*     \cgalParamType{a class model of `WritablePropertyMap` with
*                    `boost::graph_traits<PolygonMesh>::%Vertex_descriptor`
*                    as key type and `%FT` as value type}
*     \cgalParamDefault{`get(dynamic_vertex_property_t<GT::FT>(), pmesh)`}
*     \cgalParamExtra{If this parameter is omitted, mean curvatures won't be computed}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{vertex_gaussian_curvature_map}
*     \cgalParamDescription{a property map associating mean curvatures to the vertices of `pmesh`}
*     \cgalParamType{a class model of `WritablePropertyMap` with
*                    `boost::graph_traits<PolygonMesh>::%Vertex_descriptor`
*                    as key type and `%FT` as value type}
*     \cgalParamDefault{`get(dynamic_vertex_property_t<GT::FT>(), pmesh)`}
*     \cgalParamExtra{If this parameter is omitted, gaussian curvatures won't be computed}
*   \cgalParamNEnd
*
*
*   \cgalParamNBegin{vertex_principal_curvatures_and_directions_map}
*     \cgalParamDescription{a property map associating mean curvatures to the vertices of `pmesh`}
*     \cgalParamType{a class model of `WritablePropertyMap` with
*                    `boost::graph_traits<PolygonMesh>::%Vertex_descriptor`
*                    as key type and `%Principal_curvatures_and_directions` as value type}
*     \cgalParamDefault{`get(dynamic_vertex_property_t<Principal_curvatures_and_directions<GT>>(), pmesh)`}
*     \cgalParamExtra{If this parameter is omitted, mean principal won't be computed}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{ball_radius}
*     \cgalParamDescription{a scalar value specifying the radius used for expanding curvature measures
*                           by summing measures of faces inside a ball of this radius centered at the
*                           vertex expanded from. The summed face measures are weighted by their
*                           inclusion ratio inside this ball.}
*     \cgalParamType{`GT::FT`}
*     \cgalParamDefault{`-1`}
*     \cgalParamExtra{If this parameter is omitted (`-1`), the epansion is then just a sum of
*                     measures on faces around the vertex}
*   \cgalParamNEnd
*
* \cgalNamedParamsEnd
*
* @see `interpolated_corrected_mean_curvature()`
* @see `interpolated_corrected_gaussian_curvature()`
* @see `interpolated_corrected_principal_curvatures_and_directions()`
*/
template<typename PolygonMesh,
  typename NamedParameters = parameters::Default_named_parameters>
  void interpolated_corrected_curvatures(const PolygonMesh& pmesh,
    const NamedParameters& np = parameters::default_values())
{
  internal::Interpolated_corrected_curvatures_computer<PolygonMesh, NamedParameters>(pmesh, np);
}

/**
* \ingroup PMP_corrected_curvatures_grp
*
* Computes the interpolated corrected mean curvature across the mesh
* and stores it in a vertex property map `vcm`.
*
* @tparam PolygonMesh a model of `FaceListGraph`.
* @tparam VertexCurvatureMap model of `WritablePropertyMap` with
* `boost::graph_traits<PolygonMesh>::%Vertex_descriptor` as key type and `GT::FT` as value type.
* @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters".
*
* @param pmesh the polygon mesh.
* @param vcm the vertex property map in which the computed mean curvatures are stored.
* @param np an optional sequence of \ref bgl_namedparameters "Named Parameters".
*
* \cgalNamedParamsBegin
*
*   \cgalParamNBegin{vertex_point_map}
*     \cgalParamDescription{a property map associating points to the vertices of `pmesh`}
*     \cgalParamType{a class model of `ReadablePropertyMap` with
*                    `boost::graph_traits<PolygonMesh>::%Vertex_descriptor`
*                    as key type and `%Point_3` as value type}
*     \cgalParamDefault{`boost::get(CGAL::vertex_point, pmesh)`}
*     \cgalParamExtra{If this parameter is omitted, an internal property map for
*                     `CGAL::vertex_point_t` must be available in `PolygonMesh`.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{vertex_normal_map}
*     \cgalParamDescription{a property map associating normal vectors to the vertices of `pmesh`}
*     \cgalParamType{a class model of `ReadablePropertyMap` with
*                    `boost::graph_traits<PolygonMesh>::%Vertex_descriptor`
*                    as key type and `%Vector_3` as value type}
*     \cgalParamDefault{`get(dynamic_vertex_property_t<GT::Vector_3>(), pmesh)`}
*     \cgalParamExtra{If this parameter is omitted, vertex normals will be
*                     computed using compute_vertex_normals()}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{ball_radius}
*     \cgalParamDescription{a scalar value specifying the radius used for expanding curvature measures
*                           by summing measures of faces inside a ball of this radius centered at the
*                           vertex expanded from. The summed face measures are weighted by their
*                           inclusion ratio inside this ball}
*     \cgalParamType{`GT::FT`}
*     \cgalParamDefault{`-1`}
*     \cgalParamExtra{If this parameter is omitted (`-1`), the expansion will just by a weightless sum of
*                     measures on faces around the vertex}
*   \cgalParamNEnd
*
* \cgalNamedParamsEnd
*
* @see `interpolated_corrected_gaussian_curvature()`
* @see `interpolated_corrected_principal_curvatures_and_directions()`
* @see `interpolated_corrected_curvatures()`
*/

template<typename PolygonMesh, typename VertexCurvatureMap,
  typename NamedParameters = parameters::Default_named_parameters>
  void interpolated_corrected_mean_curvature(const PolygonMesh& pmesh,
    VertexCurvatureMap& vcm,
    const NamedParameters& np = parameters::default_values())
{
  interpolated_corrected_curvatures(pmesh, np.vertex_mean_curvature_map(vcm));
}

/**
* \ingroup PMP_corrected_curvatures_grp
*
* Computes the interpolated corrected gaussian curvature across the mesh
* and stores it in a vertex property map `vcm`.
*
* @tparam PolygonMesh a model of `FaceListGraph`.
* @tparam VertexCurvatureMap model of `WritablePropertyMap` with
* `boost::graph_traits<PolygonMesh>::%Vertex_descriptor` as key type and `GT::FT` as value type.
* @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters".
*
* @param pmesh the polygon mesh.
* @param vcm the vertex property map in which the computed gaussian curvatures are stored.
* @param np an optional sequence of \ref bgl_namedparameters "Named Parameters".
*
* \cgalNamedParamsBegin
*
*   \cgalParamNBegin{vertex_point_map}
*     \cgalParamDescription{a property map associating points to the vertices of `pmesh`}
*     \cgalParamType{a class model of `ReadablePropertyMap` with
*                    `boost::graph_traits<PolygonMesh>::%Vertex_descriptor`
*                    as key type and `%Point_3` as value type}
*     \cgalParamDefault{`boost::get(CGAL::vertex_point, pmesh)`}
*     \cgalParamExtra{If this parameter is omitted, an internal property map for
*                     `CGAL::vertex_point_t` must be available in `PolygonMesh`.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{vertex_normal_map}
*     \cgalParamDescription{a property map associating normal vectors to the vertices of `pmesh`}
*     \cgalParamType{a class model of `ReadablePropertyMap` with
*                    `boost::graph_traits<PolygonMesh>::%Vertex_descriptor`
*                    as key type and `%Vector_3` as value type}
*     \cgalParamDefault{`get(dynamic_vertex_property_t<GT::Vector_3>(), pmesh)`}
*     \cgalParamExtra{If this parameter is omitted, vertex normals will be
*                     computed using compute_vertex_normals()}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{ball_radius}
*     \cgalParamDescription{a scalar value specifying the radius used for expanding curvature measures
*                           by summing measures of faces inside a ball of this radius centered at the
*                           vertex expanded from. The summed face measures are weighted by their
*                           inclusion ratio inside this ball}
*     \cgalParamType{`GT::FT`}
*     \cgalParamDefault{`-1`}
*     \cgalParamExtra{If this parameter is omitted (`-1`), the expansion will just by a weightless sum of
*                     measures on faces around the vertex}
*   \cgalParamNEnd
*
* \cgalNamedParamsEnd
*
* @see `interpolated_corrected_mean_curvature()`
* @see `interpolated_corrected_principal_curvatures_and_directions()`
* @see `interpolated_corrected_curvatures()`
*/
template<typename PolygonMesh, typename VertexCurvatureMap,
  typename NamedParameters = parameters::Default_named_parameters>
  void interpolated_corrected_gaussian_curvature(const PolygonMesh& pmesh,
    VertexCurvatureMap& vcm,
    const NamedParameters& np = parameters::default_values())
{
  interpolated_corrected_curvatures(pmesh, np.vertex_gaussian_curvature_map(vcm));
}

/**
* \ingroup PMP_corrected_curvatures_grp
*
* Computes the interpolated corrected principal curvatures across the mesh
* and stores it in a vertex property map `vcm`.
*
* @tparam PolygonMesh a model of `FaceListGraph`.
* @tparam VertexCurvatureMap model of `WritablePropertyMap` with
* `boost::graph_traits<PolygonMesh>::%Vertex_descriptor` as key type and
* `%Principal_curvatures_and_directions` as value type.
* @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters".
*
* @param pmesh the polygon mesh.
* @param vcm the vertex property map in which the computed principal curvatures are stored.
* @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*
* \cgalNamedParamsBegin
*
*   \cgalParamNBegin{vertex_point_map}
*     \cgalParamDescription{a property map associating points to the vertices of `pmesh`}
*     \cgalParamType{a class model of `ReadablePropertyMap` with
*                    `boost::graph_traits<PolygonMesh>::%Vertex_descriptor`
*                    as key type and `%Point_3` as value type}
*     \cgalParamDefault{`boost::get(CGAL::vertex_point, pmesh)`}
*     \cgalParamExtra{If this parameter is omitted, an internal property map for
*                     `CGAL::vertex_point_t` must be available in `PolygonMesh`.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{vertex_normal_map}
*     \cgalParamDescription{a property map associating normal vectors to the vertices of `pmesh`}
*     \cgalParamType{a class model of `ReadablePropertyMap` with
*                    `boost::graph_traits<PolygonMesh>::%Vertex_descriptor`
*                    as key type and `%Vector_3` as value type}
*     \cgalParamDefault{`get(dynamic_vertex_property_t<GT::Vector_3>(), pmesh)`}
*     \cgalParamExtra{If this parameter is omitted, vertex normals will be
*                     computed using compute_vertex_normals()}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{ball_radius}
*     \cgalParamDescription{a scalar value specifying the radius used for expanding curvature measures
*                           by summing measures of faces inside a ball of this radius centered at the
*                           vertex expanded from. The summed face measures are weighted by their
*                           inclusion ratio inside this ball}
*     \cgalParamType{`GT::FT`}
*     \cgalParamDefault{`-1`}
*     \cgalParamExtra{If this parameter is omitted (`-1`), the expansion will just by a weightless sum of
*                     measures on faces around the vertex}
*   \cgalParamNEnd
*
* \cgalNamedParamsEnd
*
* @see `interpolated_corrected_mean_curvature()`
* @see `interpolated_corrected_gaussian_curvature()`
* @see `interpolated_corrected_curvatures()`
*/
template<typename PolygonMesh, typename VertexCurvatureMap,
  typename NamedParameters = parameters::Default_named_parameters>
  void interpolated_corrected_principal_curvatures_and_directions(const PolygonMesh& pmesh,
    VertexCurvatureMap& vcm,
    const NamedParameters& np = parameters::default_values())
{
  interpolated_corrected_curvatures(pmesh, np.vertex_principal_curvatures_and_directions_map(vcm));
}


/**
  * \ingroup PMP_corrected_curvatures_grp
  * Computes the interpolated corrected curvatures at a certain vertex, based on the provided pointers.
  * By providing mean, gaussian and/or principal curvature pointers as named parameters, the user
  * can choose which curvatures to compute.
  * The pointers are used to store the computed curvatures.
  * The user is responsible for the memory management of the pointers.
  *
  * @tparam PolygonMesh a model of `FaceListGraph`
  * @tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
  *
  * @param pmesh the polygon mesh
  * @param v the vertex of `pmesh` to compute the curvatures at
  * @param np optional sequence of \ref pmp_namedparameters "Named Parameters" among the ones listed below
  *
  * \cgalNamedParamsBegin
  *   \cgalParamNBegin{vertex_point_map}
  *     \cgalParamDescription{a property map associating points to the vertices of `pmesh`}
  *     \cgalParamType{a class model of `ReadablePropertyMap` with
  *                    `boost::graph_traits<PolygonMesh>::%Vertex_descriptor`
  *                    as key type and `%Point_3` as value type}
  *     \cgalParamDefault{`boost::get(CGAL::vertex_point, pmesh)`}
  *     \cgalParamExtra{If this parameter is omitted, an internal property map for
  *                     `CGAL::vertex_point_t` must be available in `PolygonMesh`.}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{vertex_normal_map}
  *     \cgalParamDescription{a property map associating normal vectors to the vertices of `pmesh`}
  *     \cgalParamType{a class model of `ReadablePropertyMap` with
  *                    `boost::graph_traits<PolygonMesh>::%Vertex_descriptor`
  *                    as key type and `%Vector_3` as value type}
  *     \cgalParamDefault{`get(dynamic_vertex_property_t<GT::Vector_3>(), pmesh)`}
  *     \cgalParamExtra{If this parameter is omitted, vertex normals will be
  *                     computed using compute_vertex_normals()}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{vertex_mean_curvature}
  *     \cgalParamDescription{a pointer to a scalar value to store the mean curvature at the vertex `v`}
  *     \cgalParamType{`GT::FT*`}
  *     \cgalParamDefault{`nullptr`}
  *     \cgalParamExtra{If this parameter is omitted, mean curvature won't be computed}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{vertex_gaussian_curvature}
  *     \cgalParamDescription{a pointer to a scalar value to store the gaussian curvature at the vertex `v`}
  *     \cgalParamType{`GT::FT*`}
  *     \cgalParamDefault{`nullptr`}
  *     \cgalParamExtra{If this parameter is omitted, gaussian curvature won't be computed}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{vertex_principal_curvatures_and_directions}
  *     \cgalParamDescription{a pointer to a Principal_curvatures_and_directions object to store the principal curvatures and directions at the vertex `v`}
  *     \cgalParamType{`Principal_curvatures_and_directions<GT>*`}
  *     \cgalParamDefault{`nullptr`}
  *     \cgalParamExtra{If this parameter is omitted, principal curvatures and directions won't be computed}
  *  \cgalParamNEnd
  *
  *   \cgalParamNBegin{ball_radius}
  *     \cgalParamDescription{a scalar value specifying the radius used for expanding curvature measures
  *                           by summing measures of faces inside a ball of this radius centered at the
  *                           vertex expanded from. The summed face measures are weighted by their
  *                           inclusion ratio inside this ball.}
  *     \cgalParamType{`GT::FT`}
  *     \cgalParamDefault{`-1`}
  *     \cgalParamExtra{If this parameter is omitted (`-1`), the epansion is then just a sum of
  *                     measures on faces around the vertex}
  *   \cgalParamNEnd
  *
  * \cgalNamedParamsEnd
  *
  * @see `CGAL::Polygon_mesh_processing::interpolated_corrected_curvatures()`
  * @see `CGAL::Polygon_mesh_processing::interpolated_corrected_mean_curvature_at_vertex()`
  * @see `CGAL::Polygon_mesh_processing::interpolated_corrected_gaussian_curvature_at_vertex()`
  * @see `CGAL::Polygon_mesh_processing::interpolated_corrected_principal_curvatures_and_directions_at_vertex()`
*/
template<typename PolygonMesh,
  typename NamedParameters = parameters::Default_named_parameters>
  void interpolated_corrected_curvatures_at_vertex(const PolygonMesh& pmesh,
    typename boost::graph_traits<PolygonMesh>::vertex_descriptor v,
    const NamedParameters& np = parameters::default_values())
{
  internal::interpolated_corrected_curvatures_one_vertex(pmesh, v, np);
}

/**
  * \ingroup PMP_corrected_curvatures_grp
  * computes the interpolated corrected mean curvature at a vertex of a triangle mesh.
  *
  * @tparam GT a geometric traits class that provides the nested type `FT`,
  * @tparam PolygonMesh a model of `FaceListGraph`
  * @tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
  *
  * @param pmesh the polygon mesh
  * @param v the vertex of `pmesh` to compute the mean curvature at
  * @param np optional sequence of \ref pmp_namedparameters "Named Parameters" among the ones listed below
  *
  * \cgalNamedParamsBegin
  *   \cgalParamNBegin{vertex_point_map}
  *     \cgalParamDescription{a property map associating points to the vertices of `pmesh`}
  *     \cgalParamType{a class model of `ReadablePropertyMap` with
  *                    `boost::graph_traits<PolygonMesh>::%Vertex_descriptor`
  *                    as key type and `%Point_3` as value type}
  *     \cgalParamDefault{`boost::get(CGAL::vertex_point, pmesh)`}
  *     \cgalParamExtra{If this parameter is omitted, an internal property map for
  *                     `CGAL::vertex_point_t` must be available in `PolygonMesh`.}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{vertex_normal_map}
  *     \cgalParamDescription{a property map associating normal vectors to the vertices of `pmesh`}
  *     \cgalParamType{a class model of `ReadablePropertyMap` with
  *                    `boost::graph_traits<PolygonMesh>::%Vertex_descriptor`
  *                    as key type and `%Vector_3` as value type}
  *     \cgalParamDefault{`get(dynamic_vertex_property_t<GT::Vector_3>(), pmesh)`}
  *     \cgalParamExtra{If this parameter is omitted, vertex normals will be
  *                     computed using compute_vertex_normals()}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{ball_radius}
  *     \cgalParamDescription{a scalar value specifying the radius used for expanding curvature measures
  *                           by summing measures of faces inside a ball of this radius centered at the
  *                           vertex expanded from. The summed face measures are weighted by their
  *                           inclusion ratio inside this ball}
  *     \cgalParamType{`GT::FT`}
  *     \cgalParamDefault{`-1`}
  *     \cgalParamExtra{If this parameter is omitted (`-1`), the expansion will just by a weightless sum of
  *                     measures on faces around the vertex}
  *   \cgalParamNEnd
  *
  * \cgalNamedParamsEnd
  *
  * @return the interpolated corrected mean curvature at the vertex `v`
  *
  * @see `interpolated_corrected_mean_curvature()`
  * @see `interpolated_corrected_gaussian_curvature_at_vertex()`
  * @see `interpolated_corrected_principal_curvatures_and_directions_at_vertex()`
  * @see `interpolated_corrected_curvatures_at_vertex()`
*/

template<typename GT, typename PolygonMesh,
  typename NamedParameters = parameters::Default_named_parameters>
  typename GT::FT interpolated_corrected_mean_curvature_at_vertex(const PolygonMesh& pmesh,
    typename boost::graph_traits<PolygonMesh>::vertex_descriptor v,
    const NamedParameters& np = parameters::default_values())
{
  // use interpolated_corrected_curvatures_at_vertex to compute mean curvature
  typename GT::FT mean_curvature;
  interpolated_corrected_curvatures_at_vertex(pmesh, v, np.vertex_mean_curvature(&mean_curvature));
  return mean_curvature;
}

/**
  * \ingroup PMP_corrected_curvatures_grp
  * computes the interpolated corrected Gaussian curvature at a vertex of a triangle mesh.
  *
  * @tparam GT a geometric traits class that provides the nested type `FT`,
  * @tparam PolygonMesh a model of `FaceListGraph`
  * @tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
  *
  * @param pmesh the polygon mesh
  * @param v the vertex of `pmesh` to compute the Gaussian curvature at
  * @param np optional sequence of \ref pmp_namedparameters "Named Parameters" among the ones listed below
  *
  * \cgalNamedParamsBegin
  *   \cgalParamNBegin{vertex_point_map}
  *     \cgalParamDescription{a property map associating points to the vertices of `pmesh`}
  *     \cgalParamType{a class model of `ReadablePropertyMap` with
  *                    `boost::graph_traits<PolygonMesh>::%Vertex_descriptor`
  *                    as key type and `%Point_3` as value type}
  *     \cgalParamDefault{`boost::get(CGAL::vertex_point, pmesh)`}
  *     \cgalParamExtra{If this parameter is omitted, an internal property map for
  *                     `CGAL::vertex_point_t` must be available in `PolygonMesh`.}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{vertex_normal_map}
  *     \cgalParamDescription{a property map associating normal vectors to the vertices of `pmesh`}
  *     \cgalParamType{a class model of `ReadablePropertyMap` with
  *                    `boost::graph_traits<PolygonMesh>::%Vertex_descriptor`
  *                    as key type and `%Vector_3` as value type}
  *     \cgalParamDefault{`get(dynamic_vertex_property_t<GT::Vector_3>(), pmesh)`}
  *     \cgalParamExtra{If this parameter is omitted, vertex normals will be
  *                     computed using compute_vertex_normals()}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{ball_radius}
  *     \cgalParamDescription{a scalar value specifying the radius used for expanding curvature measures
  *                           by summing measures of faces inside a ball of this radius centered at the
  *                           vertex expanded from. The summed face measures are weighted by their
  *                           inclusion ratio inside this ball}
  *     \cgalParamType{`GT::FT`}
  *     \cgalParamDefault{`-1`}
  *     \cgalParamExtra{If this parameter is omitted (`-1`), the expansion will just by a weightless sum of
  *                     measures on faces around the vertex}
  *   \cgalParamNEnd
  *
  * \cgalNamedParamsEnd
  *
  * @return the interpolated corrected Gaussian curvature at the vertex `v`
  *
  * @see `interpolated_corrected_gaussian_curvature()`
  * @see `interpolated_corrected_mean_curvature_at_vertex()`
  * @see `interpolated_corrected_principal_curvatures_and_directions_at_vertex()`
  * @see `interpolated_corrected_curvatures_at_vertex()`
*/

template<typename GT, typename PolygonMesh,
  typename NamedParameters = parameters::Default_named_parameters>
  typename GT::FT interpolated_corrected_gaussian_curvature_at_vertex(const PolygonMesh& pmesh,
    typename boost::graph_traits<PolygonMesh>::vertex_descriptor v,
    const NamedParameters& np = parameters::default_values())
{
  // use interpolated_corrected_curvatures_at_vertex to compute gaussian curvature
  typename GT::FT gc;
  interpolated_corrected_curvatures_at_vertex(pmesh, v, np.vertex_gaussian_curvature(&gc));
  return gc;
}

/**
  * \ingroup PMP_corrected_curvatures_grp
  * computes the interpolated corrected principal curvatures and directions at a vertex of a triangle mesh.
  *
  * @tparam GT the geometric traits class,
  * @tparam PolygonMesh a model of `FaceListGraph`
  * @tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
  *
  * @param pmesh the polygon mesh
  * @param v the vertex of `pmesh` to compute the principal curvatures and directions at
  * @param np optional sequence of \ref pmp_namedparameters "Named Parameters" among the ones listed below
  *
  * \cgalNamedParamsBegin
  *   \cgalParamNBegin{vertex_point_map}
  *     \cgalParamDescription{a property map associating points to the vertices of `pmesh`}
  *     \cgalParamType{a class model of `ReadablePropertyMap` with
  *                    `boost::graph_traits<PolygonMesh>::%Vertex_descriptor`
  *                    as key type and `%Point_3` as value type}
  *     \cgalParamDefault{`boost::get(CGAL::vertex_point, pmesh)`}
  *     \cgalParamExtra{If this parameter is omitted, an internal property map for
  *                     `CGAL::vertex_point_t` must be available in `PolygonMesh`.}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{vertex_normal_map}
  *     \cgalParamDescription{a property map associating normal vectors to the vertices of `pmesh`}
  *     \cgalParamType{a class model of `ReadablePropertyMap` with
  *                    `boost::graph_traits<PolygonMesh>::%Vertex_descriptor`
  *                    as key type and `%Vector_3` as value type}
  *     \cgalParamDefault{`get(dynamic_vertex_property_t<GT::Vector_3>(), pmesh)`}
  *     \cgalParamExtra{If this parameter is omitted, vertex normals will be
  *                     computed using compute_vertex_normals()}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{ball_radius}
  *     \cgalParamDescription{a scalar value specifying the radius used for expanding curvature measures
  *                           by summing measures of faces inside a ball of this radius centered at the
  *                           vertex expanded from. The summed face measures are weighted by their
  *                           inclusion ratio inside this ball}
  *     \cgalParamType{`GT::FT`}
  *     \cgalParamDefault{`-1`}
  *     \cgalParamExtra{If this parameter is omitted (`-1`), the expansion will just by a weightless sum of
  *                     measures on faces around the vertex}
  *   \cgalParamNEnd
  * \cgalNamedParamsEnd
  *
  * @return the interpolated corrected principal curvatures and directions at the vertex `v`
  *
  * @see `interpolated_corrected_principal_curvatures_and_directions()`
  * @see `interpolated_corrected_mean_curvature_at_vertex()`
  * @see `interpolated_corrected_gaussian_curvature_at_vertex()`
  * @see `interpolated_corrected_curvatures_at_vertex()`
*/

template<typename GT, typename PolygonMesh,
  typename NamedParameters = parameters::Default_named_parameters>
  Principal_curvatures_and_directions<GT> interpolated_corrected_principal_curvatures_and_directions_at_vertex(const PolygonMesh& pmesh,
    typename boost::graph_traits<PolygonMesh>::vertex_descriptor v,
    const NamedParameters& np = parameters::default_values())
{
  // use interpolated_corrected_curvatures_at_vertex to compute principal curvatures
  Principal_curvatures_and_directions<GT> pcd;
  interpolated_corrected_curvatures_at_vertex(pmesh, v, np.vertex_principal_curvatures_and_directions(&pcd));
  return pcd;
}

} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_INTERPOLATED_CORRECTED_CURVATURE_MEASURES_H