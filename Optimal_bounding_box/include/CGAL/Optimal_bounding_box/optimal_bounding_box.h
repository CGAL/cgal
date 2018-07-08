// Copyright (c) 2018 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Konstantinos Katrioplas

#ifndef CGAL_OPTIMAL_BOUNDING_BOX_OBB_H
#define CGAL_OPTIMAL_BOUNDING_BOX_OBB_H

#include <CGAL/Optimal_bounding_box/population.h>
#include <CGAL/Optimal_bounding_box/evolution.h>
#include <CGAL/Optimal_bounding_box/helper.h>

#include <CGAL/assertions.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Iso_cuboid_3.h>
#include <CGAL/Simple_cartesian.h>

#ifdef CGAL_OPTIMAL_BOUNDING_BOX_BENCHMARKS
  #include <CGAL/Timer.h>
#endif

#if defined(CGAL_EIGEN3_ENABLED)
#include <CGAL/Eigen_linear_algebra_traits.h>
#endif

#include <iostream>
#include <iterator>
#include <vector>

namespace CGAL {

namespace Optimal_bounding_box {

// works on matrices only
/// \cond SKIP_IN_MANUAL
template <typename Linear_algebra_traits, typename Vertex, typename Matrix>
void post_processing(const Matrix& points, Vertex& R, Matrix& obb)
{
  CGAL_assertion(points.cols() == 3);
  CGAL_assertion(R.rows() == 3);
  CGAL_assertion(R.cols() == 3);
  CGAL_assertion(obb.rows() == 8);
  CGAL_assertion(obb.cols() == 3);

  // 1) rotate points with R
  Matrix rotated_points(points.rows(), points.cols());
  rotated_points = points * Linear_algebra_traits::transpose(R);

  // 2) get AABB from rotated points
  typedef CGAL::Simple_cartesian<double> K;
  typedef K::Point_3 Point;
  typedef typename Linear_algebra_traits::Index index;

  // Simplex -> std::vector
  std::vector<Point> v_points;
  for(index i = 0; i < static_cast<index>(rotated_points.rows()); ++i)
  {
    Point p(rotated_points(i, 0), rotated_points(i, 1), rotated_points(i, 2));
    v_points.push_back(p);
  }
  CGAL::Bbox_3 bbox;
  bbox = bbox_3(v_points.begin(), v_points.end());
  K::Iso_cuboid_3 ic(bbox);

  Matrix aabb(8, 3);
  for(std::size_t i = 0; i < 8; ++i)
  {
    aabb.set_coef(i, 0, ic[i].x());
    aabb.set_coef(i, 1, ic[i].y());
    aabb.set_coef(i, 2, ic[i].z());
  }

  // 3) apply inverse rotation to rotated AABB
  obb = aabb * R;
}
/// \endcond

/// \ingroup OBB_grp
/// calculates the optimal bounding box.
///
/// @tparam Point the point type
/// @tparam LinearAlgebraTraits a model of `LinearAlgebraTraits`. If no instance of `LinearAlgebraTraits`
/// is provided, then `CGAL::Eigen_linear_algebra_traits` is used.
///
/// @param points the input points that are included in the optimal bounding box.
/// @param obb_points the eight points of the optimal bounding box to be calculated.
/// @param use_ch a bool flag to indicating whether to use the convex hull of the input points
/// as an optimization step.
template <typename Point, typename LinearAlgebraTraits>
void compute_optimal_bounding_box(const std::vector<Point>& points,
                                  std::vector<Point>& obb_points,
                                  LinearAlgebraTraits&,
                                  bool use_ch)
{
  CGAL_assertion(points.size() >= 3);

  if(obb_points.size() != 8)
    obb_points.resize(8);
  CGAL_assertion(obb_points.size() == 8);

  // eigen linear algebra traits
  typedef typename LinearAlgebraTraits::MatrixXd MatrixXd;
  typedef typename LinearAlgebraTraits::Matrix3d Matrix3d;
  MatrixXd points_mat;

  if(use_ch) // get the ch3
  {
    std::vector<Point> ch_points;
    CGAL::extreme_points_3(points, std::back_inserter(ch_points));
    CGAL::Optimal_bounding_box::fill_matrix(ch_points, points_mat);
  }
  else
  {
    CGAL::Optimal_bounding_box::fill_matrix(points, points_mat);
  }

#ifdef CGAL_OPTIMAL_BOUNDING_BOX_BENCHMARKS
  CGAL::Timer timer;
#endif

  std::size_t max_generations = 100;
  Population<LinearAlgebraTraits> pop(50);

#ifdef CGAL_OPTIMAL_BOUNDING_BOX_BENCHMARKS
  timer.start();
#endif

  CGAL::Optimal_bounding_box::Evolution<LinearAlgebraTraits> search_solution(pop, points_mat);

#ifdef CGAL_OPTIMAL_BOUNDING_BOX_BENCHMARKS
  timer.stop();
  std::cout << "constructor: " << timer.time() << std::endl;
  timer.reset();
  timer.start();
#endif

  search_solution.evolve(max_generations);

#ifdef CGAL_OPTIMAL_BOUNDING_BOX_BENCHMARKS
  timer.stop();
  std::cout << "evolve: " << timer.time() << std::endl;
  timer.reset();
  timer.start();
#endif

  Matrix3d rotation = search_solution.get_best();

#ifdef CGAL_OPTIMAL_BOUNDING_BOX_BENCHMARKS
  timer.stop();
  std::cout << "get best: " << timer.time() << std::endl;
#endif

  MatrixXd obb; // may be preallocated at compile time
  obb.resize(8, 3);

#ifdef CGAL_OPTIMAL_BOUNDING_BOX_BENCHMARKS
  timer.reset();
  timer.start();
#endif

  post_processing<LinearAlgebraTraits>(points_mat, rotation, obb);

#ifdef CGAL_OPTIMAL_BOUNDING_BOX_BENCHMARKS
  timer.stop();
  std::cout << "post porcessing: " << timer.time() << std::endl;
#endif

  // matrix -> vector
  for(std::size_t i = 0; i < 8; ++i)
  {
    Point p(obb(i, 0), obb(i, 1), obb(i, 2));
    obb_points[i] = p;
  }
}

template <typename Point>
void compute_optimal_bounding_box(const std::vector<Point>& points,
                                  std::vector<Point>& obb_points,
                                  bool use_ch)
{
#if defined(CGAL_EIGEN3_ENABLED)
  typedef CGAL::Eigen_linear_algebra_traits Linear_algebra_traits;
#else
  #pragma message("Error: You must either provide linear traits or link CGAL with the Eigen library")
  Linear_algebra_traits; // no parameter provided, and Eigen is not enabled --> don't compile!
#endif

  Linear_algebra_traits la_traits;
  compute_optimal_bounding_box(points, obb_points, la_traits, use_ch);
}

/// \ingroup OBB_grp
/// calculates the optimal bounding box.
///
/// @tparam PolygonMesh a model of `FaceListGraph` 
/// @tparam LinearAlgebraTraits a model of `LinearAlgebraTraits`. If no instance of `LinearAlgebraTraits`
/// is provided, then `CGAL::Eigen_linear_algebra_traits` is used.
///
/// @param pmesh the input mesh.
/// @param obbmesh the hexaedron mesh to be built out of the optimal bounding box.
/// @param la_traits an instance of the linear algebra traits.
/// @param use_ch a bool flag to indicating whether to use the convex hull of the input points
/// as an optimization step.
template <typename PolygonMesh, typename LinearAlgebraTraits>
void compute_optimal_bounding_box(const PolygonMesh& pmesh,
                                  PolygonMesh& obbmesh,
                                  LinearAlgebraTraits& la_traits,
                                  bool use_ch)
{
  CGAL_assertion(vertices(pmesh).size() >= 3);

  if(vertices(pmesh).size() <= 3)
  {
    std::cerr << "Not enough points in the mesh!\n";
    return;
  }

  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
  typedef typename boost::property_map<PolygonMesh, CGAL::vertex_point_t>::type Vpm;
  typedef typename boost::property_traits<Vpm>::value_type Point;

  std::vector<Point> points;
  Vpm pmap = get(boost::vertex_point, pmesh);
  BOOST_FOREACH(vertex_descriptor v, vertices(pmesh))
    points.push_back(get(pmap, v));


  std::vector<Point> obb_points;
  compute_optimal_bounding_box(points, obb_points, la_traits, use_ch);

  CGAL::make_hexahedron(obb_points[0], obb_points[1], obb_points[2], obb_points[3],
                        obb_points[4], obb_points[5], obb_points[6], obb_points[7], obbmesh);
}

template <typename PolygonMesh>
void compute_optimal_bounding_box(const PolygonMesh& pmesh,
                                  PolygonMesh& obbmesh,
                                  bool use_ch)
{
#if defined(CGAL_EIGEN3_ENABLED)
  typedef CGAL::Eigen_linear_algebra_traits Linear_algebra_traits;
#else
  #pragma message("Error: You must either provide linear traits or link CGAL with the Eigen library")
  Linear_algebra_traits; // no parameter provided, and Eigen is not enabled --> don't compile!
#endif

  Linear_algebra_traits la_traits;
  compute_optimal_bounding_box(pmesh, obbmesh, la_traits, use_ch);
}

} // end namespace Optimal_bounding_box
} // end namespace CGAL

#endif // CGAL_OPTIMAL_BOUNDING_BOX_OBB_H
