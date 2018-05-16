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

#ifndef CGAL_OBB_H
#define CGAL_OBB_H

#include <CGAL/Optimal_bounding_box/optimization_algorithms.h>
#include <CGAL/Optimal_bounding_box/population.h>
#include <CGAL/Optimal_bounding_box/evolution.h>
#include <vector>
#include <fstream>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Iso_cuboid_3.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#if defined(CGAL_EIGEN3_ENABLED)
#include <CGAL/Eigen_linear_algebra_traits.h>
#endif

namespace CGAL {
namespace Optimal_bounding_box {


#if defined(CGAL_EIGEN3_ENABLED)
typedef CGAL::Eigen_linear_algebra_traits Linear_algebra_traits;
#endif


template <typename Vertex, typename Matrix>
void evolution(Vertex& R, Matrix& points, std::size_t max_generations) // todo: points is const
{
  CGAL_assertion(points.rows() >= 3);
  CGAL_assertion(points.cols() == 3);
  CGAL_assertion(R.rows() == 3);
  CGAL_assertion(R.cols() == 3);

  Population<Vertex> pop(50);
  std::size_t nelder_mead_iterations = 20;

  double prev_fit_value = 0;
  double new_fit_value = 0;
  double tolerance = 1e-2;
  int stale = 0;

  for(std::size_t t = 0; t < max_generations; ++t)
  {

#ifdef OBB_DEBUG
    std::cout << "generation= " << t << "\n";
#endif

    genetic_algorithm(pop, points);

#ifdef OBB_DEBUG
    //std::cout << "pop after genetic" << std::endl;
    //pop.show_population();
    //std::cout << std::endl;
#endif

    for(std::size_t s = 0; s < pop.size(); ++s)
      nelder_mead(pop[s], points, nelder_mead_iterations);

#ifdef OBB_DEBUG
    //std::cout << "pop after nelder mead: " << std::endl;
    //pop.show_population();
    //std::cout << std::endl;

    // debugging
    Fitness_map<Matrix> fitness_map(pop, points);
    Matrix R_now = fitness_map.get_best();
    //std::cout << "det= " << R_now.determinant() << std::endl;
    std::cout << "det= " << determinant(R_now) << std::endl;
#endif

    // stopping criteria
    Fitness_map<Vertex, Matrix> fitness_map(pop, points);
    new_fit_value = fitness_map.get_best_fitness_value(points);

    double difference = new_fit_value - prev_fit_value;
    if(abs(difference) < tolerance * new_fit_value)
      stale++;

    if(stale == 5)
      break;

    prev_fit_value = new_fit_value;
  }

  Fitness_map<Vertex, Matrix> fitness_map(pop, points);
  R = fitness_map.get_best();
}

// works on matrices only
template <typename Vertex, typename Matrix_dynamic, typename Matrix_fixed>
void post_processing(const Matrix_dynamic& points, Vertex& R, Matrix_fixed& obb)
{
  CGAL_assertion(points.cols() == 3);
  CGAL_assertion(R.rows() == 3);
  CGAL_assertion(R.cols() == 3);
  CGAL_assertion(obb.rows() == 8);
  CGAL_assertion(obb.cols() == 3);

  // 1) rotate points with R
  Matrix_dynamic rotated_points(points.rows(), points.cols());
  rotated_points = points * Linear_algebra_traits::transpose(R);

  // 2) get AABB from rotated points
  typedef CGAL::Simple_cartesian<double> K;
  typedef K::Point_3 Point;
  std::vector<Point> v_points; // Simplex -> std::vector
  for(std::size_t i = 0; i < rotated_points.rows(); ++i)
  {
    //Point p(rotated_points.coeff(i, 0), rotated_points.coeff(i, 1), rotated_points.coeff(i, 2));
    Point p(rotated_points(i, 0), rotated_points(i, 1), rotated_points(i, 2));
    v_points.push_back(p);
  }
  CGAL::Bbox_3 bbox;
  bbox = bbox_3(v_points.begin(), v_points.end());
  K::Iso_cuboid_3 ic(bbox);

  Matrix_fixed aabb;
  // preallocate sanity: if Matrix is not preallocated in compile time
  if(aabb.cols() != 3 && aabb.rows() != 8)
    aabb.resize(8, 3);

  for(std::size_t i = 0; i < 8; ++i)
  {
    /*
    aabb.coeffRef(i, 0) = ic[i].x();
    aabb.coeffRef(i, 1) = ic[i].y();
    aabb.coeffRef(i, 2) = ic[i].z();
    */

    aabb.set_coef(i, 0, ic[i].x());
    aabb.set_coef(i, 1, ic[i].y());
    aabb.set_coef(i, 2, ic[i].z());


  }

  // 3) apply inverse rotation to rotated AABB
  obb = aabb * R;
}


// to be moved into a helper
template <typename Matrix, typename Point>
void fill_matrix(std::vector<Point>& v_points, Matrix& points_mat)
{
  points_mat.resize(v_points.size(), 3);
  for(std::size_t i = 0; i < v_points.size(); ++i)
  {
    Point p = v_points[i];
    /*
    points_mat.coeffRef(i, 0) = p.x();
    points_mat.coeffRef(i, 1) = p.y();
    points_mat.coeffRef(i, 2) = p.z();
    */

    points_mat.set_coef(i, 0, p.x());
    points_mat.set_coef(i, 1, p.y());
    points_mat.set_coef(i, 2, p.z());

  }
}


// with linear algebra traits
template <typename Point>
void find_obb(std::vector<Point>& points, std::vector<Point>& obb_points, bool use_ch)
{
  CGAL_assertion(points.size() >= 3);

  if(obb_points.size() != 8) // temp sanity until the API is decided
    obb_points.resize(8);
  CGAL_assertion(obb_points.size() == 8);

  // could be used from a t. parameter
  typedef CGAL::Eigen_linear_algebra_traits Linear_algebra_traits;
  typedef typename Linear_algebra_traits::MatrixXd MatrixXd;
  typedef typename Linear_algebra_traits::Matrix3d Matrix3d;

  MatrixXd points_mat;

  // get the ch3
  if(use_ch)
  {
    // find the ch - todo: template kernel
    typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
    CGAL::Polyhedron_3<Kernel> poly;
    convex_hull_3(points.begin(), points.end(), poly);
    std::vector<Kernel::Point_3> ch_points(poly.points_begin(), poly.points_end());

    // points: vector -> matrix
    fill_matrix(ch_points, points_mat);
  }
  else
  {
    // points: vector -> matrix
    fill_matrix(points, points_mat);
  }

  std::size_t max_generations = 100;
  Population<Matrix3d> pop(50);
  CGAL::Optimal_bounding_box::Evolution<Linear_algebra_traits>
      search_solution(pop, points_mat);
  search_solution.evolve(max_generations);
  Matrix3d rotation = search_solution.get_best();

  MatrixXd obb; // could be preallocated at compile time
  obb.resize(8, 3);

  post_processing(points_mat, rotation, obb);

  // matrix -> vector
  for(std::size_t i = 0; i < 8; ++i)
  {
    Point p(obb(i, 0), obb(i, 1), obb(i, 2));
    obb_points[i] = p;
  }

}


/*
  //USES ::EVOLUTION

/// @param points point coordinates of the input mesh.
/// @param obb_points the 8 points of the obb.
/// @param use convex hull or not.
///
/// todo named parameters: max iterations, population size, tolerance.
template <typename Point>
void find_obb(std::vector<Point>& points, std::vector<Point>& obb_points, bool use_ch, bool use_eigen)
{
  CGAL_assertion(points.size() >= 3);

  if(obb_points.size() != 8) // temp sanity until the API is decided
    obb_points.resize(8);
  CGAL_assertion(obb_points.size() == 8);

  // using eigen internally
  typedef Eigen::MatrixXd MatrixXd; // for point data
  typedef Eigen::Matrix3d Matrix3d; // for matrices in simplices

  MatrixXd points_mat;

  // get the ch3
  if(use_ch)
  {
    // find the ch - todo: template kernel
    typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
    CGAL::Polyhedron_3<Kernel> poly;
    convex_hull_3(points.begin(), points.end(), poly);
    std::vector<Kernel::Point_3> ch_points(poly.points_begin(), poly.points_end());

    // points: vector -> matrix
    fill_matrix(ch_points, points_mat);
  }
  else
  {
    // points: vector -> matrix
    fill_matrix(points, points_mat);
  }

  //Matrix3d R; // todo: sort out preallocation
  Matrix3d R(3, 3);
  std::size_t max_generations = 100;
  CGAL::Optimal_bounding_box::evolution(R, points_mat, max_generations);

  //Eigen::Matrix<double, 8, 3> obb;

  MatrixXd obb; // compile time preallocation??
  obb.resize(8, 3);

  CGAL::Optimal_bounding_box::post_processing(points_mat, R, obb);

  // matrix -> vector
  for(std::size_t i = 0; i < 8; ++i)
  {
    //Point p(obb.coeff(i, 0), obb.coeff(i, 1), obb.coeff(i, 2));
    Point p(obb(i, 0), obb(i, 1), obb(i, 2));
    obb_points[i] = p;
  }
}

*/

/// @param pmesh the input mesh.
/// @param obbmesh the obb in a hexahedron pmesh.
/// @param use convex hull or not.
///
/// todo named parameters: max iterations, population size, tolerance.
template <typename PolygonMesh>
void find_obb(PolygonMesh& pmesh, PolygonMesh& obbmesh, bool use_ch)
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
  //typedef typename Kernel_traits<Point>::Kernel Kernel;

  std::vector<Point> points;
  Vpm pmap = get(boost::vertex_point, pmesh);
  BOOST_FOREACH(vertex_descriptor v, vertices(pmesh))
    points.push_back(get(pmap, v));

  std::vector<Point> obb_points;
  find_obb(points, obb_points, use_ch);

  CGAL::make_hexahedron(obb_points[0], obb_points[1], obb_points[2], obb_points[3],
      obb_points[4], obb_points[5],
      obb_points[6], obb_points[7], obbmesh);
}


}} // end namespaces






#endif //CGAL_OBB_H



