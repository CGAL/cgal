// Copyright (c) 2014 INRIA Sophia-Antipolis (France).
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
//
//
// Author(s)     : Clement Jamin
//
//******************************************************************************
// File Description :
//
//******************************************************************************

#ifndef CGAL_TC_TEST_TEST_UTILITIES_H
#define CGAL_TC_TEST_TEST_UTILITIES_H

#include <CGAL/Random.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/Tangential_complex/Point_cloud.h>

template <typename Kernel, typename Point_container>
std::vector<typename Point_container::value_type>
sparsify_point_set(
  const Kernel &k, Point_container const& input_pts,
  typename Kernel::FT min_squared_dist)
{
  typedef typename Point_container::value_type Point;
  typedef typename CGAL::Tangential_complex_::Point_cloud_data_structure<Kernel,
                                                    Point_container>  Points_ds;

  typename Kernel::Squared_distance_d sqdist = k.squared_distance_d_object();

  // Create the output container and push the first point into it
  std::vector<typename Point_container::value_type> output;
  typename Point_container::const_iterator it_pt = input_pts.begin();
  output.push_back(*it_pt);
  ++it_pt;

  // Parse the following points, and add them if they are not too close to
  // the other points
  std::size_t c = 1;
  for ( ;
       it_pt != input_pts.end();
       ++it_pt)
  {
    Points_ds points_ds(output, 0, c);
    if (points_ds.query_ANN(*it_pt, 1).begin()->second >= min_squared_dist)
    {
      output.push_back(*it_pt);
      ++c;
    }
  }

  return output;
}

template <typename Kernel>
std::vector<typename Kernel::Point_d> generate_points_on_plane(std::size_t num_points)
{
  typedef typename Kernel::Point_d Point;
  typedef typename Kernel::FT FT;
  CGAL::Random rng;
  std::vector<Point> points;
  points.reserve(num_points);
  for (std::size_t i = 0 ; i != num_points ; ++i)
  {
    FT x = rng.get_double(0, 5);
    FT y = rng.get_double(0, 5);
    points.push_back(Kernel().construct_point_d_object()(x, y, 0));
  }
  return points;
}

template <typename Kernel>
std::vector<typename Kernel::Point_d> generate_points_on_moment_curve(
  std::size_t num_points, int dim,
  typename Kernel::FT min_x , typename Kernel::FT max_x)
{
  typedef typename Kernel::Point_d Point;
  typedef typename Kernel::FT FT;
  CGAL::Random rng;
  std::vector<Point> points;
  points.reserve(num_points);
  for (std::size_t i = 0 ; i != num_points ; ++i)
  {
    FT x = rng.get_double(min_x, max_x);
    std::vector<FT> coords;
    coords.reserve(dim);
    for (int p = 1 ; p <= dim ; ++p)
      coords.push_back(std::pow(CGAL::to_double(x), p));
    points.push_back(
      Kernel().construct_point_d_object()(dim, coords.begin(), coords.end()));
  }
  return points;
}

template <typename Kernel>
std::vector<typename Kernel::Point_d> generate_points_on_circle_2(
  std::size_t num_points, double radius)
{
  typedef typename Kernel::Point_d Point;
  CGAL::Random_points_on_circle_2<Point> generator(radius);
  std::vector<Point> points;
  points.reserve(num_points);
  for (std::size_t i = 0 ; i != num_points ; ++i)
    points.push_back(*generator++);
  return points;
}

template <typename Kernel>
std::vector<typename Kernel::Point_d> generate_points_on_sphere_3(
  std::size_t num_points, double radius)
{
  typedef typename Kernel::Point_d Point;
  CGAL::Random_points_on_sphere_3<Point> generator(radius);
  std::vector<Point> points;
  points.reserve(num_points);
  for (std::size_t i = 0 ; i != num_points ; ++i)
    points.push_back(*generator++);
  return points;
}

template <typename Kernel>
std::vector<typename Kernel::Point_d> generate_points_on_sphere_d(
  std::size_t num_points, int dim, double radius)
{
  typedef typename Kernel::Point_d Point;
  CGAL::Random_points_on_sphere_d<Point> generator(dim, radius);
  std::vector<Point> points;
  points.reserve(num_points);
  for (std::size_t i = 0 ; i != num_points ; ++i)
    points.push_back(*generator++);
  return points;
}

// a = big radius, b = small radius
template <typename Kernel>
std::vector<typename Kernel::Point_d> generate_points_on_klein_bottle_3D(
  std::size_t num_points, double a, double b, bool uniform = false)
{
  typedef typename Kernel::Point_d Point;
  typedef typename Kernel::FT FT;
  CGAL::Random rng;

  // if uniform
  std::size_t num_lines = (std::size_t)sqrt(num_points);
  std::size_t num_cols = num_points/num_lines + 1;

  std::vector<Point> points;
  points.reserve(num_points);
  for (std::size_t i = 0 ; i != num_points ; ++i)
  {
    FT u, v;
    if (uniform)
    {
      std::size_t k1 = i / num_lines;
      std::size_t k2 = i % num_lines;
      u = 6.2832 * k1 / num_lines;
      v = 6.2832 * k2 / num_lines;
    }
    else
    {
      u = rng.get_double(0, 6.2832);
      v = rng.get_double(0, 6.2832);
    }
    double tmp = cos(u/2)*sin(v) - sin(u/2)*sin(2.*v);
    points.push_back(Kernel().construct_point_d_object()(
      (a + b*tmp)*cos(u),
      (a + b*tmp)*sin(u),
      b*(sin(u/2)*sin(v) + cos(u/2)*sin(2.*v))));
  }
  return points;
}

// a = big radius, b = small radius
template <typename Kernel>
std::vector<typename Kernel::Point_d> generate_points_on_klein_bottle_4D(
  std::size_t num_points, double a, double b, bool uniform = false)
{
  typedef typename Kernel::Point_d Point;
  typedef typename Kernel::FT FT;
  CGAL::Random rng;

  // if uniform
  std::size_t num_lines = (std::size_t)sqrt(num_points);
  std::size_t num_cols = num_points/num_lines + 1;

  std::vector<Point> points;
  points.reserve(num_points);
  for (std::size_t i = 0 ; i != num_points ; ++i)
  {
    FT u, v;
    if (uniform)
    {
      std::size_t k1 = i / num_lines;
      std::size_t k2 = i % num_lines;
      u = 6.2832 * k1 / num_lines;
      v = 6.2832 * k2 / num_lines;
    }
    else
    {
      u = rng.get_double(0, 6.2832);
      v = rng.get_double(0, 6.2832);
    }
    points.push_back(Kernel().construct_point_d_object()(
      (a + b*cos(v))*cos(u) /*+ rng.get_double(0, 0.01)*/,
      (a + b*cos(v))*sin(u) /*+ rng.get_double(0, 0.01)*/,
      b*sin(v)*cos(u/2)     /*+ rng.get_double(0, 0.01)*/,
      b*sin(v)*sin(u/2)     /*+ rng.get_double(0, 0.01)*/) );
  }
  return points;
}


// a = big radius, b = small radius
template <typename Kernel>
std::vector<typename Kernel::Point_d>
generate_points_on_klein_bottle_variant_5D(
  std::size_t num_points, double a, double b, bool uniform = false)
{
  typedef typename Kernel::Point_d Point;
  typedef typename Kernel::FT FT;
  CGAL::Random rng;

  // if uniform
  std::size_t num_lines = (std::size_t)sqrt(num_points);
  std::size_t num_cols = num_points/num_lines + 1;

  std::vector<Point> points;
  points.reserve(num_points);
  for (std::size_t i = 0 ; i != num_points ; ++i)
  {
    FT u, v;
    if (uniform)
    {
      std::size_t k1 = i / num_lines;
      std::size_t k2 = i % num_lines;
      u = 6.2832 * k1 / num_lines;
      v = 6.2832 * k2 / num_lines;
    }
    else
    {
      u = rng.get_double(0, 6.2832);
      v = rng.get_double(0, 6.2832);
    }
    FT x1 = (a + b*cos(v))*cos(u);
    FT x2 = (a + b*cos(v))*sin(u);
    FT x3 = b*sin(v)*cos(u/2);
    FT x4 = b*sin(v)*sin(u/2);
    FT x5 = x1 + x2 + x3 + x4;

    points.push_back(Kernel().construct_point_d_object()(x1, x2, x3, x4, x5) );
  }
  return points;
}

#endif // CGAL_MESH_3_TEST_TEST_UTILITIES_H
