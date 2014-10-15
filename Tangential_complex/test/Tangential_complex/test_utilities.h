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

template <typename Kernel>
std::vector<typename Kernel::Point_d> generate_points_on_plane(std::size_t num_points)
{
  typedef typename Kernel::Point_d Point;
  typedef typename Kernel::FT FT;
  CGAL::Random rng;
  std::vector<Point> points;
  points.reserve(NUM_POINTS);
  for (int i = 0 ; i != NUM_POINTS ; ++i)
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
  points.reserve(NUM_POINTS);
  for (int i = 0 ; i != NUM_POINTS ; ++i)
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
  CGAL::Random_points_on_circle_2<Point> generator(radius);
  std::vector<Point> points;
  points.reserve(NUM_POINTS);
  for (int i = 0 ; i != NUM_POINTS ; ++i)
    points.push_back(*generator++);
  return points;
}

template <typename Kernel>
std::vector<typename Kernel::Point_d> generate_points_on_sphere_3(
  std::size_t num_points, double radius)
{
  CGAL::Random_points_on_sphere_3<Point> generator(radius);
  std::vector<Point> points;
  points.reserve(NUM_POINTS);
  for (int i = 0 ; i != NUM_POINTS ; ++i)
    points.push_back(*generator++);
  return points;
}

template <typename Kernel>
std::vector<typename Kernel::Point_d> generate_points_on_sphere_d(
  std::size_t num_points, int dim, double radius)
{
  CGAL::Random_points_on_sphere_d<Point> generator(dim, radius);
  std::vector<Point> points;
  points.reserve(NUM_POINTS);
  for (int i = 0 ; i != NUM_POINTS ; ++i)
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
  int num_lines = (int)sqrt(NUM_POINTS);
  int num_cols = NUM_POINTS/num_lines + 1;

  std::vector<Point> points;
  points.reserve(NUM_POINTS);
  for (int i = 0 ; i != NUM_POINTS ; ++i)
  {
    FT u, v;
    if (uniform)
    {
      int k1 = i / num_lines;
      int k2 = i % num_lines;
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
  int num_lines = (int)sqrt(NUM_POINTS);
  int num_cols = NUM_POINTS/num_lines + 1;

  std::vector<Point> points;
  points.reserve(NUM_POINTS);
  for (int i = 0 ; i != NUM_POINTS ; ++i)
  {
    FT u, v;
    if (uniform)
    {
      int k1 = i / num_lines;
      int k2 = i % num_lines;
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

#endif // CGAL_MESH_3_TEST_TEST_UTILITIES_H
