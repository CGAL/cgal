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

// Actually, this is very slow because the "m_points_ds->insert"
// cleans the tree, which is thus built at each query_ANN call
#ifdef CGAL_TC_USE_SLOW_BUT_ACCURATE_SPARSIFIER
template <typename Kernel, typename Point_container>
class Point_sparsifier
{
public:
  typedef typename Kernel::FT                         FT;
  typedef typename Point_container::value_type        Point;
  typedef typename CGAL::Tangential_complex_::Point_cloud_data_structure<
    Kernel, Point_container>                          Points_ds;
  typedef typename Points_ds::KNS_range               KNS_range;

  // We can't instantiate m_points_ds right now since it requires that
  // points is not empty (which be the case here)
  Point_sparsifier(Point_container &points, 
                   FT sparsity = FT(INPUT_SPARSITY*INPUT_SPARSITY))
    : m_points(points), m_sparsity(sparsity), m_points_ds(NULL)
  {}

  bool try_to_insert_point(const Point &p)
  {
    if (m_points_ds == NULL)
    {
      m_points.push_back(p);
      m_points_ds = new Points_ds(m_points);
      m_points_ds->insert(0);
      return true;
    }
    else
    {
      KNS_range kns_range = m_points_ds->query_ANN(p, 1, false);
      if (kns_range.begin()->second >= m_sparsity)
      {
        m_points.push_back(p);
        m_points_ds->insert(m_points.size() - 1);
        return true;
      }
    }
    
    return false;
  }

private:
  Point_container     & m_points;
  Points_ds           * m_points_ds;
  FT                    m_sparsity;
};
#endif

template <typename Kernel, typename Point_container>
std::vector<typename Point_container::value_type>
sparsify_point_set(
  const Kernel &k, Point_container const& input_pts,
  typename Kernel::FT min_squared_dist)
{
  typedef typename Point_container::value_type  Point;
  typedef typename CGAL::Tangential_complex_::Point_cloud_data_structure<
    Kernel, Point_container>                    Points_ds;
  typedef typename Points_ds::INS_iterator      INS_iterator;
  typedef typename Points_ds::INS_range         INS_range;

  typename Kernel::Squared_distance_d sqdist = k.squared_distance_d_object();
  
#ifdef CGAL_TC_PROFILING
    Wall_clock_timer t;
#endif

  // Create the output container
  std::vector<typename Point_container::value_type> output;
  
  Points_ds points_ds(input_pts);

  // Parse the following points, and add them if they are not too close to
  // the other points
  std::size_t pt_idx = 0;
  for (typename Point_container::const_iterator it_pt = input_pts.begin() ;
       it_pt != input_pts.end();
       ++it_pt, ++pt_idx)
  {
    INS_range ins_range = points_ds.query_incremental_ANN(*it_pt);

    // Drop it if there is another point that:
    // - is closer that min_squared_dist
    // - and has a higher index
    for (INS_iterator nn_it = ins_range.begin() ;
        nn_it != ins_range.end() ;
        ++nn_it)
    {
      std::size_t neighbor_point_idx = nn_it->first;
      typename Kernel::FT sq_dist = nn_it->second;
      // The neighbor is further, we keep the point
      if (sq_dist >= min_squared_dist)
      {
        output.push_back(*it_pt);
        break;
      }
      // The neighbor is close and it has a higher index
      else if (neighbor_point_idx > pt_idx)
        break; // We drop the point
      // Otherwise, we go the next closest point
    }
  }
  
#ifdef CGAL_TC_PROFILING
    std::cerr << "Point set sparsified in " << t.elapsed()
              << " seconds." << std::endl;
#endif

  return output;
}

template<typename Point, typename OutputIterator>
bool load_points_from_file(
  const std::string &filename, OutputIterator points)
{
  std::ifstream in(filename);
  if (!in.is_open())
  {
    std::cerr << "Could not open '" << filename << "'" << std::endl;
    return false;
  }

  Point p;
  int dim_from_file;
  in >> dim_from_file;

  while(in >> p)
    *points++ = p;

  return true;
}

template <typename Kernel>
std::vector<typename Kernel::Point_d> generate_points_on_plane(std::size_t num_points)
{
  typedef typename Kernel::Point_d Point;
  typedef typename Kernel::FT FT;
  CGAL::Random rng;
  std::vector<Point> points;
  points.reserve(num_points);
#ifdef CGAL_TC_USE_SLOW_BUT_ACCURATE_SPARSIFIER
  Point_sparsifier<Kernel, std::vector<Point> > sparsifier(points);
#endif
  for (std::size_t i = 0 ; i < num_points ; )
  {
    FT x = rng.get_double(0, 5);
    FT y = rng.get_double(0, 5);
    Point p = Kernel().construct_point_d_object()(x, y, 0);
#ifdef CGAL_TC_USE_SLOW_BUT_ACCURATE_SPARSIFIER
    if (sparsifier.try_to_insert_point(p))
      ++i;
#else
    points.push_back(p);
    ++i;
#endif
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
#ifdef CGAL_TC_USE_SLOW_BUT_ACCURATE_SPARSIFIER
  Point_sparsifier<Kernel, std::vector<Point> > sparsifier(points);
#endif
  for (std::size_t i = 0 ; i < num_points ; )
  {
    FT x = rng.get_double(min_x, max_x);
    std::vector<FT> coords;
    coords.reserve(dim);
    for (int p = 1 ; p <= dim ; ++p)
      coords.push_back(std::pow(CGAL::to_double(x), p));
    Point p = Kernel().construct_point_d_object()(
      dim, coords.begin(), coords.end());
#ifdef CGAL_TC_USE_SLOW_BUT_ACCURATE_SPARSIFIER
    if (sparsifier.try_to_insert_point(p))
      ++i;
#else
    points.push_back(p);    
    ++i;
#endif
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
#ifdef CGAL_TC_USE_SLOW_BUT_ACCURATE_SPARSIFIER
  Point_sparsifier<Kernel, std::vector<Point> > sparsifier(points);
#endif
  for (std::size_t i = 0 ; i < num_points ; )
  {
    Point p = *generator++;
#ifdef CGAL_TC_USE_SLOW_BUT_ACCURATE_SPARSIFIER
    if (sparsifier.try_to_insert_point(p))
      ++i;
#else
    points.push_back(p);
    ++i;
#endif
  }
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
#ifdef CGAL_TC_USE_SLOW_BUT_ACCURATE_SPARSIFIER
  Point_sparsifier<Kernel, std::vector<Point> > sparsifier(points);
#endif
  for (std::size_t i = 0 ; i < num_points ; )
  {
    Point p = *generator++;
#ifdef CGAL_TC_USE_SLOW_BUT_ACCURATE_SPARSIFIER
    if (sparsifier.try_to_insert_point(p))
      ++i;
#else
    points.push_back(p);
    ++i;
#endif
  }
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
#ifdef CGAL_TC_USE_SLOW_BUT_ACCURATE_SPARSIFIER
  Point_sparsifier<Kernel, std::vector<Point> > sparsifier(points);
#endif
  for (std::size_t i = 0 ; i < num_points ; )
  {
    Point p = *generator++;
#ifdef CGAL_TC_USE_SLOW_BUT_ACCURATE_SPARSIFIER
    if (sparsifier.try_to_insert_point(p))
      ++i;
#else
    points.push_back(p);
    ++i;
#endif
  }
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
#ifdef CGAL_TC_USE_SLOW_BUT_ACCURATE_SPARSIFIER
  Point_sparsifier<Kernel, std::vector<Point> > sparsifier(points);
#endif
  for (std::size_t i = 0 ; i < num_points ; )
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
    Point p = Kernel().construct_point_d_object()(
      (a + b*tmp)*cos(u),
      (a + b*tmp)*sin(u),
      b*(sin(u/2)*sin(v) + cos(u/2)*sin(2.*v)));
#ifdef CGAL_TC_USE_SLOW_BUT_ACCURATE_SPARSIFIER
    if (sparsifier.try_to_insert_point(p))
      ++i;
#else
    points.push_back(p);
    ++i;
#endif
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
#ifdef CGAL_TC_USE_SLOW_BUT_ACCURATE_SPARSIFIER
  Point_sparsifier<Kernel, std::vector<Point> > sparsifier(points);
#endif
  for (std::size_t i = 0 ; i < num_points ; )
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
    Point p = Kernel().construct_point_d_object()(
      Kernel().construct_point_d_object()(
        (a + b*cos(v))*cos(u) /*+ rng.get_double(0, 0.01)*/,
        (a + b*cos(v))*sin(u) /*+ rng.get_double(0, 0.01)*/,
        b*sin(v)*cos(u/2)     /*+ rng.get_double(0, 0.01)*/,
        b*sin(v)*sin(u/2)     /*+ rng.get_double(0, 0.01)*/) );
#ifdef CGAL_TC_USE_SLOW_BUT_ACCURATE_SPARSIFIER
    if (sparsifier.try_to_insert_point(p))
      ++i;
#else
    points.push_back(p);
    ++i;
#endif
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
#ifdef CGAL_TC_USE_SLOW_BUT_ACCURATE_SPARSIFIER
  Point_sparsifier<Kernel, std::vector<Point> > sparsifier(points);
#endif
  for (std::size_t i = 0 ; i < num_points ; )
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
    
    Point p = Kernel().construct_point_d_object()(x1, x2, x3, x4, x5);
#ifdef CGAL_TC_USE_SLOW_BUT_ACCURATE_SPARSIFIER
    if (sparsifier.try_to_insert_point(p))
      ++i;
#else
    points.push_back(p);
    ++i;
#endif
  }
  return points;
}

#endif // CGAL_MESH_3_TEST_TEST_UTILITIES_H
