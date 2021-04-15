// Copyright (c) 2014  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Clement Jamin

#include <CGAL/point_generators_3.h>

#include <cassert>
#include <iostream>
#include <vector>

#include <boost/mpl/if.hpp>

template <typename Parallel_triangulation,
          typename WeightedTag = typename Parallel_triangulation::Weighted_tag>
struct PTR_random_pts_generator
{
  typedef typename Parallel_triangulation::Point       Point;

  void operator()(const int num, CGAL::Random& rnd, std::vector<Point>& points) const
  {
    CGAL::Random_points_in_cube_3<Point> gen(1., rnd);

    points.reserve(num);
    for(int i=0; i!=num; ++i)
      points.push_back(*gen++);
  }
};

template <typename Parallel_triangulation>
struct PTR_random_pts_generator<Parallel_triangulation, CGAL::Tag_true /*weighted triangulation*/>
{
  typedef typename Parallel_triangulation::Bare_point      Bare_point;
  typedef typename Parallel_triangulation::Weighted_point  Weighted_point;

  void operator()(const int num, CGAL::Random& rnd, std::vector<Weighted_point>& points) const
  {
    CGAL::Random_points_in_cube_3<Bare_point> gen(1., rnd);

    points.reserve(num);
    for(int i=0; i!=num; ++i)
      points.push_back(Weighted_point(*gen++, rnd.get_double(-1., 1.)));
  }
};

template <class Parallel_triangulation>
void
_test_cls_parallel_triangulation_3(const Parallel_triangulation &)
{
  typedef Parallel_triangulation                                Cls;

  typedef typename Cls::Vertex_handle                           Vertex_handle;
  typedef typename Cls::Point                                   Point;

  typedef std::vector<Point>                                    Point_container;

  CGAL::Random rnd;
  std::cout << "Seed: " << rnd.get_seed() << std::endl;

  const int num_insert = 5000;
  Point_container points;
  PTR_random_pts_generator<Parallel_triangulation> points_gen;
  points_gen(num_insert, rnd, points);

  // Construct the locking data-structure, using the bounding-box of the points
  typename Cls::Lock_data_structure locking_ds(CGAL::Bbox_3(-1., -1., -1., 1., 1., 1.), 50);

  // Contruct the triangulation in parallel
  std::cout << "Construction and parallel insertion" << std::endl;
  Cls tr(points.begin(), points.end(), &locking_ds);
  std::cout << "Triangulation has " << tr.number_of_vertices() << " vertices" << std::endl;
  assert(tr.is_valid());

  std::cout << "Parallel removal" << std::endl;
  std::vector<Vertex_handle> vertices_to_remove;
  typename Cls::Finite_vertices_iterator vit = tr.finite_vertices_begin();

  const std::size_t num_remove = tr.number_of_vertices() / 10;
  std::cout << "Removing " << num_remove << " from " << tr.number_of_vertices() << " vertices" << std::endl;

  for(std::size_t i=0 ; i<num_remove; ++i) {
    vertices_to_remove.push_back(vit++);
  }

  // Parallel remove
  tr.remove(vertices_to_remove.begin(), vertices_to_remove.end());
  std::cout << "Now, " << tr.number_of_vertices() << " vertices are left" << std::endl;
  assert(tr.is_valid());

  tr.clear();
  assert(tr.is_valid());
  assert(tr.dimension()==-1);
  assert(tr.number_of_vertices()==0);
}
