// Copyright (c) 2014  INRIA Sophia-Antipolis (France).
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

#include <CGAL/point_generators_3.h>

#include <cassert>
#include <iostream>
#include <vector>

#include <boost/mpl/if.hpp>

template <class Parallel_triangulation>
void
_test_cls_parallel_triangulation_3(const Parallel_triangulation &)
{
  const int NUM_INSERTED_POINTS = 5000;

  typedef Parallel_triangulation                                Cls;
  typedef typename Cls::Vertex_handle                           Vertex_handle;
  typedef typename boost::mpl::if_<typename Cls::Weighted_tag,
                                   typename Cls::Point, Cls>::type::Point 
                                                                Point;
  
  CGAL::Random_points_in_cube_3<Point> rnd(1.);

  // Construction from a vector of points
  std::vector<Point> points;
  points.reserve(NUM_INSERTED_POINTS);
  for (int i = 0; i != NUM_INSERTED_POINTS; ++i)
    points.push_back(*rnd++);
  
  // Construct the locking data-structure, using the bounding-box of the points
  typename Cls::Lock_data_structure locking_ds(
    CGAL::Bbox_3(-1., -1., -1., 1., 1., 1.), 50);
  // Contruct the triangulation in parallel
  std::cout << "Construction and parallel insertion" << std::endl;
  Cls tr(points.begin(), points.end(), &locking_ds);

  std::cout << "Parallel removal" << std::endl;
  // Remove the first 100,000 vertices
  std::vector<Vertex_handle> vertices_to_remove;
  typename Cls::Finite_vertices_iterator vit = tr.finite_vertices_begin();
  for (int i = 0 ; i < NUM_INSERTED_POINTS/10 ; ++i)
    vertices_to_remove.push_back(vit++);
  // Parallel remove
  tr.remove(vertices_to_remove.begin(), vertices_to_remove.end());
  
  assert(tr.is_valid());
  tr.clear();
  assert(tr.is_valid());
  assert(tr.dimension()==-1);
  assert(tr.number_of_vertices()==0);
}
