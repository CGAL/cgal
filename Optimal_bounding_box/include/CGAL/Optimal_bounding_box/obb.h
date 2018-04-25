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
#include <vector>
#include <fstream>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Iso_cuboid_3.h>


namespace CGAL {
namespace Optimal_bounding_box {



template <typename Simplex>
void find_global_minimum(Population<Simplex>& pop, Simplex points)
{

  genetic_algorithm(pop, points);

  std::size_t nelder_mead_iterations = 1;

  for(Simplex s : pop)
    nelder_mead(s, points, nelder_mead_iterations);

  //optional random mutations
}


template <typename Simplex>
void visualize_obb(Simplex data_points)
{


  typedef CGAL::Simple_cartesian<double> K;
  typedef K::Point_3 Point;
  typedef CGAL::Surface_mesh<Point> Mesh;

  // Simplex -> std::vector
  std::vector<Point> points;

  for(int i = 0; i < data_points.rows(); ++i)
  {
    Point p(data_points(i, 0), data_points(i, 1), data_points(i, 2));
    points.push_back(p);
  }


  CGAL::Bbox_3 bbox;
  bbox = bbox_3(points.begin(), points.end());
  K::Iso_cuboid_3 ic(bbox);

  Mesh mesh;
  CGAL::make_hexahedron(ic[0], ic[1], ic[2], ic[3], ic[4], ic[5], ic[6], ic[7], mesh);

  std::ofstream out("/tmp/box.off");
  out << mesh;
  out.close();


}





}} // end namespaces






#endif //CGAL_OBB_H



