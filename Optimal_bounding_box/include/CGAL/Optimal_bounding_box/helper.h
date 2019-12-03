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

#ifndef CGAL_OPTIMAL_BOUNDING_BOX_HELPER_H
#define CGAL_OPTIMAL_BOUNDING_BOX_HELPER_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Bbox_3.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/number_utils.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <fstream>
#include <string>
#include <vector>

namespace CGAL {
namespace Optimal_bounding_box {

template <typename Matrix, typename Point>
void fill_matrix(const std::vector<Point>& v_points, Matrix& points_mat)
{
  points_mat.resize(v_points.size(), 3);
  for(std::size_t i = 0; i < v_points.size(); ++i)
  {
    Point p = v_points[i];
    points_mat.set_coef(i, 0, CGAL::to_double(p.x()));
    points_mat.set_coef(i, 1, CGAL::to_double(p.y()));
    points_mat.set_coef(i, 2, CGAL::to_double(p.z()));
  }
}

template <typename SurfaceMesh, typename Matrix>
void sm_to_matrix(SurfaceMesh& sm, Matrix& mat)
{
  typedef typename boost::property_map<SurfaceMesh, boost::vertex_point_t>::const_type Vpm;
  typedef typename boost::property_traits<Vpm>::reference Point_ref;
  typedef typename boost::graph_traits<SurfaceMesh>::vertex_descriptor vertex_descriptor;
  Vpm vpm = get(boost::vertex_point, sm);

  mat.resize(vertices(sm).size(), 3);
  std::size_t i = 0;
  for(vertex_descriptor v : vertices(sm))
  {
    Point_ref p = get(vpm, v);
    mat.set_coef(i, 0, CGAL::to_double(p.x()));
    mat.set_coef(i, 1, CGAL::to_double(p.y()));
    mat.set_coef(i, 2, CGAL::to_double(p.z()));
    ++i;
  }
}

template <typename Point>
double calculate_volume(const std::vector<Point>& points)
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

  CGAL::Bbox_3 bbox = bbox_3(points.begin(), points.end());
  K::Iso_cuboid_3 ic(bbox);
  return ic.volume();
}

// it is called after post processing in debug only
template <typename Matrix>
void matrix_to_mesh_and_draw(Matrix& data_points, std::string filename)
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

  Mesh mesh;
  CGAL::make_hexahedron(points[0], points[1], points[2], points[3],
                        points[4], points[5], points[6], points[7], mesh);

  std::ofstream out(filename);
  out << mesh;
  out.close();
}

} // end namespace Optimal_bounding_box
} // end namespace CGAL

#endif // CGAL_OPTIMAL_BOUNDING_BOX_HELPER_H
