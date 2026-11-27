// Copyright (c) 2026 GeometryFactory
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Ange Clement

#ifndef CGAL_FEATURE_GRAPH_INTERNAL_IMAGE_KHALIMSKY_COORDINATES_H
#define CGAL_FEATURE_GRAPH_INTERNAL_IMAGE_KHALIMSKY_COORDINATES_H

#include <CGAL/license/Feature_graph.h>

#include <array>
#include <CGAL/assertions.h>

namespace CGAL
{
namespace Feature_graph
{
namespace internal
{

template <typename Integer, typename IntegerPoint>
Integer get_Khalimsky_index_from_Khalimsky_coordinates(
  const Integer& max_x, const Integer& max_y,
  const IntegerPoint& integer_coordinates)
{
  CGAL_assertion(integer_coordinates[0] >= 0);
  CGAL_assertion(integer_coordinates[1] >= 0);
  CGAL_assertion(integer_coordinates[2] >= 0);
  return
    integer_coordinates[2]*max_x*max_y*4 + integer_coordinates[1]*max_x*2 + integer_coordinates[0];
}

template<typename Integer, typename Point>
Integer get_Khalimsky_index_from_image_coordinates(
  const Point& origin, const Point& voxel_size, const Integer& max_x, const Integer& max_y,
  const Point& point)
{
  CGAL_assertion(voxel_size[0] != 0);
  CGAL_assertion(voxel_size[1] != 0);
  CGAL_assertion(voxel_size[2] != 0);
  double image_space_x = (point[0] - static_cast<double>(origin[0])) /
    static_cast<double>(voxel_size[0]);
  double image_space_y = (point[1] - static_cast<double>(origin[1])) /
    static_cast<double>(voxel_size[1]);
  double image_space_z = (point[2] - static_cast<double>(origin[2])) /
    static_cast<double>(voxel_size[2]);
  // round to nearest coords (+.5) in Khalimsky coordinates (*2)
  std::array<Integer, 3> khalimsky_coordinates = {
    static_cast<Integer>(2*image_space_x +.5),
    static_cast<Integer>(2*image_space_y +.5),
    static_cast<Integer>(2*image_space_z +.5)
  };
  // get index
  return get_Khalimsky_index_from_Khalimsky_coordinates(max_x, max_y, khalimsky_coordinates);
}

} /* namespace internal */

} /* namespace Feature_graph */

} /* namespace CGAL */


#endif // CGAL_FEATURE_GRAPH_INTERNAL_IMAGE_KHALIMSKY_COORDINATES_H