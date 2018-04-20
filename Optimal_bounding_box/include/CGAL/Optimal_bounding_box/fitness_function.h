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

#ifndef CGAL_FITNESS_FUNCTION_H
#define CGAL_FITNESS_FUNCTION_H

#include <CGAL/Bbox_3.h>
#include <vector>


namespace CGAL {
namespace Optimal_bounding_box {

template<typename Matrix>
const double compute_fitness(Matrix& R, Matrix& data)
{

  // R: rotation matrix
  CGAL_assertion(R.cols() == 3);
  CGAL_assertion(R.rows() == 3);
  // data: points
  CGAL_assertion(data.cols() == 3);
  CGAL_assertion(data.rows() >= 3);

  CGAL_assertion(R.rows() == data.cols());

  // rotate points
  Matrix RT = R.transpose();
  Matrix rotated_data = data * RT;
  CGAL_assertion(rotated_data.cols() == data.cols());
  CGAL_assertion(rotated_data.rows() == data.rows());

  // AABB: take mins and maxs
  double xmin = rotated_data.col(0).minCoeff();
  double xmax = rotated_data.col(0).maxCoeff();
  double ymin = rotated_data.col(1).minCoeff();
  double ymax = rotated_data.col(1).maxCoeff();
  double zmin = rotated_data.col(2).minCoeff();
  double zmax = rotated_data.col(2).maxCoeff();

  double x_dim = abs(xmax - xmin);
  double y_dim = abs(ymax - ymin);
  double z_dim = abs(zmax - zmin);

  // volume
  return (x_dim * y_dim * z_dim);

}




}} // end namespaces






#endif //CGAL_FITNESS_FUNCTION_H



