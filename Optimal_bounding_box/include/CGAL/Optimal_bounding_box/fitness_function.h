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
#include <CGAL/Optimal_bounding_box/population.h>
#include <limits>

#include <CGAL/Eigen_linear_algebra_traits.h>

namespace CGAL {
namespace Optimal_bounding_box {

template <typename Linear_algebra_traits>
double compute_fitness(const typename Linear_algebra_traits::Matrix3d& R,
                             typename Linear_algebra_traits::MatrixXd& data)
{
  // R: rotation matrix
  CGAL_assertion(R.cols() == 3);
  CGAL_assertion(R.rows() == 3);
  // data: points
  CGAL_assertion(data.cols() == 3);
  CGAL_assertion(data.rows() >= 3);

  typedef typename Linear_algebra_traits::Vector3d Vector3d;

  double xmin, xmax, ymin, ymax, zmin, zmax;
  for(std::size_t i = 0; i < data.rows(); i++){

    Vector3d vec = Linear_algebra_traits::row(data, i);
    vec = R * vec;

    if(i == 0){
      xmin = xmax = vec.coeff(0);
      ymin = ymax = vec.coeff(1);
      zmin = zmax = vec.coeff(2);
    }else {
      if(vec.coeff(0) < xmin) xmin = vec.coeff(0);
      if(vec.coeff(1) < ymin) ymin = vec.coeff(1);
      if(vec.coeff(2) < zmin) zmin = vec.coeff(2);
      if(vec.coeff(0) > xmax) xmax = vec.coeff(0);
      if(vec.coeff(1) > ymax) ymax = vec.coeff(1);
      if(vec.coeff(2) > zmax) zmax = vec.coeff(2);
    }
  }

  CGAL_assertion(xmax > xmin);
  CGAL_assertion(ymax > ymin);
  CGAL_assertion(zmax > zmin);

  // volume
  return ((xmax - xmin) * (ymax - ymin) * (zmax - zmin));
}

template <typename Vertex, typename Matrix>
struct Fitness_map
{
  typedef CGAL::Eigen_linear_algebra_traits Linear_algebra_traits; // to be added as a parameter

  Fitness_map(Population<Vertex>& p, Matrix& points) : pop(p), points(points)
  {}

  const Vertex get_best()
  {
    std::size_t simplex_id, vertex_id;
    double best_fitness = std::numeric_limits<int>::max();
    for(std::size_t i = 0; i < pop.size(); ++i)
    {
      for(std::size_t j =0; j < 4; ++j)
      {
        const Vertex vertex = pop[i][j];
        const double fitness = compute_fitness<Linear_algebra_traits>(vertex, points);
        if (fitness < best_fitness)
        {
          simplex_id = i;
          vertex_id = j;
          best_fitness = fitness;
        }
      }
    }

    return pop[simplex_id][vertex_id];
  }

  double get_best_fitness_value(Matrix& data)
  {
    const Vertex best_mat = get_best();
    return compute_fitness<Linear_algebra_traits>(best_mat, data);
  }

  Population<Vertex>& pop;
  Matrix& points;
};

}} // end namespaces






#endif //CGAL_FITNESS_FUNCTION_H



