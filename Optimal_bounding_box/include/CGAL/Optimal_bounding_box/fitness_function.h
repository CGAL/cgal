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

namespace CGAL {
namespace Optimal_bounding_box {

template<typename Vertex, typename Matrix>
const double compute_fitness(const Vertex& R, const Matrix& data)
{

  // R: rotation matrix
  CGAL_assertion(R.cols() == 3);
  CGAL_assertion(R.rows() == 3);
  // data: points
  CGAL_assertion(data.cols() == 3);
  CGAL_assertion(data.rows() >= 3);

  // rotate points
  Vertex RT = R.transpose();
  Matrix rotated_data;
  rotated_data = data * RT;
  CGAL_assertion(rotated_data.cols() == data.cols());
  CGAL_assertion(rotated_data.rows() == data.rows());

  // AABB: take mins and maxs
  double xmin = rotated_data.col(0).minCoeff();
  double xmax = rotated_data.col(0).maxCoeff();
  double ymin = rotated_data.col(1).minCoeff();
  double ymax = rotated_data.col(1).maxCoeff();
  double zmin = rotated_data.col(2).minCoeff();
  double zmax = rotated_data.col(2).maxCoeff();

  double x_dim = abs(xmax - xmin); // abs needed?
  double y_dim = abs(ymax - ymin);
  double z_dim = abs(zmax - zmin);

  // volume
  return (x_dim * y_dim * z_dim);

}

template <typename Vertex, typename Matrix>
struct Fitness_map // -> a free function
{
  Fitness_map(Population<Vertex>& p, Matrix& points) : pop(p), points(points)
  {}

  Vertex get_best()
  {
    std::size_t count_vertices = 0;

    std::size_t simplex_id;
    std::size_t vertex_id;
    double best_fitness = std::numeric_limits<int>::max();
    for(std::size_t i = 0; i < pop.size(); ++i)
    {
      for(std::size_t j =0; j < 4; ++j)
      {
        const Vertex vertex = pop[i][j];
        //std::cout << "i= "<< i << " j=" << j<<"\n vertex= " << vertex << std::endl;
        ++count_vertices;
        //std::cout << "vertex = " << vertex << std::endl;

        const double fitness = compute_fitness(vertex, points);
        //std::cout << "fitness = " << fitness << std::endl;
        if (fitness < best_fitness)
        {
          simplex_id = i;
          vertex_id = j;
          best_fitness = fitness;
        }
      }
    }

    std::vector<Vertex> best_simplex = pop[simplex_id];
    //Matrix temp = best_simplex[vertex_id];

    return best_simplex[vertex_id];
  }


  double get_best_fitness_value(const Matrix& data)
  {
    Vertex best_mat = get_best();
    return compute_fitness(best_mat, data);
  }


  const Matrix points; // or just a reference?
  Population<Vertex> pop;
};







}} // end namespaces






#endif //CGAL_FITNESS_FUNCTION_H



