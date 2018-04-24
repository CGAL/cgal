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

#ifndef CGAL_OPTIMAL_BOUNDING_BOX_OPTIMIZATION_ALGORITHMS_H
#define CGAL_OPTIMAL_BOUNDING_BOX_OPTIMIZATION_ALGORITHMS_H

#include <vector>
#include <CGAL/Random.h>
#include <boost/iterator/counting_iterator.hpp>

#include <CGAL/Optimal_bounding_box/fitness_function.h>
#include <CGAL/Optimal_bounding_box/linear_algebra.h>


namespace CGAL {

namespace Optimal_bounding_box {


template<typename Vertex>
struct Comparator
{
  Comparator(std::vector<double> in) : fitness(in) {}

  inline bool operator() (std::size_t& i, std::size_t& j)
  {
    return fitness[i] < fitness[j];
  }

  std::vector<double> fitness;
};



// points: point coords
// simplex: 4 roation matrices are its vertices
template<typename Matrix>
void nelder_mead(std::vector<Matrix>& simplex, Matrix& points, std::size_t nb_iterations)
{

  CGAL_assertion(simplex.size() == 4); // tetrahedron

  std::vector<double> fitness(4);
  std::vector<std::size_t> indices( boost::counting_iterator<std::size_t>( 0 ),
                                    boost::counting_iterator<std::size_t>( simplex.size() ) );


  for(std::size_t t = 0; t < nb_iterations; ++t)
  {

    for(int i = 0; i < 4; ++i)
    {
      fitness[i] = compute_fitness(simplex[i], points);
    }

    /*
    for(const Matrix& v : simplex)
      fitness.push_back(compute_fitness(v, points));
    */
    CGAL_assertion(fitness.size() == 4);
    CGAL_assertion(indices.size() == 4);


    // get indices of sorted sequence
    Comparator<Matrix> compare_indices(fitness);
    std::sort(indices.begin(), indices.end(), compare_indices);


    /*
    for(int i =0 ; i < 4; ++i)
    {
      std::cout << simplex[i] << "\n";
      std::cout << "fitness= " << fitness[i] << "\n";
      std::cout << "index= " << indices[i] << "\n\n";
    }
    std::cout << std::endl;
    */


    // new sorted simplex & fitness
    std::vector<Matrix> s_simplex(4);
    std::vector<double> s_fitness(4);
    for(int i = 0; i < 4; ++i)
    {
      s_simplex[i] = simplex[indices[i]];
      s_fitness[i] = fitness[indices[i]];
    }

    simplex = s_simplex; // swap?
    fitness = s_fitness;

    // centroid
    const Matrix v_centroid = centroid(simplex[0], simplex[1], simplex[2]);

    // find worst's vertex reflection
    const Matrix v_worst = simplex[3];
    const Matrix v_refl = reflection(v_centroid, v_worst);
    const double f_refl = compute_fitness(v_refl, points);

    if(f_refl < fitness[2])
    {
      if(f_refl >= fitness[0]) // if reflected point is not better than the best
      {
        // do reflection
        simplex[3] = v_refl;
      }
      else
      {
        // expansion
        const Matrix v_expand = expansion(v_centroid, v_worst, v_refl);
        const double f_expand = compute_fitness(v_expand, points);
        if(f_expand < f_refl)
          simplex[3] = v_expand;
        else
          simplex[3] = v_refl;
      }
    }
    else // if reflected vertex is not better
    {
      const Matrix v_mean = mean(v_centroid, v_worst);
      const double f_mean = compute_fitness(v_mean, points);
      if(f_mean <= fitness[3])
        // contraction of worst
        simplex[3] = v_mean;
      else
      {
        // reduction: move all vertices towards the best
        for(int i=1; i < 4; ++i)
        {
          simplex[i] = mean(simplex[i], simplex[0]); // todo: test order of addition
        }
      }
    }

    CGAL_assertion(simplex.size() == 4); // tetrahedron


  } // iterations




}










} } // end namespaces







#endif //CGAL_OPTIMAL_BOUNDING_BOX_OPTIMIZATION_ALGORITHMS_H
