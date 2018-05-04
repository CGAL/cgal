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
#include <algorithm>
#include <CGAL/Random.h>
#include <boost/iterator/counting_iterator.hpp>

#include <CGAL/Optimal_bounding_box/fitness_function.h>
#include <CGAL/Optimal_bounding_box/linear_algebra.h>
#include <CGAL/Optimal_bounding_box/population.h>


namespace CGAL {

namespace Optimal_bounding_box {


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
// simplex: 4 rotation matrices are its vertices
template<typename Vertex, typename Matrix>
void nelder_mead(std::vector<Vertex>& simplex, const Matrix& points, std::size_t nb_iterations)
{

  CGAL_assertion(simplex.size() == 4); // tetrahedron

  std::vector<double> fitness(4);
  std::vector<std::size_t> indices( boost::counting_iterator<std::size_t>( 0 ),
                                    boost::counting_iterator<std::size_t>( simplex.size() ) );

  for(std::size_t t = 0; t < nb_iterations; ++t)
  {
    for(std::size_t i = 0; i < 4; ++i)
    {
      fitness[i] = compute_fitness(simplex[i], points);
    }

    CGAL_assertion(fitness.size() == 4);
    CGAL_assertion(indices.size() == 4);

    // get indices of sorted sequence
    Comparator compare_indices(fitness);
    std::sort(indices.begin(), indices.end(), compare_indices);

    // new sorted simplex & fitness
    std::vector<Vertex> s_simplex(4);
    std::vector<double> s_fitness(4);
    for(int i = 0; i < 4; ++i)
    {
      s_simplex[i] = simplex[indices[i]];
      s_fitness[i] = fitness[indices[i]];
    }

    simplex = s_simplex;
    fitness = s_fitness;

    // centroid
    const Vertex v_centroid = centroid(simplex[0], simplex[1], simplex[2]);

    // find worst's vertex reflection
    const Vertex v_worst = simplex[3];
    const Vertex v_refl = reflection(v_centroid, v_worst);
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
        const Vertex v_expand = expansion(v_centroid, v_worst, v_refl);
        const double f_expand = compute_fitness(v_expand, points);
        if(f_expand < f_refl)
          simplex[3] = v_expand;
        else
          simplex[3] = v_refl;
      }
    }
    else // if reflected vertex is not better
    {
      const Vertex v_mean = mean(v_centroid, v_worst);
      const double f_mean = compute_fitness(v_mean, points);
      if(f_mean <= fitness[3])
        // contraction of worst
        simplex[3] = v_mean;
      else
      {
        // reduction: move all vertices towards the best
        for(std::size_t i=1; i < 4; ++i)
        {
          simplex[i] = mean(simplex[i], simplex[0]);
        }
      }
    }

    CGAL_assertion(simplex.size() == 4); // tetrahedron

  } // iterations
}

struct Random_int_generator
{
  Random_int_generator(int l, int h) : low(l), high(h) {}

  int operator() ()
  {
    return random_gen.get_int(low, high);
  }

  CGAL::Random random_gen;
  int low;
  int high;
};


template <typename Simplex, typename Matrix>
void genetic_algorithm(Population<Simplex>& pop, const Matrix& points)
{
  // random permutations
  std::size_t m = pop.size();

  //groups 1,2 : size m/2  groups 3,4 : size (m - m/2).   m/2 is floored
  std::size_t size_first_group = m/2;
  std::size_t size_second_group = m - m/2;

  std::vector<std::size_t> ids1(m/2), ids2(m/2);
  std::vector<std::size_t> ids3(m - m/2), ids4(m - m/2);

  CGAL::Random rng;

  //Random_int_generator rgen(0, m);
  //std::generate(indices.begin(), indices.end(), rgen);

  std::generate(ids1.begin(), ids1.end(),
                [&rng, &m] ()
                { return rng.get_int(0, m); });

  std::generate(ids2.begin(), ids2.end(),
                [&rng, &m] ()
                { return rng.get_int(0, m); });

  std::generate(ids3.begin(), ids3.end(),
                [&rng, &m] ()
                { return rng.get_int(0, m); });

  std::generate(ids4.begin(), ids4.end(),
                [&rng, &m] ()
                { return rng.get_int(0, m); });

  Population<Simplex> group1(m/2), group2(m/2);
  Population<Simplex> group3(m - m/2), group4(m - m/2);

  for(std::size_t i = 0; i < ids1.size(); ++i)
    group1[i] = pop[ids1[i]];

  for(std::size_t i = 0; i < ids2.size(); ++i)
    group2[i] = pop[ids2[i]];

  for(std::size_t i = 0; i < ids3.size(); ++i)
    group3[i] = pop[ids3[i]];

  for(std::size_t i = 0; i < ids4.size(); ++i)
    group4[i] = pop[ids4[i]];

#ifdef OBB_DEBUG
  check_det(group1);
  check_det(group2);
  check_det(group3);
  check_det(group4);
#endif

  // crossover I
  Population<Simplex> offspringsA(size_first_group);
  double bias = 0.1;

  for(int i = 0; i < size_first_group; ++i)
  {
    std::vector<Simplex> offspring(4);
    for(int j = 0; j < 4; ++j)
    {
      double r = rng.get_double();
      double fitnessA = compute_fitness(group1[i][j], points);
      double fitnessB = compute_fitness(group2[i][j], points);
      double threshold;
      if(fitnessA < fitnessB)
        threshold = 0.5 + bias;
      else
        threshold = 0.5 - bias;
      if(r < threshold) // choose A
        offspring[j] = group1[i][j];
      else // choose B
        offspring[j] = group2[i][j];
    }
    offspringsA[i] = offspring;
  }

#ifdef OBB_DEBUG
  std::cout << "offspringsA: \n" ;
  check_det(offspringsA);
  std::cin.get();
#endif

  // crossover II
  Population<Simplex> offspringsB(size_second_group);
  bias = 0.1;

  for(int i = 0; i < size_second_group; ++i)
  {
    std::vector<Simplex> offspring(4);
    for(int j = 0; j < 4; ++j)
    {
      double fitnessA = compute_fitness(group3[i][j], points);
      double fitnessB = compute_fitness(group4[i][j], points);
      double lambda;
      if(fitnessA < fitnessB)
        lambda = 0.5 + bias;
      else
        lambda = 0.5 - bias;
      // combine information from A and B
      offspring[j] = lambda * group3[i][j] + lambda * group4[i][j];

    }

    // qr factorization of the offspring
    qr_factorization(offspring);
    offspringsB[i] = offspring;
  }

#ifdef OBB_DEBUG
  std::cout << "offspringsB: \n" ;
  check_det(offspringsB);
#endif

  CGAL_assertion(offspringsA.size() == size_first_group);
  CGAL_assertion(offspringsB.size() == size_second_group);
  CGAL_assertion(offspringsA.size() + offspringsB.size() == pop.size());

  // next generatrion
  for(std::size_t i = 0; i < size_first_group; ++i)
  {
    pop[i] = offspringsA[i];
  }

  for(std::size_t i = 0; i < size_second_group; ++i)
  {
    pop[size_first_group + i] = offspringsB[i];
  }

}

#ifdef OBB_DEBUG
template <typename Simplex>
void check_det(Population<Simplex>& pop)
{
  for(int i = 0; i < pop.size(); ++i)
  {
    for(int j = 0; j < 4; ++j)
    {
      auto A = pop[i][j]; // Simplex
      std::cout << A.determinant() << std::endl;
    }
  }
}
#endif

} } // end namespaces







#endif //CGAL_OPTIMAL_BOUNDING_BOX_OPTIMIZATION_ALGORITHMS_H
