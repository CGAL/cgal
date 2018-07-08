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
// Author(s)     : Konstantinos Katrioplas

#ifndef CGAL_OPTIMAL_BOUNDING_BOX_EVOLUTION_H
#define CGAL_OPTIMAL_BOUNDING_BOX_EVOLUTION_H

#include <CGAL/Optimal_bounding_box/population.h>
#include <CGAL/Optimal_bounding_box/nelder_mead_functions.h>
#include <CGAL/Optimal_bounding_box/fitness_function.h>

#include <CGAL/Random.h>
#include <CGAL/number_utils.h>

#include <algorithm>
#include <iostream>
#include <vector>

namespace CGAL {

namespace Optimal_bounding_box {

template <typename Linear_algebra_traits>
class Evolution
{
  typedef typename Linear_algebra_traits::MatrixXd MatrixXd;
  typedef typename Linear_algebra_traits::Matrix3d Matrix3d;
  typedef typename Linear_algebra_traits::Vector3d Vector3d;

public:
  Evolution(Population<Linear_algebra_traits>& pop, MatrixXd& points)
    : population(pop), point_data(points)
  {}

  void genetic_algorithm()
  {
    // random permutations
    std::size_t m = population.size();

    //groups 1,2 : size m/2  groups 3,4 : size (m - m/2).   m/2 is floored
    std::size_t size_first_group = m/2;
    std::size_t size_second_group = m - m/2;

    std::vector<std::size_t> ids1(m/2), ids2(m/2);
    std::vector<std::size_t> ids3(m - m/2), ids4(m - m/2);

    CGAL::Random rng;
    int im = static_cast<int>(m);
    std::generate(ids1.begin(), ids1.end(),
                  [&rng, &im] () { return rng.get_int(0, im); });
    std::generate(ids2.begin(), ids2.end(),
                  [&rng, &im] () { return rng.get_int(0, im); });
    std::generate(ids3.begin(), ids3.end(),
                  [&rng, &im] () { return rng.get_int(0, im); });
    std::generate(ids4.begin(), ids4.end(),
                  [&rng, &im] () { return rng.get_int(0, im); });

    Population<Linear_algebra_traits> group1(m/2), group2(m/2);
    Population<Linear_algebra_traits> group3(m - m/2), group4(m - m/2);

    for(std::size_t i = 0; i < ids1.size(); ++i)
      group1[i] = population[ids1[i]];

    for(std::size_t i = 0; i < ids2.size(); ++i)
      group2[i] = population[ids2[i]];

    for(std::size_t i = 0; i < ids3.size(); ++i)
      group3[i] = population[ids3[i]];

    for(std::size_t i = 0; i < ids4.size(); ++i)
      group4[i] = population[ids4[i]];

#ifdef CGAL_OPTIMAL_BOUNDING_BOX_DEBUG
    check_det(group1);
    check_det(group2);
    check_det(group3);
    check_det(group4);
#endif

    // crossover I
    Population<Linear_algebra_traits> offspringsA(size_first_group);
    double bias = 0.1;

    for(std::size_t i = 0; i < size_first_group; ++i)
    {
      std::vector<Matrix3d> offspring(4);
      for(int j = 0; j < 4; ++j)
      {
        double r = rng.get_double();
        double fitnessA = compute_fitness<Linear_algebra_traits>(group1[i][j], point_data);
        double fitnessB = compute_fitness<Linear_algebra_traits>(group2[i][j], point_data);
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

#ifdef CGAL_OPTIMAL_BOUNDING_BOX_DEBUG
    std::cout << "offspringsA: \n" ;
    check_det(offspringsA);
#endif

    // crossover II
    Population<Linear_algebra_traits> offspringsB(size_second_group);
    bias = 0.1;

    for(std::size_t i = 0; i < size_second_group; ++i)
    {
      std::vector<Matrix3d> offspring(4);
      for(int j = 0; j < 4; ++j)
      {
        double fitnessA = compute_fitness<Linear_algebra_traits>(group3[i][j], point_data);
        double fitnessB = compute_fitness<Linear_algebra_traits>(group4[i][j], point_data);
        double lambda;
        if(fitnessA < fitnessB)
          lambda = 0.5 + bias;
        else
          lambda = 0.5 - bias;
        // combine information from A and B
        offspring[j] = lambda * group3[i][j] + lambda * group4[i][j];
      }

      // qr factorization of the offspring
      Linear_algebra_traits::qr_factorization(offspring);
      offspringsB[i] = offspring;
    }

#ifdef CGAL_OPTIMAL_BOUNDING_BOX_DEBUG
    std::cout << "offspringsB: \n" ;
    check_det(offspringsB);
#endif

    CGAL_assertion(offspringsA.size() == size_first_group);
    CGAL_assertion(offspringsB.size() == size_second_group);
    CGAL_assertion(offspringsA.size() + offspringsB.size() == population.size());

    // next generatrion
    for(std::size_t i = 0; i < size_first_group; ++i)
      population[i] = offspringsA[i];

    for(std::size_t i = 0; i < size_second_group; ++i)
      population[size_first_group + i] = offspringsB[i];
  }

  void evolve(std::size_t generations)
  {
    // hardcoded nelder_mead_iterations
    std::size_t nelder_mead_iterations = 20;

    // stopping criteria prameters
    double prev_fit_value = 0;
    double new_fit_value = 0;
    double tolerance = 1e-2;
    int stale = 0;

    for(std::size_t t = 0; t < generations; ++t)
    {
#ifdef CGAL_OPTIMAL_BOUNDING_BOX_DEBUG
      std::cout << "generation= " << t << "\n";
#endif

      genetic_algorithm();

#ifdef CGAL_OPTIMAL_BOUNDING_BOX_DEBUG
      //std::cout << "pop after genetic" << std::endl;
      //pop.show_population();
      //std::cout << std::endl;
#endif

      for(std::size_t s = 0; s < population.size(); ++s)
        nelder_mead<Linear_algebra_traits>(population[s], point_data, nelder_mead_iterations);

#ifdef CGAL_OPTIMAL_BOUNDING_BOX_DEBUG
      //std::cout << "pop after nelder mead: " << std::endl;
      //pop.show_population();
      //std::cout << std::endl;
      Fitness_map<Linear_algebra_traits, Matrix3d, MatrixXd> fitness_map_debug(population, point_data);
      Matrix3d R_now = fitness_map_debug.get_best();
      std::cout << "det= " << Linear_algebra_traits::determinant(R_now) << std::endl;
#endif

      // stopping criteria
      Fitness_map<Linear_algebra_traits, Matrix3d, MatrixXd> fitness_map(population, point_data);
      new_fit_value = fitness_map.get_best_fitness_value();
      double difference = new_fit_value - prev_fit_value;

      if(CGAL::abs(difference) < tolerance * new_fit_value)
        stale++;

      if(stale == 5)
        break;

      prev_fit_value = new_fit_value;
    }
  }

  const Matrix3d get_best()
  {
    Fitness_map<Linear_algebra_traits, Matrix3d, MatrixXd> fitness_map(population, point_data);
    return fitness_map.get_best();
  }

private:
  // data
  Population<Linear_algebra_traits> population;
  MatrixXd point_data;
};

#ifdef CGAL_OPTIMAL_BOUNDING_BOX_DEBUG
template <typename Simplex>
void check_det(Population<Simplex>& pop)
{
  for(std::size_t i = 0; i < pop.size(); ++i)
  {
    for(std::size_t j = 0; j < 4; ++j)
    {
      auto A = pop[i][j]; // Simplex
      std::cout << Linear_algebra_traits::determinant(A) << std::endl;
    }
  }
}
#endif

} // end namespace Optimal_bounding_box
} // end namespace CGAL

#endif // CGAL_OPTIMAL_BOUNDING_BOX_EVOLUTION_H
