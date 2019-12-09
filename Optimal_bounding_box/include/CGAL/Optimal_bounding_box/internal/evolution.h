// Copyright (c) 2018-2019 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Konstantinos Katrioplas
//                 Mael Rouxel-Labbé
//
#ifndef CGAL_OPTIMAL_BOUNDING_BOX_EVOLUTION_H
#define CGAL_OPTIMAL_BOUNDING_BOX_EVOLUTION_H

#include <CGAL/license/Optimal_bounding_box.h>

#include <CGAL/Optimal_bounding_box/internal/fitness_function.h>
#include <CGAL/Optimal_bounding_box/internal/nelder_mead_functions.h>
#include <CGAL/Optimal_bounding_box/internal/population.h>

#include <CGAL/Random.h>
#include <CGAL/number_utils.h>

#include <algorithm>
#include <array>
#include <iostream>
#include <vector>

namespace CGAL {
namespace Optimal_bounding_box {
namespace internal {

template <typename PointRange, typename Traits>
class Evolution
{
public:
  typedef typename Traits::FT                                               FT;
  typedef typename Traits::Matrix                                           Matrix;

  typedef Population<Traits>                                                Population;
  typedef typename Population::Simplex                                      Simplex;

  typedef Fitness_map<Population, PointRange, Traits>                       Fitness_map;

  Evolution(const PointRange& points,
            CGAL::Random& rng,
            const Traits& traits)
    :
      m_population(traits),
      m_rng(rng), // @todo just a parameter of genetic_algorithm() ?
      m_points(points),
      m_traits(traits)
  { }

  void genetic_algorithm(const std::size_t population_size = 50)
  {
    // random permutations
    m_population.initialize(population_size, m_rng);

    //groups 1,2 : size m/2  groups 3,4 : size (m - m/2).   m/2 is floored
    const std::size_t m = population_size;
    const std::size_t first_group_size = m / 2;
    const std::size_t second_group_size = m - first_group_size;

    std::vector<std::size_t> group1(first_group_size), group2(first_group_size);
    std::vector<std::size_t> group3(second_group_size), group4(second_group_size);

    int im = static_cast<int>(m);
    std::generate(group1.begin(), group1.end(), [&]{ return m_rng.get_int(0, im); });
    std::generate(group2.begin(), group2.end(), [&]{ return m_rng.get_int(0, im); });
    std::generate(group3.begin(), group3.end(), [&]{ return m_rng.get_int(0, im); });
    std::generate(group4.begin(), group4.end(), [&]{ return m_rng.get_int(0, im); });

    // crossover I
    double bias = 0.1;

    std::vector<Simplex> new_simplices(m);

    for(std::size_t i=0; i<first_group_size; ++i)
    {
      std::array<Matrix, 4> offspring;
      for(int j=0; j<4; ++j)
      {
        const double r = m_rng.get_double();
        const double fitnessA = compute_fitness<Traits>(m_population[group1[i]][j], m_points);
        const double fitnessB = compute_fitness<Traits>(m_population[group2[i]][j], m_points);
        const double threshold = (fitnessA < fitnessB) ? (0.5 + bias) : (0.5 - bias);

        if(r < threshold)
          offspring[j] = m_population[group1[i]][j];
        else
          offspring[j] = m_population[group2[i]][j];
      }

      new_simplices[i] = std::move(offspring);
    }

    // crossover II
    bias = 0.1; // @fixme should the bias change? What should be the initial value?

    for(std::size_t i=0; i<second_group_size; ++i)
    {
      std::array<Matrix, 4> offspring;
      for(int j=0; j<4; ++j)
      {
        const double fitnessA = compute_fitness<Traits>(m_population[group3[i]][j], m_points);
        const double fitnessB = compute_fitness<Traits>(m_population[group4[i]][j], m_points);
        const double lambda = (fitnessA < fitnessB) ? (0.5 + bias) : (0.5 - bias);
        const double rambda = 1 - lambda; // the 'l' in 'lambda' stands for left

        // combine information from A and B
        Matrix new_vertex(3, 3);

        const Matrix& lm = m_population[group3[i]][j];
        const Matrix& rm = m_population[group4[i]][j];

        // just avoiding having to add matrix sums and scalar multiplications to the concept
        new_vertex.set(0, 0, lambda * lm(0, 0) + rambda * rm(0, 0));
        new_vertex.set(0, 1, lambda * lm(0, 1) + rambda * rm(0, 1));
        new_vertex.set(0, 2, lambda * lm(0, 2) + rambda * rm(0, 2));

        new_vertex.set(1, 0, lambda * lm(1, 0) + rambda * rm(1, 0));
        new_vertex.set(1, 1, lambda * lm(1, 1) + rambda * rm(1, 1));
        new_vertex.set(1, 2, lambda * lm(1, 2) + rambda * rm(1, 2));

        new_vertex.set(2, 0, lambda * lm(2, 0) + rambda * rm(2, 0));
        new_vertex.set(2, 1, lambda * lm(2, 1) + rambda * rm(2, 1));
        new_vertex.set(2, 2, lambda * lm(2, 2) + rambda * rm(2, 2));

        offspring[j] = m_traits.qr_factorization(new_vertex);
      }

      new_simplices[first_group_size + i] = std::move(offspring);
    }

    m_population.simplices() = std::move(new_simplices);
  }

  void evolve(const std::size_t generations)
  {
    // hardcoded population size
    const std::size_t population_size = 50;

    // hardcoded nelder_mead_iterations
    const std::size_t nelder_mead_iterations = 20;

    // stopping criteria prameters
    double prev_fit_value = 0.;
    double new_fit_value = 0.;
    const double tolerance = 1e-2;
    int stale = 0;

    for(std::size_t t=0; t<generations; ++t)
    {
#ifdef CGAL_OPTIMAL_BOUNDING_BOX_DEBUG
      std::cout << "generation = " << t << "\n";
#endif

      genetic_algorithm(population_size);

#ifdef CGAL_OPTIMAL_BOUNDING_BOX_DEBUG
      //std::cout << "pop after genetic" << std::endl;
      //pop.show_population();
      //std::cout << std::endl;
#endif

      for(std::size_t s=0; s<population_size; ++s)
        nelder_mead(m_population[s], m_points, nelder_mead_iterations, m_traits);

      // stopping criteria
      Fitness_map fitness_map(m_population, m_points);

#ifdef CGAL_OPTIMAL_BOUNDING_BOX_DEBUG
      //std::cout << "pop after nelder mead: " << std::endl;
      //pop.show_population();
      //std::cout << std::endl;
      const Matrix& R_now = fitness_map.get_best();
      std::cout << "det = " << m_traits.determinant(R_now) << std::endl;
#endif

      new_fit_value = fitness_map.get_best_fitness_value();
      const double difference = new_fit_value - prev_fit_value;

      if(CGAL::abs(difference) < tolerance * new_fit_value)
        ++stale;

      if(stale == 5)
        break;

      prev_fit_value = new_fit_value;
    }
  }

  const Matrix& get_best()
  {
    Fitness_map fitness_map(m_population, m_points);
    return fitness_map.get_best();
  }

private:
  Population m_population;
  CGAL::Random m_rng;
  const PointRange& m_points;
  const Traits& m_traits;
};

} // namespace internal
} // namespace Optimal_bounding_box
} // namespace CGAL

#endif // CGAL_OPTIMAL_BOUNDING_BOX_EVOLUTION_H
