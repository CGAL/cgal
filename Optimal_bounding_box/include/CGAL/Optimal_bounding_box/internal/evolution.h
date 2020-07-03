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
// Author(s)     : Mael Rouxel-Labbé
//                 Konstantinos Katrioplas
//
#ifndef CGAL_OPTIMAL_BOUNDING_BOX_EVOLUTION_H
#define CGAL_OPTIMAL_BOUNDING_BOX_EVOLUTION_H

#include <CGAL/license/Optimal_bounding_box.h>

#include <CGAL/Optimal_bounding_box/internal/fitness_function.h>
#include <CGAL/Optimal_bounding_box/internal/helper.h>
#include <CGAL/Optimal_bounding_box/internal/nelder_mead_functions.h>
#include <CGAL/Optimal_bounding_box/internal/optimize_2.h>
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

  typedef internal::Population<Traits>                                      Population;
  typedef typename Population::Simplex                                      Simplex;
  typedef typename Population::Vertex                                       Vertex;

  Evolution(const PointRange& points,
            CGAL::Random& rng,
            const Traits& traits)
    :
      m_best_v(nullptr),
      m_population(traits),
      m_rng(rng),
      m_points(points),
      m_traits(traits)
  { }

  void genetic_algorithm()
  {
    // This evolves an existing population
    CGAL_precondition(m_population.size() != 0);

    //groups 1,2 : size = floor(m/2) groups 3,4 : size = ceil(m/2).
    const std::size_t m = m_population.size();
    const std::size_t first_group_size = m / 2;
    const std::size_t second_group_size = m - first_group_size;

    std::vector<std::size_t> group1(first_group_size), group2(first_group_size);
    std::vector<std::size_t> group3(second_group_size), group4(second_group_size);

    int im = static_cast<int>(m);
    std::generate(group1.begin(), group1.end(), [&]{ return m_rng.get_int(0, im); });
    std::generate(group2.begin(), group2.end(), [&]{ return m_rng.get_int(0, im); });
    std::generate(group3.begin(), group3.end(), [&]{ return m_rng.get_int(0, im); });
    std::generate(group4.begin(), group4.end(), [&]{ return m_rng.get_int(0, im); });

    // crossover I, pick A or B
    const FT lweight = 0.4, uweight = 0.6;

    std::vector<Simplex> new_simplices(m);

    for(std::size_t i=0; i<first_group_size; ++i)
    {
      Simplex offspring;
      for(int j=0; j<4; ++j)
      {
        const FT r{m_rng.get_double()};
        const FT fitnessA = m_population[group1[i]][j].fitness();
        const FT fitnessB = m_population[group2[i]][j].fitness();
        const FT threshold = (fitnessA < fitnessB) ? uweight : lweight;

        if(r < threshold)
          offspring[j] = m_population[group1[i]][j];
        else
          offspring[j] = m_population[group2[i]][j];
      }

      new_simplices[i] = std::move(offspring);
    }

    // crossover II, combine information from A and B
    for(std::size_t i=0; i<second_group_size; ++i)
    {
      Simplex offspring;
      for(int j=0; j<4; ++j)
      {
        const FT fitnessA = m_population[group3[i]][j].fitness();
        const FT fitnessB = m_population[group4[i]][j].fitness();
        const FT lambda = (fitnessA < fitnessB) ? uweight : lweight;
        const FT rambda = FT(1) - lambda; // because the 'l' in 'lambda' stands for left

        const Matrix& lm = m_population[group3[i]][j].matrix();
        const Matrix& rm = m_population[group4[i]][j].matrix();

        offspring[j] = Vertex{m_traits.get_Q(lambda*lm + rambda*rm), m_points, m_traits};
      }

      new_simplices[first_group_size + i] = std::move(offspring);
    }

    m_population.simplices() = std::move(new_simplices);
  }

  // @todo re-enable random mutations as it is -on theory- useful, but don't allow them to override
  // the best (current) candidates
  void evolve(const std::size_t max_generations,
              const std::size_t population_size,
              const std::size_t nelder_mead_iterations,
              const std::size_t max_random_mutations = 0)
  {
    // stopping criteria prameters
    FT prev_fit_value = 0;
    const FT tolerance = 1e-10;
    int stale = 0;

    m_population.initialize(population_size, m_points, m_rng);

    std::size_t gen_iter = 0;
    for(;;)
    {
#ifdef CGAL_OPTIMAL_BOUNDING_BOX_DEBUG_PP
      std::cout << "- - - - generation #" << gen_iter << "\n";
#endif

      genetic_algorithm();

#ifdef CGAL_OPTIMAL_BOUNDING_BOX_DEBUG_PP
      std::cout << "population after genetic" << std::endl;
      pop.show_population();
      std::cout << std::endl;
#endif

      for(std::size_t s=0; s<population_size; ++s)
        nelder_mead(m_population[s], nelder_mead_iterations, m_points, m_traits);

#ifdef CGAL_OPTIMAL_BOUNDING_BOX_DEBUG_PP
      std::cout << "population after nelder mead: " << std::endl;
      pop.show_population();
      std::cout << std::endl;
#endif

      m_best_v = &(m_population.get_best_vertex());
      Matrix& best_m = m_best_v->matrix();

#ifdef CGAL_OPTIMAL_BOUNDING_BOX_DEBUG_PP
      std::cout << "new best matrix: " << std::endl << best_m << std::endl;
      std::cout << "fitness: " << m_best_v->fitness() << std::endl;
#endif

      // optimize the current best rotation by using the exact OBB 2D algorithm
      // along the axes of the current best OBB
      Optimizer_along_axes<Traits> optimizer_2D;
      optimizer_2D(best_m, m_points, m_traits);
      m_best_v->fitness() = compute_fitness(best_m, m_points, m_traits);

      // stopping criteria
      const FT new_fit_value = m_best_v->fitness();
      const FT difference = new_fit_value - prev_fit_value;

#ifdef CGAL_OPTIMAL_BOUNDING_BOX_DEBUG_PP
      std::cout << "post 2D optimization matrix: " << std::endl << best_m << std::endl;
      std::cout << "new fit value: " << new_fit_value << std::endl;
      std::cout << "difference: " << difference << std::endl;
#endif

      if(CGAL::abs(difference) < tolerance * new_fit_value)
        ++stale;

      if(stale == 5 || gen_iter++ >= max_generations)
        break;

      prev_fit_value = new_fit_value;

      // random mutations, swap #random_mutations random simplices with a new, random simplex
      if(max_random_mutations <= 0)
        continue;

      CGAL_warning(max_random_mutations <= population_size);
      const int random_mutations = m_rng.get_int(0, static_cast<int>(max_random_mutations+1));
      for(int i=0; i<random_mutations; ++i)
      {
        const int random_pos = m_rng.get_int(0, static_cast<int>(population_size));
        m_population[random_pos] = m_population.create_simplex(m_points, m_rng);
      }
    }
  }

  const Vertex& get_best_vertex() const
  {
    CGAL_assertion(m_best_v != nullptr);
    return *m_best_v;
  }

private:
  Vertex* m_best_v;
  Population m_population;

  CGAL::Random& m_rng;
  const PointRange& m_points;
  const Traits& m_traits;
};

} // namespace internal
} // namespace Optimal_bounding_box
} // namespace CGAL

#endif // CGAL_OPTIMAL_BOUNDING_BOX_EVOLUTION_H
