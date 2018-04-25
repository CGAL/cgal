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

#ifndef CGAL_OPTIMAL_BOUNDING_BOX_POPULATION_H
#define CGAL_OPTIMAL_BOUNDING_BOX_POPULATION_H

#include <vector>
#include <CGAL/Random.h>

#include <CGAL/Optimal_bounding_box/fitness_function.h>
#include <CGAL/Optimal_bounding_box/linear_algebra.h>


namespace CGAL {

namespace Optimal_bounding_box {



template<typename Matrix>
class Population
{

  typedef std::vector<Matrix> Individual;

public:
  Population(std::size_t size) : n(size), random_generator(CGAL::Random(1))
  {
    //std::cout << "seed= " << random_generator.get_seed() << std::endl;

    // reserve pop space
    pop.reserve(n);

    // create individuals
    for(std::size_t i = 0 ; i < n; ++i)
    {
      Individual member(4);
      create_individual(member);
      CGAL_assertion(member.size() == 4);
      pop.push_back(member);
    }
  }

  void show_population()
  {
    std::size_t id = 0;
    for(const Individual i : pop)
    {
      CGAL_assertion(i.size() == 4);
      std:: cout << "Individual: "<< id++ << std::endl;
      for(const Matrix R : i)
      {
        std::cout << R; // eigen out
        std::cout << "\n\n";
      }

      std:: cout << std:: endl;
    }
  }

  std::size_t size()
  {
    return n;
  }


  // access individual
  Individual& operator[](std::size_t i)
  {
    CGAL_assertion(i < n);
    return pop[i];
  }

  const Individual& operator[](std::size_t i) const
  {
    CGAL_assertion(i < n);
    return pop[i];
  }


private:

  // create random population
  void create_individual(Individual& member)
  {
    CGAL_assertion(member.size() == 4);

    for(std::size_t i = 0; i < 4; ++i)
    {
      // create matrix
      Matrix R(3,3);
      create_simplex(R);
      Matrix Q;
      qr_factorization(R, Q);
      member[i] = Q;
    }

    CGAL_assertion(member.size() == 4);
  }

  void create_simplex(Matrix& R)
  {
    CGAL_assertion(R.rows() == 3);
    CGAL_assertion(R.cols() == 3);

    for(std::size_t i = 0; i < 3; ++i)
    {
      for(std::size_t j = 0; j < 3; ++j)
      {
        R(i, j) = random_generator.get_double();
      }
    }
  }




  CGAL::Random random_generator;
  std::size_t n;
  std::vector<Individual> pop;

};



} } // end namespaces







#endif //CGAL_OPTIMAL_BOUNDING_BOX_POPULATION_H
