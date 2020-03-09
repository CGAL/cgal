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
#ifndef CGAL_OPTIMAL_BOUNDING_BOX_POPULATION_H
#define CGAL_OPTIMAL_BOUNDING_BOX_POPULATION_H

#include <CGAL/license/Optimal_bounding_box.h>

#include <CGAL/assertions.h>
#include <CGAL/Random.h>

#include <utility>
#include <vector>

namespace CGAL {
namespace Optimal_bounding_box {
namespace internal {

template<typename Traits>
class Population
{
public:
  typedef typename Traits::FT                                 FT;
  typedef typename Traits::Matrix                             Matrix;
  typedef Matrix                                              Vertex;

  typedef std::array<Matrix, 4>                               Simplex;
  typedef std::vector<Simplex>                                Simplex_container;

private:
  Matrix create_vertex(CGAL::Random& rng) const
  {
    Matrix R;

    for(std::size_t i=0; i<3; ++i)
      for(std::size_t j=0; j<3; ++j)
        R.set(i, j, FT(rng.get_double()));

    return R;
  }

  // create random population
  Simplex create_simplex(CGAL::Random& rng) const
  {
    Simplex simplex;
    for(std::size_t i=0; i<4; ++i)
      simplex[i] = m_traits.get_Q(create_vertex(rng));

    return simplex;
  }

public:
  Population(const Traits& traits) : m_traits(traits) { }

  void initialize(std::size_t population_size,
                  CGAL::Random& rng)
  {
    m_pop.clear();
    m_pop.reserve(population_size);
    for(std::size_t i=0; i<population_size; ++i)
      m_pop.emplace_back(create_simplex(rng));
  }

  // Access
  std::size_t size() const { return m_pop.size(); }
  Simplex& operator[](const std::size_t i) { CGAL_assertion(i < m_pop.size()); return m_pop[i]; }
  const Simplex& operator[](const std::size_t i) const { CGAL_assertion(i < m_pop.size()); return m_pop[i]; }
  Simplex_container& simplices() { return m_pop; }

  // Debug
#ifdef CGAL_OPTIMAL_BOUNDING_BOX_DEBUG
  void show_population();
#endif

private:
  std::vector<Simplex> m_pop;
  const Traits& m_traits;
};

#ifdef CGAL_OPTIMAL_BOUNDING_BOX_DEBUG
template <typename Matrix>
void Population<Matrix>::show_population()
{
  std::size_t id = 0;
  for(const Simplex& simplex : m_pop)
  {
    std::cout << "Simplex: " << id++ << std::endl;
    for(const Matrix& R : simplex)
      std::cout << R << "\n\n";
    std:: cout << std:: endl;
  }
}
#endif

} // namespace internal
} // namespace Optimal_bounding_box
} // namespace CGAL

#endif // CGAL_OPTIMAL_BOUNDING_BOX_POPULATION_H
