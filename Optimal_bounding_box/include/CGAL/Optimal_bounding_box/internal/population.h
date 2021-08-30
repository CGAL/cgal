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
#ifndef CGAL_OPTIMAL_BOUNDING_BOX_POPULATION_H
#define CGAL_OPTIMAL_BOUNDING_BOX_POPULATION_H

#include <CGAL/license/Optimal_bounding_box.h>

#include <CGAL/Optimal_bounding_box/internal/fitness_function.h>

#include <CGAL/assertions.h>
#include <CGAL/Random.h>

#include <utility>
#include <vector>

namespace CGAL {
namespace Optimal_bounding_box {
namespace internal {

template<typename Traits>
struct Vertex_with_fitness
{
  typedef typename Traits::FT                                 FT;
  typedef typename Traits::Matrix                             Matrix;

  Vertex_with_fitness() { CGAL_assertion_code(m_is_val_initialized = false;) }
  Vertex_with_fitness(const Matrix m, const FT v) : m_mat(std::move(m)), m_val(v)
  {
    CGAL_assertion_code(m_is_val_initialized = true;)
  }

  template <typename PointRange>
  Vertex_with_fitness(const Matrix m,
                      const PointRange& points,
                      const Traits& traits)
    :
      m_mat(std::move(m))
  {
    m_val = compute_fitness(m, points, traits);
    CGAL_assertion_code(m_is_val_initialized = true;)
  }

  Matrix& matrix() { return m_mat; }
  const Matrix& matrix() const { return m_mat; }
  FT& fitness() { return m_val; }
  FT fitness() const { CGAL_assertion(m_is_val_initialized); return m_val; }

private:
  Matrix m_mat;
  FT m_val;
  CGAL_assertion_code(bool m_is_val_initialized;)
};

template<typename Traits>
class Population
{
public:
  typedef typename Traits::FT                                 FT;
  typedef typename Traits::Matrix                             Matrix;

  typedef Vertex_with_fitness<Traits>                         Vertex;
  typedef std::array<Vertex, 4>                               Simplex;
  typedef std::vector<Simplex>                                Simplex_container;

public:
  Population(const Traits& traits) : m_traits(traits) { }

  // Access
  std::size_t size() const { return m_simplices.size(); }
  Simplex& operator[](const std::size_t i) { CGAL_assertion(i < m_simplices.size()); return m_simplices[i]; }
  const Simplex& operator[](const std::size_t i) const { CGAL_assertion(i < m_simplices.size()); return m_simplices[i]; }
  Simplex_container& simplices() { return m_simplices; }

private:
  Matrix create_random_matrix(CGAL::Random& rng) const
  {
    Matrix m;

    for(std::size_t i=0; i<3; ++i)
      for(std::size_t j=0; j<3; ++j)
        m.set(i, j, FT(rng.get_double()));

    return m;
  }

public:
  template <typename PointRange>
  Simplex create_simplex(const PointRange& points,
                         CGAL::Random& rng) const
  {
    Simplex s;
    for(std::size_t i=0; i<4; ++i)
      s[i] = Vertex{m_traits.get_Q(create_random_matrix(rng)), points, m_traits};

    return s;
  }

  // create random population
  template <typename PointRange>
  void initialize(const std::size_t population_size,
                  const PointRange& points,
                  CGAL::Random& rng)
  {
    m_simplices.clear();
    m_simplices.reserve(population_size);
    for(std::size_t i=0; i<population_size; ++i)
      m_simplices.emplace_back(create_simplex(points, rng));
  }

  Vertex& get_best_vertex()
  {
    std::size_t simplex_id = static_cast<std::size_t>(-1), vertex_id = static_cast<std::size_t>(-1);
    FT best_fitness = FT{(std::numeric_limits<double>::max)()};
    for(std::size_t i=0, ps=m_simplices.size(); i<ps; ++i)
    {
      for(std::size_t j=0; j<4; ++j)
      {
        const Vertex& vertex = m_simplices[i][j];
        const FT fitness = vertex.fitness();
        if(fitness < best_fitness)
        {
          simplex_id = i;
          vertex_id = j;
          best_fitness = fitness;
        }
      }
    }

    return m_simplices[simplex_id][vertex_id];
  }

  // Debug
#ifdef CGAL_OPTIMAL_BOUNDING_BOX_DEBUG
  void show_population() const
  {
    std::size_t id = 0;
    for(const Simplex& s : m_simplices)
    {
      std::cout << "Simplex: " << id++ << std::endl;
      for(const Matrix& m : s)
        std::cout << m << "\n\n";
      std:: cout << std:: endl;
    }
  }
#endif

private:
  std::vector<Simplex> m_simplices;

  const Traits& m_traits;
};

} // namespace internal
} // namespace Optimal_bounding_box
} // namespace CGAL

#endif // CGAL_OPTIMAL_BOUNDING_BOX_POPULATION_H
