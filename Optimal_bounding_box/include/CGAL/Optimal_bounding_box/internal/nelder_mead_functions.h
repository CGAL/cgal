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
#ifndef CGAL_OPTIMAL_BOUNDING_BOX_NEALDER_MEAD_FUNCTIONS_H
#define CGAL_OPTIMAL_BOUNDING_BOX_NEALDER_MEAD_FUNCTIONS_H

#include <CGAL/license/Optimal_bounding_box.h>

#include <CGAL/Optimal_bounding_box/internal/fitness_function.h>
#include <CGAL/Optimal_bounding_box/internal/helper.h>

#include <CGAL/assertions.h>

#include <algorithm>
#include <array>

namespace CGAL {
namespace Optimal_bounding_box {
namespace internal {

template <typename Matrix>
Matrix reflection(const Matrix& S_centroid,
                  const Matrix& S_worst)
{
  return S_centroid * transpose(S_worst) * S_centroid;
}

template <typename Matrix>
Matrix expansion(const Matrix& S_centroid,
                 const Matrix& S_worst,
                 const Matrix& S_reflection)
{
  return S_centroid * transpose(S_worst) * S_reflection;
}

template <typename Matrix, typename Traits>
Matrix mean(const Matrix& m1,
            const Matrix& m2,
            const Traits& traits)
{
  typedef typename Traits::FT                                 FT;

  constexpr half = FT(1) / FT(2);
  const Matrix reduction = half * (m1 + m2);

  return traits.get_Q(reduction);
}

template <typename Matrix, typename Traits>
const Matrix nm_centroid(const Matrix& S1,
                         const Matrix& S2,
                         const Matrix& S3,
                         const Traits& traits)
{
  typedef typename Traits::FT                                 FT;

  constexpr FT third = FT(1) / FT(3);
  const Matrix mean = third * (S1 + S2 + S3);

  return traits.get_Q(mean);
}

// It's a 3D simplex with 4 rotation matrices as vertices
template <typename Simplex, typename PointRange, typename Traits>
void nelder_mead(Simplex& simplex,
                 const std::size_t nelder_mead_iterations,
                 const PointRange& points,
                 const Traits& traits)
{
  typedef typename Simplex::value_type                        Vertex;
  typedef typename Traits::FT                                 FT;
  typedef typename Traits::Matrix                             Matrix;

  for(std::size_t t=0; t<nelder_mead_iterations; ++t)
  {
    std::sort(simplex.begin(), simplex.end(),
              [](const Vertex& vi, const Vertex& vj) -> bool
              { return vi.fitness_value() < vj.fitness_value(); });

    // centroid
    const Matrix centroid_m = nm_centroid(simplex[0].matrix(),
                                          simplex[1].matrix(),
                                          simplex[2].matrix(), traits);

    // find worst's vertex reflection
    const Matrix& worst_m = simplex[3].matrix();
    const Matrix refl_m = reflection(centroid_m, worst_m);
    const FT refl_f = compute_fitness(refl_m, points, traits);

    // if reflected point is better than the second worst
    if(refl_f < simplex[2].fitness_value())
    {
      // if reflected point is not better than the best
      if(refl_f >= simplex[0].fitness_value())
      {
        // reflection
        simplex[3] = Vertex(refl_m, refl_f);
      }
      else
      {
        // expansion
        const Matrix expand_m = expansion(centroid_m, worst_m, refl_m);
        const FT expand_f = compute_fitness(expand_m, points, traits);
        if(expand_f < refl_f)
          simplex[3] = Vertex(expand_m, expand_f);
        else
          simplex[3] = Vertex(refl_m, refl_f);
      }
    }
    else // reflected vertex is worse
    {
      const Matrix mean_m = mean(centroid_m, worst_m, traits);
      const FT mean_f = compute_fitness(mean_m, points, traits);

      if(mean_f <= simplex[3].fitness_value())
      {
        // contraction of worst
        simplex[3] = Vertex(mean_m, mean_f);
      }
      else
      {
        // reduction: move all vertices towards the best
        for(std::size_t i=1; i<4; ++i)
          simplex[i] = Vertex(mean(simplex[i].matrix(), simplex[0].matrix(), traits), points, traits);
      }
    }
  } // nelder mead iterations
}

} // namespace internal
} // namespace Optimal_bounding_box
} // namespace CGAL

#endif
