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

  const Matrix reduction = FT(0.5) * (m1 + m2);

  return traits.get_Q(reduction);
}

template <typename Matrix, typename Traits>
const Matrix nm_centroid(const Matrix& S1,
                         const Matrix& S2,
                         const Matrix& S3,
                         const Traits& traits)
{
  typedef typename Traits::FT                                 FT;

  const Matrix mean = (FT(1) / FT(3)) * (S1 + S2 + S3);

  return traits.get_Q(mean);
}

// It's a 3D simplex with 4 rotation matrices as vertices
template <typename Simplex, typename PointRange, typename Traits>
void nelder_mead(Simplex& simplex,
                 const PointRange& points,
                 const std::size_t nelder_mead_iterations,
                 const Traits& traits)
{
  typedef typename Traits::FT                                 FT;
  typedef typename Traits::Matrix                             Matrix;

  std::array<FT, 4> fitness;
  std::array<std::size_t, 4> indices = {{ 0, 1, 2, 3 }};

  for(std::size_t t=0; t<nelder_mead_iterations; ++t)
  {
    for(std::size_t i=0; i<4; ++i)
      fitness[i] = compute_fitness<Traits>(simplex[i], points);

    // get indices of sorted sequence
    std::sort(indices.begin(), indices.end(),
              [&fitness](const std::size_t i, const std::size_t j) -> bool
                        { return fitness[i] < fitness[j]; });

    // new sorted simplex & fitness
    Simplex s_simplex;
    std::array<FT, 4> s_fitness;
    for(int i=0; i<4; ++i)
    {
      s_simplex[i] = simplex[indices[i]];
      s_fitness[i] = fitness[indices[i]];
    }

    simplex = std::move(s_simplex);
    fitness = std::move(s_fitness);

    // centroid
    const Matrix v_centroid = nm_centroid(simplex[0], simplex[1], simplex[2], traits);

    // find worst's vertex reflection
    const Matrix& v_worst = simplex[3];
    const Matrix v_refl = reflection(v_centroid, v_worst);
    const FT f_refl = compute_fitness<Traits>(v_refl, points);

    if(f_refl < fitness[2])
    {
      if(f_refl >= fitness[0]) // if reflected point is not better than the best
      {
        // reflection
        simplex[3] = std::move(v_refl);
      }
      else
      {
        // expansion
        const Matrix v_expand = expansion(v_centroid, v_worst, v_refl);
        const FT f_expand = compute_fitness<Traits>(v_expand, points);
        if(f_expand < f_refl)
          simplex[3] = std::move(v_expand);
        else
          simplex[3] = std::move(v_refl);
      }
    }
    else // if reflected vertex is not better
    {
      const Matrix v_mean = mean(v_centroid, v_worst, traits);
      const FT f_mean = compute_fitness<Traits>(v_mean, points);
      if(f_mean <= fitness[3])
      {
        // contraction of worst
        simplex[3] = std::move(v_mean);
      }
      else
      {
        // reduction: move all vertices towards the best
        for(std::size_t i=1; i<4; ++i)
          simplex[i] = mean(simplex[i], simplex[0], traits);
      }
    }
  } // nelder mead iterations
}

} // namespace internal
} // namespace Optimal_bounding_box
} // namespace CGAL

#endif
