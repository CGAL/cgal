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

#include <CGAL/Optimal_bounding_box/internal/fitness_function.h>

#include <CGAL/assertions.h>

#include <algorithm>
#include <array>

namespace CGAL {
namespace Optimal_bounding_box {

template <typename Matrix, typename Traits>
Matrix reflection(const Matrix& S_centroid,
                  const Matrix& S_worst,
                  const Traits& traits)
{
  CGAL_assertion(S_centroid.number_of_rows() == 3 && S_centroid.number_of_columns() == 3);
  CGAL_assertion(S_worst.number_of_rows() == 3 && S_worst.number_of_columns() == 3);

  return S_centroid * traits.transpose(S_worst) * S_centroid;
}

template <typename Matrix, typename Traits>
Matrix expansion(const Matrix& S_centroid,
                 const Matrix& S_worst,
                 const Matrix& S_reflection,
                 const Traits& traits)
{
  CGAL_assertion(S_centroid.number_of_rows() == 3 && S_centroid.number_of_columns() == 3);
  CGAL_assertion(S_worst.number_of_rows() == 3 && S_worst.number_of_columns() == 3);
  CGAL_assertion(S_reflection.number_of_rows() == 3 && S_reflection.number_of_columns() == 3);

  return S_centroid * traits.transpose(S_worst) * S_reflection;
}

template <typename Matrix, typename Traits>
Matrix mean(const Matrix& m1,
            const Matrix& m2,
            const Traits& traits)
{
  // same API for reduction
  CGAL_assertion(m1.number_of_rows() == 3 && m1.number_of_columns() == 3);
  CGAL_assertion(m2.number_of_rows() == 3 && m2.number_of_columns() == 3);

  const Matrix reduction = 0.5 * m1 + 0.5 * m2;
  const Matrix Q = traits.qr_factorization(reduction);
  const typename Traits::FT det = traits.determinant(Q);

  return (1. / det) * Q;
}

template <typename Matrix, typename Traits>
const Matrix nm_centroid(const Matrix& S1,
                         const Matrix& S2,
                         const Matrix& S3,
                         const Traits& traits)
{
  const Matrix mean = (1./3.) * (S1 + S2 + S3);
  const Matrix Q = traits.qr_factorization(mean);
  const typename Traits::FT det = traits.determinant(Q);

  return (1. / det) * Q;
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
    const Matrix v_refl = reflection(v_centroid, v_worst, traits);
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
        const Matrix v_expand = expansion(v_centroid, v_worst, v_refl, traits);
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

} // end namespace Optimal_bounding_box
} // end namespace CGAL

#endif
