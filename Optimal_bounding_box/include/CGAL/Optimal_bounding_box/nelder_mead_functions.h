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

#ifndef CGAL_OPTIMAL_BOUNDING_BOX_NEALDER_MEAD_FUNCTIONS_H
#define CGAL_OPTIMAL_BOUNDING_BOX_NEALDER_MEAD_FUNCTIONS_H

#include <CGAL/assertions.h>
#include <boost/iterator/counting_iterator.hpp>
#include <CGAL/Optimal_bounding_box/fitness_function.h>
#include <vector>

namespace CGAL {

namespace Optimal_bounding_box {

template <typename Linear_algebra_traits, typename Matrix>
const Matrix reflection(const Matrix& S_centroid, const Matrix& S_worst)
{
  CGAL_assertion(S_centroid.rows() == 3);
  CGAL_assertion(S_centroid.rows() == 3);
  CGAL_assertion(S_worst.cols() == 3);
  CGAL_assertion(S_worst.cols() == 3);

  return S_centroid * Linear_algebra_traits::transpose(S_worst) * S_centroid;
}

template <typename Linear_algebra_traits, typename Matrix>
const Matrix expansion(const Matrix& S_centroid, const Matrix& S_worst, const Matrix& S_reflection)
{
  CGAL_assertion(S_centroid.rows() == 3);
  CGAL_assertion(S_centroid.rows() == 3);
  CGAL_assertion(S_worst.cols() == 3);
  CGAL_assertion(S_worst.cols() == 3);
  CGAL_assertion(S_reflection.cols() == 3);
  CGAL_assertion(S_reflection.cols() == 3);

  return S_centroid * Linear_algebra_traits::transpose(S_worst) * S_reflection;
}

template <typename Linear_algebra_traits, typename Matrix>
Matrix mean(const Matrix& m1, const Matrix& m2)
{
  // same API for reduction
  CGAL_assertion(m1.rows() == 3);
  CGAL_assertion(m1.rows() == 3);
  CGAL_assertion(m2.cols() == 3);
  CGAL_assertion(m2.cols() == 3);

  Matrix reduction = 0.5 * m1 + 0.5 * m2;
  Matrix Q = Linear_algebra_traits::qr_factorization(reduction);
  double det = Linear_algebra_traits::determinant(Q);
  return Q / det;
}

template <typename Linear_algebra_traits, typename Matrix>
const Matrix nm_centroid(const Matrix& S1, const Matrix& S2, const Matrix& S3)
{
  Matrix mean = (S1 + S2 + S3) / 3.0;
  Matrix Q = Linear_algebra_traits::qr_factorization(mean);
  double det = Linear_algebra_traits::determinant(Q);
  return Q / det;
}

// needed in nelder mead algorithm
struct Comparator
{
  Comparator(const std::vector<double>& in) : fitness(in) {}

  inline bool operator() (std::size_t& i, std::size_t& j) {
    return fitness[i] < fitness[j];
  }

  const std::vector<double>& fitness;
};

// simplex: 4 rotation matrices are its vertices
template <typename Linear_algebra_traits>
void nelder_mead(std::vector<typename Linear_algebra_traits::Matrix3d>& simplex,
                 const typename Linear_algebra_traits::MatrixXd& point_data,
                 std::size_t nelder_mead_iterations)
{
  CGAL_assertion(simplex.size() == 4); // tetrahedron


  typedef typename Linear_algebra_traits::Matrix3d Matrix3d;

  std::vector<double> fitness(4);
  std::vector<std::size_t> indices(boost::counting_iterator<std::size_t>(0),
                                   boost::counting_iterator<std::size_t>(simplex.size()));

  for(std::size_t t = 0; t < nelder_mead_iterations; ++t)
  {
    for(std::size_t i = 0; i < 4; ++i)
    {
      fitness[i] = compute_fitness<Linear_algebra_traits>(simplex[i], point_data);
    }

    CGAL_assertion(fitness.size() == 4);
    CGAL_assertion(indices.size() == 4);

    // get indices of sorted sequence
    Comparator compare_indices(fitness);
    std::sort(indices.begin(), indices.end(), compare_indices);

    // new sorted simplex & fitness
    std::vector<Matrix3d> s_simplex(4);
    std::vector<double> s_fitness(4);
    for(int i = 0; i < 4; ++i)
    {
      s_simplex[i] = simplex[indices[i]];
      s_fitness[i] = fitness[indices[i]];
    }

    simplex = s_simplex;
    fitness = s_fitness;

    // centroid
    const Matrix3d v_centroid = nm_centroid<Linear_algebra_traits>(simplex[0], simplex[1], simplex[2]);

    // find worst's vertex reflection
    const Matrix3d v_worst = simplex[3];
    const Matrix3d v_refl = reflection<Linear_algebra_traits>(v_centroid, v_worst);
    const double f_refl = compute_fitness<Linear_algebra_traits>(v_refl, point_data);

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
        const Matrix3d v_expand = expansion<Linear_algebra_traits>(v_centroid, v_worst, v_refl);
        const double f_expand = compute_fitness<Linear_algebra_traits>(v_expand, point_data);
        if(f_expand < f_refl)
          simplex[3] = v_expand;
        else
          simplex[3] = v_refl;
      }
    }
    else // if reflected vertex is not better
    {
      const Matrix3d v_mean = mean<Linear_algebra_traits>(v_centroid, v_worst);
      const double f_mean = compute_fitness<Linear_algebra_traits>(v_mean, point_data);
      if(f_mean <= fitness[3])
        // contraction of worst
        simplex[3] = v_mean;
      else
      {
        // reduction: move all vertices towards the best
        for(std::size_t i=1; i < 4; ++i)
        {
          simplex[i] = mean<Linear_algebra_traits>(simplex[i], simplex[0]);
        }
      }
    }

    CGAL_assertion(simplex.size() == 4); // tetrahedron
  } // iterations
}

} // end namespace Optimal_bounding_box
} // end namespace CGAL

#endif
