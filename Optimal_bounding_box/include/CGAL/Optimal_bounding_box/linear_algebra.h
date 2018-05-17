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

#ifndef CGAL_OPTIMAL_BOUNDING_BOX_H
#define CGAL_OPTIMAL_BOUNDING_BOX_H

#include <CGAL/assertions.h>
#include <vector>

#if defined(CGAL_EIGEN3_ENABLED)
#include <CGAL/Eigen_linear_algebra_traits.h>
#endif



namespace CGAL {
namespace Optimal_bounding_box {

#if defined(CGAL_EIGEN3_ENABLED)
typedef CGAL::Eigen_linear_algebra_traits Linear_algebra_traits;
#endif

template <typename Matrix>
const Matrix reflection(const Matrix& S_centroid, const Matrix& S_worst)
{
  CGAL_assertion(S_centroid.rows() == 3);
  CGAL_assertion(S_centroid.rows() == 3);
  CGAL_assertion(S_worst.cols() == 3);
  CGAL_assertion(S_worst.cols() == 3);

  return S_centroid * Linear_algebra_traits::transpose(S_worst) * S_centroid;
}

template <typename Matrix>
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

template <typename Matrix>
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

template <typename Matrix>
const Matrix centroid(const Matrix& S1, const Matrix& S2, const Matrix& S3)
{
  Matrix mean = (S1 + S2 + S3) / 3.0;
  Matrix Q = Linear_algebra_traits::qr_factorization(mean);
  double det = Linear_algebra_traits::determinant(Q);
  return Q / det;
}


}} // end namespaces






#endif //CGAL_OPTIMAL_BOUNDING_BOX_H


