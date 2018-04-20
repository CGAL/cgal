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
#include <Eigen/QR>


namespace CGAL {
namespace Optimal_bounding_box {


template<typename Matrix>
void qr_factorization(Matrix& A, Matrix& Q)
{
  Eigen::HouseholderQR<Matrix> qr(A);
  Q = qr.householderQ();
}

template <typename Matrix>
const Matrix reflection(Matrix& S_centroid, Matrix& S_worst)
{
  CGAL_assertion(S_centroid.rows() == 3);
  CGAL_assertion(S_centroid.rows() == 3);
  CGAL_assertion(S_worst.cols() == 3);
  CGAL_assertion(S_worst.cols() == 3);

  return S_centroid * S_worst.transpose() * S_centroid;
}

template <typename Matrix>
const Matrix expansion(Matrix& S_centroid, Matrix& S_worst, Matrix& S_reflection)
{
  CGAL_assertion(S_centroid.rows() == 3);
  CGAL_assertion(S_centroid.rows() == 3);
  CGAL_assertion(S_worst.cols() == 3);
  CGAL_assertion(S_worst.cols() == 3);
  CGAL_assertion(S_reflection.cols() == 3);
  CGAL_assertion(S_reflection.cols() == 3);

  return S_centroid * S_worst.transpose() * S_reflection;
}

template <typename Matrix>
const Matrix constraction(Matrix& S_centroid, Matrix& S_worst)
{
  // same API for reduction
  CGAL_assertion(S_centroid.rows() == 3);
  CGAL_assertion(S_centroid.rows() == 3);
  CGAL_assertion(S_worst.cols() == 3);
  CGAL_assertion(S_worst.cols() == 3);

  Matrix contract = 0.5 * S_centroid + 0.5 * S_worst;
  Matrix Q;
  qr_factorization(contract, Q);
  double det = Q.determinant();
  return Q / det;
}








}} // end namespaces






#endif //CGAL_OPTIMAL_BOUNDING_BOX_H


