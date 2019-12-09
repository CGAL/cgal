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
#ifndef CGAL_OPTIMAL_BOUNDING_BOX_HELPER_H
#define CGAL_OPTIMAL_BOUNDING_BOX_HELPER_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Bbox_3.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/number_utils.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#ifdef CGAL_EIGEN3_ENABLED
#include <CGAL/Eigen_matrix.h>
#include <Eigen/QR>
#endif

#include <fstream>
#include <string>
#include <vector>

namespace CGAL {
namespace Optimal_bounding_box {

#ifdef CGAL_EIGEN3_ENABLED
// @tmp move the header include too
template <typename K>
class Optimal_bounding_box_traits
  : public K
{
public:
  typedef typename K::FT                               FT;
  typedef CGAL::Eigen_matrix<FT, 3, 3>                 Matrix;

private:
  typedef typename Matrix::EigenType                   EigenType;

public:
  explicit Optimal_bounding_box_traits(const K& k = K()) : K(k) { }

  /// Get the transpose of a `Matrix<NT>` matrix
  Matrix transpose(const Matrix& mat) const
  {
    return Matrix(mat.eigen_object().transpose());
  }

  /// Get the determinant of a `Matrix<NT>` matrix
  FT determinant(const Matrix& matrix) const
  {
    return matrix.eigen_object().determinant();
  }

  /// Performs QR decomposition of matrix A to a unitary matrix and an upper triagonal
  /// and returns the unitary matrix.
  Matrix qr_factorization(const Matrix& A) const
  {
    Eigen::HouseholderQR<EigenType> qr(A.eigen_object());
    CGAL_assertion(CGAL::abs(determinant(Matrix(EigenType(qr.householderQ()))) - 1.) < 0.000001);

    return Matrix(EigenType(qr.householderQ()));
  }
};
#endif

} // end namespace Optimal_bounding_box
} // end namespace CGAL

#endif // CGAL_OPTIMAL_BOUNDING_BOX_HELPER_H
