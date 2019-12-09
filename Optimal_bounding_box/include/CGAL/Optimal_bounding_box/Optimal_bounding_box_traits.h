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
//
#ifndef CGAL_OPTIMAL_BOUNDING_BOX_BOUNDING_BOX_TRAITS_H
#define CGAL_OPTIMAL_BOUNDING_BOX_BOUNDING_BOX_TRAITS_H

#include <CGAL/license/Optimal_bounding_box.h>

#ifdef CGAL_EIGEN3_ENABLED
#include <CGAL/Eigen_matrix.h>
#include <Eigen/QR>
#endif

namespace CGAL {
namespace Optimal_bounding_box {

#if defined(CGAL_EIGEN3_ENABLED) || defined(DOXYGEN_RUNNING)

/// \ingroup PkgOptimalBoundingBoxClasses
///
/// The class `CGAL::Optimal_bounding_box::Optimal_bounding_box_traits` is a traits type
/// to be used with the functions `CGAL::optimal_bounding_box()`.
/// It uses the third party library \ref thirdpartyEigen "Eigen", which must therefore
/// be available on the system for this class to be used.
///
/// \tparam K must be a model of `Kernel`
///
/// \cgalModels `OptimalBoundingBoxTraits`
///
template <typename K>
class Optimal_bounding_box_traits
  : public K
{
public:
  /// The field number type
  typedef typename K::FT                               FT;

  /// The matrix type
  typedef CGAL::Eigen_matrix<FT, 3, 3>                 Matrix;

private:
  typedef typename Matrix::EigenType                   EigenType;

public:
  /// Constructor from the base kernel
  explicit Optimal_bounding_box_traits(const K& k = K()) : K(k) { }

  /// Get the transpose of a matrix
  Matrix transpose(const Matrix& mat) const
  {
    return Matrix(mat.eigen_object().transpose());
  }

  /// Get the determinant of a matrix
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
#endif // defined(CGAL_EIGEN3_ENABLED) || defined(DOXYGEN_RUNNING)

} // namespace Optimal_bounding_box
} // namespace CGAL

#endif // CGAL_OPTIMAL_BOUNDING_BOX_BOUNDING_BOX_TRAITS_H
