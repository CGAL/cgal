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

#if defined(CGAL_EIGEN3_ENABLED) || defined(DOXYGEN_RUNNING)

/// \ingroup PkgOptimalBoundingBoxClasses
///
/// The class `CGAL::Oriented_bounding_box_traits` is a traits type
/// to be used with the functions `CGAL::oriented_bounding_box()`.
/// It uses the third party library \ref thirdpartyEigen "Eigen", which must therefore
/// be available on the system for this class to be used.
///
/// \tparam K must be a model of `Kernel`
///
/// \cgalModels `OrientedBoundingBoxTraits`
///
template <typename K>
class Oriented_bounding_box_traits
{
public:
  /// The field number type
  typedef typename K::FT                               FT;

  /// The point type
  typedef typename K::Point_3                          Point_3;

  /// The affine transformation type
  typedef typename K::Aff_transformation_3             Aff_transformation_3;

  /// The axis-aligned bounding box construction object
  typedef typename K::Construct_bbox_3                 Construct_bbox_3;

  /// The matrix type
  typedef CGAL::Eigen_matrix<FT, 3, 3>                 Matrix;

private:
  typedef typename Matrix::EigenType                   EigenType;

public:
  /// Offers `construct_bbox_3_object()(const Point_3&)`
  Construct_bbox_3 construct_bbox_3_object() const { return Construct_bbox_3(); }

  /// Returns the transpose of a matrix
  Matrix transpose(const Matrix& mat) const
  {
    return Matrix(mat.eigen_object().transpose());
  }

  /// Returns the determinant of a matrix
  FT compute_determinant(const Matrix& matrix) const
  {
    return matrix.eigen_object().determinant();
  }

  /// Performs QR decomposition of matrix A to a unitary matrix and an upper triagonal
  /// and returns the unitary matrix.
  Matrix get_Q(const Matrix& A) const
  {
    Eigen::HouseholderQR<EigenType> qr(A.eigen_object());
    CGAL_assertion(CGAL::abs(compute_determinant(Matrix(EigenType(qr.householderQ()))) - 1.) < 0.000001);

    return Matrix(EigenType(qr.householderQ()));
  }
};
#endif // defined(CGAL_EIGEN3_ENABLED) || defined(DOXYGEN_RUNNING)

} // namespace CGAL

#endif // CGAL_OPTIMAL_BOUNDING_BOX_BOUNDING_BOX_TRAITS_H
