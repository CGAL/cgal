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

#include <CGAL/Aff_transformation_3.h>

#ifdef CGAL_EIGEN3_ENABLED
#include <CGAL/Eigen_matrix.h>
#include <CGAL/Eigen_vector.h>
#include <Eigen/QR>
#endif

namespace CGAL {

#if defined(CGAL_EIGEN3_ENABLED) || defined(DOXYGEN_RUNNING)

/// \ingroup PkgOptimalBoundingBoxClasses
///
/// The class `CGAL::Oriented_bounding_box_traits_3` is a traits type
/// to be used with the overloads of the function `CGAL::oriented_bounding_box()`.
///
/// \attention This class requires the \ref thirdpartyEigen "Eigen" library.
///
/// \tparam K must be a model of `Kernel`
///
/// \cgalModels `OrientedBoundingBoxTraits_3`
///
template <typename K>
class Oriented_bounding_box_traits_3
  : public K
{
public:
  /// The field number type
  typedef typename K::FT                               FT;

  /// The affine transformation type
  typedef typename CGAL::Aff_transformation_3<K>       Aff_transformation_3;

  /// The matrix type
  typedef CGAL::Eigen_matrix<FT, 3, 3>                 Matrix;

  /// The matrix type
  typedef CGAL::Eigen_vector<FT, 3>                    Vector;

private:
  typedef typename Matrix::EigenType                   EigenType;

public:
  /// Performs the QR-decomposition of the matrix `m` to a unitary matrix and an upper triagonal
  /// and returns the unitary matrix
  static Matrix get_Q(const Matrix& m)
  {
    Eigen::HouseholderQR<EigenType> qr(m.eigen_object());
    return Matrix(EigenType(qr.householderQ()));
  }
};
#endif // defined(CGAL_EIGEN3_ENABLED) || defined(DOXYGEN_RUNNING)

} // namespace CGAL

#endif // CGAL_OPTIMAL_BOUNDING_BOX_BOUNDING_BOX_TRAITS_H
