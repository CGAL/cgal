// Copyright (c) 2015 GeometryFactory (France), All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_DEFAULT_DIAGONALIZE_TRAITS_H
#define CGAL_DEFAULT_DIAGONALIZE_TRAITS_H

#ifdef CGAL_EIGEN3_ENABLED
#include <CGAL/Eigen_diagonalize_traits.h>
#else
#include <CGAL/Diagonalize_traits.h>
#endif

/// \cond SKIP_IN_MANUAL

namespace CGAL {

/// \ingroup PkgSolverInterfaceRef
///
/// The class `Default_diagonalize_traits` is a wrapper designed to automatically
/// use `Eigen_diagonalize_traits` if Eigen is available and otherwise use
/// the fallback `Diagonalize_traits` class of %CGAL.
///
/// \tparam FT Number type
/// \tparam dim Dimension of the matrices and vectors
///
/// \cgalModels `DiagonalizeTraits`
template <typename FT, unsigned int dim = 3>
class Default_diagonalize_traits
{
#ifdef CGAL_EIGEN3_ENABLED
  typedef Eigen_diagonalize_traits<FT, dim>             Base;
#else
  typedef Diagonalize_traits<FT, dim>                   Base;
#endif

public:
  typedef std::array<FT, dim>                         Vector;
  typedef std::array<FT, dim*dim>                     Matrix;
  typedef std::array<FT, (dim * (dim+1) / 2)>         Covariance_matrix;

  /// Fill `eigenvalues` with the eigenvalues of the covariance matrix represented by `cov`.
  /// Eigenvalues are sorted by increasing order.
  /// \return `true` if the operation was successful and `false` otherwise.
  static bool diagonalize_selfadjoint_covariance_matrix(const Covariance_matrix& cov,
                                                        Vector& eigenvalues)
  {
    return Base::diagonalize_selfadjoint_covariance_matrix(cov, eigenvalues);
  }

  /// Fill `eigenvalues` with the eigenvalues and `eigenvectors` with
  /// the eigenvectors of the covariance matrix represented by `cov`.
  /// Eigenvalues are sorted by increasing order.
  /// \return `true` if the operation was successful and `false` otherwise.
  static bool diagonalize_selfadjoint_covariance_matrix(const Covariance_matrix& cov,
                                                        Vector& eigenvalues,
                                                        Matrix& eigenvectors)
  {
    return Base::diagonalize_selfadjoint_covariance_matrix(cov, eigenvalues, eigenvectors);
  }

  /// Extract the eigenvector associated to the largest eigenvalue
  /// of the covariance matrix represented by `cov`.
  /// \return `true` if the operation was successful and `false` otherwise.
  static bool extract_largest_eigenvector_of_covariance_matrix(const Covariance_matrix& cov,
                                                               Vector& normal)
  {
    return Base::extract_largest_eigenvector_of_covariance_matrix(cov, normal);
  }
};

} // namespace CGAL

/// \endcond

#endif // CGAL_DEFAULT_DIAGONALIZE_TRAITS_H
