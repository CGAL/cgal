// Copyright (c) 2014  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Jocelyn Meyron and Quentin Mérigot
//

#ifndef CGAL_EIGEN_DIAGONALIZE_TRAITS_H
#define CGAL_EIGEN_DIAGONALIZE_TRAITS_H

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

// If the matrix to diagonalize is of dimension 2x2 or 3x3, Eigen
// provides a faster implementation using a closed-form
// algorithm. However, it offers less precision. See:
// https://eigen.tuxfamily.org/dox/classEigen_1_1SelfAdjointEigenSolver.html
// This is usually acceptable for CGAL algorithms but one might want
// to use the slower but more accurate version. In that case, just
// uncomment the following line:

//#define DO_NOT_USE_EIGEN_COMPUTEDIRECT_FOR_DIAGONALIZATION

#include <CGAL/array.h>

namespace CGAL {

/// \ingroup PkgSolverInterfaceRef
///
/// The class `Eigen_diagonalize_traits` provides an interface to the
/// diagonalization of covariance matrices of \ref thirdpartyEigen
/// "Eigen".
///
/// \ref thirdpartyEigen "Eigen" version 3.1 (or later) must be available on the system.
///
/// \tparam FT Number type
/// \tparam dim Dimension of the matrices and vectors
///
/// \cgalModels `DiagonalizeTraits`
///
/// \sa http://eigen.tuxfamily.org/index.php?title=Main_Page
template <typename FT, unsigned int dim = 3>
class Eigen_diagonalize_traits
{
public:
  typedef std::array<FT, dim>                  Vector;
  typedef std::array<FT, dim*dim>              Matrix;
  typedef std::array<FT, (dim * (dim+1) / 2)>  Covariance_matrix;

private:
  typedef Eigen::Matrix<FT, dim, dim>            EigenMatrix;
  typedef Eigen::Matrix<FT, dim, 1>              EigenVector;

  /// Construct the covariance matrix
  static EigenMatrix construct_covariance_matrix(const Covariance_matrix& cov)
  {
    EigenMatrix m;

    for(std::size_t i=0; i<dim; ++i)
    {
      for(std::size_t j=i; j<dim; ++j)
      {
        m(i,j) = static_cast<float>(cov[(dim * i) + j - ((i * (i+1)) / 2)]);

        if(i != j)
          m(j,i) = m(i,j);
      }
    }

    return m;
  }

  /// Fill `eigenvalues` with the eigenvalues and `eigenvectors` with
  /// the eigenvectors of the selfadjoint matrix represented by `m`.
  /// Eigenvalues are sorted by increasing order.
  /// \return `true` if the operation was successful and `false` otherwise.
  static bool diagonalize_selfadjoint_matrix(EigenMatrix& m,
                                             EigenMatrix& eigenvectors,
                                             EigenVector& eigenvalues)
  {
    Eigen::SelfAdjointEigenSolver<EigenMatrix> eigensolver;

#ifndef DO_NOT_USE_EIGEN_COMPUTEDIRECT_FOR_DIAGONALIZATION
    if(dim == 2 || dim == 3)
      eigensolver.computeDirect(m);
    else
#endif
      eigensolver.compute(m);

    if(eigensolver.info() != Eigen::Success)
      return false;

    eigenvalues = eigensolver.eigenvalues();
    eigenvectors = eigensolver.eigenvectors();

    return true;
  }

public:
  /// Fill `eigenvalues` with the eigenvalues of the covariance matrix represented by `cov`.
  /// Eigenvalues are sorted by increasing order.
  /// \return `true` if the operation was successful and `false` otherwise.
  static bool diagonalize_selfadjoint_covariance_matrix(const Covariance_matrix& cov,
                                                        Vector& eigenvalues)
  {
    EigenMatrix m = construct_covariance_matrix(cov);

    // Diagonalizing the matrix
    EigenVector eigenvalues_;
    EigenMatrix eigenvectors_;
    bool res = diagonalize_selfadjoint_matrix(m, eigenvectors_, eigenvalues_);

    if(res)
    {
      for(std::size_t i=0; i<dim; ++i)
        eigenvalues[i] = static_cast<FT>(eigenvalues_[i]);
    }

    return res;
  }

  /// Fill `eigenvalues` with the eigenvalues and `eigenvectors` with
  /// the eigenvectors of the covariance matrix represented by `cov`.
  /// Eigenvalues are sorted by increasing order.
  /// \return `true` if the operation was successful and `false` otherwise.
  static bool diagonalize_selfadjoint_covariance_matrix(const Covariance_matrix& cov,
                                                        Vector& eigenvalues,
                                                        Matrix& eigenvectors)
  {
    EigenMatrix m = construct_covariance_matrix(cov);

    // Diagonalizing the matrix
    EigenVector eigenvalues_;
    EigenMatrix eigenvectors_;
    bool res = diagonalize_selfadjoint_matrix(m, eigenvectors_, eigenvalues_);

    if(res)
    {
      for(std::size_t i=0; i<dim; ++i)
      {
        eigenvalues[i] = static_cast<FT>(eigenvalues_[i]);

        for(std::size_t j=0; j<dim; ++j)
          eigenvectors[dim*i + j] = static_cast<FT>(eigenvectors_(j,i));
      }
    }
    else{
      for(std::size_t i=0; i<dim; ++i)
      {
        eigenvalues[i] = static_cast<FT>(0.);

        for(std::size_t j=0; j<dim; ++j)
          eigenvectors[dim*i + j] = static_cast<FT>(0.);
      }
    }

    return res;
  }

  /// Extract the eigenvector associated to the largest eigenvalue
  /// of the covariance matrix represented by `cov`.
  /// \return `true` if the operation was successful and `false` otherwise.
  static bool extract_largest_eigenvector_of_covariance_matrix(const Covariance_matrix& cov,
                                                               Vector& normal)
  {
    // Construct covariance matrix
    EigenMatrix m = construct_covariance_matrix(cov);

    // Diagonalizing the matrix
    EigenVector eigenvalues;
    EigenMatrix eigenvectors;
    if(! diagonalize_selfadjoint_matrix(m, eigenvectors, eigenvalues))
      return false;

    // Eigenvalues are sorted by increasing order
    for(unsigned int i=0; i<dim; ++i)
      normal[i] = static_cast<FT> (eigenvectors(i, dim-1));

    return true;
  }
};

} // namespace CGAL

#endif // CGAL_EIGEN_DIAGONALIZE_TRAITS_H
