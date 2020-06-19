// Copyright (c) 2012  INRIA Bordeaux Sud-Ouest (France), All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Gael Guennebaud

#ifndef CGAL_SOLVER_INTERFACE_EIGEN_MATRIX_H
#define CGAL_SOLVER_INTERFACE_EIGEN_MATRIX_H

#include <CGAL/basic.h> // include basic.h before testing #defines

#include <CGAL/Eigen_sparse_matrix.h> // for backward compatibility

#include <Eigen/Dense>

namespace CGAL {

/*!
\ingroup PkgSolverInterfaceRef

The class `Eigen_matrix` is a wrapper around `Eigen` matrix type
<a href="http://eigen.tuxfamily.org/dox/classEigen_1_1Matrix.html">`Eigen::Matrix`</a>.

\cgalModels `SvdTraits::Matrix`

\tparam T Number type.
\tparam D1 Number of rows, or Dynamic
\tparam D2 Number of columns, or Dynamic

\sa `CGAL::Eigen_vector<T>`
\sa `CGAL::Eigen_sparse_matrix<T>`
\sa `CGAL::Eigen_sparse_symmetric_matrix<T>`
*/
template <class FT, int D1 = ::Eigen::Dynamic, int D2 = ::Eigen::Dynamic>
struct Eigen_matrix
  : public ::Eigen::Matrix<FT, D1, D2>
{
  /// The internal matrix type from \ref thirdpartyEigen "Eigen".
  typedef ::Eigen::Matrix<FT, D1, D2> EigenType;

  /// Constructs a null matrix.
  Eigen_matrix() { }

  /// Constructs an uninitialized matrix with `nr` rows and `nc` columns.
  /// This is useful for dynamic-size matrices.
  /// For fixed-size matrices, it is redundant to pass these parameters.
  Eigen_matrix(std::size_t nr, std::size_t nc) : EigenType(nr, nc) { }

  /// Constructs a matrix from an Eigen matrix.
  Eigen_matrix(const EigenType& b) : EigenType(b) { }

  /// Returns the matrix number of rows.
  std::size_t number_of_rows() const { return this->rows(); }
  /// Returns the matrix number of columns.
  std::size_t number_of_columns() const { return this->cols(); }

  /// Returns the value of the matrix at position (i,j).
  FT operator()( std::size_t i , std::size_t j ) const { return EigenType::operator()(i,j); }

  /// Writes access to a matrix coefficient: `a_ij` <- `val`.
  void set(std::size_t i, std::size_t j, FT value) { this->coeffRef(i,j) = value; }

  /// Returns the internal matrix, with type `EigenType`.
  const EigenType& eigen_object() const { return static_cast<const EigenType&>(*this); }

public:
  /// \cond SKIP_IN_MANUAL
  friend Eigen_matrix operator*(const FT c, const Eigen_matrix& M)
  {
    return Eigen_matrix(c * M.eigen_object());
  }

  friend Eigen_matrix operator*(const Eigen_matrix& M0, const Eigen_matrix& M1)
  {
    return Eigen_matrix(M0.eigen_object() * M1.eigen_object());
  }

  friend Eigen_matrix operator+(const Eigen_matrix& M0, const Eigen_matrix& M1)
  {
    return Eigen_matrix(M0.eigen_object() + M1.eigen_object());
  }
  /// \endcond
};

} //namespace CGAL

#endif // CGAL_SOLVER_INTERFACE_EIGEN_MATRIX_H
