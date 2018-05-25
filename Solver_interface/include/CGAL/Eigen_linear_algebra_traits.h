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

#ifndef CGAL_EIGEN_LINEAR_ALGEBRA_TRAITS_H
#define CGAL_EIGEN_LINEAR_ALGEBRA_TRAITS_H

#include <CGAL/basic.h>
#include <Eigen/Dense>
#include <Eigen/QR>


namespace CGAL {


template<typename T, int D = Eigen::Dynamic>
class Eigen_dense_vector;

/*!
\ingroup PkgSolver

The class `Eigen_dense_matrix` is a wrapper around \ref thirdpartyEigen "Eigen" matrix type
<a href="http://eigen.tuxfamily.org/dox/classEigen_1_1SparseMatrix.html">`Eigen::DenseMatrix`</a>.

\tparam T Number type.
\tparam D1 Number of rows, or Dynamic.
\tparam D1 Number of cols, or Dynamic.

\sa `CGAL::Eigen_dense_vector<T, D>`
\sa `CGAL::Eigen_linear_algebra_traits`
*/
template<typename T, int D1 = Eigen::Dynamic, int D2 = Eigen::Dynamic>
class Eigen_dense_matrix
{
public:

  /// The internal matrix type from \ref thirdpartyEigen "Eigen".
  typedef Eigen::Matrix<T, D1, D2> EigenType;

  /// Create a dense matrix
  Eigen_dense_matrix(std::size_t nrows, std::size_t ncols)
    : m_matrix(static_cast<int>(nrows), static_cast<int>(ncols))
  {
    CGAL_assertion(m_matrix.rows() > 0);
    CGAL_assertion(m_matrix.cols() > 0);
  }

  /// Create a dense matrix
  Eigen_dense_matrix(int nrows, int ncols)
    : m_matrix(nrows, ncols)
  {
    CGAL_assertion(m_matrix.rows() > 0);
    CGAL_assertion(m_matrix.cols() > 0);
  }

  /// Create a dense matrix out of a \ref thirdpartyEigen "Eigen" matrix type
  Eigen_dense_matrix(const EigenType& eigen_mat)
    : m_matrix(eigen_mat) {}

  Eigen_dense_matrix() : m_matrix() {}

  /// Read access to a matrix coefficient.
  ///
  /// \pre 0 <= i < row_dimension().
  /// \pre 0 <= j < column_dimension().
  T& operator() (int i_, int j_)
  {
    return m_matrix(i_, j_);
  }

  /// Write access to a matrix coefficient: a_ij <- val
  ///
  /// \pre 0 <= i < row_dimension().
  /// \pre 0 <= j < column_dimension().
  void set_coef(std::size_t i_, std::size_t j_, T val)
  {
    int i = static_cast<int>(i_);
    int j = static_cast<int>(j_);
    CGAL_precondition(i < m_matrix.rows());
    CGAL_precondition(j < m_matrix.cols());

    m_matrix.coeffRef(i,j) = val;
  }

  /// Return the matrix number of rows
  std::size_t rows() const {return m_matrix.rows();}
  /// Return the matrix number of cols
  std::size_t cols() const {return m_matrix.cols();}

  /// Resize to i rows and j cols
  void resize(int i_, int j_) { m_matrix.resize(i_, j_);}

  const T& coeff(int i_) const
  {
    return m_matrix.coeff(i_);
  }

  EigenType m_matrix;
};

/*!
\ingroup PkgSolver

The class `Eigen_vector` is a wrapper around \ref thirdpartyEigen "Eigen" dense vector
type <a href="http://eigen.tuxfamily.org/dox/classEigen_1_1Matrix.html"> </a>,
which is a simple array of numbers.

\cgalModels `SvdTraits::Vector`
\cgalModels `SparseLinearAlgebraTraits_d::Vector`.

\tparam T Number type.
\tparam D Number of colums, or Dynamic.

\sa `CGAL::Eigen_dense_matrix<T>`
\sa `CGAL::Eigen_linear_algebra_traits<T>`
*/
template <typename T, int D>
class Eigen_dense_vector
{
private:

  /// The internal vector type from \ref thirdpartyEigen "Eigen".
  typedef Eigen::Matrix<T, D, 1> EigenType;

public:

  /// Create a dense vector out of a \ref thirdpartyEigen "Eigen" vector type
  Eigen_dense_vector(const EigenType&  vec) : m_vector(vec) {}

  /// Read and write and access to a vector coefficient: `a_i`
  const T& coeff(std::size_t i)
  {
    CGAL_assertion(i >= 0);
    CGAL_assertion(i < static_cast<std::size_t>(D));
    return m_vector.coeff(i);
  }

  EigenType m_vector;
};


/*!
\ingroup PkgSolver

The class `Eigen_linear_algebra_traits` provides an interface to linear algebra functionalities of \ref thirdpartyEigen "Eigen".
\ref thirdpartyEigen "Eigen" version 3.1 (or later) must be available on the system.

\sa `CGAL::Eigen_dense_matrix<T, D1, D2>`
\sa `CGAL::Eigen_dense_vector<T, D>`
\sa http://eigen.tuxfamily.org

Example
-------------- 

\code{.cpp}

typedef CGAL::Eigen_linear_algebra_traits Linear_algebra_traits;

// dynamic matrix at run time to store a large amount of data
typedef Linear_algebra_traits::MatrixXd MatrixXd;

// preallocated 3x3 matrix at compile time
typedef Linear_algebra_traits::Matrix3d Matrix3d;

// preallocated 3-cols vector at compile time
typedef Linear_algebra_traits::Vector3d Vector3d;

\endcode
*/

class Eigen_linear_algebra_traits
{
public:
  typedef double NT;
  typedef int Index;

  // dynamic size at run time
  typedef CGAL::Eigen_dense_matrix<NT> MatrixXd;

  // dynamic rows in run time, fixed cols in compile time
  typedef CGAL::Eigen_dense_matrix<NT, Eigen::Dynamic, 3> MatrixX3d;

  // fixed at compile time
  typedef CGAL::Eigen_dense_matrix<NT, 3, 3> Matrix3d;

  // fixed at compile time
  typedef CGAL::Eigen_dense_vector<NT, 3> Vector3d;

  /// Get the transpose of a `CGAL::Eigen_dense_matrix<T, D1, D2>` matrix
  template <class Matrix>
  static Matrix transpose(const Matrix& mat)
  {
    return Matrix(mat.m_matrix.transpose());
  }

  /// Get the determinant of a `CGAL::Eigen_dense_matrix<T, D1, D2>` matrix
  template <class Matrix>
  static NT determinant(const Matrix& mat)
  {
    return mat.m_matrix.determinant();
  }

  /// Performs QR decomposition of matrix A to a unitary matrix and an upper triagonal 
  /// and returns the unitary matrix.
  template <class NT, int D1, int D2>
  static CGAL::Eigen_dense_matrix<NT, D1, D2> qr_factorization(const CGAL::Eigen_dense_matrix<NT, D1, D2>& A)
  {
    Eigen::HouseholderQR<Eigen::Matrix<NT, D1, D2> > qr(A.m_matrix);
    return CGAL::Eigen_dense_matrix<NT, D1, D2>(qr.householderQ());
  }

  template <class Matrix>
  static void qr_factorization(std::vector<Matrix>& simplex)
  {
    for(std::size_t i = 0; i < simplex.size(); ++i)
    {
      Matrix mat = simplex[i].m_matrix;
      simplex[i] = qr_factorization(mat);
    }
  }

  // CGAL::Eigen_dense_vector<NT, D2> : the returned type with D2 may be -1 in compile time,
  // and may not be equal to the expected type.
  // Eigen manages to return a precompiled row out of a dynamic matrix but I don't know how.

  /// Get the row vector out of a `CGAL::Eigen_dense_matrix<T, D1, D2>`. The result is stored in a 
  /// preallocated at compile time 3-column `CGAL::Eigen_dense_vector<T, 3>`
  template <class NT, int D1, int D2>
  static CGAL::Eigen_dense_vector<NT, 3> row3(const CGAL::Eigen_dense_matrix<NT, D1, D2>& A,
                                              int i)
  {
    return CGAL::Eigen_dense_vector<NT, 3>(A.m_matrix.row(i));
  }

};


/// Matrix multiplication. If the columns of A and the rows of B are equal at compile time, 
/// the product is stored at a preallocated at compile time `CGAL::Eigen_dense_matrix`. Otherwise, 
/// the product is stored in a dynamic at run time matrix.
template <class NT, int D1, int D2, int D3>
const CGAL::Eigen_dense_matrix<NT, D1, D3> operator* (const CGAL::Eigen_dense_matrix<NT, D1, D2 >& A,
                                                      const CGAL::Eigen_dense_matrix<NT, D2, D3 >& B)
{
  return CGAL::Eigen_dense_matrix<NT, D1, D3>(A.m_matrix * B.m_matrix);
}

// D2 and D3 may not be equal at compile time, but equal at run time!
// This overload returns a dynamic matrix.
template <class NT, int D1, int D2, int D3, int D4>
const CGAL::Eigen_dense_matrix<NT> operator* (const CGAL::Eigen_dense_matrix<NT, D1, D2 >& A,
                                              const CGAL::Eigen_dense_matrix<NT, D3, D4 >& B)
{
  return CGAL::Eigen_dense_matrix<NT>(A.m_matrix * B.m_matrix);
}

// scalar - matrix multiplication
template <class NT, int D1, int D2>
const CGAL::Eigen_dense_matrix<NT, D1, D2> operator* (const NT& scalar,
                                                      const CGAL::Eigen_dense_matrix<NT, D1, D2>& B)
{
  return CGAL::Eigen_dense_matrix<NT, D1, D2>(scalar * B.m_matrix);
}

template <class NT, int D1, int D2>
const CGAL::Eigen_dense_matrix<NT, D1, D2> operator* (const CGAL::Eigen_dense_matrix<NT, D1, D2>& A,
                                                      const NT& scalar)
{
  return CGAL::Eigen_dense_matrix<NT, D1, D2>(A.m_matrix * scalar);
}

template <class NT, int D1, int D2>
const CGAL::Eigen_dense_matrix<NT, D1, D2> operator/ (const CGAL::Eigen_dense_matrix<NT, D1, D2>& A,
                                                      const double& scalar)
{
  return CGAL::Eigen_dense_matrix<NT, D1, D2>(A.m_matrix / scalar);
}

template <class NT, int D1, int D2>
const CGAL::Eigen_dense_matrix<NT, D1, D2> operator/ (const double& scalar,
                                                      const CGAL::Eigen_dense_matrix<NT, D1, D2> & A)
{
  return CGAL::Eigen_dense_matrix<NT, D1, D2> (scalar / A.m_matrix);
}

// addition
template <class NT, int D1, int D2>
const CGAL::Eigen_dense_matrix<NT, D1, D2> operator+ (const CGAL::Eigen_dense_matrix<NT, D1, D2> & A,
                                                      const CGAL::Eigen_dense_matrix<NT, D1, D2> & B)
{
  return CGAL::Eigen_dense_matrix<NT, D1, D2> (A.m_matrix + B.m_matrix);
}

// vector - matrix multiplication
template <class NT, int D1, int D2>
const Eigen_dense_vector<NT, D1> operator* (const CGAL::Eigen_dense_matrix<NT, D1, D2>& A,
                                            const CGAL::Eigen_dense_vector<NT, D2>& V)
{
  return Eigen_dense_vector<NT, D1>(A.m_matrix * V.m_vector);
}



}
// end namespace



#endif // CGAL_EIGEN_LINEAR_ALGEBRA_TRAITS_H
