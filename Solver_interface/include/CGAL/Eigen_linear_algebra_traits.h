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

template<typename T, int D1 = Eigen::Dynamic, int D2 = Eigen::Dynamic>
class Eigen_dense_matrix
{
public:
  typedef Eigen::Matrix<T, D1, D2> EigenType;

  Eigen_dense_matrix(std::size_t nrows, std::size_t ncols)
    : m_matrix(static_cast<int>(nrows), static_cast<int>(ncols))
  {
    CGAL_assertion(m_matrix.rows() > 0);
    CGAL_assertion(m_matrix.cols() > 0);
  }

  Eigen_dense_matrix(int nrows, int ncols)
    : m_matrix(nrows, ncols)
  {
    CGAL_assertion(m_matrix.rows() > 0);
    CGAL_assertion(m_matrix.cols() > 0);
  }

  Eigen_dense_matrix(const EigenType& eigen_mat)
    : m_matrix(eigen_mat) {}

  Eigen_dense_matrix() : m_matrix() {}

  T& operator() (int i_, int j_)
  {
    return m_matrix(i_, j_);
  }

  void set_coef(std::size_t i_, std::size_t j_, T val)
  {
    int i = static_cast<int>(i_);
    int j = static_cast<int>(j_);
    CGAL_precondition(i < m_matrix.rows());
    CGAL_precondition(j < m_matrix.cols());

    m_matrix.coeffRef(i,j) = val;
  }

  std::size_t rows() const {return m_matrix.rows();}
  std::size_t cols() const {return m_matrix.cols();}

  void resize(int i_, int j_) { m_matrix.resize(i_, j_);}

  const T& coeff(int i_) const
  {
    return m_matrix.coeff(i_);
  }

  mutable EigenType m_matrix;
};


template <typename T, int D>
class Eigen_dense_vector
{
private:
  typedef Eigen::Matrix<T, D, 1> EigenType;

public:
  Eigen_dense_vector(const EigenType&  vec) : m_vector(vec) {}

  const T& coeff(std::size_t i)
  {
    CGAL_assertion(i >= 0);
    CGAL_assertion(i < static_cast<std::size_t>(D));
    return m_vector.coeff(i);
  }

  mutable EigenType m_vector;
};


class Eigen_linear_algebra_traits
{
public:
  typedef double NT;

  // dynamic size at run time
  typedef CGAL::Eigen_dense_matrix<NT> MatrixXd;

  // dynamic rows in run time, fixed cols in compile time
  typedef CGAL::Eigen_dense_matrix<NT, Eigen::Dynamic, 3> MatrixX3d;

  // fixed at compile time
  typedef CGAL::Eigen_dense_matrix<NT, 3, 3> Matrix3d;

  // fixed at compile time
  typedef CGAL::Eigen_dense_vector<NT, 3> Vector3d;

  template <class Matrix>
  static Matrix transpose(const Matrix& mat)
  {
    return Matrix(mat.m_matrix.transpose());
  }

  template <class Matrix>
  static NT determinant(const Matrix& mat)
  {
    return mat.m_matrix.determinant();
  }

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

  template <class NT, int D1, int D2>
  static CGAL::Eigen_dense_vector<NT, D2> row(const CGAL::Eigen_dense_matrix<NT, D1, D2>& A,
                                              int i)
  {
    return CGAL::Eigen_dense_vector<NT, D2>(A.m_matrix.row(i));
  }
};


// matrix multiplication
/*
template <class NT, int D1, int D2>
const CGAL::Eigen_dense_matrix<NT> operator* (const CGAL::Eigen_dense_matrix<NT>& A,
                                              const CGAL::Eigen_dense_matrix<NT, D1, D2 >& B)
{
  return CGAL::Eigen_dense_matrix<NT>(A.m_matrix * B.m_matrix);
}

template <class NT, int D1, int D2>
const CGAL::Eigen_dense_matrix<NT> operator* (const CGAL::Eigen_dense_matrix<NT, D1, D2 >& A,
                                              const CGAL::Eigen_dense_matrix<NT>& B)
{
  return CGAL::Eigen_dense_matrix<NT>(A.m_matrix * B.m_matrix);
}
*/

template <class NT, int D1, int D2, int D3>
const CGAL::Eigen_dense_matrix<NT, D1, D3> operator* (const CGAL::Eigen_dense_matrix<NT, D1, D2 >& A,
                                                      const CGAL::Eigen_dense_matrix<NT, D2, D3 >& B)
{
  return CGAL::Eigen_dense_matrix<NT, D1, D3>(A.m_matrix * B.m_matrix);
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
