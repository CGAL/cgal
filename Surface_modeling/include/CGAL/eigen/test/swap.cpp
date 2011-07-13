// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2009 Benoit Jacob <jacob.benoit.1@gmail.com>
//
// Eigen is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3 of the License, or (at your option) any later version.
//
// Alternatively, you can redistribute it and/or
// modify it under the terms of the GNU General Public License as
// published by the Free Software Foundation; either version 2 of
// the License, or (at your option) any later version.
//
// Eigen is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License or the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License and a copy of the GNU General Public License along with
// Eigen. If not, see <http://www.gnu.org/licenses/>.

#define EIGEN_NO_STATIC_ASSERT
#include "main.h"

template<typename T>
struct other_matrix_type
{
  typedef int type;
};

template<typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
struct other_matrix_type<Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> >
{
  typedef Matrix<_Scalar, _Rows, _Cols, _Options^RowMajor, _MaxRows, _MaxCols> type;
};

template<typename MatrixType> void swap(const MatrixType& m)
{
  typedef typename other_matrix_type<MatrixType>::type OtherMatrixType;
  typedef typename MatrixType::Scalar Scalar;

  eigen_assert((!internal::is_same<MatrixType,OtherMatrixType>::value));
  typename MatrixType::Index rows = m.rows();
  typename MatrixType::Index cols = m.cols();
  
  // construct 3 matrix guaranteed to be distinct
  MatrixType m1 = MatrixType::Random(rows,cols);
  MatrixType m2 = MatrixType::Random(rows,cols) + Scalar(100) * MatrixType::Identity(rows,cols);
  OtherMatrixType m3 = OtherMatrixType::Random(rows,cols) + Scalar(200) * OtherMatrixType::Identity(rows,cols);
  
  MatrixType m1_copy = m1;
  MatrixType m2_copy = m2;
  OtherMatrixType m3_copy = m3;
  
  // test swapping 2 matrices of same type
  m1.swap(m2);
  VERIFY_IS_APPROX(m1,m2_copy);
  VERIFY_IS_APPROX(m2,m1_copy);
  m1 = m1_copy;
  m2 = m2_copy;
  
  // test swapping 2 matrices of different types
  m1.swap(m3);
  VERIFY_IS_APPROX(m1,m3_copy);
  VERIFY_IS_APPROX(m3,m1_copy);
  m1 = m1_copy;
  m3 = m3_copy;
  
  // test swapping matrix with expression
  m1.swap(m2.block(0,0,rows,cols));
  VERIFY_IS_APPROX(m1,m2_copy);
  VERIFY_IS_APPROX(m2,m1_copy);
  m1 = m1_copy;
  m2 = m2_copy;

  // test swapping two expressions of different types
  m1.transpose().swap(m3.transpose());
  VERIFY_IS_APPROX(m1,m3_copy);
  VERIFY_IS_APPROX(m3,m1_copy);
  m1 = m1_copy;
  m3 = m3_copy;
  
  // test assertion on mismatching size -- matrix case
  VERIFY_RAISES_ASSERT(m1.swap(m1.row(0)));
  // test assertion on mismatching size -- xpr case
  VERIFY_RAISES_ASSERT(m1.row(0).swap(m1));
}

void test_swap()
{
  CALL_SUBTEST_1( swap(Matrix3f()) ); // fixed size, no vectorization 
  CALL_SUBTEST_2( swap(Matrix4d()) ); // fixed size, possible vectorization 
  CALL_SUBTEST_3( swap(MatrixXd(3,3)) ); // dyn size, no vectorization 
  CALL_SUBTEST_4( swap(MatrixXf(30,30)) ); // dyn size, possible vectorization 
}
