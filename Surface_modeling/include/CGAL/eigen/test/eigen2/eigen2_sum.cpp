// This file is part of Eigen, a lightweight C++ template library
// for linear algebra. Eigen itself is part of the KDE project.
//
// Copyright (C) 2008 Benoit Jacob <jacob.benoit.1@gmail.com>
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

#include "main.h"

template<typename MatrixType> void matrixSum(const MatrixType& m)
{
  typedef typename MatrixType::Scalar Scalar;

  int rows = m.rows();
  int cols = m.cols();

  MatrixType m1 = MatrixType::Random(rows, cols);

  VERIFY_IS_MUCH_SMALLER_THAN(MatrixType::Zero(rows, cols).sum(), Scalar(1));
  VERIFY_IS_APPROX(MatrixType::Ones(rows, cols).sum(), Scalar(float(rows*cols))); // the float() here to shut up excessive MSVC warning about int->complex conversion being lossy
  Scalar x = Scalar(0);
  for(int i = 0; i < rows; i++) for(int j = 0; j < cols; j++) x += m1(i,j);
  VERIFY_IS_APPROX(m1.sum(), x);
}

template<typename VectorType> void vectorSum(const VectorType& w)
{
  typedef typename VectorType::Scalar Scalar;
  int size = w.size();

  VectorType v = VectorType::Random(size);
  for(int i = 1; i < size; i++)
  {
    Scalar s = Scalar(0);
    for(int j = 0; j < i; j++) s += v[j];
    VERIFY_IS_APPROX(s, v.start(i).sum());
  }

  for(int i = 0; i < size-1; i++)
  {
    Scalar s = Scalar(0);
    for(int j = i; j < size; j++) s += v[j];
    VERIFY_IS_APPROX(s, v.end(size-i).sum());
  }

  for(int i = 0; i < size/2; i++)
  {
    Scalar s = Scalar(0);
    for(int j = i; j < size-i; j++) s += v[j];
    VERIFY_IS_APPROX(s, v.segment(i, size-2*i).sum());
  }
}

void test_eigen2_sum()
{
  for(int i = 0; i < g_repeat; i++) {
    CALL_SUBTEST_1( matrixSum(Matrix<float, 1, 1>()) );
    CALL_SUBTEST_2( matrixSum(Matrix2f()) );
    CALL_SUBTEST_3( matrixSum(Matrix4d()) );
    CALL_SUBTEST_4( matrixSum(MatrixXcf(3, 3)) );
    CALL_SUBTEST_5( matrixSum(MatrixXf(8, 12)) );
    CALL_SUBTEST_6( matrixSum(MatrixXi(8, 12)) );
  }
  for(int i = 0; i < g_repeat; i++) {
    CALL_SUBTEST_5( vectorSum(VectorXf(5)) );
    CALL_SUBTEST_7( vectorSum(VectorXd(10)) );
    CALL_SUBTEST_5( vectorSum(VectorXf(33)) );
  }
}
