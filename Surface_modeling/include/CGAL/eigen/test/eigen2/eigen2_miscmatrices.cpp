// This file is part of Eigen, a lightweight C++ template library
// for linear algebra. Eigen itself is part of the KDE project.
//
// Copyright (C) 2006-2008 Benoit Jacob <jacob.benoit.1@gmail.com>
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

template<typename MatrixType> void miscMatrices(const MatrixType& m)
{
  /* this test covers the following files:
     DiagonalMatrix.h Ones.h
  */

  typedef typename MatrixType::Scalar Scalar;
  typedef Matrix<Scalar, MatrixType::RowsAtCompileTime, 1> VectorType;
  typedef Matrix<Scalar, 1, MatrixType::ColsAtCompileTime> RowVectorType;
  int rows = m.rows();
  int cols = m.cols();

  int r = ei_random<int>(0, rows-1), r2 = ei_random<int>(0, rows-1), c = ei_random<int>(0, cols-1);
  VERIFY_IS_APPROX(MatrixType::Ones(rows,cols)(r,c), static_cast<Scalar>(1));
  MatrixType m1 = MatrixType::Ones(rows,cols);
  VERIFY_IS_APPROX(m1(r,c), static_cast<Scalar>(1));
  VectorType v1 = VectorType::Random(rows);
  v1[0];
  Matrix<Scalar, MatrixType::RowsAtCompileTime, MatrixType::RowsAtCompileTime>
  square = v1.asDiagonal();
  if(r==r2) VERIFY_IS_APPROX(square(r,r2), v1[r]);
  else VERIFY_IS_MUCH_SMALLER_THAN(square(r,r2), static_cast<Scalar>(1));
  square = MatrixType::Zero(rows, rows);
  square.diagonal() = VectorType::Ones(rows);
  VERIFY_IS_APPROX(square, MatrixType::Identity(rows, rows));
}

void test_eigen2_miscmatrices()
{
  for(int i = 0; i < g_repeat; i++) {
    CALL_SUBTEST_1( miscMatrices(Matrix<float, 1, 1>()) );
    CALL_SUBTEST_2( miscMatrices(Matrix4d()) );
    CALL_SUBTEST_3( miscMatrices(MatrixXcf(3, 3)) );
    CALL_SUBTEST_4( miscMatrices(MatrixXi(8, 12)) );
    CALL_SUBTEST_5( miscMatrices(MatrixXcd(20, 20)) );
  }
}
