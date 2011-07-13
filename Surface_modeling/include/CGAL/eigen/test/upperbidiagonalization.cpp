// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2010 Benoit Jacob <jacob.benoit.1@gmail.com>
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
#include <Eigen/SVD>

template<typename MatrixType> void upperbidiag(const MatrixType& m)
{
  const typename MatrixType::Index rows = m.rows();
  const typename MatrixType::Index cols = m.cols();

  typedef typename MatrixType::Scalar Scalar;
  typedef Matrix<typename MatrixType::RealScalar, MatrixType::RowsAtCompileTime,  MatrixType::ColsAtCompileTime> RealMatrixType;

  MatrixType a = MatrixType::Random(rows,cols);
  internal::UpperBidiagonalization<MatrixType> ubd(a);
  RealMatrixType b(rows, cols);
  b.setZero();
  b.block(0,0,cols,cols) = ubd.bidiagonal();
  MatrixType c = ubd.householderU() * b * ubd.householderV().adjoint();
  VERIFY_IS_APPROX(a,c);
}

void test_upperbidiagonalization()
{
  for(int i = 0; i < g_repeat; i++) {
   CALL_SUBTEST_1( upperbidiag(MatrixXf(3,3)) );
   CALL_SUBTEST_2( upperbidiag(MatrixXd(17,12)) );
   CALL_SUBTEST_3( upperbidiag(MatrixXcf(20,20)) );
   CALL_SUBTEST_4( upperbidiag(MatrixXcd(16,15)) );
   CALL_SUBTEST_5( upperbidiag(Matrix<float,6,4>()) );
   CALL_SUBTEST_6( upperbidiag(Matrix<float,5,5>()) );
   CALL_SUBTEST_7( upperbidiag(Matrix<double,4,3>()) );
  }
}
