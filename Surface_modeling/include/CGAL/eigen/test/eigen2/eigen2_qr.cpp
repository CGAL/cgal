// This file is part of Eigen, a lightweight C++ template library
// for linear algebra. Eigen itself is part of the KDE project.
//
// Copyright (C) 2008 Gael Guennebaud <g.gael@free.fr>
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
#include <Eigen/QR>

template<typename MatrixType> void qr(const MatrixType& m)
{
  /* this test covers the following files:
     QR.h
  */
  int rows = m.rows();
  int cols = m.cols();

  typedef typename MatrixType::Scalar Scalar;
  typedef Matrix<Scalar, MatrixType::ColsAtCompileTime, MatrixType::ColsAtCompileTime> SquareMatrixType;
  typedef Matrix<Scalar, MatrixType::ColsAtCompileTime, 1> VectorType;

  MatrixType a = MatrixType::Random(rows,cols);
  QR<MatrixType> qrOfA(a);
  VERIFY_IS_APPROX(a, qrOfA.matrixQ() * qrOfA.matrixR());
  VERIFY_IS_NOT_APPROX(a+MatrixType::Identity(rows, cols), qrOfA.matrixQ() * qrOfA.matrixR());

  #if 0 // eigenvalues module not yet ready
  SquareMatrixType b = a.adjoint() * a;

  // check tridiagonalization
  Tridiagonalization<SquareMatrixType> tridiag(b);
  VERIFY_IS_APPROX(b, tridiag.matrixQ() * tridiag.matrixT() * tridiag.matrixQ().adjoint());

  // check hessenberg decomposition
  HessenbergDecomposition<SquareMatrixType> hess(b);
  VERIFY_IS_APPROX(b, hess.matrixQ() * hess.matrixH() * hess.matrixQ().adjoint());
  VERIFY_IS_APPROX(tridiag.matrixT(), hess.matrixH());
  b = SquareMatrixType::Random(cols,cols);
  hess.compute(b);
  VERIFY_IS_APPROX(b, hess.matrixQ() * hess.matrixH() * hess.matrixQ().adjoint());
  #endif
}

void test_eigen2_qr()
{
  for(int i = 0; i < 1; i++) {
    CALL_SUBTEST_1( qr(Matrix2f()) );
    CALL_SUBTEST_2( qr(Matrix4d()) );
    CALL_SUBTEST_3( qr(MatrixXf(12,8)) );
    CALL_SUBTEST_4( qr(MatrixXcd(5,5)) );
    CALL_SUBTEST_4( qr(MatrixXcd(7,3)) );
  }

#ifdef EIGEN_TEST_PART_5
  // small isFullRank test
  {
    Matrix3d mat;
    mat << 1, 45, 1, 2, 2, 2, 1, 2, 3;
    VERIFY(mat.qr().isFullRank());
    mat << 1, 1, 1, 2, 2, 2, 1, 2, 3;
    //always returns true in eigen2support
    //VERIFY(!mat.qr().isFullRank());
  }

#endif
}
