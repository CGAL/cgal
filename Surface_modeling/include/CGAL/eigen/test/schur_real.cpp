// This file is part of Eigen, a lightweight C++ template library
// for linear algebra. Eigen itself is part of the KDE project.
//
// Copyright (C) 2010 Jitse Niesen <jitse@maths.leeds.ac.uk>
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
#include <limits>
#include <Eigen/Eigenvalues>

template<typename MatrixType> void verifyIsQuasiTriangular(const MatrixType& T)
{
  typedef typename MatrixType::Index Index;

  const Index size = T.cols();
  typedef typename MatrixType::Scalar Scalar;

  // Check T is lower Hessenberg
  for(int row = 2; row < size; ++row) {
    for(int col = 0; col < row - 1; ++col) {
      VERIFY(T(row,col) == Scalar(0));
    }
  }

  // Check that any non-zero on the subdiagonal is followed by a zero and is
  // part of a 2x2 diagonal block with imaginary eigenvalues.
  for(int row = 1; row < size; ++row) {
    if (T(row,row-1) != Scalar(0)) {
      VERIFY(row == size-1 || T(row+1,row) == 0);
      Scalar tr = T(row-1,row-1) + T(row,row);
      Scalar det = T(row-1,row-1) * T(row,row) - T(row-1,row) * T(row,row-1);
      VERIFY(4 * det > tr * tr);
    }
  }
}

template<typename MatrixType> void schur(int size = MatrixType::ColsAtCompileTime)
{
  // Test basic functionality: T is quasi-triangular and A = U T U*
  for(int counter = 0; counter < g_repeat; ++counter) {
    MatrixType A = MatrixType::Random(size, size);
    RealSchur<MatrixType> schurOfA(A);
    VERIFY_IS_EQUAL(schurOfA.info(), Success);
    MatrixType U = schurOfA.matrixU();
    MatrixType T = schurOfA.matrixT();
    verifyIsQuasiTriangular(T);
    VERIFY_IS_APPROX(A, U * T * U.transpose());
  }

  // Test asserts when not initialized
  RealSchur<MatrixType> rsUninitialized;
  VERIFY_RAISES_ASSERT(rsUninitialized.matrixT());
  VERIFY_RAISES_ASSERT(rsUninitialized.matrixU());
  VERIFY_RAISES_ASSERT(rsUninitialized.info());
  
  // Test whether compute() and constructor returns same result
  MatrixType A = MatrixType::Random(size, size);
  RealSchur<MatrixType> rs1;
  rs1.compute(A);
  RealSchur<MatrixType> rs2(A);
  VERIFY_IS_EQUAL(rs1.info(), Success);
  VERIFY_IS_EQUAL(rs2.info(), Success);
  VERIFY_IS_EQUAL(rs1.matrixT(), rs2.matrixT());
  VERIFY_IS_EQUAL(rs1.matrixU(), rs2.matrixU());

  // Test computation of only T, not U
  RealSchur<MatrixType> rsOnlyT(A, false);
  VERIFY_IS_EQUAL(rsOnlyT.info(), Success);
  VERIFY_IS_EQUAL(rs1.matrixT(), rsOnlyT.matrixT());
  VERIFY_RAISES_ASSERT(rsOnlyT.matrixU());

  if (size > 2)
  {
    // Test matrix with NaN
    A(0,0) = std::numeric_limits<typename MatrixType::Scalar>::quiet_NaN();
    RealSchur<MatrixType> rsNaN(A);
    VERIFY_IS_EQUAL(rsNaN.info(), NoConvergence);
  }
}

void test_schur_real()
{
  CALL_SUBTEST_1(( schur<Matrix4f>() ));
  CALL_SUBTEST_2(( schur<MatrixXd>(internal::random<int>(1,50)) ));
  CALL_SUBTEST_3(( schur<Matrix<float, 1, 1> >() ));
  CALL_SUBTEST_4(( schur<Matrix<double, 3, 3, Eigen::RowMajor> >() ));

  // Test problem size constructors
  CALL_SUBTEST_5(RealSchur<MatrixXf>(10));
}
