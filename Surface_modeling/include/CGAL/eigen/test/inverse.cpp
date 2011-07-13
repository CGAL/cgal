// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2008 Gael Guennebaud <gael.guennebaud@inria.fr>
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
#include <Eigen/LU>

template<typename MatrixType> void inverse(const MatrixType& m)
{
  typedef typename MatrixType::Index Index;
  /* this test covers the following files:
     Inverse.h
  */
  Index rows = m.rows();
  Index cols = m.cols();

  typedef typename MatrixType::Scalar Scalar;
  typedef typename NumTraits<Scalar>::Real RealScalar;
  typedef Matrix<Scalar, MatrixType::ColsAtCompileTime, 1> VectorType;

  MatrixType m1(rows, cols),
             m2(rows, cols),
             mzero = MatrixType::Zero(rows, cols),
             identity = MatrixType::Identity(rows, rows);
  createRandomPIMatrixOfRank(rows,rows,rows,m1);
  m2 = m1.inverse();
  VERIFY_IS_APPROX(m1, m2.inverse() );

  VERIFY_IS_APPROX((Scalar(2)*m2).inverse(), m2.inverse()*Scalar(0.5));

  VERIFY_IS_APPROX(identity, m1.inverse() * m1 );
  VERIFY_IS_APPROX(identity, m1 * m1.inverse() );

  VERIFY_IS_APPROX(m1, m1.inverse().inverse() );

  // since for the general case we implement separately row-major and col-major, test that
  VERIFY_IS_APPROX(MatrixType(m1.transpose().inverse()), MatrixType(m1.inverse().transpose()));

#if !defined(EIGEN_TEST_PART_5) && !defined(EIGEN_TEST_PART_6)
  //computeInverseAndDetWithCheck tests
  //First: an invertible matrix
  bool invertible;
  RealScalar det;

  m2.setZero();
  m1.computeInverseAndDetWithCheck(m2, det, invertible);
  VERIFY(invertible);
  VERIFY_IS_APPROX(identity, m1*m2);
  VERIFY_IS_APPROX(det, m1.determinant());

  m2.setZero();
  m1.computeInverseWithCheck(m2, invertible);
  VERIFY(invertible);
  VERIFY_IS_APPROX(identity, m1*m2);

  //Second: a rank one matrix (not invertible, except for 1x1 matrices)
  VectorType v3 = VectorType::Random(rows);
  MatrixType m3 = v3*v3.transpose(), m4(rows,cols);
  m3.computeInverseAndDetWithCheck(m4, det, invertible);
  VERIFY( rows==1 ? invertible : !invertible );
  VERIFY_IS_MUCH_SMALLER_THAN(internal::abs(det-m3.determinant()), RealScalar(1));
  m3.computeInverseWithCheck(m4, invertible);
  VERIFY( rows==1 ? invertible : !invertible );
#endif

  // check in-place inversion
  if(MatrixType::RowsAtCompileTime>=2 && MatrixType::RowsAtCompileTime<=4)
  {
    // in-place is forbidden
    VERIFY_RAISES_ASSERT(m1 = m1.inverse());
  }
  else
  {
    m2 = m1.inverse();
    m1 = m1.inverse();
    VERIFY_IS_APPROX(m1,m2);
  }
}

void test_inverse()
{
  int s;
  for(int i = 0; i < g_repeat; i++) {
    CALL_SUBTEST_1( inverse(Matrix<double,1,1>()) );
    CALL_SUBTEST_2( inverse(Matrix2d()) );
    CALL_SUBTEST_3( inverse(Matrix3f()) );
    CALL_SUBTEST_4( inverse(Matrix4f()) );
    CALL_SUBTEST_4( inverse(Matrix<float,4,4,DontAlign>()) );
    s = internal::random<int>(50,320);
    CALL_SUBTEST_5( inverse(MatrixXf(s,s)) );
    s = internal::random<int>(25,100);
    CALL_SUBTEST_6( inverse(MatrixXcd(s,s)) );
    CALL_SUBTEST_7( inverse(Matrix4d()) );
    CALL_SUBTEST_7( inverse(Matrix<double,4,4,DontAlign>()) );
  }
}
