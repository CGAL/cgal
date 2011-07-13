// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2008 Gael Guennebaud <gael.guennebaud@inria.fr>
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
  typedef typename MatrixType::Index Index;

  Index rows = m.rows();
  Index cols = m.cols();

  typedef typename MatrixType::Scalar Scalar;
  typedef Matrix<Scalar, MatrixType::RowsAtCompileTime, MatrixType::RowsAtCompileTime> MatrixQType;
  typedef Matrix<Scalar, MatrixType::ColsAtCompileTime, 1> VectorType;

  MatrixType a = MatrixType::Random(rows,cols);
  HouseholderQR<MatrixType> qrOfA(a);

  MatrixQType q = qrOfA.householderQ();
  VERIFY_IS_UNITARY(q);

  MatrixType r = qrOfA.matrixQR().template triangularView<Upper>();
  VERIFY_IS_APPROX(a, qrOfA.householderQ() * r);
}

template<typename MatrixType, int Cols2> void qr_fixedsize()
{
  enum { Rows = MatrixType::RowsAtCompileTime, Cols = MatrixType::ColsAtCompileTime };
  typedef typename MatrixType::Scalar Scalar;
  Matrix<Scalar,Rows,Cols> m1 = Matrix<Scalar,Rows,Cols>::Random();
  HouseholderQR<Matrix<Scalar,Rows,Cols> > qr(m1);

  Matrix<Scalar,Rows,Cols> r = qr.matrixQR();
  // FIXME need better way to construct trapezoid
  for(int i = 0; i < Rows; i++) for(int j = 0; j < Cols; j++) if(i>j) r(i,j) = Scalar(0);

  VERIFY_IS_APPROX(m1, qr.householderQ() * r);

  Matrix<Scalar,Cols,Cols2> m2 = Matrix<Scalar,Cols,Cols2>::Random(Cols,Cols2);
  Matrix<Scalar,Rows,Cols2> m3 = m1*m2;
  m2 = Matrix<Scalar,Cols,Cols2>::Random(Cols,Cols2);
  m2 = qr.solve(m3);
  VERIFY_IS_APPROX(m3, m1*m2);
}

template<typename MatrixType> void qr_invertible()
{
  typedef typename NumTraits<typename MatrixType::Scalar>::Real RealScalar;
  typedef typename MatrixType::Scalar Scalar;

  int size = internal::random<int>(10,50);

  MatrixType m1(size, size), m2(size, size), m3(size, size);
  m1 = MatrixType::Random(size,size);

  if (internal::is_same<RealScalar,float>::value)
  {
    // let's build a matrix more stable to inverse
    MatrixType a = MatrixType::Random(size,size*2);
    m1 += a * a.adjoint();
  }

  HouseholderQR<MatrixType> qr(m1);
  m3 = MatrixType::Random(size,size);
  m2 = qr.solve(m3);
  VERIFY_IS_APPROX(m3, m1*m2);

  // now construct a matrix with prescribed determinant
  m1.setZero();
  for(int i = 0; i < size; i++) m1(i,i) = internal::random<Scalar>();
  RealScalar absdet = internal::abs(m1.diagonal().prod());
  m3 = qr.householderQ(); // get a unitary
  m1 = m3 * m1 * m3;
  qr.compute(m1);
  VERIFY_IS_APPROX(absdet, qr.absDeterminant());
  VERIFY_IS_APPROX(internal::log(absdet), qr.logAbsDeterminant());
}

template<typename MatrixType> void qr_verify_assert()
{
  MatrixType tmp;

  HouseholderQR<MatrixType> qr;
  VERIFY_RAISES_ASSERT(qr.matrixQR())
  VERIFY_RAISES_ASSERT(qr.solve(tmp))
  VERIFY_RAISES_ASSERT(qr.householderQ())
  VERIFY_RAISES_ASSERT(qr.absDeterminant())
  VERIFY_RAISES_ASSERT(qr.logAbsDeterminant())
}

void test_qr()
{
  for(int i = 0; i < g_repeat; i++) {
   CALL_SUBTEST_1( qr(MatrixXf(internal::random<int>(1,200),internal::random<int>(1,200))) );
   CALL_SUBTEST_2( qr(MatrixXcd(internal::random<int>(1,200),internal::random<int>(1,200))) );
   CALL_SUBTEST_3(( qr_fixedsize<Matrix<float,3,4>, 2 >() ));
   CALL_SUBTEST_4(( qr_fixedsize<Matrix<double,6,2>, 4 >() ));
   CALL_SUBTEST_5(( qr_fixedsize<Matrix<double,2,5>, 7 >() ));
   CALL_SUBTEST_11( qr(Matrix<float,1,1>()) );
  }

  for(int i = 0; i < g_repeat; i++) {
    CALL_SUBTEST_1( qr_invertible<MatrixXf>() );
    CALL_SUBTEST_6( qr_invertible<MatrixXd>() );
    CALL_SUBTEST_7( qr_invertible<MatrixXcf>() );
    CALL_SUBTEST_8( qr_invertible<MatrixXcd>() );
  }

  CALL_SUBTEST_9(qr_verify_assert<Matrix3f>());
  CALL_SUBTEST_10(qr_verify_assert<Matrix3d>());
  CALL_SUBTEST_1(qr_verify_assert<MatrixXf>());
  CALL_SUBTEST_6(qr_verify_assert<MatrixXd>());
  CALL_SUBTEST_7(qr_verify_assert<MatrixXcf>());
  CALL_SUBTEST_8(qr_verify_assert<MatrixXcd>());

  // Test problem size constructors
  CALL_SUBTEST_12(HouseholderQR<MatrixXf>(10, 20));
}
