// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
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

template<typename MatrixType> void linearStructure(const MatrixType& m)
{
  /* this test covers the following files:
     CwiseUnaryOp.h, CwiseBinaryOp.h, SelfCwiseBinaryOp.h 
  */
  typedef typename MatrixType::Index Index;
  typedef typename MatrixType::Scalar Scalar;

  Index rows = m.rows();
  Index cols = m.cols();

  // this test relies a lot on Random.h, and there's not much more that we can do
  // to test it, hence I consider that we will have tested Random.h
  MatrixType m1 = MatrixType::Random(rows, cols),
             m2 = MatrixType::Random(rows, cols),
             m3(rows, cols),
             mzero = MatrixType::Zero(rows, cols);

  Scalar s1 = internal::random<Scalar>();
  while (internal::abs(s1)<1e-3) s1 = internal::random<Scalar>();

  Index r = internal::random<Index>(0, rows-1),
        c = internal::random<Index>(0, cols-1);

  VERIFY_IS_APPROX(-(-m1),                  m1);
  VERIFY_IS_APPROX(m1+m1,                   2*m1);
  VERIFY_IS_APPROX(m1+m2-m1,                m2);
  VERIFY_IS_APPROX(-m2+m1+m2,               m1);
  VERIFY_IS_APPROX(m1*s1,                   s1*m1);
  VERIFY_IS_APPROX((m1+m2)*s1,              s1*m1+s1*m2);
  VERIFY_IS_APPROX((-m1+m2)*s1,             -s1*m1+s1*m2);
  m3 = m2; m3 += m1;
  VERIFY_IS_APPROX(m3,                      m1+m2);
  m3 = m2; m3 -= m1;
  VERIFY_IS_APPROX(m3,                      m2-m1);
  m3 = m2; m3 *= s1;
  VERIFY_IS_APPROX(m3,                      s1*m2);
  if(!NumTraits<Scalar>::IsInteger)
  {
    m3 = m2; m3 /= s1;
    VERIFY_IS_APPROX(m3,                    m2/s1);
  }

  // again, test operator() to check const-qualification
  VERIFY_IS_APPROX((-m1)(r,c), -(m1(r,c)));
  VERIFY_IS_APPROX((m1-m2)(r,c), (m1(r,c))-(m2(r,c)));
  VERIFY_IS_APPROX((m1+m2)(r,c), (m1(r,c))+(m2(r,c)));
  VERIFY_IS_APPROX((s1*m1)(r,c), s1*(m1(r,c)));
  VERIFY_IS_APPROX((m1*s1)(r,c), (m1(r,c))*s1);
  if(!NumTraits<Scalar>::IsInteger)
    VERIFY_IS_APPROX((m1/s1)(r,c), (m1(r,c))/s1);

  // use .block to disable vectorization and compare to the vectorized version
  VERIFY_IS_APPROX(m1+m1.block(0,0,rows,cols), m1+m1);
  VERIFY_IS_APPROX(m1.cwiseProduct(m1.block(0,0,rows,cols)), m1.cwiseProduct(m1));
  VERIFY_IS_APPROX(m1 - m1.block(0,0,rows,cols), m1 - m1);
  VERIFY_IS_APPROX(m1.block(0,0,rows,cols) * s1, m1 * s1);
}

void test_linearstructure()
{
  for(int i = 0; i < g_repeat; i++) {
    CALL_SUBTEST_1( linearStructure(Matrix<float, 1, 1>()) );
    CALL_SUBTEST_2( linearStructure(Matrix2f()) );
    CALL_SUBTEST_3( linearStructure(Vector3d()) );
    CALL_SUBTEST_4( linearStructure(Matrix4d()) );
    CALL_SUBTEST_5( linearStructure(MatrixXcf(3, 3)) );
    CALL_SUBTEST_6( linearStructure(MatrixXf(8, 12)) );
    CALL_SUBTEST_7( linearStructure(MatrixXi(8, 12)) );
    CALL_SUBTEST_8( linearStructure(MatrixXcd(20, 20)) );
    CALL_SUBTEST_9( linearStructure(ArrayXXf(12, 8)) );
  }
}
