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
#include <Eigen/Array>

template<typename MatrixType> void array(const MatrixType& m)
{
  /* this test covers the following files:
     Array.cpp
  */

  typedef typename MatrixType::Scalar Scalar;
  typedef typename NumTraits<Scalar>::Real RealScalar;
  typedef Matrix<Scalar, MatrixType::RowsAtCompileTime, 1> VectorType;

  int rows = m.rows();
  int cols = m.cols();

  MatrixType m1 = MatrixType::Random(rows, cols),
             m2 = MatrixType::Random(rows, cols),
             m3(rows, cols);

  Scalar  s1 = ei_random<Scalar>(),
          s2 = ei_random<Scalar>();

  // scalar addition
  VERIFY_IS_APPROX(m1.cwise() + s1, s1 + m1.cwise());
  VERIFY_IS_APPROX(m1.cwise() + s1, MatrixType::Constant(rows,cols,s1) + m1);
  VERIFY_IS_APPROX((m1*Scalar(2)).cwise() - s2, (m1+m1) - MatrixType::Constant(rows,cols,s2) );
  m3 = m1;
  m3.cwise() += s2;
  VERIFY_IS_APPROX(m3, m1.cwise() + s2);
  m3 = m1;
  m3.cwise() -= s1;
  VERIFY_IS_APPROX(m3, m1.cwise() - s1);

  // reductions
  VERIFY_IS_APPROX(m1.colwise().sum().sum(), m1.sum());
  VERIFY_IS_APPROX(m1.rowwise().sum().sum(), m1.sum());
  if (!ei_isApprox(m1.sum(), (m1+m2).sum()))
    VERIFY_IS_NOT_APPROX(((m1+m2).rowwise().sum()).sum(), m1.sum());
  VERIFY_IS_APPROX(m1.colwise().sum(), m1.colwise().redux(internal::scalar_sum_op<Scalar>()));
}

template<typename MatrixType> void comparisons(const MatrixType& m)
{
  typedef typename MatrixType::Scalar Scalar;
  typedef typename NumTraits<Scalar>::Real RealScalar;
  typedef Matrix<Scalar, MatrixType::RowsAtCompileTime, 1> VectorType;

  int rows = m.rows();
  int cols = m.cols();

  int r = ei_random<int>(0, rows-1),
      c = ei_random<int>(0, cols-1);

  MatrixType m1 = MatrixType::Random(rows, cols),
             m2 = MatrixType::Random(rows, cols),
             m3(rows, cols);

  VERIFY(((m1.cwise() + Scalar(1)).cwise() > m1).all());
  VERIFY(((m1.cwise() - Scalar(1)).cwise() < m1).all());
  if (rows*cols>1)
  {
    m3 = m1;
    m3(r,c) += 1;
    VERIFY(! (m1.cwise() < m3).all() );
    VERIFY(! (m1.cwise() > m3).all() );
  }
  
  // comparisons to scalar
  VERIFY( (m1.cwise() != (m1(r,c)+1) ).any() );
  VERIFY( (m1.cwise() > (m1(r,c)-1) ).any() );
  VERIFY( (m1.cwise() < (m1(r,c)+1) ).any() );
  VERIFY( (m1.cwise() == m1(r,c) ).any() );
  
  // test Select
  VERIFY_IS_APPROX( (m1.cwise()<m2).select(m1,m2), m1.cwise().min(m2) );
  VERIFY_IS_APPROX( (m1.cwise()>m2).select(m1,m2), m1.cwise().max(m2) );
  Scalar mid = (m1.cwise().abs().minCoeff() + m1.cwise().abs().maxCoeff())/Scalar(2);
  for (int j=0; j<cols; ++j)
  for (int i=0; i<rows; ++i)
    m3(i,j) = ei_abs(m1(i,j))<mid ? 0 : m1(i,j);
  VERIFY_IS_APPROX( (m1.cwise().abs().cwise()<MatrixType::Constant(rows,cols,mid))
                        .select(MatrixType::Zero(rows,cols),m1), m3);
  // shorter versions:
  VERIFY_IS_APPROX( (m1.cwise().abs().cwise()<MatrixType::Constant(rows,cols,mid))
                        .select(0,m1), m3);
  VERIFY_IS_APPROX( (m1.cwise().abs().cwise()>=MatrixType::Constant(rows,cols,mid))
                        .select(m1,0), m3);
  // even shorter version:
  VERIFY_IS_APPROX( (m1.cwise().abs().cwise()<mid).select(0,m1), m3);
  
  // count
  VERIFY(((m1.cwise().abs().cwise()+1).cwise()>RealScalar(0.1)).count() == rows*cols);
  VERIFY_IS_APPROX(((m1.cwise().abs().cwise()+1).cwise()>RealScalar(0.1)).colwise().count().template cast<int>(), RowVectorXi::Constant(cols,rows));
  VERIFY_IS_APPROX(((m1.cwise().abs().cwise()+1).cwise()>RealScalar(0.1)).rowwise().count().template cast<int>(), VectorXi::Constant(rows, cols));
}

template<typename VectorType> void lpNorm(const VectorType& v)
{
  VectorType u = VectorType::Random(v.size());

  VERIFY_IS_APPROX(u.template lpNorm<Infinity>(), u.cwise().abs().maxCoeff());
  VERIFY_IS_APPROX(u.template lpNorm<1>(), u.cwise().abs().sum());
  VERIFY_IS_APPROX(u.template lpNorm<2>(), ei_sqrt(u.cwise().abs().cwise().square().sum()));
  VERIFY_IS_APPROX(ei_pow(u.template lpNorm<5>(), typename VectorType::RealScalar(5)), u.cwise().abs().cwise().pow(5).sum());
}

void test_eigen2_array()
{
  for(int i = 0; i < g_repeat; i++) {
    CALL_SUBTEST_1( array(Matrix<float, 1, 1>()) );
    CALL_SUBTEST_2( array(Matrix2f()) );
    CALL_SUBTEST_3( array(Matrix4d()) );
    CALL_SUBTEST_4( array(MatrixXcf(3, 3)) );
    CALL_SUBTEST_5( array(MatrixXf(8, 12)) );
    CALL_SUBTEST_6( array(MatrixXi(8, 12)) );
  }
  for(int i = 0; i < g_repeat; i++) {
    CALL_SUBTEST_1( comparisons(Matrix<float, 1, 1>()) );
    CALL_SUBTEST_2( comparisons(Matrix2f()) );
    CALL_SUBTEST_3( comparisons(Matrix4d()) );
    CALL_SUBTEST_5( comparisons(MatrixXf(8, 12)) );
    CALL_SUBTEST_6( comparisons(MatrixXi(8, 12)) );
  }
  for(int i = 0; i < g_repeat; i++) {
    CALL_SUBTEST_1( lpNorm(Matrix<float, 1, 1>()) );
    CALL_SUBTEST_2( lpNorm(Vector2f()) );
    CALL_SUBTEST_3( lpNorm(Vector3d()) );
    CALL_SUBTEST_4( lpNorm(Vector4f()) );
    CALL_SUBTEST_5( lpNorm(VectorXf(16)) );
    CALL_SUBTEST_7( lpNorm(VectorXcd(10)) );
  }
}
