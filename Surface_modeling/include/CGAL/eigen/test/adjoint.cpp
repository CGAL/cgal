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

#define EIGEN_NO_STATIC_ASSERT

#include "main.h"

template<typename MatrixType> void adjoint(const MatrixType& m)
{
  /* this test covers the following files:
     Transpose.h Conjugate.h Dot.h
  */
  typedef typename MatrixType::Index Index;
  typedef typename MatrixType::Scalar Scalar;
  typedef typename NumTraits<Scalar>::Real RealScalar;
  typedef Matrix<Scalar, MatrixType::RowsAtCompileTime, 1> VectorType;
  typedef Matrix<Scalar, MatrixType::RowsAtCompileTime, MatrixType::RowsAtCompileTime> SquareMatrixType;
  
  Index rows = m.rows();
  Index cols = m.cols();

  MatrixType m1 = MatrixType::Random(rows, cols),
             m2 = MatrixType::Random(rows, cols),
             m3(rows, cols),
             mzero = MatrixType::Zero(rows, cols),
             identity = SquareMatrixType::Identity(rows, rows),
             square = SquareMatrixType::Random(rows, rows);
  VectorType v1 = VectorType::Random(rows),
             v2 = VectorType::Random(rows),
             v3 = VectorType::Random(rows),
             vzero = VectorType::Zero(rows);

  Scalar s1 = internal::random<Scalar>(),
         s2 = internal::random<Scalar>();

  // check basic compatibility of adjoint, transpose, conjugate
  VERIFY_IS_APPROX(m1.transpose().conjugate().adjoint(),    m1);
  VERIFY_IS_APPROX(m1.adjoint().conjugate().transpose(),    m1);

  // check multiplicative behavior
  VERIFY_IS_APPROX((m1.adjoint() * m2).adjoint(),           m2.adjoint() * m1);
  VERIFY_IS_APPROX((s1 * m1).adjoint(),                     internal::conj(s1) * m1.adjoint());

  // check basic properties of dot, norm, norm2
  typedef typename NumTraits<Scalar>::Real RealScalar;
  
  RealScalar ref = NumTraits<Scalar>::IsInteger ? 0 : std::max((s1 * v1 + s2 * v2).norm(),v3.norm());
  VERIFY(test_isApproxWithRef((s1 * v1 + s2 * v2).dot(v3),     internal::conj(s1) * v1.dot(v3) + internal::conj(s2) * v2.dot(v3), ref));
  VERIFY(test_isApproxWithRef(v3.dot(s1 * v1 + s2 * v2),       s1*v3.dot(v1)+s2*v3.dot(v2), ref));
  VERIFY_IS_APPROX(internal::conj(v1.dot(v2)),               v2.dot(v1));
  VERIFY_IS_APPROX(internal::real(v1.dot(v1)),                v1.squaredNorm());
  if(!NumTraits<Scalar>::IsInteger)
    VERIFY_IS_APPROX(v1.squaredNorm(),                v1.norm() * v1.norm());
  VERIFY_IS_MUCH_SMALLER_THAN(internal::abs(vzero.dot(v1)),  static_cast<RealScalar>(1));

  // check compatibility of dot and adjoint
  
  ref = NumTraits<Scalar>::IsInteger ? 0 : std::max(std::max(v1.norm(),v2.norm()),std::max((square * v2).norm(),(square.adjoint() * v1).norm()));
  VERIFY(test_isApproxWithRef(v1.dot(square * v2), (square.adjoint() * v1).dot(v2), ref));

  // like in testBasicStuff, test operator() to check const-qualification
  Index r = internal::random<Index>(0, rows-1),
      c = internal::random<Index>(0, cols-1);
  VERIFY_IS_APPROX(m1.conjugate()(r,c), internal::conj(m1(r,c)));
  VERIFY_IS_APPROX(m1.adjoint()(c,r), internal::conj(m1(r,c)));

  if(!NumTraits<Scalar>::IsInteger)
  {
    // check that Random().normalized() works: tricky as the random xpr must be evaluated by
    // normalized() in order to produce a consistent result.
    VERIFY_IS_APPROX(VectorType::Random(rows).normalized().norm(), RealScalar(1));
  }

  // check inplace transpose
  m3 = m1;
  m3.transposeInPlace();
  VERIFY_IS_APPROX(m3,m1.transpose());
  m3.transposeInPlace();
  VERIFY_IS_APPROX(m3,m1);

  // check inplace adjoint
  m3 = m1;
  m3.adjointInPlace();
  VERIFY_IS_APPROX(m3,m1.adjoint());
  m3.transposeInPlace();
  VERIFY_IS_APPROX(m3,m1.conjugate());

  // check mixed dot product
  typedef Matrix<RealScalar, MatrixType::RowsAtCompileTime, 1> RealVectorType;
  RealVectorType rv1 = RealVectorType::Random(rows);
  VERIFY_IS_APPROX(v1.dot(rv1.template cast<Scalar>()), v1.dot(rv1));
  VERIFY_IS_APPROX(rv1.template cast<Scalar>().dot(v1), rv1.dot(v1));
}

void test_adjoint()
{
  for(int i = 0; i < g_repeat; i++) {
    CALL_SUBTEST_1( adjoint(Matrix<float, 1, 1>()) );
    CALL_SUBTEST_2( adjoint(Matrix3d()) );
    CALL_SUBTEST_3( adjoint(Matrix4f()) );
    CALL_SUBTEST_4( adjoint(MatrixXcf(internal::random<int>(1,50), internal::random<int>(1,50))) );
    CALL_SUBTEST_5( adjoint(MatrixXi(internal::random<int>(1,50), internal::random<int>(1,50))) );
    CALL_SUBTEST_6( adjoint(MatrixXf(internal::random<int>(1,50), internal::random<int>(1,50))) );
  }
  // test a large matrix only once
  CALL_SUBTEST_7( adjoint(Matrix<float, 100, 100>()) );

#ifdef EIGEN_TEST_PART_4
  {
    MatrixXcf a(10,10), b(10,10);
    VERIFY_RAISES_ASSERT(a = a.transpose());
    VERIFY_RAISES_ASSERT(a = a.transpose() + b);
    VERIFY_RAISES_ASSERT(a = b + a.transpose());
    VERIFY_RAISES_ASSERT(a = a.conjugate().transpose());
    VERIFY_RAISES_ASSERT(a = a.adjoint());
    VERIFY_RAISES_ASSERT(a = a.adjoint() + b);
    VERIFY_RAISES_ASSERT(a = b + a.adjoint());

    // no assertion should be triggered for these cases:
    a.transpose() = a.transpose();
    a.transpose() += a.transpose();
    a.transpose() += a.transpose() + b;
    a.transpose() = a.adjoint();
    a.transpose() += a.adjoint();
    a.transpose() += a.adjoint() + b;
  }
#endif
}

