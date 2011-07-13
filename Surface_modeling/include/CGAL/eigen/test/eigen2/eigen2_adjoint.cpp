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

template<typename MatrixType> void adjoint(const MatrixType& m)
{
  /* this test covers the following files:
     Transpose.h Conjugate.h Dot.h
  */

  typedef typename MatrixType::Scalar Scalar;
  typedef typename NumTraits<Scalar>::Real RealScalar;
  typedef Matrix<Scalar, MatrixType::RowsAtCompileTime, 1> VectorType;
  typedef Matrix<Scalar, MatrixType::RowsAtCompileTime, MatrixType::RowsAtCompileTime> SquareMatrixType;
  int rows = m.rows();
  int cols = m.cols();

  RealScalar largerEps = test_precision<RealScalar>();
  if (ei_is_same_type<RealScalar,float>::ret)
    largerEps = RealScalar(1e-3f);

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

  Scalar s1 = ei_random<Scalar>(),
         s2 = ei_random<Scalar>();

  // check basic compatibility of adjoint, transpose, conjugate
  VERIFY_IS_APPROX(m1.transpose().conjugate().adjoint(),    m1);
  VERIFY_IS_APPROX(m1.adjoint().conjugate().transpose(),    m1);

  // check multiplicative behavior
  VERIFY_IS_APPROX((m1.adjoint() * m2).adjoint(),           m2.adjoint() * m1);
  VERIFY_IS_APPROX((s1 * m1).adjoint(),                     ei_conj(s1) * m1.adjoint());

  // check basic properties of dot, norm, norm2
  typedef typename NumTraits<Scalar>::Real RealScalar;
  VERIFY(ei_isApprox((s1 * v1 + s2 * v2).eigen2_dot(v3),   s1 * v1.eigen2_dot(v3) + s2 * v2.eigen2_dot(v3), largerEps));
  VERIFY(ei_isApprox(v3.eigen2_dot(s1 * v1 + s2 * v2),     ei_conj(s1)*v3.eigen2_dot(v1)+ei_conj(s2)*v3.eigen2_dot(v2), largerEps));
  VERIFY_IS_APPROX(ei_conj(v1.eigen2_dot(v2)),               v2.eigen2_dot(v1));
  VERIFY_IS_APPROX(ei_real(v1.eigen2_dot(v1)),               v1.squaredNorm());
  if(NumTraits<Scalar>::HasFloatingPoint)
    VERIFY_IS_APPROX(v1.squaredNorm(),                      v1.norm() * v1.norm());
  VERIFY_IS_MUCH_SMALLER_THAN(ei_abs(vzero.eigen2_dot(v1)),  static_cast<RealScalar>(1));
  if(NumTraits<Scalar>::HasFloatingPoint)
    VERIFY_IS_MUCH_SMALLER_THAN(vzero.norm(),         static_cast<RealScalar>(1));

  // check compatibility of dot and adjoint
  VERIFY(ei_isApprox(v1.eigen2_dot(square * v2), (square.adjoint() * v1).eigen2_dot(v2), largerEps));

  // like in testBasicStuff, test operator() to check const-qualification
  int r = ei_random<int>(0, rows-1),
      c = ei_random<int>(0, cols-1);
  VERIFY_IS_APPROX(m1.conjugate()(r,c), ei_conj(m1(r,c)));
  VERIFY_IS_APPROX(m1.adjoint()(c,r), ei_conj(m1(r,c)));

  if(NumTraits<Scalar>::HasFloatingPoint)
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
  
}

void test_eigen2_adjoint()
{
  for(int i = 0; i < g_repeat; i++) {
    CALL_SUBTEST_1( adjoint(Matrix<float, 1, 1>()) );
    CALL_SUBTEST_2( adjoint(Matrix3d()) );
    CALL_SUBTEST_3( adjoint(Matrix4f()) );
    CALL_SUBTEST_4( adjoint(MatrixXcf(4, 4)) );
    CALL_SUBTEST_5( adjoint(MatrixXi(8, 12)) );
    CALL_SUBTEST_6( adjoint(MatrixXf(21, 21)) );
  }
  // test a large matrix only once
  CALL_SUBTEST_7( adjoint(Matrix<float, 100, 100>()) );
}

