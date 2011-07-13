// This file is part of Eigen, a lightweight C++ template library
// for linear algebra. Eigen itself is part of the KDE project.
//
// Copyright (C) 2008-2009 Gael Guennebaud <gael.guennebaud@inria.fr>
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
#include <Eigen/LU>

/* Check that two column vectors are approximately equal upto permutations,
   by checking that the k-th power sums are equal for k = 1, ..., vec1.rows() */
template<typename VectorType>
void verify_is_approx_upto_permutation(const VectorType& vec1, const VectorType& vec2)
{
  typedef typename NumTraits<typename VectorType::Scalar>::Real RealScalar;

  VERIFY(vec1.cols() == 1);
  VERIFY(vec2.cols() == 1);
  VERIFY(vec1.rows() == vec2.rows());
  for (int k = 1; k <= vec1.rows(); ++k)
  {
    VERIFY_IS_APPROX(vec1.array().pow(RealScalar(k)).sum(), vec2.array().pow(RealScalar(k)).sum());
  }
}


template<typename MatrixType> void eigensolver(const MatrixType& m)
{
  typedef typename MatrixType::Index Index;
  /* this test covers the following files:
     ComplexEigenSolver.h, and indirectly ComplexSchur.h
  */
  Index rows = m.rows();
  Index cols = m.cols();

  typedef typename MatrixType::Scalar Scalar;
  typedef typename NumTraits<Scalar>::Real RealScalar;
  typedef Matrix<Scalar, MatrixType::RowsAtCompileTime, 1> VectorType;
  typedef Matrix<RealScalar, MatrixType::RowsAtCompileTime, 1> RealVectorType;
  typedef typename std::complex<typename NumTraits<typename MatrixType::Scalar>::Real> Complex;

  MatrixType a = MatrixType::Random(rows,cols);
  MatrixType symmA =  a.adjoint() * a;

  ComplexEigenSolver<MatrixType> ei0(symmA);
  VERIFY_IS_EQUAL(ei0.info(), Success);
  VERIFY_IS_APPROX(symmA * ei0.eigenvectors(), ei0.eigenvectors() * ei0.eigenvalues().asDiagonal());

  ComplexEigenSolver<MatrixType> ei1(a);
  VERIFY_IS_EQUAL(ei1.info(), Success);
  VERIFY_IS_APPROX(a * ei1.eigenvectors(), ei1.eigenvectors() * ei1.eigenvalues().asDiagonal());
  // Note: If MatrixType is real then a.eigenvalues() uses EigenSolver and thus
  // another algorithm so results may differ slightly
  verify_is_approx_upto_permutation(a.eigenvalues(), ei1.eigenvalues());

  ComplexEigenSolver<MatrixType> eiNoEivecs(a, false);
  VERIFY_IS_EQUAL(eiNoEivecs.info(), Success);
  VERIFY_IS_APPROX(ei1.eigenvalues(), eiNoEivecs.eigenvalues());

  // Regression test for issue #66
  MatrixType z = MatrixType::Zero(rows,cols);
  ComplexEigenSolver<MatrixType> eiz(z);
  VERIFY((eiz.eigenvalues().cwiseEqual(0)).all());

  MatrixType id = MatrixType::Identity(rows, cols);
  VERIFY_IS_APPROX(id.operatorNorm(), RealScalar(1));

  if (rows > 1)
  {
    // Test matrix with NaN
    a(0,0) = std::numeric_limits<typename MatrixType::RealScalar>::quiet_NaN();
    ComplexEigenSolver<MatrixType> eiNaN(a);
    VERIFY_IS_EQUAL(eiNaN.info(), NoConvergence);
  }
}

template<typename MatrixType> void eigensolver_verify_assert(const MatrixType& m)
{
  ComplexEigenSolver<MatrixType> eig;
  VERIFY_RAISES_ASSERT(eig.eigenvectors());
  VERIFY_RAISES_ASSERT(eig.eigenvalues());

  MatrixType a = MatrixType::Random(m.rows(),m.cols());
  eig.compute(a, false);
  VERIFY_RAISES_ASSERT(eig.eigenvectors());
}

void test_eigensolver_complex()
{
  for(int i = 0; i < g_repeat; i++) {
    CALL_SUBTEST_1( eigensolver(Matrix4cf()) );
    CALL_SUBTEST_2( eigensolver(MatrixXcd(14,14)) );
    CALL_SUBTEST_3( eigensolver(Matrix<std::complex<float>, 1, 1>()) );
    CALL_SUBTEST_4( eigensolver(Matrix3f()) );
  }

  CALL_SUBTEST_1( eigensolver_verify_assert(Matrix4cf()) );
  CALL_SUBTEST_2( eigensolver_verify_assert(MatrixXcd(14,14)) );
  CALL_SUBTEST_3( eigensolver_verify_assert(Matrix<std::complex<float>, 1, 1>()) );
  CALL_SUBTEST_4( eigensolver_verify_assert(Matrix3f()) );

  // Test problem size constructors
  CALL_SUBTEST_5(ComplexEigenSolver<MatrixXf>(10));
}
