// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2008 Gael Guennebaud <gael.guennebaud@inria.fr>
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

#ifdef HAS_GSL
#include "gsl_helper.h"
#endif

template<typename MatrixType> void eigensolver(const MatrixType& m)
{
  typedef typename MatrixType::Index Index;
  /* this test covers the following files:
     EigenSolver.h
  */
  Index rows = m.rows();
  Index cols = m.cols();

  typedef typename MatrixType::Scalar Scalar;
  typedef typename NumTraits<Scalar>::Real RealScalar;
  typedef Matrix<Scalar, MatrixType::RowsAtCompileTime, 1> VectorType;
  typedef Matrix<RealScalar, MatrixType::RowsAtCompileTime, 1> RealVectorType;
  typedef typename std::complex<typename NumTraits<typename MatrixType::Scalar>::Real> Complex;

  MatrixType a = MatrixType::Random(rows,cols);
  MatrixType a1 = MatrixType::Random(rows,cols);
  MatrixType symmA =  a.adjoint() * a + a1.adjoint() * a1;

  EigenSolver<MatrixType> ei0(symmA);
  VERIFY_IS_EQUAL(ei0.info(), Success);
  VERIFY_IS_APPROX(symmA * ei0.pseudoEigenvectors(), ei0.pseudoEigenvectors() * ei0.pseudoEigenvalueMatrix());
  VERIFY_IS_APPROX((symmA.template cast<Complex>()) * (ei0.pseudoEigenvectors().template cast<Complex>()),
    (ei0.pseudoEigenvectors().template cast<Complex>()) * (ei0.eigenvalues().asDiagonal()));

  EigenSolver<MatrixType> ei1(a);
  VERIFY_IS_EQUAL(ei1.info(), Success);
  VERIFY_IS_APPROX(a * ei1.pseudoEigenvectors(), ei1.pseudoEigenvectors() * ei1.pseudoEigenvalueMatrix());
  VERIFY_IS_APPROX(a.template cast<Complex>() * ei1.eigenvectors(),
                   ei1.eigenvectors() * ei1.eigenvalues().asDiagonal());
  VERIFY_IS_APPROX(a.eigenvalues(), ei1.eigenvalues());

  EigenSolver<MatrixType> eiNoEivecs(a, false);
  VERIFY_IS_EQUAL(eiNoEivecs.info(), Success);
  VERIFY_IS_APPROX(ei1.eigenvalues(), eiNoEivecs.eigenvalues());
  VERIFY_IS_APPROX(ei1.pseudoEigenvalueMatrix(), eiNoEivecs.pseudoEigenvalueMatrix());

  MatrixType id = MatrixType::Identity(rows, cols);
  VERIFY_IS_APPROX(id.operatorNorm(), RealScalar(1));

  if (rows > 2)
  {
    // Test matrix with NaN
    a(0,0) = std::numeric_limits<typename MatrixType::RealScalar>::quiet_NaN();
    EigenSolver<MatrixType> eiNaN(a);
    VERIFY_IS_EQUAL(eiNaN.info(), NoConvergence);
  }
}

template<typename MatrixType> void eigensolver_verify_assert(const MatrixType& m)
{
  EigenSolver<MatrixType> eig;
  VERIFY_RAISES_ASSERT(eig.eigenvectors());
  VERIFY_RAISES_ASSERT(eig.pseudoEigenvectors());
  VERIFY_RAISES_ASSERT(eig.pseudoEigenvalueMatrix());
  VERIFY_RAISES_ASSERT(eig.eigenvalues());

  MatrixType a = MatrixType::Random(m.rows(),m.cols());
  eig.compute(a, false);
  VERIFY_RAISES_ASSERT(eig.eigenvectors());
  VERIFY_RAISES_ASSERT(eig.pseudoEigenvectors());
}

void test_eigensolver_generic()
{
  for(int i = 0; i < g_repeat; i++) {
    CALL_SUBTEST_1( eigensolver(Matrix4f()) );
    CALL_SUBTEST_2( eigensolver(MatrixXd(17,17)) );

    // some trivial but implementation-wise tricky cases
    CALL_SUBTEST_2( eigensolver(MatrixXd(1,1)) );
    CALL_SUBTEST_2( eigensolver(MatrixXd(2,2)) );
    CALL_SUBTEST_3( eigensolver(Matrix<double,1,1>()) );
    CALL_SUBTEST_4( eigensolver(Matrix2d()) );
  }

  CALL_SUBTEST_1( eigensolver_verify_assert(Matrix4f()) );
  CALL_SUBTEST_2( eigensolver_verify_assert(MatrixXd(17,17)) );
  CALL_SUBTEST_3( eigensolver_verify_assert(Matrix<double,1,1>()) );
  CALL_SUBTEST_4( eigensolver_verify_assert(Matrix2d()) );

  // Test problem size constructors
  CALL_SUBTEST_5(EigenSolver<MatrixXf>(10));
}
