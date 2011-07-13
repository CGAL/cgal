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

template<typename MatrixType> void selfadjointeigensolver(const MatrixType& m)
{
  typedef typename MatrixType::Index Index;
  /* this test covers the following files:
     EigenSolver.h, SelfAdjointEigenSolver.h (and indirectly: Tridiagonalization.h)
  */
  Index rows = m.rows();
  Index cols = m.cols();

  typedef typename MatrixType::Scalar Scalar;
  typedef typename NumTraits<Scalar>::Real RealScalar;
  typedef Matrix<Scalar, MatrixType::RowsAtCompileTime, 1> VectorType;
  typedef Matrix<RealScalar, MatrixType::RowsAtCompileTime, 1> RealVectorType;
  typedef typename std::complex<typename NumTraits<typename MatrixType::Scalar>::Real> Complex;

  RealScalar largerEps = 10*test_precision<RealScalar>();

  MatrixType a = MatrixType::Random(rows,cols);
  MatrixType a1 = MatrixType::Random(rows,cols);
  MatrixType symmA =  a.adjoint() * a + a1.adjoint() * a1;
  symmA.template triangularView<StrictlyUpper>().setZero();

  MatrixType b = MatrixType::Random(rows,cols);
  MatrixType b1 = MatrixType::Random(rows,cols);
  MatrixType symmB = b.adjoint() * b + b1.adjoint() * b1;
  symmB.template triangularView<StrictlyUpper>().setZero();

  SelfAdjointEigenSolver<MatrixType> eiSymm(symmA);
  // generalized eigen pb
  GeneralizedSelfAdjointEigenSolver<MatrixType> eiSymmGen(symmA, symmB);

  #ifdef HAS_GSL
  if (internal::is_same<RealScalar,double>::value)
  {
    // restore symmA and symmB.
    symmA = MatrixType(symmA.template selfadjointView<Lower>());
    symmB = MatrixType(symmB.template selfadjointView<Lower>());
    typedef GslTraits<Scalar> Gsl;
    typename Gsl::Matrix gEvec=0, gSymmA=0, gSymmB=0;
    typename GslTraits<RealScalar>::Vector gEval=0;
    RealVectorType _eval;
    MatrixType _evec;
    convert<MatrixType>(symmA, gSymmA);
    convert<MatrixType>(symmB, gSymmB);
    convert<MatrixType>(symmA, gEvec);
    gEval = GslTraits<RealScalar>::createVector(rows);

    Gsl::eigen_symm(gSymmA, gEval, gEvec);
    convert(gEval, _eval);
    convert(gEvec, _evec);

    // test gsl itself !
    VERIFY((symmA * _evec).isApprox(_evec * _eval.asDiagonal(), largerEps));

    // compare with eigen
    VERIFY_IS_APPROX(_eval, eiSymm.eigenvalues());
    VERIFY_IS_APPROX(_evec.cwiseAbs(), eiSymm.eigenvectors().cwiseAbs());

    // generalized pb
    Gsl::eigen_symm_gen(gSymmA, gSymmB, gEval, gEvec);
    convert(gEval, _eval);
    convert(gEvec, _evec);
    // test GSL itself:
    VERIFY((symmA * _evec).isApprox(symmB * (_evec * _eval.asDiagonal()), largerEps));

    // compare with eigen
    MatrixType normalized_eivec = eiSymmGen.eigenvectors()*eiSymmGen.eigenvectors().colwise().norm().asDiagonal().inverse();
    VERIFY_IS_APPROX(_eval, eiSymmGen.eigenvalues());
    VERIFY_IS_APPROX(_evec.cwiseAbs(), normalized_eivec.cwiseAbs());

    Gsl::free(gSymmA);
    Gsl::free(gSymmB);
    GslTraits<RealScalar>::free(gEval);
    Gsl::free(gEvec);
  }
  #endif

  VERIFY_IS_EQUAL(eiSymm.info(), Success);
  VERIFY((symmA.template selfadjointView<Lower>() * eiSymm.eigenvectors()).isApprox(
          eiSymm.eigenvectors() * eiSymm.eigenvalues().asDiagonal(), largerEps));
  VERIFY_IS_APPROX(symmA.template selfadjointView<Lower>().eigenvalues(), eiSymm.eigenvalues());

  SelfAdjointEigenSolver<MatrixType> eiSymmNoEivecs(symmA, false);
  VERIFY_IS_EQUAL(eiSymmNoEivecs.info(), Success);
  VERIFY_IS_APPROX(eiSymm.eigenvalues(), eiSymmNoEivecs.eigenvalues());

  // generalized eigen problem Ax = lBx
  eiSymmGen.compute(symmA, symmB,Ax_lBx);
  VERIFY_IS_EQUAL(eiSymmGen.info(), Success);
  VERIFY((symmA.template selfadjointView<Lower>() * eiSymmGen.eigenvectors()).isApprox(
          symmB.template selfadjointView<Lower>() * (eiSymmGen.eigenvectors() * eiSymmGen.eigenvalues().asDiagonal()), largerEps));

  // generalized eigen problem BAx = lx
  eiSymmGen.compute(symmA, symmB,BAx_lx);
  VERIFY_IS_EQUAL(eiSymmGen.info(), Success);
  VERIFY((symmB.template selfadjointView<Lower>() * (symmA.template selfadjointView<Lower>() * eiSymmGen.eigenvectors())).isApprox(
         (eiSymmGen.eigenvectors() * eiSymmGen.eigenvalues().asDiagonal()), largerEps));

  // generalized eigen problem ABx = lx
  eiSymmGen.compute(symmA, symmB,ABx_lx);
  VERIFY_IS_EQUAL(eiSymmGen.info(), Success);
  VERIFY((symmA.template selfadjointView<Lower>() * (symmB.template selfadjointView<Lower>() * eiSymmGen.eigenvectors())).isApprox(
         (eiSymmGen.eigenvectors() * eiSymmGen.eigenvalues().asDiagonal()), largerEps));


  MatrixType sqrtSymmA = eiSymm.operatorSqrt();
  VERIFY_IS_APPROX(MatrixType(symmA.template selfadjointView<Lower>()), sqrtSymmA*sqrtSymmA);
  VERIFY_IS_APPROX(sqrtSymmA, symmA.template selfadjointView<Lower>()*eiSymm.operatorInverseSqrt());

  MatrixType id = MatrixType::Identity(rows, cols);
  VERIFY_IS_APPROX(id.template selfadjointView<Lower>().operatorNorm(), RealScalar(1));

  SelfAdjointEigenSolver<MatrixType> eiSymmUninitialized;
  VERIFY_RAISES_ASSERT(eiSymmUninitialized.info());
  VERIFY_RAISES_ASSERT(eiSymmUninitialized.eigenvalues());
  VERIFY_RAISES_ASSERT(eiSymmUninitialized.eigenvectors());
  VERIFY_RAISES_ASSERT(eiSymmUninitialized.operatorSqrt());
  VERIFY_RAISES_ASSERT(eiSymmUninitialized.operatorInverseSqrt());

  eiSymmUninitialized.compute(symmA, false);
  VERIFY_RAISES_ASSERT(eiSymmUninitialized.eigenvectors());
  VERIFY_RAISES_ASSERT(eiSymmUninitialized.operatorSqrt());
  VERIFY_RAISES_ASSERT(eiSymmUninitialized.operatorInverseSqrt());

  // test Tridiagonalization's methods
  Tridiagonalization<MatrixType> tridiag(symmA);
  // FIXME tridiag.matrixQ().adjoint() does not work
  VERIFY_IS_APPROX(MatrixType(symmA.template selfadjointView<Lower>()), tridiag.matrixQ() * tridiag.matrixT().eval() * MatrixType(tridiag.matrixQ()).adjoint());
  
  if (rows > 1)
  {
    // Test matrix with NaN
    symmA(0,0) = std::numeric_limits<typename MatrixType::RealScalar>::quiet_NaN();
    SelfAdjointEigenSolver<MatrixType> eiSymmNaN(symmA);
    VERIFY_IS_EQUAL(eiSymmNaN.info(), NoConvergence);
  }
}

void test_eigensolver_selfadjoint()
{
  for(int i = 0; i < g_repeat; i++) {
    // very important to test a 3x3 matrix since we provide a special path for it
    CALL_SUBTEST_1( selfadjointeigensolver(Matrix3f()) );
    CALL_SUBTEST_2( selfadjointeigensolver(Matrix4d()) );
    CALL_SUBTEST_3( selfadjointeigensolver(MatrixXf(10,10)) );
    CALL_SUBTEST_4( selfadjointeigensolver(MatrixXd(19,19)) );
    CALL_SUBTEST_5( selfadjointeigensolver(MatrixXcd(17,17)) );
    
    CALL_SUBTEST_9( selfadjointeigensolver(Matrix<std::complex<double>,Dynamic,Dynamic,RowMajor>(17,17)) );

    // some trivial but implementation-wise tricky cases
    CALL_SUBTEST_4( selfadjointeigensolver(MatrixXd(1,1)) );
    CALL_SUBTEST_4( selfadjointeigensolver(MatrixXd(2,2)) );
    CALL_SUBTEST_6( selfadjointeigensolver(Matrix<double,1,1>()) );
    CALL_SUBTEST_7( selfadjointeigensolver(Matrix<double,2,2>()) );
  }

  // Test problem size constructors
  CALL_SUBTEST_8(SelfAdjointEigenSolver<MatrixXf>(10));
  CALL_SUBTEST_8(Tridiagonalization<MatrixXf>(10));
}

