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

#define EIGEN_NO_ASSERTION_CHECKING
#include "main.h"
#include <Eigen/Cholesky>
#include <Eigen/LU>

#ifdef HAS_GSL
#include "gsl_helper.h"
#endif

template<typename MatrixType> void cholesky(const MatrixType& m)
{
  /* this test covers the following files:
     LLT.h LDLT.h
  */
  int rows = m.rows();
  int cols = m.cols();

  typedef typename MatrixType::Scalar Scalar;
  typedef typename NumTraits<Scalar>::Real RealScalar;
  typedef Matrix<Scalar, MatrixType::RowsAtCompileTime, MatrixType::RowsAtCompileTime> SquareMatrixType;
  typedef Matrix<Scalar, MatrixType::RowsAtCompileTime, 1> VectorType;

  MatrixType a0 = MatrixType::Random(rows,cols);
  VectorType vecB = VectorType::Random(rows), vecX(rows);
  MatrixType matB = MatrixType::Random(rows,cols), matX(rows,cols);
  SquareMatrixType symm =  a0 * a0.adjoint();
  // let's make sure the matrix is not singular or near singular
  MatrixType a1 = MatrixType::Random(rows,cols);
  symm += a1 * a1.adjoint();

  #ifdef HAS_GSL
  if (ei_is_same_type<RealScalar,double>::ret)
  {
    typedef GslTraits<Scalar> Gsl;
    typename Gsl::Matrix gMatA=0, gSymm=0;
    typename Gsl::Vector gVecB=0, gVecX=0;
    convert<MatrixType>(symm, gSymm);
    convert<MatrixType>(symm, gMatA);
    convert<VectorType>(vecB, gVecB);
    convert<VectorType>(vecB, gVecX);
    Gsl::cholesky(gMatA);
    Gsl::cholesky_solve(gMatA, gVecB, gVecX);
    VectorType vecX(rows), _vecX, _vecB;
    convert(gVecX, _vecX);
    symm.llt().solve(vecB, &vecX);
    Gsl::prod(gSymm, gVecX, gVecB);
    convert(gVecB, _vecB);
    // test gsl itself !
    VERIFY_IS_APPROX(vecB, _vecB);
    VERIFY_IS_APPROX(vecX, _vecX);

    Gsl::free(gMatA);
    Gsl::free(gSymm);
    Gsl::free(gVecB);
    Gsl::free(gVecX);
  }
  #endif

  {
    LDLT<SquareMatrixType> ldlt(symm);
    VERIFY(ldlt.isPositiveDefinite());
    // in eigen3, LDLT is pivoting
    //VERIFY_IS_APPROX(symm, ldlt.matrixL() * ldlt.vectorD().asDiagonal() * ldlt.matrixL().adjoint());
    ldlt.solve(vecB, &vecX);
    VERIFY_IS_APPROX(symm * vecX, vecB);
    ldlt.solve(matB, &matX);
    VERIFY_IS_APPROX(symm * matX, matB);
  }

  {
    LLT<SquareMatrixType> chol(symm);
    VERIFY(chol.isPositiveDefinite());
    VERIFY_IS_APPROX(symm, chol.matrixL() * chol.matrixL().adjoint());
    chol.solve(vecB, &vecX);
    VERIFY_IS_APPROX(symm * vecX, vecB);
    chol.solve(matB, &matX);
    VERIFY_IS_APPROX(symm * matX, matB);
  }

#if 0 // cholesky is not rank-revealing anyway
  // test isPositiveDefinite on non definite matrix
  if (rows>4)
  {
    SquareMatrixType symm =  a0.block(0,0,rows,cols-4) * a0.block(0,0,rows,cols-4).adjoint();
    LLT<SquareMatrixType> chol(symm);
    VERIFY(!chol.isPositiveDefinite());
    LDLT<SquareMatrixType> cholnosqrt(symm);
    VERIFY(!cholnosqrt.isPositiveDefinite());
  }
#endif
}

void test_eigen2_cholesky()
{
  for(int i = 0; i < g_repeat; i++) {
    CALL_SUBTEST_1( cholesky(Matrix<double,1,1>()) );
    CALL_SUBTEST_2( cholesky(Matrix2d()) );
    CALL_SUBTEST_3( cholesky(Matrix3f()) );
    CALL_SUBTEST_4( cholesky(Matrix4d()) );
    CALL_SUBTEST_5( cholesky(MatrixXcd(7,7)) );
    CALL_SUBTEST_6( cholesky(MatrixXf(17,17)) );
    CALL_SUBTEST_7( cholesky(MatrixXd(33,33)) );
  }
}
