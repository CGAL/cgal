// This file is part of Eigen, a lightweight C++ template library
// for linear algebra. Eigen itself is part of the KDE project.
//
// Copyright (C) 2008 Gael Guennebaud <gael.guennebaud@gmail.com>
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

template<typename MatrixType> void triangular(const MatrixType& m)
{
  typedef typename MatrixType::Scalar Scalar;
  typedef typename NumTraits<Scalar>::Real RealScalar;
  typedef Matrix<Scalar, MatrixType::RowsAtCompileTime, 1> VectorType;

  RealScalar largerEps = 10*test_precision<RealScalar>();

  int rows = m.rows();
  int cols = m.cols();

  MatrixType m1 = MatrixType::Random(rows, cols),
             m2 = MatrixType::Random(rows, cols),
             m3(rows, cols),
             m4(rows, cols),
             r1(rows, cols),
             r2(rows, cols),
             mzero = MatrixType::Zero(rows, cols),
             mones = MatrixType::Ones(rows, cols),
             identity = Matrix<Scalar, MatrixType::RowsAtCompileTime, MatrixType::RowsAtCompileTime>
                              ::Identity(rows, rows),
             square = Matrix<Scalar, MatrixType::RowsAtCompileTime, MatrixType::RowsAtCompileTime>
                              ::Random(rows, rows);
  VectorType v1 = VectorType::Random(rows),
             v2 = VectorType::Random(rows),
             vzero = VectorType::Zero(rows);

  MatrixType m1up = m1.template part<Eigen::UpperTriangular>();
  MatrixType m2up = m2.template part<Eigen::UpperTriangular>();

  if (rows*cols>1)
  {
    VERIFY(m1up.isUpperTriangular());
    VERIFY(m2up.transpose().isLowerTriangular());
    VERIFY(!m2.isLowerTriangular());
  }

//   VERIFY_IS_APPROX(m1up.transpose() * m2, m1.upper().transpose().lower() * m2);

  // test overloaded operator+=
  r1.setZero();
  r2.setZero();
  r1.template part<Eigen::UpperTriangular>() +=  m1;
  r2 += m1up;
  VERIFY_IS_APPROX(r1,r2);

  // test overloaded operator=
  m1.setZero();
  m1.template part<Eigen::UpperTriangular>() = (m2.transpose() * m2).lazy();
  m3 = m2.transpose() * m2;
  VERIFY_IS_APPROX(m3.template part<Eigen::LowerTriangular>().transpose(), m1);

  // test overloaded operator=
  m1.setZero();
  m1.template part<Eigen::LowerTriangular>() = (m2.transpose() * m2).lazy();
  VERIFY_IS_APPROX(m3.template part<Eigen::LowerTriangular>(), m1);

  VERIFY_IS_APPROX(m3.template part<Diagonal>(), m3.diagonal().asDiagonal());

  m1 = MatrixType::Random(rows, cols);
  for (int i=0; i<rows; ++i)
    while (ei_abs2(m1(i,i))<1e-3) m1(i,i) = ei_random<Scalar>();

  Transpose<MatrixType> trm4(m4);
  // test back and forward subsitution
  m3 = m1.template part<Eigen::LowerTriangular>();
  VERIFY(m3.template marked<Eigen::LowerTriangular>().solveTriangular(m3).cwise().abs().isIdentity(test_precision<RealScalar>()));
  VERIFY(m3.transpose().template marked<Eigen::UpperTriangular>()
    .solveTriangular(m3.transpose()).cwise().abs().isIdentity(test_precision<RealScalar>()));
  // check M * inv(L) using in place API
  m4 = m3;
  m3.transpose().template marked<Eigen::UpperTriangular>().solveTriangularInPlace(trm4);
  VERIFY(m4.cwise().abs().isIdentity(test_precision<RealScalar>()));

  m3 = m1.template part<Eigen::UpperTriangular>();
  VERIFY(m3.template marked<Eigen::UpperTriangular>().solveTriangular(m3).cwise().abs().isIdentity(test_precision<RealScalar>()));
  VERIFY(m3.transpose().template marked<Eigen::LowerTriangular>()
    .solveTriangular(m3.transpose()).cwise().abs().isIdentity(test_precision<RealScalar>()));
  // check M * inv(U) using in place API
  m4 = m3;
  m3.transpose().template marked<Eigen::LowerTriangular>().solveTriangularInPlace(trm4);
  VERIFY(m4.cwise().abs().isIdentity(test_precision<RealScalar>()));

  m3 = m1.template part<Eigen::UpperTriangular>();
  VERIFY(m2.isApprox(m3 * (m3.template marked<Eigen::UpperTriangular>().solveTriangular(m2)), largerEps));
  m3 = m1.template part<Eigen::LowerTriangular>();
  VERIFY(m2.isApprox(m3 * (m3.template marked<Eigen::LowerTriangular>().solveTriangular(m2)), largerEps));

  VERIFY((m1.template part<Eigen::UpperTriangular>() * m2.template part<Eigen::UpperTriangular>()).isUpperTriangular());

  // test swap
  m1.setOnes();
  m2.setZero();
  m2.template part<Eigen::UpperTriangular>().swap(m1);
  m3.setZero();
  m3.template part<Eigen::UpperTriangular>().setOnes();
  VERIFY_IS_APPROX(m2,m3);

}

void selfadjoint()
{
  Matrix2i m;
  m << 1, 2,
       3, 4;

  Matrix2i m1 = Matrix2i::Zero();
  m1.part<SelfAdjoint>() = m;
  Matrix2i ref1;
  ref1 << 1, 2,
          2, 4;
  VERIFY(m1 == ref1);
  
  Matrix2i m2 = Matrix2i::Zero();
  m2.part<SelfAdjoint>() = m.part<UpperTriangular>();
  Matrix2i ref2;
  ref2 << 1, 2,
          2, 4;
  VERIFY(m2 == ref2);
 
  Matrix2i m3 = Matrix2i::Zero();
  m3.part<SelfAdjoint>() = m.part<LowerTriangular>();
  Matrix2i ref3;
  ref3 << 1, 0,
          0, 4;
  VERIFY(m3 == ref3);
  
  // example inspired from bug 159
  int array[] = {1, 2, 3, 4};
  Matrix2i::Map(array).part<SelfAdjoint>() = Matrix2i::Random().part<LowerTriangular>();
  
  std::cout << "hello\n" << array << std::endl;
}

void test_eigen2_triangular()
{
  CALL_SUBTEST_8( selfadjoint() );
  for(int i = 0; i < g_repeat ; i++) {
    CALL_SUBTEST_1( triangular(Matrix<float, 1, 1>()) );
    CALL_SUBTEST_2( triangular(Matrix<float, 2, 2>()) );
    CALL_SUBTEST_3( triangular(Matrix3d()) );
    CALL_SUBTEST_4( triangular(MatrixXcf(4, 4)) );
    CALL_SUBTEST_5( triangular(Matrix<std::complex<float>,8, 8>()) );
    CALL_SUBTEST_6( triangular(MatrixXd(17,17)) );
    CALL_SUBTEST_7( triangular(Matrix<float,Dynamic,Dynamic,RowMajor>(5, 5)) );
  }
}
