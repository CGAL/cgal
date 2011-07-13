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

template<typename MatrixType> void basicStuff(const MatrixType& m)
{
  typedef typename MatrixType::Scalar Scalar;
  typedef Matrix<Scalar, MatrixType::RowsAtCompileTime, 1> VectorType;

  int rows = m.rows();
  int cols = m.cols();

  // this test relies a lot on Random.h, and there's not much more that we can do
  // to test it, hence I consider that we will have tested Random.h
  MatrixType m1 = MatrixType::Random(rows, cols),
             m2 = MatrixType::Random(rows, cols),
             m3(rows, cols),
             mzero = MatrixType::Zero(rows, cols),
             identity = Matrix<Scalar, MatrixType::RowsAtCompileTime, MatrixType::RowsAtCompileTime>
                              ::Identity(rows, rows),
             square = Matrix<Scalar, MatrixType::RowsAtCompileTime, MatrixType::RowsAtCompileTime>::Random(rows, rows);
  VectorType v1 = VectorType::Random(rows),
             v2 = VectorType::Random(rows),
             vzero = VectorType::Zero(rows);

  Scalar x = ei_random<Scalar>();

  int r = ei_random<int>(0, rows-1),
      c = ei_random<int>(0, cols-1);

  m1.coeffRef(r,c) = x;
  VERIFY_IS_APPROX(x, m1.coeff(r,c));
  m1(r,c) = x;
  VERIFY_IS_APPROX(x, m1(r,c));
  v1.coeffRef(r) = x;
  VERIFY_IS_APPROX(x, v1.coeff(r));
  v1(r) = x;
  VERIFY_IS_APPROX(x, v1(r));
  v1[r] = x;
  VERIFY_IS_APPROX(x, v1[r]);

  VERIFY_IS_APPROX(               v1,    v1);
  VERIFY_IS_NOT_APPROX(           v1,    2*v1);
  VERIFY_IS_MUCH_SMALLER_THAN(    vzero, v1);
  if(NumTraits<Scalar>::HasFloatingPoint)
    VERIFY_IS_MUCH_SMALLER_THAN(  vzero, v1.norm());
  VERIFY_IS_NOT_MUCH_SMALLER_THAN(v1,    v1);
  VERIFY_IS_APPROX(               vzero, v1-v1);
  VERIFY_IS_APPROX(               m1,    m1);
  VERIFY_IS_NOT_APPROX(           m1,    2*m1);
  VERIFY_IS_MUCH_SMALLER_THAN(    mzero, m1);
  VERIFY_IS_NOT_MUCH_SMALLER_THAN(m1,    m1);
  VERIFY_IS_APPROX(               mzero, m1-m1);

  // always test operator() on each read-only expression class,
  // in order to check const-qualifiers.
  // indeed, if an expression class (here Zero) is meant to be read-only,
  // hence has no _write() method, the corresponding MatrixBase method (here zero())
  // should return a const-qualified object so that it is the const-qualified
  // operator() that gets called, which in turn calls _read().
  VERIFY_IS_MUCH_SMALLER_THAN(MatrixType::Zero(rows,cols)(r,c), static_cast<Scalar>(1));

  // now test copying a row-vector into a (column-)vector and conversely.
  square.col(r) = square.row(r).eval();
  Matrix<Scalar, 1, MatrixType::RowsAtCompileTime> rv(rows);
  Matrix<Scalar, MatrixType::RowsAtCompileTime, 1> cv(rows);
  rv = square.row(r);
  cv = square.col(r);
  VERIFY_IS_APPROX(rv, cv.transpose());

  if(cols!=1 && rows!=1 && MatrixType::SizeAtCompileTime!=Dynamic)
  {
    VERIFY_RAISES_ASSERT(m1 = (m2.block(0,0, rows-1, cols-1)));
  }

  VERIFY_IS_APPROX(m3 = m1,m1);
  MatrixType m4;
  VERIFY_IS_APPROX(m4 = m1,m1);

  // test swap
  m3 = m1;
  m1.swap(m2);
  VERIFY_IS_APPROX(m3, m2);
  if(rows*cols>=3)
  {
    VERIFY_IS_NOT_APPROX(m3, m1);
  }
}

void test_eigen2_basicstuff()
{
  for(int i = 0; i < g_repeat; i++) {
    CALL_SUBTEST_1( basicStuff(Matrix<float, 1, 1>()) );
    CALL_SUBTEST_2( basicStuff(Matrix4d()) );
    CALL_SUBTEST_3( basicStuff(MatrixXcf(3, 3)) );
    CALL_SUBTEST_4( basicStuff(MatrixXi(8, 12)) );
    CALL_SUBTEST_5( basicStuff(MatrixXcd(20, 20)) );
    CALL_SUBTEST_6( basicStuff(Matrix<float, 100, 100>()) );
    CALL_SUBTEST_7( basicStuff(Matrix<long double,Dynamic,Dynamic>(10,10)) );
  }
}
