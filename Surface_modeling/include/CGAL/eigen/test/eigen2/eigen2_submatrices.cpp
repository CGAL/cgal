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

// check minor separately in order to avoid the possible creation of a zero-sized
// array. Comes from a compilation error with gcc-3.4 or gcc-4 with -ansi -pedantic.
// Another solution would be to declare the array like this: T m_data[Size==0?1:Size]; in ei_matrix_storage
// but this is probably not bad to raise such an error at compile time...
template<typename Scalar, int _Rows, int _Cols> struct CheckMinor
{
    typedef Matrix<Scalar, _Rows, _Cols> MatrixType;
    CheckMinor(MatrixType& m1, int r1, int c1)
    {
        int rows = m1.rows();
        int cols = m1.cols();

        Matrix<Scalar, Dynamic, Dynamic> mi = m1.minor(0,0).eval();
        VERIFY_IS_APPROX(mi, m1.block(1,1,rows-1,cols-1));
        mi = m1.minor(r1,c1);
        VERIFY_IS_APPROX(mi.transpose(), m1.transpose().minor(c1,r1));
        //check operator(), both constant and non-constant, on minor()
        m1.minor(r1,c1)(0,0) = m1.minor(0,0)(0,0);
    }
};

template<typename Scalar> struct CheckMinor<Scalar,1,1>
{
    typedef Matrix<Scalar, 1, 1> MatrixType;
    CheckMinor(MatrixType&, int, int) {}
};

template<typename MatrixType> void submatrices(const MatrixType& m)
{
  /* this test covers the following files:
     Row.h Column.h Block.h Minor.h DiagonalCoeffs.h
  */
  typedef typename MatrixType::Scalar Scalar;
  typedef typename MatrixType::RealScalar RealScalar;
  typedef Matrix<Scalar, MatrixType::RowsAtCompileTime, 1> VectorType;
  typedef Matrix<Scalar, 1, MatrixType::ColsAtCompileTime> RowVectorType;
  int rows = m.rows();
  int cols = m.cols();

  MatrixType m1 = MatrixType::Random(rows, cols),
             m2 = MatrixType::Random(rows, cols),
             m3(rows, cols),
             mzero = MatrixType::Zero(rows, cols),
             ones = MatrixType::Ones(rows, cols),
             identity = Matrix<Scalar, MatrixType::RowsAtCompileTime, MatrixType::RowsAtCompileTime>
                              ::Identity(rows, rows),
             square = Matrix<Scalar, MatrixType::RowsAtCompileTime, MatrixType::RowsAtCompileTime>
                              ::Random(rows, rows);
  VectorType v1 = VectorType::Random(rows),
             v2 = VectorType::Random(rows),
             v3 = VectorType::Random(rows),
             vzero = VectorType::Zero(rows);

  Scalar s1 = ei_random<Scalar>();

  int r1 = ei_random<int>(0,rows-1);
  int r2 = ei_random<int>(r1,rows-1);
  int c1 = ei_random<int>(0,cols-1);
  int c2 = ei_random<int>(c1,cols-1);

  //check row() and col()
  VERIFY_IS_APPROX(m1.col(c1).transpose(), m1.transpose().row(c1));
  VERIFY_IS_APPROX(square.row(r1).eigen2_dot(m1.col(c1)), (square.lazy() * m1.conjugate())(r1,c1));
  //check operator(), both constant and non-constant, on row() and col()
  m1.row(r1) += s1 * m1.row(r2);
  m1.col(c1) += s1 * m1.col(c2);

  //check block()
  Matrix<Scalar,Dynamic,Dynamic> b1(1,1); b1(0,0) = m1(r1,c1);
  RowVectorType br1(m1.block(r1,0,1,cols));
  VectorType bc1(m1.block(0,c1,rows,1));
  VERIFY_IS_APPROX(b1, m1.block(r1,c1,1,1));
  VERIFY_IS_APPROX(m1.row(r1), br1);
  VERIFY_IS_APPROX(m1.col(c1), bc1);
  //check operator(), both constant and non-constant, on block()
  m1.block(r1,c1,r2-r1+1,c2-c1+1) = s1 * m2.block(0, 0, r2-r1+1,c2-c1+1);
  m1.block(r1,c1,r2-r1+1,c2-c1+1)(r2-r1,c2-c1) = m2.block(0, 0, r2-r1+1,c2-c1+1)(0,0);

  //check minor()
  CheckMinor<Scalar, MatrixType::RowsAtCompileTime, MatrixType::ColsAtCompileTime> checkminor(m1,r1,c1);

  //check diagonal()
  VERIFY_IS_APPROX(m1.diagonal(), m1.transpose().diagonal());
  m2.diagonal() = 2 * m1.diagonal();
  m2.diagonal()[0] *= 3;
  VERIFY_IS_APPROX(m2.diagonal()[0], static_cast<Scalar>(6) * m1.diagonal()[0]);

  enum {
    BlockRows = EIGEN_SIZE_MIN_PREFER_FIXED(MatrixType::RowsAtCompileTime,2),
    BlockCols = EIGEN_SIZE_MIN_PREFER_FIXED(MatrixType::ColsAtCompileTime,5)
  };
  if (rows>=5 && cols>=8)
  {
    // test fixed block() as lvalue
    m1.template block<BlockRows,BlockCols>(1,1) *= s1;
    // test operator() on fixed block() both as constant and non-constant
    m1.template block<BlockRows,BlockCols>(1,1)(0, 3) = m1.template block<2,5>(1,1)(1,2);
    // check that fixed block() and block() agree
    Matrix<Scalar,Dynamic,Dynamic> b = m1.template block<BlockRows,BlockCols>(3,3);
    VERIFY_IS_APPROX(b, m1.block(3,3,BlockRows,BlockCols));
  }

  if (rows>2)
  {
    // test sub vectors
    VERIFY_IS_APPROX(v1.template start<2>(), v1.block(0,0,2,1));
    VERIFY_IS_APPROX(v1.template start<2>(), v1.start(2));
    VERIFY_IS_APPROX(v1.template start<2>(), v1.segment(0,2));
    VERIFY_IS_APPROX(v1.template start<2>(), v1.template segment<2>(0));
    int i = rows-2;
    VERIFY_IS_APPROX(v1.template end<2>(), v1.block(i,0,2,1));
    VERIFY_IS_APPROX(v1.template end<2>(), v1.end(2));
    VERIFY_IS_APPROX(v1.template end<2>(), v1.segment(i,2));
    VERIFY_IS_APPROX(v1.template end<2>(), v1.template segment<2>(i));
    i = ei_random(0,rows-2);
    VERIFY_IS_APPROX(v1.segment(i,2), v1.template segment<2>(i));
  }

  // stress some basic stuffs with block matrices
  VERIFY(ei_real(ones.col(c1).sum()) == RealScalar(rows));
  VERIFY(ei_real(ones.row(r1).sum()) == RealScalar(cols));

  VERIFY(ei_real(ones.col(c1).eigen2_dot(ones.col(c2))) == RealScalar(rows));
  VERIFY(ei_real(ones.row(r1).eigen2_dot(ones.row(r2))) == RealScalar(cols));
}

void test_eigen2_submatrices()
{
  for(int i = 0; i < g_repeat; i++) {
    CALL_SUBTEST_1( submatrices(Matrix<float, 1, 1>()) );
    CALL_SUBTEST_2( submatrices(Matrix4d()) );
    CALL_SUBTEST_3( submatrices(MatrixXcf(3, 3)) );
    CALL_SUBTEST_4( submatrices(MatrixXi(8, 12)) );
    CALL_SUBTEST_5( submatrices(MatrixXcd(20, 20)) );
    CALL_SUBTEST_6( submatrices(MatrixXf(20, 20)) );
  }
}
