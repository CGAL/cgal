// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2006-2010 Benoit Jacob <jacob.benoit.1@gmail.com>
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

template<typename MatrixType> void diagonal(const MatrixType& m)
{
  typedef typename MatrixType::Index Index;
  typedef typename MatrixType::Scalar Scalar;
  typedef typename MatrixType::RealScalar RealScalar;
  typedef Matrix<Scalar, MatrixType::RowsAtCompileTime, 1> VectorType;
  typedef Matrix<Scalar, 1, MatrixType::ColsAtCompileTime> RowVectorType;

  Index rows = m.rows();
  Index cols = m.cols();

  MatrixType m1 = MatrixType::Random(rows, cols),
             m2 = MatrixType::Random(rows, cols);

  //check diagonal()
  VERIFY_IS_APPROX(m1.diagonal(), m1.transpose().diagonal());
  m2.diagonal() = 2 * m1.diagonal();
  m2.diagonal()[0] *= 3;

  if (rows>2)
  {
    enum {
      N1 = MatrixType::RowsAtCompileTime>1 ?  1 : 0,
      N2 = MatrixType::RowsAtCompileTime>2 ? -2 : 0
    };

    // check sub/super diagonal
    m2.template diagonal<N1>() = 2 * m1.template diagonal<N1>();
    m2.template diagonal<N1>()[0] *= 3;
    VERIFY_IS_APPROX(m2.template diagonal<N1>()[0], static_cast<Scalar>(6) * m1.template diagonal<N1>()[0]);

    m2.template diagonal<N2>() = 2 * m1.template diagonal<N2>();
    m2.template diagonal<N2>()[0] *= 3;
    VERIFY_IS_APPROX(m2.template diagonal<N2>()[0], static_cast<Scalar>(6) * m1.template diagonal<N2>()[0]);

    m2.diagonal(N1) = 2 * m1.diagonal(N1);
    m2.diagonal(N1)[0] *= 3;
    VERIFY_IS_APPROX(m2.diagonal(N1)[0], static_cast<Scalar>(6) * m1.diagonal(N1)[0]);

    m2.diagonal(N2) = 2 * m1.diagonal(N2);
    m2.diagonal(N2)[0] *= 3;
    VERIFY_IS_APPROX(m2.diagonal(N2)[0], static_cast<Scalar>(6) * m1.diagonal(N2)[0]);
  }
}

void test_diagonal()
{
  for(int i = 0; i < g_repeat; i++) {
    CALL_SUBTEST_1( diagonal(Matrix<float, 1, 1>()) );
    CALL_SUBTEST_2( diagonal(Matrix4d()) );
    CALL_SUBTEST_2( diagonal(MatrixXcf(3, 3)) );
    CALL_SUBTEST_2( diagonal(MatrixXi(8, 12)) );
    CALL_SUBTEST_2( diagonal(MatrixXcd(20, 20)) );
    CALL_SUBTEST_1( diagonal(MatrixXf(21, 19)) );
    CALL_SUBTEST_1( diagonal(Matrix<float,Dynamic,4>(3, 4)) );
  }
}
