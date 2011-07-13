// This file is part of Eigen, a lightweight C++ template library
// for linear algebra. Eigen itself is part of the KDE project.
//
// Copyright (C) 2008 Gael Guennebaud <g.gael@free.fr>
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

// this hack is needed to make this file compiles with -pedantic (gcc)
#ifdef __GNUC__
#define throw(X)
#endif
// discard stack allocation as that too bypasses malloc
#define EIGEN_STACK_ALLOCATION_LIMIT 0
// any heap allocation will raise an assert
#define EIGEN_NO_MALLOC

#include "main.h"

template<typename MatrixType> void nomalloc(const MatrixType& m)
{
  /* this test check no dynamic memory allocation are issued with fixed-size matrices
  */

  typedef typename MatrixType::Scalar Scalar;
  typedef Matrix<Scalar, MatrixType::RowsAtCompileTime, 1> VectorType;

  int rows = m.rows();
  int cols = m.cols();

  MatrixType m1 = MatrixType::Random(rows, cols),
             m2 = MatrixType::Random(rows, cols),
             m3(rows, cols),
             mzero = MatrixType::Zero(rows, cols),
             identity = Matrix<Scalar, MatrixType::RowsAtCompileTime, MatrixType::RowsAtCompileTime>
                              ::Identity(rows, rows),
             square = Matrix<Scalar, MatrixType::RowsAtCompileTime, MatrixType::RowsAtCompileTime>
                              ::Random(rows, rows);
  VectorType v1 = VectorType::Random(rows),
             v2 = VectorType::Random(rows),
             vzero = VectorType::Zero(rows);

  Scalar s1 = ei_random<Scalar>();

  int r = ei_random<int>(0, rows-1),
      c = ei_random<int>(0, cols-1);

  VERIFY_IS_APPROX((m1+m2)*s1,              s1*m1+s1*m2);
  VERIFY_IS_APPROX((m1+m2)(r,c), (m1(r,c))+(m2(r,c)));
  VERIFY_IS_APPROX(m1.cwise() * m1.block(0,0,rows,cols), m1.cwise() * m1);
  VERIFY_IS_APPROX((m1*m1.transpose())*m2,  m1*(m1.transpose()*m2));
}

void test_eigen2_nomalloc()
{
  // check that our operator new is indeed called:
  VERIFY_RAISES_ASSERT(MatrixXd dummy = MatrixXd::Random(3,3));
  CALL_SUBTEST_1( nomalloc(Matrix<float, 1, 1>()) );
  CALL_SUBTEST_2( nomalloc(Matrix4d()) );
  CALL_SUBTEST_3( nomalloc(Matrix<float,32,32>()) );
}
