// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2009 Keir Mierle <mierle@gmail.com>
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

template<DenseIndex rows, DenseIndex cols>
void resizeLikeTest()
{
  MatrixXf A(rows, cols);
  MatrixXf B;
  Matrix<double, rows, cols> C;
  B.resizeLike(A);
  C.resizeLike(B);  // Shouldn't crash.
  VERIFY(B.rows() == rows && B.cols() == cols);

  VectorXf x(rows);
  RowVectorXf y;
  y.resizeLike(x);
  VERIFY(y.rows() == 1 && y.cols() == rows);

  y.resize(cols);
  x.resizeLike(y);
  VERIFY(x.rows() == cols && x.cols() == 1);
}

void resizeLikeTest12() { resizeLikeTest<1,2>(); }
void resizeLikeTest1020() { resizeLikeTest<10,20>(); }
void resizeLikeTest31() { resizeLikeTest<3,1>(); }

void test_resize()
{
  CALL_SUBTEST(resizeLikeTest12() );
  CALL_SUBTEST(resizeLikeTest1020() );
  CALL_SUBTEST(resizeLikeTest31() );
}
