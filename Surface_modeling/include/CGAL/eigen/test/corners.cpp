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

#define COMPARE_CORNER(A,B) \
  VERIFY_IS_EQUAL(matrix.A, matrix.B); \
  VERIFY_IS_EQUAL(const_matrix.A, const_matrix.B);

template<typename MatrixType> void corners(const MatrixType& m)
{
  typedef typename MatrixType::Index Index;
  Index rows = m.rows();
  Index cols = m.cols();

  Index r = internal::random<Index>(1,rows);
  Index c = internal::random<Index>(1,cols);

  MatrixType matrix = MatrixType::Random(rows,cols);
  const MatrixType const_matrix = MatrixType::Random(rows,cols);

  COMPARE_CORNER(topLeftCorner(r,c), block(0,0,r,c));
  COMPARE_CORNER(topRightCorner(r,c), block(0,cols-c,r,c));
  COMPARE_CORNER(bottomLeftCorner(r,c), block(rows-r,0,r,c));
  COMPARE_CORNER(bottomRightCorner(r,c), block(rows-r,cols-c,r,c));

  Index sr = internal::random<Index>(1,rows) - 1;
  Index nr = internal::random<Index>(1,rows-sr);
  Index sc = internal::random<Index>(1,cols) - 1;
  Index nc = internal::random<Index>(1,cols-sc);

  COMPARE_CORNER(topRows(r), block(0,0,r,cols));
  COMPARE_CORNER(middleRows(sr,nr), block(sr,0,nr,cols));
  COMPARE_CORNER(bottomRows(r), block(rows-r,0,r,cols));
  COMPARE_CORNER(leftCols(c), block(0,0,rows,c));
  COMPARE_CORNER(middleCols(sc,nc), block(0,sc,rows,nc));
  COMPARE_CORNER(rightCols(c), block(0,cols-c,rows,c));
}

template<typename MatrixType, int CRows, int CCols, int SRows, int SCols> void corners_fixedsize()
{
  MatrixType matrix = MatrixType::Random();
  const MatrixType const_matrix = MatrixType::Random();

  enum {
    rows = MatrixType::RowsAtCompileTime,
    cols = MatrixType::ColsAtCompileTime,
    r = CRows,
    c = CCols,
	sr = SRows,
	sc = SCols
  };

  VERIFY_IS_EQUAL((matrix.template topLeftCorner<r,c>()), (matrix.template block<r,c>(0,0)));
  VERIFY_IS_EQUAL((matrix.template topRightCorner<r,c>()), (matrix.template block<r,c>(0,cols-c)));
  VERIFY_IS_EQUAL((matrix.template bottomLeftCorner<r,c>()), (matrix.template block<r,c>(rows-r,0)));
  VERIFY_IS_EQUAL((matrix.template bottomRightCorner<r,c>()), (matrix.template block<r,c>(rows-r,cols-c)));

  VERIFY_IS_EQUAL((matrix.template topRows<r>()), (matrix.template block<r,cols>(0,0)));
  VERIFY_IS_EQUAL((matrix.template middleRows<r>(sr)), (matrix.template block<r,cols>(sr,0)));
  VERIFY_IS_EQUAL((matrix.template bottomRows<r>()), (matrix.template block<r,cols>(rows-r,0)));
  VERIFY_IS_EQUAL((matrix.template leftCols<c>()), (matrix.template block<rows,c>(0,0)));
  VERIFY_IS_EQUAL((matrix.template middleCols<c>(sc)), (matrix.template block<rows,c>(0,sc)));
  VERIFY_IS_EQUAL((matrix.template rightCols<c>()), (matrix.template block<rows,c>(0,cols-c)));

  VERIFY_IS_EQUAL((const_matrix.template topLeftCorner<r,c>()), (const_matrix.template block<r,c>(0,0)));
  VERIFY_IS_EQUAL((const_matrix.template topRightCorner<r,c>()), (const_matrix.template block<r,c>(0,cols-c)));
  VERIFY_IS_EQUAL((const_matrix.template bottomLeftCorner<r,c>()), (const_matrix.template block<r,c>(rows-r,0)));
  VERIFY_IS_EQUAL((const_matrix.template bottomRightCorner<r,c>()), (const_matrix.template block<r,c>(rows-r,cols-c)));

  VERIFY_IS_EQUAL((const_matrix.template topRows<r>()), (const_matrix.template block<r,cols>(0,0)));
  VERIFY_IS_EQUAL((const_matrix.template middleRows<r>(sr)), (const_matrix.template block<r,cols>(sr,0)));
  VERIFY_IS_EQUAL((const_matrix.template bottomRows<r>()), (const_matrix.template block<r,cols>(rows-r,0)));
  VERIFY_IS_EQUAL((const_matrix.template leftCols<c>()), (const_matrix.template block<rows,c>(0,0)));
  VERIFY_IS_EQUAL((const_matrix.template middleCols<c>(sc)), (const_matrix.template block<rows,c>(0,sc)));
  VERIFY_IS_EQUAL((const_matrix.template rightCols<c>()), (const_matrix.template block<rows,c>(0,cols-c)));
}

void test_corners()
{
  for(int i = 0; i < g_repeat; i++) {
    CALL_SUBTEST_1( corners(Matrix<float, 1, 1>()) );
    CALL_SUBTEST_2( corners(Matrix4d()) );
    CALL_SUBTEST_3( corners(Matrix<int,10,12>()) );
    CALL_SUBTEST_4( corners(MatrixXcf(5, 7)) );
    CALL_SUBTEST_5( corners(MatrixXf(21, 20)) );

    CALL_SUBTEST_1(( corners_fixedsize<Matrix<float, 1, 1>, 1, 1, 0, 0>() ));
    CALL_SUBTEST_2(( corners_fixedsize<Matrix4d,2,2,1,1>() ));
    CALL_SUBTEST_3(( corners_fixedsize<Matrix<int,10,12>,4,7,5,2>() ));
  }
}
