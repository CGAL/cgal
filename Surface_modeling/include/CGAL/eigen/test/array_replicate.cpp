// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2009 Gael Guennebaud <gael.guennebaud@inria.fr>
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

template<typename MatrixType> void replicate(const MatrixType& m)
{
  /* this test covers the following files:
     Replicate.cpp
  */
  typedef typename MatrixType::Index Index;
  typedef typename MatrixType::Scalar Scalar;
  typedef typename NumTraits<Scalar>::Real RealScalar;
  typedef Matrix<Scalar, MatrixType::RowsAtCompileTime, 1> VectorType;
  typedef Matrix<Scalar, Dynamic, Dynamic> MatrixX;
  typedef Matrix<Scalar, Dynamic, 1> VectorX;

  Index rows = m.rows();
  Index cols = m.cols();

  MatrixType m1 = MatrixType::Random(rows, cols),
             m2 = MatrixType::Random(rows, cols);

  VectorType v1 = VectorType::Random(rows);

  MatrixX x1, x2;
  VectorX vx1;

  int  f1 = internal::random<int>(1,10),
       f2 = internal::random<int>(1,10);

  x1.resize(rows*f1,cols*f2);
  for(int j=0; j<f2; j++)
  for(int i=0; i<f1; i++)
    x1.block(i*rows,j*cols,rows,cols) = m1;
  VERIFY_IS_APPROX(x1, m1.replicate(f1,f2));

  x2.resize(2*rows,3*cols);
  x2 << m2, m2, m2,
        m2, m2, m2;
  VERIFY_IS_APPROX(x2, (m2.template replicate<2,3>()));

  x2.resize(rows,f1);
  for (int j=0; j<f1; ++j)
    x2.col(j) = v1;
  VERIFY_IS_APPROX(x2, v1.rowwise().replicate(f1));

  vx1.resize(rows*f2);
  for (int j=0; j<f2; ++j)
    vx1.segment(j*rows,rows) = v1;
  VERIFY_IS_APPROX(vx1, v1.colwise().replicate(f2));
}

void test_array_replicate()
{
  for(int i = 0; i < g_repeat; i++) {
    CALL_SUBTEST_1( replicate(Matrix<float, 1, 1>()) );
    CALL_SUBTEST_2( replicate(Vector2f()) );
    CALL_SUBTEST_3( replicate(Vector3d()) );
    CALL_SUBTEST_4( replicate(Vector4f()) );
    CALL_SUBTEST_5( replicate(VectorXf(16)) );
    CALL_SUBTEST_6( replicate(VectorXcd(10)) );
  }
}
