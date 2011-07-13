// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2011 Benoit Jacob <jacob.benoit.1@gmail.com>
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

template<typename MatrixType> void zeroSizedMatrix()
{
  MatrixType t1;

  if (MatrixType::SizeAtCompileTime == Dynamic)
  {
    if (MatrixType::RowsAtCompileTime == Dynamic)
      VERIFY(t1.rows() == 0);
    if (MatrixType::ColsAtCompileTime == Dynamic)
      VERIFY(t1.cols() == 0);

    if (MatrixType::RowsAtCompileTime == Dynamic && MatrixType::ColsAtCompileTime == Dynamic)
    {
      MatrixType t2(0, 0);
      VERIFY(t2.rows() == 0);
      VERIFY(t2.cols() == 0);
    }
  }
}

template<typename VectorType> void zeroSizedVector()
{
  VectorType t1;

  if (VectorType::SizeAtCompileTime == Dynamic)
  {
    VERIFY(t1.size() == 0);
    VectorType t2(DenseIndex(0)); // DenseIndex disambiguates with 0-the-null-pointer (error with gcc 4.4 and MSVC8)
    VERIFY(t2.size() == 0);
  }
}

void test_zerosized()
{
  zeroSizedMatrix<Matrix2d>();
  zeroSizedMatrix<Matrix3i>();
  zeroSizedMatrix<Matrix<float, 2, Dynamic> >();
  zeroSizedMatrix<MatrixXf>();
  
  zeroSizedVector<Vector2d>();
  zeroSizedVector<Vector3i>();
  zeroSizedVector<VectorXf>();
}
