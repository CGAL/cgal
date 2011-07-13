// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2008 Gael Guennebaud <gael.guennebaud@inria.fr>
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

template<typename MatrixType> void verifySizeOf(const MatrixType&)
{
  typedef typename MatrixType::Scalar Scalar;
  if (MatrixType::RowsAtCompileTime!=Dynamic && MatrixType::ColsAtCompileTime!=Dynamic)
    VERIFY(sizeof(MatrixType)==sizeof(Scalar)*size_t(MatrixType::SizeAtCompileTime));
  else
    VERIFY(sizeof(MatrixType)==sizeof(Scalar*) + 2 * sizeof(typename MatrixType::Index));
}

void test_sizeof()
{
  CALL_SUBTEST(verifySizeOf(Matrix<float, 1, 1>()) );
  CALL_SUBTEST(verifySizeOf(Matrix4d()) );
  CALL_SUBTEST(verifySizeOf(Matrix<double, 4, 2>()) );
  CALL_SUBTEST(verifySizeOf(Matrix<bool, 7, 5>()) );
  CALL_SUBTEST(verifySizeOf(MatrixXcf(3, 3)) );
  CALL_SUBTEST(verifySizeOf(MatrixXi(8, 12)) );
  CALL_SUBTEST(verifySizeOf(MatrixXcd(20, 20)) );
  CALL_SUBTEST(verifySizeOf(Matrix<float, 100, 100>()) );
  
  VERIFY(sizeof(std::complex<float>) == 2*sizeof(float));
  VERIFY(sizeof(std::complex<double>) == 2*sizeof(double));
}
