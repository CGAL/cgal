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

void test_commainitializer()
{
  Matrix3d m3;
  Matrix4d m4;

  VERIFY_RAISES_ASSERT( (m3 << 1, 2, 3, 4, 5, 6, 7, 8) );
  
  #ifndef _MSC_VER
  VERIFY_RAISES_ASSERT( (m3 << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10) );
  #endif

  double data[] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  Matrix3d ref = Map<Matrix<double,3,3,RowMajor> >(data);

  m3 = Matrix3d::Random();
  m3 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  VERIFY_IS_APPROX(m3, ref );

  Vector3d vec[3];
  vec[0] << 1, 4, 7;
  vec[1] << 2, 5, 8;
  vec[2] << 3, 6, 9;
  m3 = Matrix3d::Random();
  m3 << vec[0], vec[1], vec[2];
  VERIFY_IS_APPROX(m3, ref);

  vec[0] << 1, 2, 3;
  vec[1] << 4, 5, 6;
  vec[2] << 7, 8, 9;
  m3 = Matrix3d::Random();
  m3 << vec[0].transpose(),
        4, 5, 6,
        vec[2].transpose();
  VERIFY_IS_APPROX(m3, ref);
}
