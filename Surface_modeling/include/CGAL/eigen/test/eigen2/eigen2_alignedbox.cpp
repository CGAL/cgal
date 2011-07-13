// This file is part of Eigen, a lightweight C++ template library
// for linear algebra. Eigen itself is part of the KDE project.
//
// Copyright (C) 2008 Gael Guennebaud <g.gael@free.fr>
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
#include <Eigen/Geometry>
#include <Eigen/LU>
#include <Eigen/QR>

template<typename BoxType> void alignedbox(const BoxType& _box)
{
  /* this test covers the following files:
     AlignedBox.h
  */

  const int dim = _box.dim();
  typedef typename BoxType::Scalar Scalar;
  typedef typename NumTraits<Scalar>::Real RealScalar;
  typedef Matrix<Scalar, BoxType::AmbientDimAtCompileTime, 1> VectorType;

  VectorType p0 = VectorType::Random(dim);
  VectorType p1 = VectorType::Random(dim);
  RealScalar s1 = ei_random<RealScalar>(0,1);

  BoxType b0(dim);
  BoxType b1(VectorType::Random(dim),VectorType::Random(dim));
  BoxType b2;

  b0.extend(p0);
  b0.extend(p1);
  VERIFY(b0.contains(p0*s1+(Scalar(1)-s1)*p1));
  VERIFY(!b0.contains(p0 + (1+s1)*(p1-p0)));

  (b2 = b0).extend(b1);
  VERIFY(b2.contains(b0));
  VERIFY(b2.contains(b1));
  VERIFY_IS_APPROX(b2.clamp(b0), b0);

  // casting
  const int Dim = BoxType::AmbientDimAtCompileTime;
  typedef typename GetDifferentType<Scalar>::type OtherScalar;
  AlignedBox<OtherScalar,Dim> hp1f = b0.template cast<OtherScalar>();
  VERIFY_IS_APPROX(hp1f.template cast<Scalar>(),b0);
  AlignedBox<Scalar,Dim> hp1d = b0.template cast<Scalar>();
  VERIFY_IS_APPROX(hp1d.template cast<Scalar>(),b0);
}

void test_eigen2_alignedbox()
{
  for(int i = 0; i < g_repeat; i++) {
    CALL_SUBTEST_1( alignedbox(AlignedBox<float,2>()) );
    CALL_SUBTEST_2( alignedbox(AlignedBox<float,3>()) );
    CALL_SUBTEST_3( alignedbox(AlignedBox<double,4>()) );
  }
}
