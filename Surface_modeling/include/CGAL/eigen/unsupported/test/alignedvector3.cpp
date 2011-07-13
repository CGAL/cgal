// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2009 Gael Guennebaud <g.gael@free.fr>
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
#include <unsupported/Eigen/AlignedVector3>

template<typename Scalar>
void alignedvector3()
{
  Scalar s1 = internal::random<Scalar>();
  Scalar s2 = internal::random<Scalar>();
  typedef Matrix<Scalar,3,1> RefType;
  typedef Matrix<Scalar,3,3> Mat33;
  typedef AlignedVector3<Scalar> FastType;
  RefType  r1(RefType::Random()), r2(RefType::Random()), r3(RefType::Random()),
           r4(RefType::Random()), r5(RefType::Random()), r6(RefType::Random());
  FastType f1(r1), f2(r2), f3(r3), f4(r4), f5(r5), f6(r6);
  Mat33 m1(Mat33::Random());
  
  VERIFY_IS_APPROX(f1,r1);
  VERIFY_IS_APPROX(f4,r4);

  VERIFY_IS_APPROX(f4+f1,r4+r1);
  VERIFY_IS_APPROX(f4-f1,r4-r1);
  VERIFY_IS_APPROX(f4+f1-f2,r4+r1-r2);
  VERIFY_IS_APPROX(f4+=f3,r4+=r3);
  VERIFY_IS_APPROX(f4-=f5,r4-=r5);
  VERIFY_IS_APPROX(f4-=f5+f1,r4-=r5+r1);
  VERIFY_IS_APPROX(f5+f1-s1*f2,r5+r1-s1*r2);
  VERIFY_IS_APPROX(f5+f1/s2-s1*f2,r5+r1/s2-s1*r2);
  
  VERIFY_IS_APPROX(m1*f4,m1*r4);
  VERIFY_IS_APPROX(f4.transpose()*m1,r4.transpose()*m1);
  
  VERIFY_IS_APPROX(f2.dot(f3),r2.dot(r3));
  VERIFY_IS_APPROX(f2.cross(f3),r2.cross(r3));
  VERIFY_IS_APPROX(f2.norm(),r2.norm());

  VERIFY_IS_APPROX(f2.normalized(),r2.normalized());

  VERIFY_IS_APPROX((f2+f1).normalized(),(r2+r1).normalized());
  
  f2.normalize();
  r2.normalize();
  VERIFY_IS_APPROX(f2,r2);
}

void test_alignedvector3()
{
  for(int i = 0; i < g_repeat; i++) {
    CALL_SUBTEST( alignedvector3<float>() );
  }
}
