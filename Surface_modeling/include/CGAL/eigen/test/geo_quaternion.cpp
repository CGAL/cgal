// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2008-2009 Gael Guennebaud <gael.guennebaud@inria.fr>
// Copyright (C) 2009 Mathieu Gautier <mathieu.gautier@cea.fr>
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
#include <Eigen/SVD>

template<typename Scalar, int Options> void quaternion(void)
{
  /* this test covers the following files:
     Quaternion.h
  */

  typedef Matrix<Scalar,3,3> Matrix3;
  typedef Matrix<Scalar,3,1> Vector3;
  typedef Quaternion<Scalar,Options> Quaternionx;
  typedef AngleAxis<Scalar> AngleAxisx;

  Scalar largeEps = test_precision<Scalar>();
  if (internal::is_same<Scalar,float>::value)
    largeEps = 1e-3f;

  Scalar eps = internal::random<Scalar>() * Scalar(1e-2);

  Vector3 v0 = Vector3::Random(),
          v1 = Vector3::Random(),
          v2 = Vector3::Random(),
          v3 = Vector3::Random();

  Scalar a = internal::random<Scalar>(-Scalar(M_PI), Scalar(M_PI));

  // Quaternion: Identity(), setIdentity();
  Quaternionx q1, q2;
  q2.setIdentity();
  VERIFY_IS_APPROX(Quaternionx(Quaternionx::Identity()).coeffs(), q2.coeffs());
  q1.coeffs().setRandom();
  VERIFY_IS_APPROX(q1.coeffs(), (q1*q2).coeffs());

  // concatenation
  q1 *= q2;

  q1 = AngleAxisx(a, v0.normalized());
  q2 = AngleAxisx(a, v1.normalized());

  // angular distance
  Scalar refangle = internal::abs(AngleAxisx(q1.inverse()*q2).angle());
  if (refangle>Scalar(M_PI))
    refangle = Scalar(2)*Scalar(M_PI) - refangle;

  if((q1.coeffs()-q2.coeffs()).norm() > 10*largeEps)
  {
    VERIFY_IS_MUCH_SMALLER_THAN(internal::abs(q1.angularDistance(q2) - refangle), Scalar(1));
  }

  // rotation matrix conversion
  VERIFY_IS_APPROX(q1 * v2, q1.toRotationMatrix() * v2);
  VERIFY_IS_APPROX(q1 * q2 * v2,
    q1.toRotationMatrix() * q2.toRotationMatrix() * v2);

  VERIFY(  (q2*q1).isApprox(q1*q2, largeEps)
        || !(q2 * q1 * v2).isApprox(q1.toRotationMatrix() * q2.toRotationMatrix() * v2));

  q2 = q1.toRotationMatrix();
  VERIFY_IS_APPROX(q1*v1,q2*v1);


  // angle-axis conversion
  AngleAxisx aa = AngleAxisx(q1);
  VERIFY_IS_APPROX(q1 * v1, Quaternionx(aa) * v1);

  // Do not execute the test if the rotation angle is almost zero, or
  // the rotation axis and v1 are almost parallel.
  if (internal::abs(aa.angle()) > 5*test_precision<Scalar>()
      && (aa.axis() - v1.normalized()).norm() < 1.99
      && (aa.axis() + v1.normalized()).norm() < 1.99) 
  {
    VERIFY_IS_NOT_APPROX(q1 * v1, Quaternionx(AngleAxisx(aa.angle()*2,aa.axis())) * v1);
  }

  // from two vector creation
  VERIFY_IS_APPROX( v2.normalized(),(q2.setFromTwoVectors(v1, v2)*v1).normalized());
  VERIFY_IS_APPROX( v1.normalized(),(q2.setFromTwoVectors(v1, v1)*v1).normalized());
  VERIFY_IS_APPROX(-v1.normalized(),(q2.setFromTwoVectors(v1,-v1)*v1).normalized());
  if (internal::is_same<Scalar,double>::value)
  {
    v3 = (v1.array()+eps).matrix();
    VERIFY_IS_APPROX( v3.normalized(),(q2.setFromTwoVectors(v1, v3)*v1).normalized());
    VERIFY_IS_APPROX(-v3.normalized(),(q2.setFromTwoVectors(v1,-v3)*v1).normalized());
  }

  // inverse and conjugate
  VERIFY_IS_APPROX(q1 * (q1.inverse() * v1), v1);
  VERIFY_IS_APPROX(q1 * (q1.conjugate() * v1), v1);

  // test casting
  Quaternion<float> q1f = q1.template cast<float>();
  VERIFY_IS_APPROX(q1f.template cast<Scalar>(),q1);
  Quaternion<double> q1d = q1.template cast<double>();
  VERIFY_IS_APPROX(q1d.template cast<Scalar>(),q1);
}

template<typename Scalar> void mapQuaternion(void){
  typedef Map<Quaternion<Scalar>, Aligned> MQuaternionA;
  typedef Map<Quaternion<Scalar> > MQuaternionUA;
  typedef Quaternion<Scalar> Quaternionx;

  EIGEN_ALIGN16 Scalar array1[4];
  EIGEN_ALIGN16 Scalar array2[4];
  EIGEN_ALIGN16 Scalar array3[4+1];
  Scalar* array3unaligned = array3+1;

  MQuaternionA(array1).coeffs().setRandom();
  (MQuaternionA(array2)) = MQuaternionA(array1);
  (MQuaternionUA(array3unaligned)) = MQuaternionA(array1);

  Quaternionx q1 = MQuaternionA(array1);
  Quaternionx q2 = MQuaternionA(array2);
  Quaternionx q3 = MQuaternionUA(array3unaligned);

  VERIFY_IS_APPROX(q1.coeffs(), q2.coeffs());
  VERIFY_IS_APPROX(q1.coeffs(), q3.coeffs());
  #ifdef EIGEN_VECTORIZE
  VERIFY_RAISES_ASSERT((MQuaternionA(array3unaligned)));
  #endif
}

template<typename Scalar> void quaternionAlignment(void){
  typedef Quaternion<Scalar,AutoAlign> QuaternionA;
  typedef Quaternion<Scalar,DontAlign> QuaternionUA;

  EIGEN_ALIGN16 Scalar array1[4];
  EIGEN_ALIGN16 Scalar array2[4];
  EIGEN_ALIGN16 Scalar array3[4+1];
  Scalar* arrayunaligned = array3+1;

  QuaternionA *q1 = ::new(reinterpret_cast<void*>(array1)) QuaternionA;
  QuaternionUA *q2 = ::new(reinterpret_cast<void*>(array2)) QuaternionUA;
  QuaternionUA *q3 = ::new(reinterpret_cast<void*>(arrayunaligned)) QuaternionUA;

  q1->coeffs().setRandom();
  *q2 = *q1;
  *q3 = *q1;

  VERIFY_IS_APPROX(q1->coeffs(), q2->coeffs());
  VERIFY_IS_APPROX(q1->coeffs(), q3->coeffs());
  #ifdef EIGEN_VECTORIZE
  VERIFY_RAISES_ASSERT((::new(reinterpret_cast<void*>(arrayunaligned)) QuaternionA));
  #endif
}

void test_geo_quaternion()
{
  for(int i = 0; i < g_repeat; i++) {
    CALL_SUBTEST_1(( quaternion<float,AutoAlign>() ));
    CALL_SUBTEST_2(( quaternion<double,AutoAlign>() ));
    CALL_SUBTEST_3(( quaternion<float,DontAlign>() ));
    CALL_SUBTEST_4(( quaternion<double,DontAlign>() ));
    CALL_SUBTEST_5(( quaternionAlignment<float>() ));
    CALL_SUBTEST_6(( quaternionAlignment<double>() ));
    CALL_SUBTEST( mapQuaternion<float>() );
    CALL_SUBTEST( mapQuaternion<double>() );
  }
}
