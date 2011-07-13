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
#include <Eigen/SVD>

template<typename Scalar> void geometry(void)
{
  /* this test covers the following files:
     Cross.h Quaternion.h, Transform.cpp
  */

  typedef Matrix<Scalar,2,2> Matrix2;
  typedef Matrix<Scalar,3,3> Matrix3;
  typedef Matrix<Scalar,4,4> Matrix4;
  typedef Matrix<Scalar,2,1> Vector2;
  typedef Matrix<Scalar,3,1> Vector3;
  typedef Matrix<Scalar,4,1> Vector4;
  typedef Quaternion<Scalar> Quaternionx;
  typedef AngleAxis<Scalar> AngleAxisx;
  typedef Transform<Scalar,2> Transform2;
  typedef Transform<Scalar,3> Transform3;
  typedef Scaling<Scalar,2> Scaling2;
  typedef Scaling<Scalar,3> Scaling3;
  typedef Translation<Scalar,2> Translation2;
  typedef Translation<Scalar,3> Translation3;

  Scalar largeEps = test_precision<Scalar>();
  if (ei_is_same_type<Scalar,float>::ret)
    largeEps = 1e-2f;

  Vector3 v0 = Vector3::Random(),
    v1 = Vector3::Random(),
    v2 = Vector3::Random();
  Vector2 u0 = Vector2::Random();
  Matrix3 matrot1;

  Scalar a = ei_random<Scalar>(-Scalar(M_PI), Scalar(M_PI));

  // cross product
  VERIFY_IS_MUCH_SMALLER_THAN(v1.cross(v2).eigen2_dot(v1), Scalar(1));
  Matrix3 m;
  m << v0.normalized(),
      (v0.cross(v1)).normalized(),
      (v0.cross(v1).cross(v0)).normalized();
  VERIFY(m.isUnitary());

  // Quaternion: Identity(), setIdentity();
  Quaternionx q1, q2;
  q2.setIdentity();
  VERIFY_IS_APPROX(Quaternionx(Quaternionx::Identity()).coeffs(), q2.coeffs());
  q1.coeffs().setRandom();
  VERIFY_IS_APPROX(q1.coeffs(), (q1*q2).coeffs());

  // unitOrthogonal
  VERIFY_IS_MUCH_SMALLER_THAN(u0.unitOrthogonal().eigen2_dot(u0), Scalar(1));
  VERIFY_IS_MUCH_SMALLER_THAN(v0.unitOrthogonal().eigen2_dot(v0), Scalar(1));
  VERIFY_IS_APPROX(u0.unitOrthogonal().norm(), Scalar(1));
  VERIFY_IS_APPROX(v0.unitOrthogonal().norm(), Scalar(1));


  VERIFY_IS_APPROX(v0, AngleAxisx(a, v0.normalized()) * v0);
  VERIFY_IS_APPROX(-v0, AngleAxisx(Scalar(M_PI), v0.unitOrthogonal()) * v0);
  VERIFY_IS_APPROX(ei_cos(a)*v0.squaredNorm(), v0.eigen2_dot(AngleAxisx(a, v0.unitOrthogonal()) * v0));
  m = AngleAxisx(a, v0.normalized()).toRotationMatrix().adjoint();
  VERIFY_IS_APPROX(Matrix3::Identity(), m * AngleAxisx(a, v0.normalized()));
  VERIFY_IS_APPROX(Matrix3::Identity(), AngleAxisx(a, v0.normalized()) * m);

  q1 = AngleAxisx(a, v0.normalized());
  q2 = AngleAxisx(a, v1.normalized());

  // angular distance
  Scalar refangle = ei_abs(AngleAxisx(q1.inverse()*q2).angle());
  if (refangle>Scalar(M_PI))
    refangle = Scalar(2)*Scalar(M_PI) - refangle;
  
  if((q1.coeffs()-q2.coeffs()).norm() > 10*largeEps)
  {
    VERIFY(ei_isApprox(q1.angularDistance(q2), refangle, largeEps));
  }

  // rotation matrix conversion
  VERIFY_IS_APPROX(q1 * v2, q1.toRotationMatrix() * v2);
  VERIFY_IS_APPROX(q1 * q2 * v2,
    q1.toRotationMatrix() * q2.toRotationMatrix() * v2);

  VERIFY( (q2*q1).isApprox(q1*q2, largeEps) || !(q2 * q1 * v2).isApprox(
    q1.toRotationMatrix() * q2.toRotationMatrix() * v2));

  q2 = q1.toRotationMatrix();
  VERIFY_IS_APPROX(q1*v1,q2*v1);

  matrot1 = AngleAxisx(Scalar(0.1), Vector3::UnitX())
          * AngleAxisx(Scalar(0.2), Vector3::UnitY())
          * AngleAxisx(Scalar(0.3), Vector3::UnitZ());
  VERIFY_IS_APPROX(matrot1 * v1,
       AngleAxisx(Scalar(0.1), Vector3(1,0,0)).toRotationMatrix()
    * (AngleAxisx(Scalar(0.2), Vector3(0,1,0)).toRotationMatrix()
    * (AngleAxisx(Scalar(0.3), Vector3(0,0,1)).toRotationMatrix() * v1)));

  // angle-axis conversion
  AngleAxisx aa = q1;
  VERIFY_IS_APPROX(q1 * v1, Quaternionx(aa) * v1);
  VERIFY_IS_NOT_APPROX(q1 * v1, Quaternionx(AngleAxisx(aa.angle()*2,aa.axis())) * v1);

  // from two vector creation
  VERIFY_IS_APPROX(v2.normalized(),(q2.setFromTwoVectors(v1,v2)*v1).normalized());
  VERIFY_IS_APPROX(v2.normalized(),(q2.setFromTwoVectors(v1,v2)*v1).normalized());

  // inverse and conjugate
  VERIFY_IS_APPROX(q1 * (q1.inverse() * v1), v1);
  VERIFY_IS_APPROX(q1 * (q1.conjugate() * v1), v1);

  // AngleAxis
  VERIFY_IS_APPROX(AngleAxisx(a,v1.normalized()).toRotationMatrix(),
    Quaternionx(AngleAxisx(a,v1.normalized())).toRotationMatrix());

  AngleAxisx aa1;
  m = q1.toRotationMatrix();
  aa1 = m;
  VERIFY_IS_APPROX(AngleAxisx(m).toRotationMatrix(),
    Quaternionx(m).toRotationMatrix());

  // Transform
  // TODO complete the tests !
  a = 0;
  while (ei_abs(a)<Scalar(0.1))
    a = ei_random<Scalar>(-Scalar(0.4)*Scalar(M_PI), Scalar(0.4)*Scalar(M_PI));
  q1 = AngleAxisx(a, v0.normalized());
  Transform3 t0, t1, t2;
  // first test setIdentity() and Identity()
  t0.setIdentity();
  VERIFY_IS_APPROX(t0.matrix(), Transform3::MatrixType::Identity());
  t0.matrix().setZero();
  t0 = Transform3::Identity();
  VERIFY_IS_APPROX(t0.matrix(), Transform3::MatrixType::Identity());

  t0.linear() = q1.toRotationMatrix();
  t1.setIdentity();
  t1.linear() = q1.toRotationMatrix();

  v0 << 50, 2, 1;//= ei_random_matrix<Vector3>().cwiseProduct(Vector3(10,2,0.5));
  t0.scale(v0);
  t1.prescale(v0);

  VERIFY_IS_APPROX( (t0 * Vector3(1,0,0)).norm(), v0.x());
  //VERIFY(!ei_isApprox((t1 * Vector3(1,0,0)).norm(), v0.x()));

  t0.setIdentity();
  t1.setIdentity();
  v1 << 1, 2, 3;
  t0.linear() = q1.toRotationMatrix();
  t0.pretranslate(v0);
  t0.scale(v1);
  t1.linear() = q1.conjugate().toRotationMatrix();
  t1.prescale(v1.cwise().inverse());
  t1.translate(-v0);

  VERIFY((t0.matrix() * t1.matrix()).isIdentity(test_precision<Scalar>()));

  t1.fromPositionOrientationScale(v0, q1, v1);
  VERIFY_IS_APPROX(t1.matrix(), t0.matrix());
  VERIFY_IS_APPROX(t1*v1, t0*v1);

  t0.setIdentity(); t0.scale(v0).rotate(q1.toRotationMatrix());
  t1.setIdentity(); t1.scale(v0).rotate(q1);
  VERIFY_IS_APPROX(t0.matrix(), t1.matrix());

  t0.setIdentity(); t0.scale(v0).rotate(AngleAxisx(q1));
  VERIFY_IS_APPROX(t0.matrix(), t1.matrix());

  VERIFY_IS_APPROX(t0.scale(a).matrix(), t1.scale(Vector3::Constant(a)).matrix());
  VERIFY_IS_APPROX(t0.prescale(a).matrix(), t1.prescale(Vector3::Constant(a)).matrix());

  // More transform constructors, operator=, operator*=

  Matrix3 mat3 = Matrix3::Random();
  Matrix4 mat4;
  mat4 << mat3 , Vector3::Zero() , Vector4::Zero().transpose();
  Transform3 tmat3(mat3), tmat4(mat4);
  tmat4.matrix()(3,3) = Scalar(1);
  VERIFY_IS_APPROX(tmat3.matrix(), tmat4.matrix());

  Scalar a3 = ei_random<Scalar>(-Scalar(M_PI), Scalar(M_PI));
  Vector3 v3 = Vector3::Random().normalized();
  AngleAxisx aa3(a3, v3);
  Transform3 t3(aa3);
  Transform3 t4;
  t4 = aa3;
  VERIFY_IS_APPROX(t3.matrix(), t4.matrix());
  t4.rotate(AngleAxisx(-a3,v3));
  VERIFY_IS_APPROX(t4.matrix(), Matrix4::Identity());
  t4 *= aa3;
  VERIFY_IS_APPROX(t3.matrix(), t4.matrix());

  v3 = Vector3::Random();
  Translation3 tv3(v3);
  Transform3 t5(tv3);
  t4 = tv3;
  VERIFY_IS_APPROX(t5.matrix(), t4.matrix());
  t4.translate(-v3);
  VERIFY_IS_APPROX(t4.matrix(), Matrix4::Identity());
  t4 *= tv3;
  VERIFY_IS_APPROX(t5.matrix(), t4.matrix());

  Scaling3 sv3(v3);
  Transform3 t6(sv3);
  t4 = sv3;
  VERIFY_IS_APPROX(t6.matrix(), t4.matrix());
  t4.scale(v3.cwise().inverse());
  VERIFY_IS_APPROX(t4.matrix(), Matrix4::Identity());
  t4 *= sv3;
  VERIFY_IS_APPROX(t6.matrix(), t4.matrix());

  // matrix * transform
  VERIFY_IS_APPROX(Transform3(t3.matrix()*t4).matrix(), Transform3(t3*t4).matrix());

  // chained Transform product
  VERIFY_IS_APPROX(((t3*t4)*t5).matrix(), (t3*(t4*t5)).matrix());

  // check that Transform product doesn't have aliasing problems
  t5 = t4;
  t5 = t5*t5;
  VERIFY_IS_APPROX(t5, t4*t4);

  // 2D transformation
  Transform2 t20, t21;
  Vector2 v20 = Vector2::Random();
  Vector2 v21 = Vector2::Random();
  for (int k=0; k<2; ++k)
    if (ei_abs(v21[k])<Scalar(1e-3)) v21[k] = Scalar(1e-3);
  t21.setIdentity();
  t21.linear() = Rotation2D<Scalar>(a).toRotationMatrix();
  VERIFY_IS_APPROX(t20.fromPositionOrientationScale(v20,a,v21).matrix(),
    t21.pretranslate(v20).scale(v21).matrix());

  t21.setIdentity();
  t21.linear() = Rotation2D<Scalar>(-a).toRotationMatrix();
  VERIFY( (t20.fromPositionOrientationScale(v20,a,v21)
        * (t21.prescale(v21.cwise().inverse()).translate(-v20))).matrix().isIdentity(test_precision<Scalar>()) );

  // Transform - new API
  // 3D
  t0.setIdentity();
  t0.rotate(q1).scale(v0).translate(v0);
  // mat * scaling and mat * translation
  t1 = (Matrix3(q1) * Scaling3(v0)) * Translation3(v0);
  VERIFY_IS_APPROX(t0.matrix(), t1.matrix());
  // mat * transformation and scaling * translation
  t1 = Matrix3(q1) * (Scaling3(v0) * Translation3(v0));
  VERIFY_IS_APPROX(t0.matrix(), t1.matrix());

  t0.setIdentity();
  t0.prerotate(q1).prescale(v0).pretranslate(v0);
  // translation * scaling and transformation * mat
  t1 = (Translation3(v0) * Scaling3(v0)) * Matrix3(q1);
  VERIFY_IS_APPROX(t0.matrix(), t1.matrix());
  // scaling * mat and translation * mat
  t1 = Translation3(v0) * (Scaling3(v0) * Matrix3(q1));
  VERIFY_IS_APPROX(t0.matrix(), t1.matrix());

  t0.setIdentity();
  t0.scale(v0).translate(v0).rotate(q1);
  // translation * mat and scaling * transformation
  t1 = Scaling3(v0) * (Translation3(v0) * Matrix3(q1));
  VERIFY_IS_APPROX(t0.matrix(), t1.matrix());
  // transformation * scaling
  t0.scale(v0);
  t1 = t1 * Scaling3(v0);
  VERIFY_IS_APPROX(t0.matrix(), t1.matrix());
  // transformation * translation
  t0.translate(v0);
  t1 = t1 * Translation3(v0);
  VERIFY_IS_APPROX(t0.matrix(), t1.matrix());
  // translation * transformation
  t0.pretranslate(v0);
  t1 = Translation3(v0) * t1;
  VERIFY_IS_APPROX(t0.matrix(), t1.matrix());

  // transform * quaternion
  t0.rotate(q1);
  t1 = t1 * q1;
  VERIFY_IS_APPROX(t0.matrix(), t1.matrix());

  // translation * quaternion
  t0.translate(v1).rotate(q1);
  t1 = t1 * (Translation3(v1) * q1);
  VERIFY_IS_APPROX(t0.matrix(), t1.matrix());

  // scaling * quaternion
  t0.scale(v1).rotate(q1);
  t1 = t1 * (Scaling3(v1) * q1);
  VERIFY_IS_APPROX(t0.matrix(), t1.matrix());

  // quaternion * transform
  t0.prerotate(q1);
  t1 = q1 * t1;
  VERIFY_IS_APPROX(t0.matrix(), t1.matrix());

  // quaternion * translation
  t0.rotate(q1).translate(v1);
  t1 = t1 * (q1 * Translation3(v1));
  VERIFY_IS_APPROX(t0.matrix(), t1.matrix());

  // quaternion * scaling
  t0.rotate(q1).scale(v1);
  t1 = t1 * (q1 * Scaling3(v1));
  VERIFY_IS_APPROX(t0.matrix(), t1.matrix());

  // translation * vector
  t0.setIdentity();
  t0.translate(v0);
  VERIFY_IS_APPROX(t0 * v1, Translation3(v0) * v1);

  // scaling * vector
  t0.setIdentity();
  t0.scale(v0);
  VERIFY_IS_APPROX(t0 * v1, Scaling3(v0) * v1);

  // test transform inversion
  t0.setIdentity();
  t0.translate(v0);
  t0.linear().setRandom();
  VERIFY_IS_APPROX(t0.inverse(Affine), t0.matrix().inverse());
  t0.setIdentity();
  t0.translate(v0).rotate(q1);
  VERIFY_IS_APPROX(t0.inverse(Isometry), t0.matrix().inverse());

  // test extract rotation and scaling
  t0.setIdentity();
  t0.translate(v0).rotate(q1).scale(v1);
  VERIFY_IS_APPROX(t0.rotation() * v1, Matrix3(q1) * v1);

  Matrix3 mat_rotation, mat_scaling;
  t0.setIdentity();
  t0.translate(v0).rotate(q1).scale(v1);
  t0.computeRotationScaling(&mat_rotation, &mat_scaling);
  VERIFY_IS_APPROX(t0.linear(), mat_rotation * mat_scaling);
  VERIFY_IS_APPROX(mat_rotation*mat_rotation.adjoint(), Matrix3::Identity());
  VERIFY_IS_APPROX(mat_rotation.determinant(), Scalar(1));
  t0.computeScalingRotation(&mat_scaling, &mat_rotation);
  VERIFY_IS_APPROX(t0.linear(), mat_scaling * mat_rotation);
  VERIFY_IS_APPROX(mat_rotation*mat_rotation.adjoint(), Matrix3::Identity());
  VERIFY_IS_APPROX(mat_rotation.determinant(), Scalar(1));

  // test casting
  Transform<float,3> t1f = t1.template cast<float>();
  VERIFY_IS_APPROX(t1f.template cast<Scalar>(),t1);
  Transform<double,3> t1d = t1.template cast<double>();
  VERIFY_IS_APPROX(t1d.template cast<Scalar>(),t1);

  Translation3 tr1(v0);
  Translation<float,3> tr1f = tr1.template cast<float>();
  VERIFY_IS_APPROX(tr1f.template cast<Scalar>(),tr1);
  Translation<double,3> tr1d = tr1.template cast<double>();
  VERIFY_IS_APPROX(tr1d.template cast<Scalar>(),tr1);

  Scaling3 sc1(v0);
  Scaling<float,3> sc1f = sc1.template cast<float>();
  VERIFY_IS_APPROX(sc1f.template cast<Scalar>(),sc1);
  Scaling<double,3> sc1d = sc1.template cast<double>();
  VERIFY_IS_APPROX(sc1d.template cast<Scalar>(),sc1);

  Quaternion<float> q1f = q1.template cast<float>();
  VERIFY_IS_APPROX(q1f.template cast<Scalar>(),q1);
  Quaternion<double> q1d = q1.template cast<double>();
  VERIFY_IS_APPROX(q1d.template cast<Scalar>(),q1);

  AngleAxis<float> aa1f = aa1.template cast<float>();
  VERIFY_IS_APPROX(aa1f.template cast<Scalar>(),aa1);
  AngleAxis<double> aa1d = aa1.template cast<double>();
  VERIFY_IS_APPROX(aa1d.template cast<Scalar>(),aa1);

  Rotation2D<Scalar> r2d1(ei_random<Scalar>());
  Rotation2D<float> r2d1f = r2d1.template cast<float>();
  VERIFY_IS_APPROX(r2d1f.template cast<Scalar>(),r2d1);
  Rotation2D<double> r2d1d = r2d1.template cast<double>();
  VERIFY_IS_APPROX(r2d1d.template cast<Scalar>(),r2d1);

  m = q1;
//   m.col(1) = Vector3(0,ei_random<Scalar>(),ei_random<Scalar>()).normalized();
//   m.col(0) = Vector3(-1,0,0).normalized();
//   m.col(2) = m.col(0).cross(m.col(1));
  #define VERIFY_EULER(I,J,K, X,Y,Z) { \
    Vector3 ea = m.eulerAngles(I,J,K); \
    Matrix3 m1 = Matrix3(AngleAxisx(ea[0], Vector3::Unit##X()) * AngleAxisx(ea[1], Vector3::Unit##Y()) * AngleAxisx(ea[2], Vector3::Unit##Z())); \
    VERIFY_IS_APPROX(m,  Matrix3(AngleAxisx(ea[0], Vector3::Unit##X()) * AngleAxisx(ea[1], Vector3::Unit##Y()) * AngleAxisx(ea[2], Vector3::Unit##Z()))); \
  }
  VERIFY_EULER(0,1,2, X,Y,Z);
  VERIFY_EULER(0,1,0, X,Y,X);
  VERIFY_EULER(0,2,1, X,Z,Y);
  VERIFY_EULER(0,2,0, X,Z,X);

  VERIFY_EULER(1,2,0, Y,Z,X);
  VERIFY_EULER(1,2,1, Y,Z,Y);
  VERIFY_EULER(1,0,2, Y,X,Z);
  VERIFY_EULER(1,0,1, Y,X,Y);

  VERIFY_EULER(2,0,1, Z,X,Y);
  VERIFY_EULER(2,0,2, Z,X,Z);
  VERIFY_EULER(2,1,0, Z,Y,X);
  VERIFY_EULER(2,1,2, Z,Y,Z);

  // colwise/rowwise cross product
  mat3.setRandom();
  Vector3 vec3 = Vector3::Random();
  Matrix3 mcross;
  int i = ei_random<int>(0,2);
  mcross = mat3.colwise().cross(vec3);
  VERIFY_IS_APPROX(mcross.col(i), mat3.col(i).cross(vec3));
  mcross = mat3.rowwise().cross(vec3);
  VERIFY_IS_APPROX(mcross.row(i), mat3.row(i).cross(vec3));


}

void test_eigen2_geometry()
{
  for(int i = 0; i < g_repeat; i++) {
    CALL_SUBTEST_1( geometry<float>() );
    CALL_SUBTEST_2( geometry<double>() );
  }
}
