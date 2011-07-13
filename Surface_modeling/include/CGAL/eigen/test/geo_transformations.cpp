// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2008-2009 Gael Guennebaud <gael.guennebaud@inria.fr>
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

template<typename Scalar, int Mode, int Options> void non_projective_only()
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
  typedef Transform<Scalar,2,Mode,Options> Transform2;
  typedef Transform<Scalar,3,Mode,Options> Transform3;
  typedef Transform<Scalar,2,Isometry,Options> Isometry2;
  typedef Transform<Scalar,3,Isometry,Options> Isometry3;
  typedef typename Transform3::MatrixType MatrixType;
  typedef DiagonalMatrix<Scalar,2> AlignedScaling2;
  typedef DiagonalMatrix<Scalar,3> AlignedScaling3;
  typedef Translation<Scalar,2> Translation2;
  typedef Translation<Scalar,3> Translation3;

  Vector3 v0 = Vector3::Random(),
          v1 = Vector3::Random();

  Transform3 t0, t1, t2;

  Scalar a = internal::random<Scalar>(-Scalar(M_PI), Scalar(M_PI));

  Quaternionx q1, q2;

  q1 = AngleAxisx(a, v0.normalized());

  t0 = Transform3::Identity();
  VERIFY_IS_APPROX(t0.matrix(), Transform3::MatrixType::Identity());

  t0.linear() = q1.toRotationMatrix();

  v0 << 50, 2, 1;
  t0.scale(v0);

  VERIFY_IS_APPROX( (t0 * Vector3(1,0,0)).template head<3>().norm(), v0.x());

  t0.setIdentity();
  t1.setIdentity();
  v1 << 1, 2, 3;
  t0.linear() = q1.toRotationMatrix();
  t0.pretranslate(v0);
  t0.scale(v1);
  t1.linear() = q1.conjugate().toRotationMatrix();
  t1.prescale(v1.cwiseInverse());
  t1.translate(-v0);

  VERIFY((t0 * t1).matrix().isIdentity(test_precision<Scalar>()));

  t1.fromPositionOrientationScale(v0, q1, v1);
  VERIFY_IS_APPROX(t1.matrix(), t0.matrix());
  VERIFY_IS_APPROX(t1*v1, t0*v1);

  // translation * vector
  t0.setIdentity();
  t0.translate(v0);
  VERIFY_IS_APPROX((t0 * v1).template head<3>(), Translation3(v0) * v1);

  // AlignedScaling * vector
  t0.setIdentity();
  t0.scale(v0);
  VERIFY_IS_APPROX((t0 * v1).template head<3>(), AlignedScaling3(v0) * v1);
}

template<typename Scalar, int Mode, int Options> void transformations()
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
  typedef Transform<Scalar,2,Mode,Options> Transform2;
  typedef Transform<Scalar,3,Mode,Options> Transform3;
  typedef Transform<Scalar,2,Isometry,Options> Isometry2;
  typedef Transform<Scalar,3,Isometry,Options> Isometry3;
  typedef typename Transform3::MatrixType MatrixType;
  typedef DiagonalMatrix<Scalar,2> AlignedScaling2;
  typedef DiagonalMatrix<Scalar,3> AlignedScaling3;
  typedef Translation<Scalar,2> Translation2;
  typedef Translation<Scalar,3> Translation3;

  Vector3 v0 = Vector3::Random(),
    v1 = Vector3::Random(),
    v2 = Vector3::Random();
  Vector2 u0 = Vector2::Random();
  Matrix3 matrot1, m;

  Scalar a = internal::random<Scalar>(-Scalar(M_PI), Scalar(M_PI));
  Scalar s0 = internal::random<Scalar>();

  VERIFY_IS_APPROX(v0, AngleAxisx(a, v0.normalized()) * v0);
  VERIFY_IS_APPROX(-v0, AngleAxisx(Scalar(M_PI), v0.unitOrthogonal()) * v0);
  VERIFY_IS_APPROX(internal::cos(a)*v0.squaredNorm(), v0.dot(AngleAxisx(a, v0.unitOrthogonal()) * v0));
  m = AngleAxisx(a, v0.normalized()).toRotationMatrix().adjoint();
  VERIFY_IS_APPROX(Matrix3::Identity(), m * AngleAxisx(a, v0.normalized()));
  VERIFY_IS_APPROX(Matrix3::Identity(), AngleAxisx(a, v0.normalized()) * m);

  Quaternionx q1, q2;
  q1 = AngleAxisx(a, v0.normalized());
  q2 = AngleAxisx(a, v1.normalized());

  // rotation matrix conversion
  matrot1 = AngleAxisx(Scalar(0.1), Vector3::UnitX())
          * AngleAxisx(Scalar(0.2), Vector3::UnitY())
          * AngleAxisx(Scalar(0.3), Vector3::UnitZ());
  VERIFY_IS_APPROX(matrot1 * v1,
       AngleAxisx(Scalar(0.1), Vector3(1,0,0)).toRotationMatrix()
    * (AngleAxisx(Scalar(0.2), Vector3(0,1,0)).toRotationMatrix()
    * (AngleAxisx(Scalar(0.3), Vector3(0,0,1)).toRotationMatrix() * v1)));

  // angle-axis conversion
  AngleAxisx aa = AngleAxisx(q1);
  VERIFY_IS_APPROX(q1 * v1, Quaternionx(aa) * v1);
  VERIFY_IS_NOT_APPROX(q1 * v1, Quaternionx(AngleAxisx(aa.angle()*2,aa.axis())) * v1);

  aa.fromRotationMatrix(aa.toRotationMatrix());
  VERIFY_IS_APPROX(q1 * v1, Quaternionx(aa) * v1);
  VERIFY_IS_NOT_APPROX(q1 * v1, Quaternionx(AngleAxisx(aa.angle()*2,aa.axis())) * v1);

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
  while (internal::abs(a)<Scalar(0.1))
    a = internal::random<Scalar>(-Scalar(0.4)*Scalar(M_PI), Scalar(0.4)*Scalar(M_PI));
  q1 = AngleAxisx(a, v0.normalized());
  Transform3 t0, t1, t2;

  // first test setIdentity() and Identity()
  t0.setIdentity();
  VERIFY_IS_APPROX(t0.matrix(), Transform3::MatrixType::Identity());
  t0.matrix().setZero();
  t0 = Transform3::Identity();
  VERIFY_IS_APPROX(t0.matrix(), Transform3::MatrixType::Identity());

  t0.setIdentity();
  t1.setIdentity();
  v1 << 1, 2, 3;
  t0.linear() = q1.toRotationMatrix();
  t0.pretranslate(v0);
  t0.scale(v1);
  t1.linear() = q1.conjugate().toRotationMatrix();
  t1.prescale(v1.cwiseInverse());
  t1.translate(-v0);

  VERIFY((t0 * t1).matrix().isIdentity(test_precision<Scalar>()));

  t1.fromPositionOrientationScale(v0, q1, v1);
  VERIFY_IS_APPROX(t1.matrix(), t0.matrix());

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
  if(Mode!=int(AffineCompact))
    tmat4.matrix()(3,3) = Scalar(1);
  VERIFY_IS_APPROX(tmat3.matrix(), tmat4.matrix());

  Scalar a3 = internal::random<Scalar>(-Scalar(M_PI), Scalar(M_PI));
  Vector3 v3 = Vector3::Random().normalized();
  AngleAxisx aa3(a3, v3);
  Transform3 t3(aa3);
  Transform3 t4;
  t4 = aa3;
  VERIFY_IS_APPROX(t3.matrix(), t4.matrix());
  t4.rotate(AngleAxisx(-a3,v3));
  VERIFY_IS_APPROX(t4.matrix(), MatrixType::Identity());
  t4 *= aa3;
  VERIFY_IS_APPROX(t3.matrix(), t4.matrix());

  v3 = Vector3::Random();
  Translation3 tv3(v3);
  Transform3 t5(tv3);
  t4 = tv3;
  VERIFY_IS_APPROX(t5.matrix(), t4.matrix());
  t4.translate(-v3);
  VERIFY_IS_APPROX(t4.matrix(), MatrixType::Identity());
  t4 *= tv3;
  VERIFY_IS_APPROX(t5.matrix(), t4.matrix());

  AlignedScaling3 sv3(v3);
  Transform3 t6(sv3);
  t4 = sv3;
  VERIFY_IS_APPROX(t6.matrix(), t4.matrix());
  t4.scale(v3.cwiseInverse());
  VERIFY_IS_APPROX(t4.matrix(), MatrixType::Identity());
  t4 *= sv3;
  VERIFY_IS_APPROX(t6.matrix(), t4.matrix());

  // matrix * transform
  VERIFY_IS_APPROX((t3.matrix()*t4).matrix(), (t3*t4).matrix());

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
    if (internal::abs(v21[k])<Scalar(1e-3)) v21[k] = Scalar(1e-3);
  t21.setIdentity();
  t21.linear() = Rotation2D<Scalar>(a).toRotationMatrix();
  VERIFY_IS_APPROX(t20.fromPositionOrientationScale(v20,a,v21).matrix(),
    t21.pretranslate(v20).scale(v21).matrix());

  t21.setIdentity();
  t21.linear() = Rotation2D<Scalar>(-a).toRotationMatrix();
  VERIFY( (t20.fromPositionOrientationScale(v20,a,v21)
        * (t21.prescale(v21.cwiseInverse()).translate(-v20))).matrix().isIdentity(test_precision<Scalar>()) );

  // Transform - new API
  // 3D
  t0.setIdentity();
  t0.rotate(q1).scale(v0).translate(v0);
  // mat * aligned scaling and mat * translation
  t1 = (Matrix3(q1) * AlignedScaling3(v0)) * Translation3(v0);
  VERIFY_IS_APPROX(t0.matrix(), t1.matrix());
  t1 = (Matrix3(q1) * Scaling(v0)) * Translation3(v0);
  VERIFY_IS_APPROX(t0.matrix(), t1.matrix());
  t1 = (q1 * Scaling(v0)) * Translation3(v0);
  VERIFY_IS_APPROX(t0.matrix(), t1.matrix());
  // mat * transformation and aligned scaling * translation
  t1 = Matrix3(q1) * (AlignedScaling3(v0) * Translation3(v0));
  VERIFY_IS_APPROX(t0.matrix(), t1.matrix());


  t0.setIdentity();
  t0.scale(s0).translate(v0);
  t1 = Scaling(s0) * Translation3(v0);
  VERIFY_IS_APPROX(t0.matrix(), t1.matrix());
  t0.prescale(s0);
  t1 = Scaling(s0) * t1;
  VERIFY_IS_APPROX(t0.matrix(), t1.matrix());
  
  t0 = t3;
  t0.scale(s0);
  t1 = t3 * Scaling(s0,s0,s0);
  VERIFY_IS_APPROX(t0.matrix(), t1.matrix());
  t0.prescale(s0);
  t1 = Scaling(s0,s0,s0) * t1;
  VERIFY_IS_APPROX(t0.matrix(), t1.matrix());


  t0.setIdentity();
  t0.prerotate(q1).prescale(v0).pretranslate(v0);
  // translation * aligned scaling and transformation * mat
  t1 = (Translation3(v0) * AlignedScaling3(v0)) * Transform3(q1);
  VERIFY_IS_APPROX(t0.matrix(), t1.matrix());
  // scaling * mat and translation * mat
  t1 = Translation3(v0) * (AlignedScaling3(v0) * Transform3(q1));
  VERIFY_IS_APPROX(t0.matrix(), t1.matrix());

  t0.setIdentity();
  t0.scale(v0).translate(v0).rotate(q1);
  // translation * mat and aligned scaling * transformation
  t1 = AlignedScaling3(v0) * (Translation3(v0) * Transform3(q1));
  VERIFY_IS_APPROX(t0.matrix(), t1.matrix());
  // transformation * aligned scaling
  t0.scale(v0);
  t1 *= AlignedScaling3(v0);
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

  // aligned scaling * quaternion
  t0.scale(v1).rotate(q1);
  t1 = t1 * (AlignedScaling3(v1) * q1);
  VERIFY_IS_APPROX(t0.matrix(), t1.matrix());

  // quaternion * transform
  t0.prerotate(q1);
  t1 = q1 * t1;
  VERIFY_IS_APPROX(t0.matrix(), t1.matrix());

  // quaternion * translation
  t0.rotate(q1).translate(v1);
  t1 = t1 * (q1 * Translation3(v1));
  VERIFY_IS_APPROX(t0.matrix(), t1.matrix());

  // quaternion * aligned scaling
  t0.rotate(q1).scale(v1);
  t1 = t1 * (q1 * AlignedScaling3(v1));
  VERIFY_IS_APPROX(t0.matrix(), t1.matrix());

  // test transform inversion
  t0.setIdentity();
  t0.translate(v0);
  t0.linear().setRandom();
  Matrix4 t044 = Matrix4::Zero();
  t044(3,3) = 1;
  t044.block(0,0,t0.matrix().rows(),4) = t0.matrix();
  VERIFY_IS_APPROX(t0.inverse(Affine).matrix(), t044.inverse().block(0,0,t0.matrix().rows(),4));
  t0.setIdentity();
  t0.translate(v0).rotate(q1);
  t044 = Matrix4::Zero();
  t044(3,3) = 1;
  t044.block(0,0,t0.matrix().rows(),4) = t0.matrix();
  VERIFY_IS_APPROX(t0.inverse(Isometry).matrix(), t044.inverse().block(0,0,t0.matrix().rows(),4));

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
  Transform<float,3,Mode> t1f = t1.template cast<float>();
  VERIFY_IS_APPROX(t1f.template cast<Scalar>(),t1);
  Transform<double,3,Mode> t1d = t1.template cast<double>();
  VERIFY_IS_APPROX(t1d.template cast<Scalar>(),t1);

  Translation3 tr1(v0);
  Translation<float,3> tr1f = tr1.template cast<float>();
  VERIFY_IS_APPROX(tr1f.template cast<Scalar>(),tr1);
  Translation<double,3> tr1d = tr1.template cast<double>();
  VERIFY_IS_APPROX(tr1d.template cast<Scalar>(),tr1);

  AngleAxis<float> aa1f = aa1.template cast<float>();
  VERIFY_IS_APPROX(aa1f.template cast<Scalar>(),aa1);
  AngleAxis<double> aa1d = aa1.template cast<double>();
  VERIFY_IS_APPROX(aa1d.template cast<Scalar>(),aa1);

  Rotation2D<Scalar> r2d1(internal::random<Scalar>());
  Rotation2D<float> r2d1f = r2d1.template cast<float>();
  VERIFY_IS_APPROX(r2d1f.template cast<Scalar>(),r2d1);
  Rotation2D<double> r2d1d = r2d1.template cast<double>();
  VERIFY_IS_APPROX(r2d1d.template cast<Scalar>(),r2d1);

}

template<typename Scalar> void transform_alignment()
{
  typedef Transform<Scalar,3,Projective,AutoAlign> Projective3a;
  typedef Transform<Scalar,3,Projective,DontAlign> Projective3u;

  EIGEN_ALIGN16 Scalar array1[16];
  EIGEN_ALIGN16 Scalar array2[16];
  EIGEN_ALIGN16 Scalar array3[16+1];
  Scalar* array3u = array3+1;

  Projective3a *p1 = ::new(reinterpret_cast<void*>(array1)) Projective3a;
  Projective3u *p2 = ::new(reinterpret_cast<void*>(array2)) Projective3u;
  Projective3u *p3 = ::new(reinterpret_cast<void*>(array3u)) Projective3u;
  
  p1->matrix().setRandom();
  *p2 = *p1;
  *p3 = *p1;

  VERIFY_IS_APPROX(p1->matrix(), p2->matrix());
  VERIFY_IS_APPROX(p1->matrix(), p3->matrix());
  
  VERIFY_IS_APPROX( (*p1) * (*p1), (*p2)*(*p3));
  
  #ifdef EIGEN_VECTORIZE
  VERIFY_RAISES_ASSERT((::new(reinterpret_cast<void*>(array3u)) Projective3a));
  #endif
}

void test_geo_transformations()
{
  for(int i = 0; i < g_repeat; i++) {
    CALL_SUBTEST_1(( transformations<double,Affine,AutoAlign>() ));
    CALL_SUBTEST_1(( non_projective_only<double,Affine,AutoAlign>() ));
    
    CALL_SUBTEST_2(( transformations<float,AffineCompact,AutoAlign>() ));
    CALL_SUBTEST_2(( non_projective_only<float,AffineCompact,AutoAlign>() ));

    CALL_SUBTEST_3(( transformations<double,Projective,AutoAlign>() ));
    CALL_SUBTEST_3(( transformations<double,Projective,DontAlign>() ));
    CALL_SUBTEST_3(( transform_alignment<double>() ));
    
    CALL_SUBTEST_4(( transformations<float,Affine,RowMajor|AutoAlign>() ));
    CALL_SUBTEST_4(( non_projective_only<float,Affine,RowMajor>() ));
    
    CALL_SUBTEST_5(( transformations<double,AffineCompact,RowMajor|AutoAlign>() ));
    CALL_SUBTEST_5(( non_projective_only<double,AffineCompact,RowMajor>() ));

    CALL_SUBTEST_6(( transformations<double,Projective,RowMajor|AutoAlign>() ));
    CALL_SUBTEST_6(( transformations<double,Projective,RowMajor|DontAlign>() ));
  }
}
