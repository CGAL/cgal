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
#include <Eigen/Geometry>

template<typename Scalar,int Size> void homogeneous(void)
{
  /* this test covers the following files:
     Homogeneous.h
  */

  typedef Matrix<Scalar,Size,Size> MatrixType;
  typedef Matrix<Scalar,Size,1, ColMajor> VectorType;

  typedef Matrix<Scalar,Size+1,Size> HMatrixType;
  typedef Matrix<Scalar,Size+1,1> HVectorType;

  typedef Matrix<Scalar,Size,Size+1>   T1MatrixType;
  typedef Matrix<Scalar,Size+1,Size+1> T2MatrixType;
  typedef Matrix<Scalar,Size+1,Size> T3MatrixType;

  VectorType v0 = VectorType::Random(),
             v1 = VectorType::Random(),
             ones = VectorType::Ones();

  HVectorType hv0 = HVectorType::Random(),
              hv1 = HVectorType::Random();

  MatrixType m0 = MatrixType::Random(),
             m1 = MatrixType::Random();

  HMatrixType hm0 = HMatrixType::Random(),
              hm1 = HMatrixType::Random();

  hv0 << v0, 1;
  VERIFY_IS_APPROX(v0.homogeneous(), hv0);
  VERIFY_IS_APPROX(v0, hv0.hnormalized());

  hm0 << m0, ones.transpose();
  VERIFY_IS_APPROX(m0.colwise().homogeneous(), hm0);
  VERIFY_IS_APPROX(m0, hm0.colwise().hnormalized());
  hm0.row(Size-1).setRandom();
  for(int j=0; j<Size; ++j)
    m0.col(j) = hm0.col(j).head(Size) / hm0(Size,j);
  VERIFY_IS_APPROX(m0, hm0.colwise().hnormalized());

  T1MatrixType t1 = T1MatrixType::Random();
  VERIFY_IS_APPROX(t1 * (v0.homogeneous().eval()), t1 * v0.homogeneous());
  VERIFY_IS_APPROX(t1 * (m0.colwise().homogeneous().eval()), t1 * m0.colwise().homogeneous());

  T2MatrixType t2 = T2MatrixType::Random();
  VERIFY_IS_APPROX(t2 * (v0.homogeneous().eval()), t2 * v0.homogeneous());
  VERIFY_IS_APPROX(t2 * (m0.colwise().homogeneous().eval()), t2 * m0.colwise().homogeneous());

  VERIFY_IS_APPROX((v0.transpose().rowwise().homogeneous().eval()) * t2,
                    v0.transpose().rowwise().homogeneous() * t2);
                    m0.transpose().rowwise().homogeneous().eval();
  VERIFY_IS_APPROX((m0.transpose().rowwise().homogeneous().eval()) * t2,
                    m0.transpose().rowwise().homogeneous() * t2);

  T3MatrixType t3 = T3MatrixType::Random();
  VERIFY_IS_APPROX((v0.transpose().rowwise().homogeneous().eval()) * t3,
                    v0.transpose().rowwise().homogeneous() * t3);
  VERIFY_IS_APPROX((m0.transpose().rowwise().homogeneous().eval()) * t3,
                    m0.transpose().rowwise().homogeneous() * t3);

  // test product with a Transform object
  Transform<Scalar, Size, Affine> aff;
  Transform<Scalar, Size, AffineCompact> caff;
  Transform<Scalar, Size, Projective> proj;
  Matrix<Scalar, Size, Dynamic>   pts;
  Matrix<Scalar, Size+1, Dynamic> pts1, pts2;

  aff.affine().setRandom();
  proj = caff = aff;
  pts.setRandom(Size,internal::random<int>(1,20));
  
  pts1 = pts.colwise().homogeneous();
  VERIFY_IS_APPROX(aff  * pts.colwise().homogeneous(), (aff  * pts1).colwise().hnormalized());
  VERIFY_IS_APPROX(caff * pts.colwise().homogeneous(), (caff * pts1).colwise().hnormalized());
  VERIFY_IS_APPROX(proj * pts.colwise().homogeneous(), (proj * pts1));
  
  VERIFY_IS_APPROX((aff  * pts1).colwise().hnormalized(),  aff  * pts);
  VERIFY_IS_APPROX((caff * pts1).colwise().hnormalized(), caff * pts);
  
  pts2 = pts1;
  pts2.row(Size).setRandom();
  VERIFY_IS_APPROX((aff  * pts2).colwise().hnormalized(), aff  * pts2.colwise().hnormalized());
  VERIFY_IS_APPROX((caff * pts2).colwise().hnormalized(), caff * pts2.colwise().hnormalized());
  VERIFY_IS_APPROX((proj * pts2).colwise().hnormalized(), (proj * pts2.colwise().hnormalized().colwise().homogeneous()).colwise().hnormalized());
}

void test_geo_homogeneous()
{
  for(int i = 0; i < g_repeat; i++) {
    CALL_SUBTEST_1(( homogeneous<float,1>() ));
    CALL_SUBTEST_2(( homogeneous<double,3>() ));
    CALL_SUBTEST_3(( homogeneous<double,8>() ));
  }
}
