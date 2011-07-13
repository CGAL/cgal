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
#include <Eigen/QR>

#include<iostream>
using namespace std;

template<typename BoxType> void alignedbox(const BoxType& _box)
{
  /* this test covers the following files:
     AlignedBox.h
  */
  typedef typename BoxType::Index Index;  
  typedef typename BoxType::Scalar Scalar;
  typedef typename NumTraits<Scalar>::Real RealScalar;
  typedef Matrix<Scalar, BoxType::AmbientDimAtCompileTime, 1> VectorType;

  const Index dim = _box.dim();

  VectorType p0 = VectorType::Random(dim);
  VectorType p1 = VectorType::Random(dim);
  while( p1 == p0 ){
      p1 =  VectorType::Random(dim); }
  RealScalar s1 = internal::random<RealScalar>(0,1);

  BoxType b0(dim);
  BoxType b1(VectorType::Random(dim),VectorType::Random(dim));
  BoxType b2;

  b0.extend(p0);
  b0.extend(p1);
  VERIFY(b0.contains(p0*s1+(Scalar(1)-s1)*p1));

  (b2 = b0).extend(b1);
  VERIFY(b2.contains(b0));
  VERIFY(b2.contains(b1));
  VERIFY_IS_APPROX(b2.clamp(b0), b0);


  // alignment -- make sure there is no memory alignment assertion
  BoxType *bp0 = new BoxType(dim);
  BoxType *bp1 = new BoxType(dim);
  bp0->extend(*bp1);
  delete bp0;
  delete bp1;

  // sampling
  for( int i=0; i<10; ++i )
  {
      VectorType r = b0.sample();
      VERIFY(b0.contains(r));
  }

}



template<typename BoxType>
void alignedboxCastTests(const BoxType& _box)
{
  // casting  
  typedef typename BoxType::Index Index;
  typedef typename BoxType::Scalar Scalar;
  typedef typename NumTraits<Scalar>::Real RealScalar;
  typedef Matrix<Scalar, BoxType::AmbientDimAtCompileTime, 1> VectorType;

  const Index dim = _box.dim();

  VectorType p0 = VectorType::Random(dim);
  VectorType p1 = VectorType::Random(dim);

  BoxType b0(dim);

  b0.extend(p0);
  b0.extend(p1);

  const int Dim = BoxType::AmbientDimAtCompileTime;
  typedef typename GetDifferentType<Scalar>::type OtherScalar;
  AlignedBox<OtherScalar,Dim> hp1f = b0.template cast<OtherScalar>();
  VERIFY_IS_APPROX(hp1f.template cast<Scalar>(),b0);
  AlignedBox<Scalar,Dim> hp1d = b0.template cast<Scalar>();
  VERIFY_IS_APPROX(hp1d.template cast<Scalar>(),b0);
}


void specificTest1()
{
    Vector2f m; m << -1.0f, -2.0f;
    Vector2f M; M <<  1.0f,  5.0f;

    typedef AlignedBox<float,2>  BoxType;
    BoxType box( m, M );

    Vector2f sides = M-m;
    VERIFY_IS_APPROX(sides, box.sizes() );
    VERIFY_IS_APPROX(sides[1], box.sizes()[1] );
    VERIFY_IS_APPROX(sides[1], box.sizes().maxCoeff() );
    VERIFY_IS_APPROX(sides[0], box.sizes().minCoeff() );

    VERIFY_IS_APPROX( 14.0f, box.volume() );
    VERIFY_IS_APPROX( 53.0f, box.diagonal().squaredNorm() );
    VERIFY_IS_APPROX( internal::sqrt( 53.0f ), box.diagonal().norm() );

    VERIFY_IS_APPROX( m, box.corner( BoxType::BottomLeft ) );
    VERIFY_IS_APPROX( M, box.corner( BoxType::TopRight ) );
    Vector2f bottomRight; bottomRight << M[0], m[1];
    Vector2f topLeft; topLeft << m[0], M[1];
    VERIFY_IS_APPROX( bottomRight, box.corner( BoxType::BottomRight ) );
    VERIFY_IS_APPROX( topLeft, box.corner( BoxType::TopLeft ) );
}


void specificTest2()
{
    Vector3i m; m << -1, -2, 0;
    Vector3i M; M <<  1,  5, 3;

    typedef AlignedBox<int,3>  BoxType;
    BoxType box( m, M );

    Vector3i sides = M-m;
    VERIFY_IS_APPROX(sides, box.sizes() );
    VERIFY_IS_APPROX(sides[1], box.sizes()[1] );
    VERIFY_IS_APPROX(sides[1], box.sizes().maxCoeff() );
    VERIFY_IS_APPROX(sides[0], box.sizes().minCoeff() );

    VERIFY_IS_APPROX( 42, box.volume() );
    VERIFY_IS_APPROX( 62, box.diagonal().squaredNorm() );

    VERIFY_IS_APPROX( m, box.corner( BoxType::BottomLeftFloor ) );
    VERIFY_IS_APPROX( M, box.corner( BoxType::TopRightCeil ) );
    Vector3i bottomRightFloor; bottomRightFloor << M[0], m[1], m[2];
    Vector3i topLeftFloor; topLeftFloor << m[0], M[1], m[2];
    VERIFY_IS_APPROX( bottomRightFloor, box.corner( BoxType::BottomRightFloor ) );
    VERIFY_IS_APPROX( topLeftFloor, box.corner( BoxType::TopLeftFloor ) );
}


void test_geo_alignedbox()
{
  for(int i = 0; i < g_repeat; i++)
  {
    CALL_SUBTEST_1( alignedbox(AlignedBox<float,2>()) );
    CALL_SUBTEST_2( alignedboxCastTests(AlignedBox<float,2>()) );

    CALL_SUBTEST_3( alignedbox(AlignedBox<float,3>()) );
    CALL_SUBTEST_4( alignedboxCastTests(AlignedBox<float,3>()) );

    CALL_SUBTEST_5( alignedbox(AlignedBox<double,4>()) );
    CALL_SUBTEST_6( alignedboxCastTests(AlignedBox<double,4>()) );

    CALL_SUBTEST_7( alignedbox(AlignedBox<double,1>()) );
    CALL_SUBTEST_8( alignedboxCastTests(AlignedBox<double,1>()) );

    CALL_SUBTEST_9( alignedbox(AlignedBox<int,1>()) );
    CALL_SUBTEST_10( alignedbox(AlignedBox<int,2>()) );
    CALL_SUBTEST_11( alignedbox(AlignedBox<int,3>()) );
  }
  CALL_SUBTEST_12( specificTest1() );
  CALL_SUBTEST_13( specificTest2() );
}
