// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2008 Gael Guennebaud <gael.guennebaud@inria.fr>
// Copyright (C) 2008 Benoit Jacob <jacob.benoit.1@gmail.com>
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

template<typename LineType> void parametrizedline(const LineType& _line)
{
  /* this test covers the following files:
     ParametrizedLine.h
  */
  typedef typename LineType::Index Index;
  const Index dim = _line.dim();
  typedef typename LineType::Scalar Scalar;
  typedef typename NumTraits<Scalar>::Real RealScalar;
  typedef Matrix<Scalar, LineType::AmbientDimAtCompileTime, 1> VectorType;
  typedef Matrix<Scalar, LineType::AmbientDimAtCompileTime,
                         LineType::AmbientDimAtCompileTime> MatrixType;

  VectorType p0 = VectorType::Random(dim);
  VectorType p1 = VectorType::Random(dim);

  VectorType d0 = VectorType::Random(dim).normalized();

  LineType l0(p0, d0);

  Scalar s0 = internal::random<Scalar>();
  Scalar s1 = internal::abs(internal::random<Scalar>());

  VERIFY_IS_MUCH_SMALLER_THAN( l0.distance(p0), RealScalar(1) );
  VERIFY_IS_MUCH_SMALLER_THAN( l0.distance(p0+s0*d0), RealScalar(1) );
  VERIFY_IS_APPROX( (l0.projection(p1)-p1).norm(), l0.distance(p1) );
  VERIFY_IS_MUCH_SMALLER_THAN( l0.distance(l0.projection(p1)), RealScalar(1) );
  VERIFY_IS_APPROX( Scalar(l0.distance((p0+s0*d0) + d0.unitOrthogonal() * s1)), s1 );

  // casting
  const int Dim = LineType::AmbientDimAtCompileTime;
  typedef typename GetDifferentType<Scalar>::type OtherScalar;
  ParametrizedLine<OtherScalar,Dim> hp1f = l0.template cast<OtherScalar>();
  VERIFY_IS_APPROX(hp1f.template cast<Scalar>(),l0);
  ParametrizedLine<Scalar,Dim> hp1d = l0.template cast<Scalar>();
  VERIFY_IS_APPROX(hp1d.template cast<Scalar>(),l0);
}

template<typename Scalar> void parametrizedline_alignment()
{
  typedef ParametrizedLine<Scalar,4,AutoAlign> Line4a;
  typedef ParametrizedLine<Scalar,4,DontAlign> Line4u;

  EIGEN_ALIGN16 Scalar array1[8];
  EIGEN_ALIGN16 Scalar array2[8];
  EIGEN_ALIGN16 Scalar array3[8+1];
  Scalar* array3u = array3+1;

  Line4a *p1 = ::new(reinterpret_cast<void*>(array1)) Line4a;
  Line4u *p2 = ::new(reinterpret_cast<void*>(array2)) Line4u;
  Line4u *p3 = ::new(reinterpret_cast<void*>(array3u)) Line4u;
  
  p1->origin().setRandom();
  p1->direction().setRandom();
  *p2 = *p1;
  *p3 = *p1;

  VERIFY_IS_APPROX(p1->origin(), p2->origin());
  VERIFY_IS_APPROX(p1->origin(), p3->origin());
  VERIFY_IS_APPROX(p1->direction(), p2->direction());
  VERIFY_IS_APPROX(p1->direction(), p3->direction());
  
  #ifdef EIGEN_VECTORIZE
  VERIFY_RAISES_ASSERT((::new(reinterpret_cast<void*>(array3u)) Line4a));
  #endif
}

void test_geo_parametrizedline()
{
  for(int i = 0; i < g_repeat; i++) {
    CALL_SUBTEST_1( parametrizedline(ParametrizedLine<float,2>()) );
    CALL_SUBTEST_2( parametrizedline(ParametrizedLine<float,3>()) );
    CALL_SUBTEST_3( parametrizedline(ParametrizedLine<double,4>()) );
    CALL_SUBTEST_3( parametrizedline_alignment<double>() );
    CALL_SUBTEST_4( parametrizedline(ParametrizedLine<std::complex<double>,5>()) );
  }
}
