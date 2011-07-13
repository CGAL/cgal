// This file is part of Eigen, a lightweight C++ template library
// for linear algebra. Eigen itself is part of the KDE project.
//
// Copyright (C) 2008 Gael Guennebaud <g.gael@free.fr>
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

  const int dim = _line.dim();
  typedef typename LineType::Scalar Scalar;
  typedef typename NumTraits<Scalar>::Real RealScalar;
  typedef Matrix<Scalar, LineType::AmbientDimAtCompileTime, 1> VectorType;
  typedef Matrix<Scalar, LineType::AmbientDimAtCompileTime,
                         LineType::AmbientDimAtCompileTime> MatrixType;

  VectorType p0 = VectorType::Random(dim);
  VectorType p1 = VectorType::Random(dim);

  VectorType d0 = VectorType::Random(dim).normalized();

  LineType l0(p0, d0);

  Scalar s0 = ei_random<Scalar>();
  Scalar s1 = ei_abs(ei_random<Scalar>());

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

void test_eigen2_parametrizedline()
{
  for(int i = 0; i < g_repeat; i++) {
    CALL_SUBTEST_1( parametrizedline(ParametrizedLine<float,2>()) );
    CALL_SUBTEST_2( parametrizedline(ParametrizedLine<float,3>()) );
    CALL_SUBTEST_3( parametrizedline(ParametrizedLine<double,4>()) );
    CALL_SUBTEST_4( parametrizedline(ParametrizedLine<std::complex<double>,5>()) );
  }
}
