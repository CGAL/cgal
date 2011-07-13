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

template<typename HyperplaneType> void hyperplane(const HyperplaneType& _plane)
{
  /* this test covers the following files:
     Hyperplane.h
  */

  const int dim = _plane.dim();
  typedef typename HyperplaneType::Scalar Scalar;
  typedef typename NumTraits<Scalar>::Real RealScalar;
  typedef Matrix<Scalar, HyperplaneType::AmbientDimAtCompileTime, 1> VectorType;
  typedef Matrix<Scalar, HyperplaneType::AmbientDimAtCompileTime,
                         HyperplaneType::AmbientDimAtCompileTime> MatrixType;

  VectorType p0 = VectorType::Random(dim);
  VectorType p1 = VectorType::Random(dim);

  VectorType n0 = VectorType::Random(dim).normalized();
  VectorType n1 = VectorType::Random(dim).normalized();

  HyperplaneType pl0(n0, p0);
  HyperplaneType pl1(n1, p1);
  HyperplaneType pl2 = pl1;

  Scalar s0 = ei_random<Scalar>();
  Scalar s1 = ei_random<Scalar>();

  VERIFY_IS_APPROX( n1.eigen2_dot(n1), Scalar(1) );

  VERIFY_IS_MUCH_SMALLER_THAN( pl0.absDistance(p0), Scalar(1) );
  VERIFY_IS_APPROX( pl1.signedDistance(p1 + n1 * s0), s0 );
  VERIFY_IS_MUCH_SMALLER_THAN( pl1.signedDistance(pl1.projection(p0)), Scalar(1) );
  VERIFY_IS_MUCH_SMALLER_THAN( pl1.absDistance(p1 +  pl1.normal().unitOrthogonal() * s1), Scalar(1) );

  // transform
  if (!NumTraits<Scalar>::IsComplex)
  {
    MatrixType rot = MatrixType::Random(dim,dim).qr().matrixQ();
    Scaling<Scalar,HyperplaneType::AmbientDimAtCompileTime> scaling(VectorType::Random());
    Translation<Scalar,HyperplaneType::AmbientDimAtCompileTime> translation(VectorType::Random());

    pl2 = pl1;
    VERIFY_IS_MUCH_SMALLER_THAN( pl2.transform(rot).absDistance(rot * p1), Scalar(1) );
    pl2 = pl1;
    VERIFY_IS_MUCH_SMALLER_THAN( pl2.transform(rot,Isometry).absDistance(rot * p1), Scalar(1) );
    pl2 = pl1;
    VERIFY_IS_MUCH_SMALLER_THAN( pl2.transform(rot*scaling).absDistance((rot*scaling) * p1), Scalar(1) );
    pl2 = pl1;
    VERIFY_IS_MUCH_SMALLER_THAN( pl2.transform(rot*scaling*translation)
                                 .absDistance((rot*scaling*translation) * p1), Scalar(1) );
    pl2 = pl1;
    VERIFY_IS_MUCH_SMALLER_THAN( pl2.transform(rot*translation,Isometry)
                                 .absDistance((rot*translation) * p1), Scalar(1) );
  }

  // casting
  const int Dim = HyperplaneType::AmbientDimAtCompileTime;
  typedef typename GetDifferentType<Scalar>::type OtherScalar;
  Hyperplane<OtherScalar,Dim> hp1f = pl1.template cast<OtherScalar>();
  VERIFY_IS_APPROX(hp1f.template cast<Scalar>(),pl1);
  Hyperplane<Scalar,Dim> hp1d = pl1.template cast<Scalar>();
  VERIFY_IS_APPROX(hp1d.template cast<Scalar>(),pl1);
}

template<typename Scalar> void lines()
{
  typedef Hyperplane<Scalar, 2> HLine;
  typedef ParametrizedLine<Scalar, 2> PLine;
  typedef Matrix<Scalar,2,1> Vector;
  typedef Matrix<Scalar,3,1> CoeffsType;

  for(int i = 0; i < 10; i++)
  {
    Vector center = Vector::Random();
    Vector u = Vector::Random();
    Vector v = Vector::Random();
    Scalar a = ei_random<Scalar>();
    while (ei_abs(a-1) < 1e-4) a = ei_random<Scalar>();
    while (u.norm() < 1e-4) u = Vector::Random();
    while (v.norm() < 1e-4) v = Vector::Random();

    HLine line_u = HLine::Through(center + u, center + a*u);
    HLine line_v = HLine::Through(center + v, center + a*v);

    // the line equations should be normalized so that a^2+b^2=1
    VERIFY_IS_APPROX(line_u.normal().norm(), Scalar(1));
    VERIFY_IS_APPROX(line_v.normal().norm(), Scalar(1));

    Vector result = line_u.intersection(line_v);

    // the lines should intersect at the point we called "center"
    VERIFY_IS_APPROX(result, center);

    // check conversions between two types of lines
    PLine pl(line_u); // gcc 3.3 will commit suicide if we don't name this variable
    CoeffsType converted_coeffs(HLine(pl).coeffs());
    converted_coeffs *= line_u.coeffs()(0)/converted_coeffs(0);
    VERIFY(line_u.coeffs().isApprox(converted_coeffs));
  }
}

void test_eigen2_hyperplane()
{
  for(int i = 0; i < g_repeat; i++) {
    CALL_SUBTEST_1( hyperplane(Hyperplane<float,2>()) );
    CALL_SUBTEST_2( hyperplane(Hyperplane<float,3>()) );
    CALL_SUBTEST_3( hyperplane(Hyperplane<double,4>()) );
    CALL_SUBTEST_4( hyperplane(Hyperplane<std::complex<double>,5>()) );
    CALL_SUBTEST_5( lines<float>() );
    CALL_SUBTEST_6( lines<double>() );
  }
}
