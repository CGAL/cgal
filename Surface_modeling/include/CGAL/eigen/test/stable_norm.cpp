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

template<typename T> bool isNotNaN(const T& x)
{
  return x==x;
}

// workaround aggressive optimization in ICC
template<typename T> EIGEN_DONT_INLINE  T sub(T a, T b) { return a - b; }

template<typename T> bool isFinite(const T& x)
{
  return isNotNaN(sub(x,x));
}

template<typename T> EIGEN_DONT_INLINE T copy(const T& x)
{
  return x;
}

template<typename MatrixType> void stable_norm(const MatrixType& m)
{
  /* this test covers the following files:
     StableNorm.h
  */
  typedef typename MatrixType::Index Index;
  typedef typename MatrixType::Scalar Scalar;
  typedef typename NumTraits<Scalar>::Real RealScalar;

  // Check the basic machine-dependent constants.
  {
    int ibeta, it, iemin, iemax;

    ibeta = std::numeric_limits<RealScalar>::radix;         // base for floating-point numbers
    it    = std::numeric_limits<RealScalar>::digits;        // number of base-beta digits in mantissa
    iemin = std::numeric_limits<RealScalar>::min_exponent;  // minimum exponent
    iemax = std::numeric_limits<RealScalar>::max_exponent;  // maximum exponent

    VERIFY( (!(iemin > 1 - 2*it || 1+it>iemax || (it==2 && ibeta<5) || (it<=4 && ibeta <= 3 ) || it<2))
           && "the stable norm algorithm cannot be guaranteed on this computer");
  }


  Index rows = m.rows();
  Index cols = m.cols();

  Scalar big = internal::random<Scalar>() * (std::numeric_limits<RealScalar>::max() * RealScalar(1e-4));
  Scalar small = internal::random<Scalar>() * (std::numeric_limits<RealScalar>::min() * RealScalar(1e4));

  MatrixType  vzero = MatrixType::Zero(rows, cols),
              vrand = MatrixType::Random(rows, cols),
              vbig(rows, cols),
              vsmall(rows,cols);

  vbig.fill(big);
  vsmall.fill(small);

  VERIFY_IS_MUCH_SMALLER_THAN(vzero.norm(), static_cast<RealScalar>(1));
  VERIFY_IS_APPROX(vrand.stableNorm(),      vrand.norm());
  VERIFY_IS_APPROX(vrand.blueNorm(),        vrand.norm());
  VERIFY_IS_APPROX(vrand.hypotNorm(),       vrand.norm());

  RealScalar size = static_cast<RealScalar>(m.size());

  // test isFinite
  VERIFY(!isFinite( std::numeric_limits<RealScalar>::infinity()));
  VERIFY(!isFinite(internal::sqrt(-internal::abs(big))));

  // test overflow
  VERIFY(isFinite(internal::sqrt(size)*internal::abs(big)));
  VERIFY_IS_NOT_APPROX(internal::sqrt(copy(vbig.squaredNorm())),   internal::abs(internal::sqrt(size)*big)); // here the default norm must fail
  VERIFY_IS_APPROX(vbig.stableNorm(), internal::sqrt(size)*internal::abs(big));
  VERIFY_IS_APPROX(vbig.blueNorm(),   internal::sqrt(size)*internal::abs(big));
  VERIFY_IS_APPROX(vbig.hypotNorm(),  internal::sqrt(size)*internal::abs(big));

  // test underflow
  VERIFY(isFinite(internal::sqrt(size)*internal::abs(small)));
  VERIFY_IS_NOT_APPROX(internal::sqrt(copy(vsmall.squaredNorm())),   internal::abs(internal::sqrt(size)*small)); // here the default norm must fail
  VERIFY_IS_APPROX(vsmall.stableNorm(), internal::sqrt(size)*internal::abs(small));
  VERIFY_IS_APPROX(vsmall.blueNorm(),   internal::sqrt(size)*internal::abs(small));
  VERIFY_IS_APPROX(vsmall.hypotNorm(),  internal::sqrt(size)*internal::abs(small));

// Test compilation of cwise() version
  VERIFY_IS_APPROX(vrand.colwise().stableNorm(),      vrand.colwise().norm());
  VERIFY_IS_APPROX(vrand.colwise().blueNorm(),        vrand.colwise().norm());
  VERIFY_IS_APPROX(vrand.colwise().hypotNorm(),       vrand.colwise().norm());
  VERIFY_IS_APPROX(vrand.rowwise().stableNorm(),      vrand.rowwise().norm());
  VERIFY_IS_APPROX(vrand.rowwise().blueNorm(),        vrand.rowwise().norm());
  VERIFY_IS_APPROX(vrand.rowwise().hypotNorm(),       vrand.rowwise().norm());
}

void test_stable_norm()
{
  for(int i = 0; i < g_repeat; i++) {
    CALL_SUBTEST_1( stable_norm(Matrix<float, 1, 1>()) );
    CALL_SUBTEST_2( stable_norm(Vector4d()) );
    CALL_SUBTEST_3( stable_norm(VectorXd(internal::random<int>(10,2000))) );
    CALL_SUBTEST_4( stable_norm(VectorXf(internal::random<int>(10,2000))) );
    CALL_SUBTEST_5( stable_norm(VectorXcd(internal::random<int>(10,2000))) );
  }
}
