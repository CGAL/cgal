// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2010 Hauke Heibel <hauke.heibel@gmail.com>
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

template <typename MatrixType> void run_nesting_ops(const MatrixType& _m)
{
  typename MatrixType::Nested m(_m);
  typedef typename MatrixType::Scalar Scalar;

#ifdef NDEBUG
  const bool is_debug = false;
#else
  const bool is_debug = true;
#endif

  // Make really sure that we are in debug mode! We don't want any type of
  // inlining for these tests to pass.
  VERIFY(is_debug);

  // The only intention of these tests is to ensure that this code does
  // not trigger any asserts or segmentation faults... more to come.
  VERIFY_IS_APPROX( (m.transpose() * m).diagonal().sum(), (m.transpose() * m).diagonal().sum() );
  VERIFY_IS_APPROX( (m.transpose() * m).diagonal().array().abs().sum(), (m.transpose() * m).diagonal().array().abs().sum() );

  VERIFY_IS_APPROX( (m.transpose() * m).array().abs().sum(), (m.transpose() * m).array().abs().sum() );
}

void test_nesting_ops()
{
  CALL_SUBTEST_1(run_nesting_ops(MatrixXf::Random(25,25)));
  CALL_SUBTEST_2(run_nesting_ops(MatrixXd::Random(25,25)));
  CALL_SUBTEST_3(run_nesting_ops(Matrix4f::Random()));
  CALL_SUBTEST_4(run_nesting_ops(Matrix4d::Random()));
}
