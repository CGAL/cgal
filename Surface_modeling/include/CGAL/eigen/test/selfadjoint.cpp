// This file is triangularView of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2010 Gael Guennebaud <gael.guennebaud@inria.fr>
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

// This file tests the basic selfadjointView API,
// the related products and decompositions are tested in specific files.

template<typename MatrixType> void selfadjoint(const MatrixType& m)
{
  typedef typename MatrixType::Index Index;
  typedef typename MatrixType::Scalar Scalar;
  typedef typename NumTraits<Scalar>::Real RealScalar;

  Index rows = m.rows();
  Index cols = m.cols();

  MatrixType m1 = MatrixType::Random(rows, cols),
             m3(rows, cols);

  m1.diagonal() = m1.diagonal().real().template cast<Scalar>();

  // check selfadjoint to dense
  m3 = m1.template selfadjointView<Upper>();
  VERIFY_IS_APPROX(MatrixType(m3.template triangularView<Upper>()), MatrixType(m1.template triangularView<Upper>()));
  VERIFY_IS_APPROX(m3, m3.adjoint());


  m3 = m1.template selfadjointView<Lower>();
  VERIFY_IS_APPROX(MatrixType(m3.template triangularView<Lower>()), MatrixType(m1.template triangularView<Lower>()));
  VERIFY_IS_APPROX(m3, m3.adjoint());
}

void bug_159()
{
  Matrix3d m = Matrix3d::Random().selfadjointView<Lower>(); 
}

void test_selfadjoint()
{
  for(int i = 0; i < g_repeat ; i++)
  {
    int s = internal::random<int>(1,20); EIGEN_UNUSED_VARIABLE(s);

    CALL_SUBTEST_1( selfadjoint(Matrix<float, 1, 1>()) );
    CALL_SUBTEST_2( selfadjoint(Matrix<float, 2, 2>()) );
    CALL_SUBTEST_3( selfadjoint(Matrix3cf()) );
    CALL_SUBTEST_4( selfadjoint(MatrixXcd(s,s)) );
    CALL_SUBTEST_5( selfadjoint(Matrix<float,Dynamic,Dynamic,RowMajor>(s, s)) );
  }
  
  CALL_SUBTEST_1( bug_159() );
}
