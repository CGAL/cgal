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

template<typename MatrixType> void syrk(const MatrixType& m)
{
  typedef typename MatrixType::Index Index;
  typedef typename MatrixType::Scalar Scalar;
  typedef typename NumTraits<Scalar>::Real RealScalar;
  typedef Matrix<Scalar, MatrixType::ColsAtCompileTime, Dynamic> Rhs1;
  typedef Matrix<Scalar, Dynamic, MatrixType::RowsAtCompileTime> Rhs2;
  typedef Matrix<Scalar, MatrixType::ColsAtCompileTime, Dynamic,RowMajor> Rhs3;

  Index rows = m.rows();
  Index cols = m.cols();

  MatrixType m1 = MatrixType::Random(rows, cols),
             m2 = MatrixType::Random(rows, cols);

  Rhs1 rhs1 = Rhs1::Random(internal::random<int>(1,320), cols);
  Rhs2 rhs2 = Rhs2::Random(rows, internal::random<int>(1,320));
  Rhs3 rhs3 = Rhs3::Random(internal::random<int>(1,320), rows);

  Scalar s1 = internal::random<Scalar>();
  
  Index c = internal::random<Index>(0,cols-1);

  m2.setZero();
  VERIFY_IS_APPROX((m2.template selfadjointView<Lower>().rankUpdate(rhs2,s1)._expression()),
                   ((s1 * rhs2 * rhs2.adjoint()).eval().template triangularView<Lower>().toDenseMatrix()));

  m2.setZero();
  VERIFY_IS_APPROX(m2.template selfadjointView<Upper>().rankUpdate(rhs2,s1)._expression(),
                   (s1 * rhs2 * rhs2.adjoint()).eval().template triangularView<Upper>().toDenseMatrix());

  m2.setZero();
  VERIFY_IS_APPROX(m2.template selfadjointView<Lower>().rankUpdate(rhs1.adjoint(),s1)._expression(),
                   (s1 * rhs1.adjoint() * rhs1).eval().template triangularView<Lower>().toDenseMatrix());

  m2.setZero();
  VERIFY_IS_APPROX(m2.template selfadjointView<Upper>().rankUpdate(rhs1.adjoint(),s1)._expression(),
                   (s1 * rhs1.adjoint() * rhs1).eval().template triangularView<Upper>().toDenseMatrix());

  m2.setZero();
  VERIFY_IS_APPROX(m2.template selfadjointView<Lower>().rankUpdate(rhs3.adjoint(),s1)._expression(),
                   (s1 * rhs3.adjoint() * rhs3).eval().template triangularView<Lower>().toDenseMatrix());

  m2.setZero();
  VERIFY_IS_APPROX(m2.template selfadjointView<Upper>().rankUpdate(rhs3.adjoint(),s1)._expression(),
                   (s1 * rhs3.adjoint() * rhs3).eval().template triangularView<Upper>().toDenseMatrix());
                   
  m2.setZero();
  VERIFY_IS_APPROX((m2.template selfadjointView<Lower>().rankUpdate(m1.col(c),s1)._expression()),
                   ((s1 * m1.col(c) * m1.col(c).adjoint()).eval().template triangularView<Lower>().toDenseMatrix()));
                   
  m2.setZero();
  VERIFY_IS_APPROX((m2.template selfadjointView<Upper>().rankUpdate(m1.col(c),s1)._expression()),
                   ((s1 * m1.col(c) * m1.col(c).adjoint()).eval().template triangularView<Upper>().toDenseMatrix()));
  
  m2.setZero();
  VERIFY_IS_APPROX((m2.template selfadjointView<Lower>().rankUpdate(m1.col(c).conjugate(),s1)._expression()),
                   ((s1 * m1.col(c).conjugate() * m1.col(c).conjugate().adjoint()).eval().template triangularView<Lower>().toDenseMatrix()));
                   
  m2.setZero();
  VERIFY_IS_APPROX((m2.template selfadjointView<Upper>().rankUpdate(m1.col(c).conjugate(),s1)._expression()),
                   ((s1 * m1.col(c).conjugate() * m1.col(c).conjugate().adjoint()).eval().template triangularView<Upper>().toDenseMatrix()));
  
  m2.setZero();
  VERIFY_IS_APPROX((m2.template selfadjointView<Lower>().rankUpdate(m1.row(c),s1)._expression()),
                   ((s1 * m1.row(c).transpose() * m1.row(c).transpose().adjoint()).eval().template triangularView<Lower>().toDenseMatrix()));
  
  m2.setZero();
  VERIFY_IS_APPROX((m2.template selfadjointView<Upper>().rankUpdate(m1.row(c).adjoint(),s1)._expression()),
                   ((s1 * m1.row(c).adjoint() * m1.row(c).adjoint().adjoint()).eval().template triangularView<Upper>().toDenseMatrix()));
}

void test_product_syrk()
{
  for(int i = 0; i < g_repeat ; i++)
  {
    int s;
    s = internal::random<int>(1,320);
    CALL_SUBTEST_1( syrk(MatrixXf(s, s)) );
    s = internal::random<int>(1,320);
    CALL_SUBTEST_2( syrk(MatrixXd(s, s)) );
    s = internal::random<int>(1,200);
    CALL_SUBTEST_3( syrk(MatrixXcf(s, s)) );
    s = internal::random<int>(1,200);
    CALL_SUBTEST_4( syrk(MatrixXcd(s, s)) );
  }
}
