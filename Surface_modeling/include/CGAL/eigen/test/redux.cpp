// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
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

template<typename MatrixType> void matrixRedux(const MatrixType& m)
{
  typedef typename MatrixType::Index Index;
  typedef typename MatrixType::Scalar Scalar;
  typedef typename MatrixType::RealScalar RealScalar;

  Index rows = m.rows();
  Index cols = m.cols();

  MatrixType m1 = MatrixType::Random(rows, cols);

  VERIFY_IS_MUCH_SMALLER_THAN(MatrixType::Zero(rows, cols).sum(), Scalar(1));
  VERIFY_IS_APPROX(MatrixType::Ones(rows, cols).sum(), Scalar(float(rows*cols))); // the float() here to shut up excessive MSVC warning about int->complex conversion being lossy
  Scalar s(0), p(1), minc(internal::real(m1.coeff(0))), maxc(internal::real(m1.coeff(0)));
  for(int j = 0; j < cols; j++)
  for(int i = 0; i < rows; i++)
  {
    s += m1(i,j);
    p *= m1(i,j);
    minc = std::min(internal::real(minc), internal::real(m1(i,j)));
    maxc = std::max(internal::real(maxc), internal::real(m1(i,j)));
  }
  const Scalar mean = s/Scalar(RealScalar(rows*cols));

  VERIFY_IS_APPROX(m1.sum(), s);
  VERIFY_IS_APPROX(m1.mean(), mean);
  VERIFY_IS_APPROX(m1.prod(), p);
  VERIFY_IS_APPROX(m1.real().minCoeff(), internal::real(minc));
  VERIFY_IS_APPROX(m1.real().maxCoeff(), internal::real(maxc));

  // test slice vectorization assuming assign is ok
  Index r0 = internal::random<Index>(0,rows-1);
  Index c0 = internal::random<Index>(0,cols-1);
  Index r1 = internal::random<Index>(r0+1,rows)-r0;
  Index c1 = internal::random<Index>(c0+1,cols)-c0;
  VERIFY_IS_APPROX(m1.block(r0,c0,r1,c1).sum(), m1.block(r0,c0,r1,c1).eval().sum());
  VERIFY_IS_APPROX(m1.block(r0,c0,r1,c1).mean(), m1.block(r0,c0,r1,c1).eval().mean());
  VERIFY_IS_APPROX(m1.block(r0,c0,r1,c1).prod(), m1.block(r0,c0,r1,c1).eval().prod());
  VERIFY_IS_APPROX(m1.block(r0,c0,r1,c1).real().minCoeff(), m1.block(r0,c0,r1,c1).real().eval().minCoeff());
  VERIFY_IS_APPROX(m1.block(r0,c0,r1,c1).real().maxCoeff(), m1.block(r0,c0,r1,c1).real().eval().maxCoeff());
  
  // test empty objects
  VERIFY_IS_APPROX(m1.block(r0,c0,0,0).sum(),   Scalar(0));
  VERIFY_IS_APPROX(m1.block(r0,c0,0,0).prod(),  Scalar(1));
}

template<typename VectorType> void vectorRedux(const VectorType& w)
{
  typedef typename VectorType::Index Index;
  typedef typename VectorType::Scalar Scalar;
  typedef typename NumTraits<Scalar>::Real RealScalar;
  Index size = w.size();

  VectorType v = VectorType::Random(size);
  for(int i = 1; i < size; i++)
  {
    Scalar s(0), p(1);
    RealScalar minc(internal::real(v.coeff(0))), maxc(internal::real(v.coeff(0)));
    for(int j = 0; j < i; j++)
    {
      s += v[j];
      p *= v[j];
      minc = std::min(minc, internal::real(v[j]));
      maxc = std::max(maxc, internal::real(v[j]));
    }
    VERIFY_IS_MUCH_SMALLER_THAN(internal::abs(s - v.head(i).sum()), Scalar(1));
    VERIFY_IS_APPROX(p, v.head(i).prod());
    VERIFY_IS_APPROX(minc, v.real().head(i).minCoeff());
    VERIFY_IS_APPROX(maxc, v.real().head(i).maxCoeff());
  }

  for(int i = 0; i < size-1; i++)
  {
    Scalar s(0), p(1);
    RealScalar minc(internal::real(v.coeff(i))), maxc(internal::real(v.coeff(i)));
    for(int j = i; j < size; j++)
    {
      s += v[j];
      p *= v[j];
      minc = std::min(minc, internal::real(v[j]));
      maxc = std::max(maxc, internal::real(v[j]));
    }
    VERIFY_IS_MUCH_SMALLER_THAN(internal::abs(s - v.tail(size-i).sum()), Scalar(1));
    VERIFY_IS_APPROX(p, v.tail(size-i).prod());
    VERIFY_IS_APPROX(minc, v.real().tail(size-i).minCoeff());
    VERIFY_IS_APPROX(maxc, v.real().tail(size-i).maxCoeff());
  }

  for(int i = 0; i < size/2; i++)
  {
    Scalar s(0), p(1);
    RealScalar minc(internal::real(v.coeff(i))), maxc(internal::real(v.coeff(i)));
    for(int j = i; j < size-i; j++)
    {
      s += v[j];
      p *= v[j];
      minc = std::min(minc, internal::real(v[j]));
      maxc = std::max(maxc, internal::real(v[j]));
    }
    VERIFY_IS_MUCH_SMALLER_THAN(internal::abs(s - v.segment(i, size-2*i).sum()), Scalar(1));
    VERIFY_IS_APPROX(p, v.segment(i, size-2*i).prod());
    VERIFY_IS_APPROX(minc, v.real().segment(i, size-2*i).minCoeff());
    VERIFY_IS_APPROX(maxc, v.real().segment(i, size-2*i).maxCoeff());
  }
  
  // test empty objects
  VERIFY_IS_APPROX(v.head(0).sum(),   Scalar(0));
  VERIFY_IS_APPROX(v.tail(0).prod(),  Scalar(1));
  VERIFY_RAISES_ASSERT(v.head(0).mean());
  VERIFY_RAISES_ASSERT(v.head(0).minCoeff());
  VERIFY_RAISES_ASSERT(v.head(0).maxCoeff());
}

void test_redux()
{
  for(int i = 0; i < g_repeat; i++) {
    CALL_SUBTEST_1( matrixRedux(Matrix<float, 1, 1>()) );
    CALL_SUBTEST_1( matrixRedux(Array<float, 1, 1>()) );
    CALL_SUBTEST_2( matrixRedux(Matrix2f()) );
    CALL_SUBTEST_2( matrixRedux(Array2f()) );
    CALL_SUBTEST_3( matrixRedux(Matrix4d()) );
    CALL_SUBTEST_3( matrixRedux(Array4d()) );
    CALL_SUBTEST_4( matrixRedux(MatrixXcf(3, 3)) );
    CALL_SUBTEST_4( matrixRedux(ArrayXXcf(3, 3)) );
    CALL_SUBTEST_5( matrixRedux(MatrixXd(8, 12)) );
    CALL_SUBTEST_5( matrixRedux(ArrayXXd(8, 12)) );
    CALL_SUBTEST_6( matrixRedux(MatrixXi(8, 12)) );
    CALL_SUBTEST_6( matrixRedux(ArrayXXi(8, 12)) );
  }
  for(int i = 0; i < g_repeat; i++) {
    CALL_SUBTEST_7( vectorRedux(Vector4f()) );
    CALL_SUBTEST_7( vectorRedux(Array4f()) );
    CALL_SUBTEST_5( vectorRedux(VectorXd(10)) );
    CALL_SUBTEST_5( vectorRedux(ArrayXd(10)) );
    CALL_SUBTEST_8( vectorRedux(VectorXf(33)) );
    CALL_SUBTEST_8( vectorRedux(ArrayXf(33)) );
  }
}
