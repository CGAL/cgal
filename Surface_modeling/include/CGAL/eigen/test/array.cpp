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

template<typename ArrayType> void array(const ArrayType& m)
{
  typedef typename ArrayType::Index Index;
  typedef typename ArrayType::Scalar Scalar;
  typedef typename NumTraits<Scalar>::Real RealScalar;
  typedef Array<Scalar, ArrayType::RowsAtCompileTime, 1> ColVectorType;
  typedef Array<Scalar, 1, ArrayType::ColsAtCompileTime> RowVectorType;

  Index rows = m.rows();
  Index cols = m.cols(); 

  ArrayType m1 = ArrayType::Random(rows, cols),
             m2 = ArrayType::Random(rows, cols),
             m3(rows, cols);

  ColVectorType cv1 = ColVectorType::Random(rows);
  RowVectorType rv1 = RowVectorType::Random(cols);

  Scalar  s1 = internal::random<Scalar>(),
          s2 = internal::random<Scalar>();          

  // scalar addition
  VERIFY_IS_APPROX(m1 + s1, s1 + m1);
  VERIFY_IS_APPROX(m1 + s1, ArrayType::Constant(rows,cols,s1) + m1);
  VERIFY_IS_APPROX(s1 - m1, (-m1)+s1 );
  VERIFY_IS_APPROX(m1 - s1, m1 - ArrayType::Constant(rows,cols,s1));
  VERIFY_IS_APPROX(s1 - m1, ArrayType::Constant(rows,cols,s1) - m1);
  VERIFY_IS_APPROX((m1*Scalar(2)) - s2, (m1+m1) - ArrayType::Constant(rows,cols,s2) );
  m3 = m1;
  m3 += s2;
  VERIFY_IS_APPROX(m3, m1 + s2);
  m3 = m1;
  m3 -= s1;
  VERIFY_IS_APPROX(m3, m1 - s1);  
  
  // scalar operators via Maps
  m3 = m1;
  ArrayType::Map(m1.data(), m1.rows(), m1.cols()) -= ArrayType::Map(m2.data(), m2.rows(), m2.cols());
  VERIFY_IS_APPROX(m1, m3 - m2);
  
  m3 = m1;
  ArrayType::Map(m1.data(), m1.rows(), m1.cols()) += ArrayType::Map(m2.data(), m2.rows(), m2.cols());
  VERIFY_IS_APPROX(m1, m3 + m2);
  
  m3 = m1;
  ArrayType::Map(m1.data(), m1.rows(), m1.cols()) *= ArrayType::Map(m2.data(), m2.rows(), m2.cols());
  VERIFY_IS_APPROX(m1, m3 * m2);
  
  m3 = m1;
  m2 = ArrayType::Random(rows,cols);
  m2 = (m2==0).select(1,m2);
  ArrayType::Map(m1.data(), m1.rows(), m1.cols()) /= ArrayType::Map(m2.data(), m2.rows(), m2.cols());  
  VERIFY_IS_APPROX(m1, m3 / m2);

  // reductions
  VERIFY_IS_APPROX(m1.colwise().sum().sum(), m1.sum());
  VERIFY_IS_APPROX(m1.rowwise().sum().sum(), m1.sum());
  if (!internal::isApprox(m1.sum(), (m1+m2).sum(), test_precision<Scalar>()))
      VERIFY_IS_NOT_APPROX(((m1+m2).rowwise().sum()).sum(), m1.sum());
  VERIFY_IS_APPROX(m1.colwise().sum(), m1.colwise().redux(internal::scalar_sum_op<Scalar>()));

  // vector-wise ops
  m3 = m1;
  VERIFY_IS_APPROX(m3.colwise() += cv1, m1.colwise() + cv1);
  m3 = m1;
  VERIFY_IS_APPROX(m3.colwise() -= cv1, m1.colwise() - cv1);
  m3 = m1;
  VERIFY_IS_APPROX(m3.rowwise() += rv1, m1.rowwise() + rv1);
  m3 = m1;
  VERIFY_IS_APPROX(m3.rowwise() -= rv1, m1.rowwise() - rv1);
}

template<typename ArrayType> void comparisons(const ArrayType& m)
{
  typedef typename ArrayType::Index Index;
  typedef typename ArrayType::Scalar Scalar;
  typedef typename NumTraits<Scalar>::Real RealScalar;
  typedef Array<Scalar, ArrayType::RowsAtCompileTime, 1> VectorType;

  Index rows = m.rows();
  Index cols = m.cols();

  Index r = internal::random<Index>(0, rows-1),
        c = internal::random<Index>(0, cols-1);

  ArrayType m1 = ArrayType::Random(rows, cols),
             m2 = ArrayType::Random(rows, cols),
             m3(rows, cols);            

  VERIFY(((m1 + Scalar(1)) > m1).all());
  VERIFY(((m1 - Scalar(1)) < m1).all());
  if (rows*cols>1)
  {
    m3 = m1;
    m3(r,c) += 1;
    VERIFY(! (m1 < m3).all() );
    VERIFY(! (m1 > m3).all() );
  }

  // comparisons to scalar
  VERIFY( (m1 != (m1(r,c)+1) ).any() );
  VERIFY( (m1 > (m1(r,c)-1) ).any() );
  VERIFY( (m1 < (m1(r,c)+1) ).any() );
  VERIFY( (m1 == m1(r,c) ).any() );

  // test Select
  VERIFY_IS_APPROX( (m1<m2).select(m1,m2), m1.cwiseMin(m2) );
  VERIFY_IS_APPROX( (m1>m2).select(m1,m2), m1.cwiseMax(m2) );
  Scalar mid = (m1.cwiseAbs().minCoeff() + m1.cwiseAbs().maxCoeff())/Scalar(2);
  for (int j=0; j<cols; ++j)
  for (int i=0; i<rows; ++i)
    m3(i,j) = internal::abs(m1(i,j))<mid ? 0 : m1(i,j);
  VERIFY_IS_APPROX( (m1.abs()<ArrayType::Constant(rows,cols,mid))
                        .select(ArrayType::Zero(rows,cols),m1), m3);
  // shorter versions:
  VERIFY_IS_APPROX( (m1.abs()<ArrayType::Constant(rows,cols,mid))
                        .select(0,m1), m3);
  VERIFY_IS_APPROX( (m1.abs()>=ArrayType::Constant(rows,cols,mid))
                        .select(m1,0), m3);
  // even shorter version:
  VERIFY_IS_APPROX( (m1.abs()<mid).select(0,m1), m3);

  // count
  VERIFY(((m1.abs()+1)>RealScalar(0.1)).count() == rows*cols);

  typedef Array<typename ArrayType::Index, Dynamic, 1> ArrayOfIndices;

  // TODO allows colwise/rowwise for array
  VERIFY_IS_APPROX(((m1.abs()+1)>RealScalar(0.1)).colwise().count(), ArrayOfIndices::Constant(cols,rows).transpose());
  VERIFY_IS_APPROX(((m1.abs()+1)>RealScalar(0.1)).rowwise().count(), ArrayOfIndices::Constant(rows, cols));
}

template<typename ArrayType> void array_real(const ArrayType& m)
{
  typedef typename ArrayType::Index Index;
  typedef typename ArrayType::Scalar Scalar;
  typedef typename NumTraits<Scalar>::Real RealScalar;

  Index rows = m.rows();
  Index cols = m.cols();

  ArrayType m1 = ArrayType::Random(rows, cols),
             m2 = ArrayType::Random(rows, cols),
             m3(rows, cols);

  // these these are mostly to check possible compilation issues.
  VERIFY_IS_APPROX(m1.sin(), std::sin(m1));
  VERIFY_IS_APPROX(m1.sin(), internal::sin(m1));
  VERIFY_IS_APPROX(m1.cos(), std::cos(m1));
  VERIFY_IS_APPROX(m1.cos(), internal::cos(m1));
  VERIFY_IS_APPROX(m1.asin(), std::asin(m1));
  VERIFY_IS_APPROX(m1.asin(), internal::asin(m1));
  VERIFY_IS_APPROX(m1.acos(), std::acos(m1));
  VERIFY_IS_APPROX(m1.acos(), internal::acos(m1));
  VERIFY_IS_APPROX(m1.tan(), std::tan(m1));
  VERIFY_IS_APPROX(m1.tan(), internal::tan(m1));

  VERIFY_IS_APPROX(internal::cos(m1+RealScalar(3)*m2), internal::cos((m1+RealScalar(3)*m2).eval()));
  VERIFY_IS_APPROX(std::cos(m1+RealScalar(3)*m2), std::cos((m1+RealScalar(3)*m2).eval()));

  VERIFY_IS_APPROX(m1.abs().sqrt(), std::sqrt(std::abs(m1)));
  VERIFY_IS_APPROX(m1.abs().sqrt(), internal::sqrt(internal::abs(m1)));
  VERIFY_IS_APPROX(m1.abs(), internal::sqrt(internal::abs2(m1)));

  VERIFY_IS_APPROX(internal::abs2(internal::real(m1)) + internal::abs2(internal::imag(m1)), internal::abs2(m1));
  VERIFY_IS_APPROX(internal::abs2(std::real(m1)) + internal::abs2(std::imag(m1)), internal::abs2(m1));
  if(!NumTraits<Scalar>::IsComplex)
    VERIFY_IS_APPROX(internal::real(m1), m1);

  VERIFY_IS_APPROX(m1.abs().log(), std::log(std::abs(m1)));
  VERIFY_IS_APPROX(m1.abs().log(), internal::log(internal::abs(m1)));

  VERIFY_IS_APPROX(m1.exp(), std::exp(m1));
  VERIFY_IS_APPROX(m1.exp() * m2.exp(), std::exp(m1+m2));
  VERIFY_IS_APPROX(m1.exp(), internal::exp(m1));
  VERIFY_IS_APPROX(m1.exp() / m2.exp(), std::exp(m1-m2));

  VERIFY_IS_APPROX(m1.pow(2), m1.square());
  VERIFY_IS_APPROX(std::pow(m1,2), m1.square());
  m3 = m1.abs();
  VERIFY_IS_APPROX(m3.pow(RealScalar(0.5)), m3.sqrt());
  VERIFY_IS_APPROX(std::pow(m3,RealScalar(0.5)), m3.sqrt());
}

void test_array()
{
  for(int i = 0; i < g_repeat; i++) {
    CALL_SUBTEST_1( array(Array<float, 1, 1>()) );
    CALL_SUBTEST_2( array(Array22f()) );
    CALL_SUBTEST_3( array(Array44d()) );
    CALL_SUBTEST_4( array(ArrayXXcf(3, 3)) );
    CALL_SUBTEST_5( array(ArrayXXf(8, 12)) );
    CALL_SUBTEST_6( array(ArrayXXi(8, 12)) );
  }
  for(int i = 0; i < g_repeat; i++) {
    CALL_SUBTEST_1( comparisons(Array<float, 1, 1>()) );
    CALL_SUBTEST_2( comparisons(Array22f()) );
    CALL_SUBTEST_3( comparisons(Array44d()) );
    CALL_SUBTEST_5( comparisons(ArrayXXf(8, 12)) );
    CALL_SUBTEST_6( comparisons(ArrayXXi(8, 12)) );
  }
  for(int i = 0; i < g_repeat; i++) {
    CALL_SUBTEST_1( array_real(Array<float, 1, 1>()) );
    CALL_SUBTEST_2( array_real(Array22f()) );
    CALL_SUBTEST_3( array_real(Array44d()) );
    CALL_SUBTEST_5( array_real(ArrayXXf(8, 12)) );
  }

  VERIFY((internal::is_same< internal::global_math_functions_filtering_base<int>::type, int >::value));
  VERIFY((internal::is_same< internal::global_math_functions_filtering_base<float>::type, float >::value));
  VERIFY((internal::is_same< internal::global_math_functions_filtering_base<Array2i>::type, ArrayBase<Array2i> >::value));
  typedef CwiseUnaryOp<internal::scalar_sum_op<double>, ArrayXd > Xpr;
  VERIFY((internal::is_same< internal::global_math_functions_filtering_base<Xpr>::type,
                           ArrayBase<Xpr>
                         >::value));
}
