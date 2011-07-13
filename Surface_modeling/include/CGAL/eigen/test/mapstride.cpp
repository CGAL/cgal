// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2010 Benoit Jacob <jacob.benoit.1@gmail.com>
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

template<int Alignment,typename VectorType> void map_class_vector(const VectorType& m)
{
  typedef typename VectorType::Index Index;
  typedef typename VectorType::Scalar Scalar;

  Index size = m.size();

  VectorType v = VectorType::Random(size);

  Index arraysize = 3*size;
  
  Scalar* a_array = internal::aligned_new<Scalar>(arraysize+1);
  Scalar* array = a_array;
  if(Alignment!=Aligned)
    array = (Scalar*)(ptrdiff_t(a_array) + (internal::packet_traits<Scalar>::AlignedOnScalar?sizeof(Scalar):sizeof(typename NumTraits<Scalar>::Real)));

  {
    Map<VectorType, Alignment, InnerStride<3> > map(array, size);
    map = v;
    for(int i = 0; i < size; ++i)
    {
      VERIFY(array[3*i] == v[i]);
      VERIFY(map[i] == v[i]);
    }
  }

  {
    Map<VectorType, Unaligned, InnerStride<Dynamic> > map(array, size, InnerStride<Dynamic>(2));
    map = v;
    for(int i = 0; i < size; ++i)
    {
      VERIFY(array[2*i] == v[i]);
      VERIFY(map[i] == v[i]);
    }
  }

  internal::aligned_delete(a_array, arraysize+1);
}

template<int Alignment,typename MatrixType> void map_class_matrix(const MatrixType& _m)
{
  typedef typename MatrixType::Index Index;
  typedef typename MatrixType::Scalar Scalar;

  Index rows = _m.rows(), cols = _m.cols();

  MatrixType m = MatrixType::Random(rows,cols);

  Index arraysize = 2*(rows+4)*(cols+4);

  Scalar* a_array = internal::aligned_new<Scalar>(arraysize+1);
  Scalar* array = a_array;
  if(Alignment!=Aligned)
    array = (Scalar*)(ptrdiff_t(a_array) + (internal::packet_traits<Scalar>::AlignedOnScalar?sizeof(Scalar):sizeof(typename NumTraits<Scalar>::Real)));

  // test no inner stride and some dynamic outer stride
  {
    Map<MatrixType, Alignment, OuterStride<Dynamic> > map(array, rows, cols, OuterStride<Dynamic>(m.innerSize()+1));
    map = m;
    VERIFY(map.outerStride() == map.innerSize()+1);
    for(int i = 0; i < m.outerSize(); ++i)
      for(int j = 0; j < m.innerSize(); ++j)
      {
        VERIFY(array[map.outerStride()*i+j] == m.coeffByOuterInner(i,j));
        VERIFY(map.coeffByOuterInner(i,j) == m.coeffByOuterInner(i,j));
      }
  }

  // test no inner stride and an outer stride of +4. This is quite important as for fixed-size matrices,
  // this allows to hit the special case where it's vectorizable.
  {
    enum {
      InnerSize = MatrixType::InnerSizeAtCompileTime,
      OuterStrideAtCompileTime = InnerSize==Dynamic ? Dynamic : InnerSize+4
    };
    Map<MatrixType, Alignment, OuterStride<OuterStrideAtCompileTime> >
      map(array, rows, cols, OuterStride<OuterStrideAtCompileTime>(m.innerSize()+4));
    map = m;
    VERIFY(map.outerStride() == map.innerSize()+4);
    for(int i = 0; i < m.outerSize(); ++i)
      for(int j = 0; j < m.innerSize(); ++j)
      {
        VERIFY(array[map.outerStride()*i+j] == m.coeffByOuterInner(i,j));
        VERIFY(map.coeffByOuterInner(i,j) == m.coeffByOuterInner(i,j));
      }
  }

  // test both inner stride and outer stride
  {
    Map<MatrixType, Alignment, Stride<Dynamic,Dynamic> > map(array, rows, cols, Stride<Dynamic,Dynamic>(2*m.innerSize()+1, 2));
    map = m;
    VERIFY(map.outerStride() == 2*map.innerSize()+1);
    VERIFY(map.innerStride() == 2);
    for(int i = 0; i < m.outerSize(); ++i)
      for(int j = 0; j < m.innerSize(); ++j)
      {
        VERIFY(array[map.outerStride()*i+map.innerStride()*j] == m.coeffByOuterInner(i,j));
        VERIFY(map.coeffByOuterInner(i,j) == m.coeffByOuterInner(i,j));
      }
  }

  internal::aligned_delete(a_array, arraysize+1);
}

void test_mapstride()
{
  for(int i = 0; i < g_repeat; i++) {
    EIGEN_UNUSED int maxn = 30;
    CALL_SUBTEST_1( map_class_vector<Aligned>(Matrix<float, 1, 1>()) );
    CALL_SUBTEST_1( map_class_vector<Unaligned>(Matrix<float, 1, 1>()) );
    CALL_SUBTEST_2( map_class_vector<Aligned>(Vector4d()) );
    CALL_SUBTEST_2( map_class_vector<Unaligned>(Vector4d()) );
    CALL_SUBTEST_3( map_class_vector<Aligned>(RowVector4f()) );
    CALL_SUBTEST_3( map_class_vector<Unaligned>(RowVector4f()) );
    CALL_SUBTEST_4( map_class_vector<Aligned>(VectorXcf(internal::random<int>(1,maxn))) );
    CALL_SUBTEST_4( map_class_vector<Unaligned>(VectorXcf(internal::random<int>(1,maxn))) );
    CALL_SUBTEST_5( map_class_vector<Aligned>(VectorXi(internal::random<int>(1,maxn))) );
    CALL_SUBTEST_5( map_class_vector<Unaligned>(VectorXi(internal::random<int>(1,maxn))) );

    CALL_SUBTEST_1( map_class_matrix<Aligned>(Matrix<float, 1, 1>()) );
    CALL_SUBTEST_1( map_class_matrix<Unaligned>(Matrix<float, 1, 1>()) );
    CALL_SUBTEST_2( map_class_matrix<Aligned>(Matrix4d()) );
    CALL_SUBTEST_2( map_class_matrix<Unaligned>(Matrix4d()) );
    CALL_SUBTEST_3( map_class_matrix<Aligned>(Matrix<float,3,5>()) );
    CALL_SUBTEST_3( map_class_matrix<Unaligned>(Matrix<float,3,5>()) );
    CALL_SUBTEST_3( map_class_matrix<Aligned>(Matrix<float,4,8>()) );
    CALL_SUBTEST_3( map_class_matrix<Unaligned>(Matrix<float,4,8>()) );
    CALL_SUBTEST_4( map_class_matrix<Aligned>(MatrixXcf(internal::random<int>(1,maxn),internal::random<int>(1,maxn))) );
    CALL_SUBTEST_4( map_class_matrix<Unaligned>(MatrixXcf(internal::random<int>(1,maxn),internal::random<int>(1,maxn))) );
    CALL_SUBTEST_5( map_class_matrix<Aligned>(MatrixXi(internal::random<int>(1,maxn),internal::random<int>(1,maxn))) );
    CALL_SUBTEST_5( map_class_matrix<Unaligned>(MatrixXi(internal::random<int>(1,maxn),internal::random<int>(1,maxn))) );
    CALL_SUBTEST_6( map_class_matrix<Aligned>(MatrixXcd(internal::random<int>(1,maxn),internal::random<int>(1,maxn))) );
    CALL_SUBTEST_6( map_class_matrix<Unaligned>(MatrixXcd(internal::random<int>(1,maxn),internal::random<int>(1,maxn))) );
  }
}
