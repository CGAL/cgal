// This file is part of Eigen, a lightweight C++ template library
// for linear algebra. Eigen itself is part of the KDE project.
//
// Copyright (C) 2007-2010 Benoit Jacob <jacob.benoit.1@gmail.com>
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

template<typename VectorType> void map_class_vector(const VectorType& m)
{
  typedef typename VectorType::Scalar Scalar;

  int size = m.size();

  // test Map.h
  Scalar* array1 = ei_aligned_new<Scalar>(size);
  Scalar* array2 = ei_aligned_new<Scalar>(size);
  Scalar* array3 = new Scalar[size+1];
  Scalar* array3unaligned = std::size_t(array3)%16 == 0 ? array3+1 : array3;
  
  Map<VectorType, Aligned>(array1, size) = VectorType::Random(size);
  Map<VectorType>(array2, size) = Map<VectorType>(array1, size);
  Map<VectorType>(array3unaligned, size) = Map<VectorType>((const Scalar*)array1, size); // test non-const-correctness support in eigen2
  VectorType ma1 = Map<VectorType>(array1, size);
  VectorType ma2 = Map<VectorType, Aligned>(array2, size);
  VectorType ma3 = Map<VectorType>(array3unaligned, size);
  VERIFY_IS_APPROX(ma1, ma2);
  VERIFY_IS_APPROX(ma1, ma3);
  
  ei_aligned_delete(array1, size);
  ei_aligned_delete(array2, size);
  delete[] array3;
}

template<typename MatrixType> void map_class_matrix(const MatrixType& m)
{
  typedef typename MatrixType::Scalar Scalar;

  int rows = m.rows(), cols = m.cols(), size = rows*cols;

  // test Map.h
  Scalar* array1 = ei_aligned_new<Scalar>(size);
  for(int i = 0; i < size; i++) array1[i] = Scalar(1);
  Scalar* array2 = ei_aligned_new<Scalar>(size);
  for(int i = 0; i < size; i++) array2[i] = Scalar(1);
  Scalar* array3 = new Scalar[size+1];
  for(int i = 0; i < size+1; i++) array3[i] = Scalar(1);
  Scalar* array3unaligned = std::size_t(array3)%16 == 0 ? array3+1 : array3;
  Map<MatrixType, Aligned>(array1, rows, cols) = MatrixType::Ones(rows,cols);
  Map<MatrixType>(array2, rows, cols) = Map<MatrixType>((const Scalar*)array1, rows, cols); // test non-const-correctness support in eigen2
  Map<MatrixType>(array3unaligned, rows, cols) = Map<MatrixType>(array1, rows, cols);
  MatrixType ma1 = Map<MatrixType>(array1, rows, cols);
  MatrixType ma2 = Map<MatrixType, Aligned>(array2, rows, cols);
  VERIFY_IS_APPROX(ma1, ma2);
  MatrixType ma3 = Map<MatrixType>(array3unaligned, rows, cols);
  VERIFY_IS_APPROX(ma1, ma3);
  
  ei_aligned_delete(array1, size);
  ei_aligned_delete(array2, size);
  delete[] array3;
}

template<typename VectorType> void map_static_methods(const VectorType& m)
{
  typedef typename VectorType::Scalar Scalar;

  int size = m.size();

  // test Map.h
  Scalar* array1 = ei_aligned_new<Scalar>(size);
  Scalar* array2 = ei_aligned_new<Scalar>(size);
  Scalar* array3 = new Scalar[size+1];
  Scalar* array3unaligned = std::size_t(array3)%16 == 0 ? array3+1 : array3;
  
  VectorType::MapAligned(array1, size) = VectorType::Random(size);
  VectorType::Map(array2, size) = VectorType::Map(array1, size);
  VectorType::Map(array3unaligned, size) = VectorType::Map(array1, size);
  VectorType ma1 = VectorType::Map(array1, size);
  VectorType ma2 = VectorType::MapAligned(array2, size);
  VectorType ma3 = VectorType::Map(array3unaligned, size);
  VERIFY_IS_APPROX(ma1, ma2);
  VERIFY_IS_APPROX(ma1, ma3);
  
  ei_aligned_delete(array1, size);
  ei_aligned_delete(array2, size);
  delete[] array3;
}


void test_eigen2_map()
{
  for(int i = 0; i < g_repeat; i++) {
    CALL_SUBTEST_1( map_class_vector(Matrix<float, 1, 1>()) );
    CALL_SUBTEST_2( map_class_vector(Vector4d()) );
    CALL_SUBTEST_3( map_class_vector(RowVector4f()) );
    CALL_SUBTEST_4( map_class_vector(VectorXcf(8)) );
    CALL_SUBTEST_5( map_class_vector(VectorXi(12)) );

    CALL_SUBTEST_1( map_class_matrix(Matrix<float, 1, 1>()) );
    CALL_SUBTEST_2( map_class_matrix(Matrix4d()) );
    CALL_SUBTEST_6( map_class_matrix(Matrix<float,3,5>()) );
    CALL_SUBTEST_4( map_class_matrix(MatrixXcf(ei_random<int>(1,10),ei_random<int>(1,10))) );
    CALL_SUBTEST_5( map_class_matrix(MatrixXi(ei_random<int>(1,10),ei_random<int>(1,10))) );

    CALL_SUBTEST_1( map_static_methods(Matrix<double, 1, 1>()) );
    CALL_SUBTEST_2( map_static_methods(Vector3f()) );
    CALL_SUBTEST_7( map_static_methods(RowVector3d()) );
    CALL_SUBTEST_4( map_static_methods(VectorXcd(8)) );
    CALL_SUBTEST_5( map_static_methods(VectorXf(12)) );
  }
}
