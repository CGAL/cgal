// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2009 Benoit Jacob <jacob.benoit.1@gmail.com>
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

template<typename Scalar>
void test_first_aligned_helper(Scalar *array, int size)
{
  const int packet_size = sizeof(Scalar) * internal::packet_traits<Scalar>::size;
  VERIFY(((size_t(array) + sizeof(Scalar) * internal::first_aligned(array, size)) % packet_size) == 0);
}

template<typename Scalar>
void test_none_aligned_helper(Scalar *array, int size)
{
  EIGEN_UNUSED_VARIABLE(array);
  EIGEN_UNUSED_VARIABLE(size);
  VERIFY(internal::packet_traits<Scalar>::size == 1 || internal::first_aligned(array, size) == size);
}

struct some_non_vectorizable_type { float x; };

void test_first_aligned()
{
  EIGEN_ALIGN16 float array_float[100];
  test_first_aligned_helper(array_float, 50);
  test_first_aligned_helper(array_float+1, 50);
  test_first_aligned_helper(array_float+2, 50);
  test_first_aligned_helper(array_float+3, 50);
  test_first_aligned_helper(array_float+4, 50);
  test_first_aligned_helper(array_float+5, 50);
  
  EIGEN_ALIGN16 double array_double[100];
  test_first_aligned_helper(array_double, 50);
  test_first_aligned_helper(array_double+1, 50);
  test_first_aligned_helper(array_double+2, 50);
  
  double *array_double_plus_4_bytes = (double*)(size_t(array_double)+4);
  test_none_aligned_helper(array_double_plus_4_bytes, 50);
  test_none_aligned_helper(array_double_plus_4_bytes+1, 50);
  
  some_non_vectorizable_type array_nonvec[100];
  test_first_aligned_helper(array_nonvec, 100);
  test_none_aligned_helper(array_nonvec, 100);
}
