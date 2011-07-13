// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2006-2008 Benoit Jacob <jacob.benoit.1@gmail.com>
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

#include "product.h"

void test_product_large()
{
  for(int i = 0; i < g_repeat; i++) {
    CALL_SUBTEST_1( product(MatrixXf(internal::random<int>(1,320), internal::random<int>(1,320))) );
    CALL_SUBTEST_2( product(MatrixXd(internal::random<int>(1,320), internal::random<int>(1,320))) );
    CALL_SUBTEST_3( product(MatrixXi(internal::random<int>(1,320), internal::random<int>(1,320))) );
    CALL_SUBTEST_4( product(MatrixXcf(internal::random<int>(1,150), internal::random<int>(1,150))) );
    CALL_SUBTEST_5( product(Matrix<float,Dynamic,Dynamic,RowMajor>(internal::random<int>(1,320), internal::random<int>(1,320))) );
  }

#if defined EIGEN_TEST_PART_6
  {
    // test a specific issue in DiagonalProduct
    int N = 1000000;
    VectorXf v = VectorXf::Ones(N);
    MatrixXf m = MatrixXf::Ones(N,3);
    m = (v+v).asDiagonal() * m;
    VERIFY_IS_APPROX(m, MatrixXf::Constant(N,3,2));
  }

  {
    // test deferred resizing in Matrix::operator=
    MatrixXf a = MatrixXf::Random(10,4), b = MatrixXf::Random(4,10), c = a;
    VERIFY_IS_APPROX((a = a * b), (c * b).eval());
  }

  {
    // check the functions to setup blocking sizes compile and do not segfault
    // FIXME check they do what they are supposed to do !!
    std::ptrdiff_t l1 = internal::random<int>(10000,20000);
    std::ptrdiff_t l2 = internal::random<int>(1000000,2000000);
    setCpuCacheSizes(l1,l2);
    VERIFY(l1==l1CacheSize());
    VERIFY(l2==l2CacheSize());
    std::ptrdiff_t k1 = internal::random<int>(10,100)*16;
    std::ptrdiff_t m1 = internal::random<int>(10,100)*16;
    std::ptrdiff_t n1 = internal::random<int>(10,100)*16;
    // only makes sure it compiles fine
    internal::computeProductBlockingSizes<float,float>(k1,m1,n1);
  }

  {
    // test regression in row-vector by matrix (bad Map type)
    MatrixXf mat1(10,32); mat1.setRandom();
    MatrixXf mat2(32,32); mat2.setRandom();
    MatrixXf r1 = mat1.row(2)*mat2.transpose();
    VERIFY_IS_APPROX(r1, (mat1.row(2)*mat2.transpose()).eval());

    MatrixXf r2 = mat1.row(2)*mat2;
    VERIFY_IS_APPROX(r2, (mat1.row(2)*mat2).eval());
  }
#endif
}
