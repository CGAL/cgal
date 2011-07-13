// This file is part of Eigen, a lightweight C++ template library
// for linear algebra. Eigen itself is part of the KDE project.
//
// Copyright (C) 2009 Hauke Heibel <hauke.heibel@gmail.com>
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
// WARRANTY; without even the implied warranty of MERCHANTABILITY or1 FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License or the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License and a copy of the GNU General Public License along with
// Eigen. If not, see <http://www.gnu.org/licenses/>.

#include "main.h"

#include <Eigen/Core>

using namespace Eigen;

template <typename Scalar, int Storage>
void run_matrix_tests()
{
  typedef Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Storage> MatrixType;
  typedef typename MatrixType::Index Index;

  MatrixType m, n;

  // boundary cases ...
  m = n = MatrixType::Random(50,50);
  m.conservativeResize(1,50);
  VERIFY_IS_APPROX(m, n.block(0,0,1,50));

  m = n = MatrixType::Random(50,50);
  m.conservativeResize(50,1);
  VERIFY_IS_APPROX(m, n.block(0,0,50,1));

  m = n = MatrixType::Random(50,50);
  m.conservativeResize(50,50);
  VERIFY_IS_APPROX(m, n.block(0,0,50,50));

  // random shrinking ...
  for (int i=0; i<25; ++i)
  {
    const Index rows = internal::random<Index>(1,50);
    const Index cols = internal::random<Index>(1,50);
    m = n = MatrixType::Random(50,50);
    m.conservativeResize(rows,cols);
    VERIFY_IS_APPROX(m, n.block(0,0,rows,cols));
  }

  // random growing with zeroing ...
  for (int i=0; i<25; ++i)
  {
    const Index rows = internal::random<Index>(50,75);
    const Index cols = internal::random<Index>(50,75);
    m = n = MatrixType::Random(50,50);
    m.conservativeResizeLike(MatrixType::Zero(rows,cols));
    VERIFY_IS_APPROX(m.block(0,0,n.rows(),n.cols()), n);
    VERIFY( rows<=50 || m.block(50,0,rows-50,cols).sum() == Scalar(0) );
    VERIFY( cols<=50 || m.block(0,50,rows,cols-50).sum() == Scalar(0) );
  }
}

template <typename Scalar>
void run_vector_tests()
{
  typedef Matrix<Scalar, 1, Eigen::Dynamic> MatrixType;

  MatrixType m, n;

  // boundary cases ...
  m = n = MatrixType::Random(50);
  m.conservativeResize(1);
  VERIFY_IS_APPROX(m, n.segment(0,1));

  m = n = MatrixType::Random(50);
  m.conservativeResize(50);
  VERIFY_IS_APPROX(m, n.segment(0,50));

  // random shrinking ...
  for (int i=0; i<50; ++i)
  {
    const int size = internal::random<int>(1,50);
    m = n = MatrixType::Random(50);
    m.conservativeResize(size);
    VERIFY_IS_APPROX(m, n.segment(0,size));
  }

  // random growing with zeroing ...
  for (int i=0; i<50; ++i)
  {
    const int size = internal::random<int>(50,100);
    m = n = MatrixType::Random(50);
    m.conservativeResizeLike(MatrixType::Zero(size));
    VERIFY_IS_APPROX(m.segment(0,50), n);
    VERIFY( size<=50 || m.segment(50,size-50).sum() == Scalar(0) );
  }
}

void test_conservative_resize()
{
  CALL_SUBTEST_1((run_matrix_tests<int, Eigen::RowMajor>()));
  CALL_SUBTEST_1((run_matrix_tests<int, Eigen::ColMajor>()));
  CALL_SUBTEST_2((run_matrix_tests<float, Eigen::RowMajor>()));
  CALL_SUBTEST_2((run_matrix_tests<float, Eigen::ColMajor>()));
  CALL_SUBTEST_3((run_matrix_tests<double, Eigen::RowMajor>()));
  CALL_SUBTEST_3((run_matrix_tests<double, Eigen::ColMajor>()));
  CALL_SUBTEST_4((run_matrix_tests<std::complex<float>, Eigen::RowMajor>()));
  CALL_SUBTEST_4((run_matrix_tests<std::complex<float>, Eigen::ColMajor>()));
  CALL_SUBTEST_5((run_matrix_tests<std::complex<double>, Eigen::RowMajor>()));
  CALL_SUBTEST_6((run_matrix_tests<std::complex<double>, Eigen::ColMajor>()));

  CALL_SUBTEST_1((run_vector_tests<int>()));
  CALL_SUBTEST_2((run_vector_tests<float>()));
  CALL_SUBTEST_3((run_vector_tests<double>()));
  CALL_SUBTEST_4((run_vector_tests<std::complex<float> >()));
  CALL_SUBTEST_5((run_vector_tests<std::complex<double> >()));
}
