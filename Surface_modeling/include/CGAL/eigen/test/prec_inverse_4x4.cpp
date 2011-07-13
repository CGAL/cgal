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
#include <Eigen/LU>
#include <algorithm>

template<typename MatrixType> void inverse_permutation_4x4()
{
  typedef typename MatrixType::Scalar Scalar;
  typedef typename MatrixType::RealScalar RealScalar;
  Vector4i indices(0,1,2,3);
  for(int i = 0; i < 24; ++i)
  {
    MatrixType m = PermutationMatrix<4>(indices);
    MatrixType inv = m.inverse();
    double error = double( (m*inv-MatrixType::Identity()).norm() / NumTraits<Scalar>::epsilon() );
    EIGEN_DEBUG_VAR(error)
    VERIFY(error == 0.0);
    std::next_permutation(indices.data(),indices.data()+4);
  }
}

template<typename MatrixType> void inverse_general_4x4(int repeat)
{
  typedef typename MatrixType::Scalar Scalar;
  typedef typename MatrixType::RealScalar RealScalar;
  double error_sum = 0., error_max = 0.;
  for(int i = 0; i < repeat; ++i)
  {
    MatrixType m;
    RealScalar absdet;
    do {
      m = MatrixType::Random();
      absdet = internal::abs(m.determinant());
    } while(absdet < NumTraits<Scalar>::epsilon());
    MatrixType inv = m.inverse();
    double error = double( (m*inv-MatrixType::Identity()).norm() * absdet / NumTraits<Scalar>::epsilon() );
    error_sum += error;
    error_max = std::max(error_max, error);
  }
  std::cerr << "inverse_general_4x4, Scalar = " << type_name<Scalar>() << std::endl;
  double error_avg = error_sum / repeat;
  EIGEN_DEBUG_VAR(error_avg);
  EIGEN_DEBUG_VAR(error_max);
   // FIXME that 1.25 used to be a 1.0 until the NumTraits changes on 28 April 2010, what's going wrong??
   // FIXME that 1.25 used to be 1.2 until we tested gcc 4.1 on 30 June 2010 and got 1.21.
  VERIFY(error_avg < (NumTraits<Scalar>::IsComplex ? 8.0 : 1.25));
  VERIFY(error_max < (NumTraits<Scalar>::IsComplex ? 64.0 : 20.0));
}

void test_prec_inverse_4x4()
{
  CALL_SUBTEST_1((inverse_permutation_4x4<Matrix4f>()));
  CALL_SUBTEST_1(( inverse_general_4x4<Matrix4f>(200000 * g_repeat) ));

  CALL_SUBTEST_2((inverse_permutation_4x4<Matrix<double,4,4,RowMajor> >()));
  CALL_SUBTEST_2(( inverse_general_4x4<Matrix<double,4,4,RowMajor> >(200000 * g_repeat) ));

  CALL_SUBTEST_3((inverse_permutation_4x4<Matrix4cf>()));
  CALL_SUBTEST_3((inverse_general_4x4<Matrix4cf>(50000 * g_repeat)));
}
