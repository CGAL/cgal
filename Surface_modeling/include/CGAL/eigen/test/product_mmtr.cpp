// This file is part of Eigen, a lightweight C++ template library
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

#define CHECK_MMTR(DEST, TRI, OP) {                   \
    ref2 = ref1 = DEST;                               \
    DEST.template triangularView<TRI>() OP;           \
    ref1 OP;                                          \
    ref2.template triangularView<TRI>() = ref1;       \
    VERIFY_IS_APPROX(DEST,ref2);                      \
  }

template<typename Scalar> void mmtr(int size)
{
  typedef typename NumTraits<Scalar>::Real RealScalar;

  typedef Matrix<Scalar,Dynamic,Dynamic,ColMajor> MatrixColMaj;
  typedef Matrix<Scalar,Dynamic,Dynamic,RowMajor> MatrixRowMaj;

  DenseIndex othersize = internal::random<DenseIndex>(1,200);
  
  MatrixColMaj matc(size, size);
  MatrixRowMaj matr(size, size);
  MatrixColMaj ref1(size, size), ref2(size, size);
  
  MatrixColMaj soc(size,othersize); soc.setRandom();
  MatrixColMaj osc(othersize,size); osc.setRandom();
  MatrixRowMaj sor(size,othersize); sor.setRandom();
  MatrixRowMaj osr(othersize,size); osr.setRandom();
  
  Scalar s = internal::random<Scalar>();
  
  CHECK_MMTR(matc, Lower, = s*soc*sor.adjoint());
  CHECK_MMTR(matc, Upper, = s*(soc*soc.adjoint()));
  CHECK_MMTR(matr, Lower, = s*soc*soc.adjoint());
  CHECK_MMTR(matr, Upper, = soc*(s*sor.adjoint()));
  
  CHECK_MMTR(matc, Lower, += s*soc*soc.adjoint());
  CHECK_MMTR(matc, Upper, += s*(soc*sor.transpose()));
  CHECK_MMTR(matr, Lower, += s*sor*soc.adjoint());
  CHECK_MMTR(matr, Upper, += soc*(s*soc.adjoint()));
  
  CHECK_MMTR(matc, Lower, -= s*soc*soc.adjoint());
  CHECK_MMTR(matc, Upper, -= s*(osc.transpose()*osc.conjugate()));
  CHECK_MMTR(matr, Lower, -= s*soc*soc.adjoint());
  CHECK_MMTR(matr, Upper, -= soc*(s*soc.adjoint()));
}

void test_product_mmtr()
{
  for(int i = 0; i < g_repeat ; i++)
  {
    CALL_SUBTEST_1((mmtr<float>(internal::random<int>(1,320))));
    CALL_SUBTEST_2((mmtr<double>(internal::random<int>(1,320))));
    CALL_SUBTEST_3((mmtr<std::complex<float> >(internal::random<int>(1,200))));
    CALL_SUBTEST_4((mmtr<std::complex<double> >(internal::random<int>(1,200))));
  }
}
