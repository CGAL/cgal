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

template<typename Scalar> void trmm(int size,int /*othersize*/)
{
  typedef typename NumTraits<Scalar>::Real RealScalar;

  typedef Matrix<Scalar,Dynamic,Dynamic,ColMajor> MatrixColMaj;
  typedef Matrix<Scalar,Dynamic,Dynamic,RowMajor> MatrixRowMaj;

  DenseIndex rows = size;
  DenseIndex cols = internal::random<DenseIndex>(1,size);

  MatrixColMaj  triV(rows,cols), triH(cols,rows), upTri(cols,rows), loTri(rows,cols),
                unitUpTri(cols,rows), unitLoTri(rows,cols), strictlyUpTri(cols,rows), strictlyLoTri(rows,cols);
  MatrixColMaj  ge1(rows,cols), ge2(cols,rows), ge3;
  MatrixRowMaj  rge3;

  Scalar s1 = internal::random<Scalar>(),
         s2 = internal::random<Scalar>();

  triV.setRandom();
  triH.setRandom();
  loTri = triV.template triangularView<Lower>();
  upTri = triH.template triangularView<Upper>();
  unitLoTri = triV.template triangularView<UnitLower>();
  unitUpTri = triH.template triangularView<UnitUpper>();
  strictlyLoTri = triV.template triangularView<StrictlyLower>();
  strictlyUpTri = triH.template triangularView<StrictlyUpper>();
  ge1.setRandom();
  ge2.setRandom();

  VERIFY_IS_APPROX( ge3 = triV.template triangularView<Lower>() * ge2, loTri * ge2);
  VERIFY_IS_APPROX( ge3 = ge2 * triV.template triangularView<Lower>(), ge2 * loTri);
  VERIFY_IS_APPROX( ge3 = triH.template triangularView<Upper>() * ge1, upTri * ge1);
  VERIFY_IS_APPROX( ge3 = ge1 * triH.template triangularView<Upper>(), ge1 * upTri);
  VERIFY_IS_APPROX( ge3 = (s1*triV.adjoint()).template triangularView<Upper>() * (s2*ge1), s1*loTri.adjoint() * (s2*ge1));
  VERIFY_IS_APPROX( ge3 = ge1 * triV.adjoint().template triangularView<Upper>(), ge1 * loTri.adjoint());
  VERIFY_IS_APPROX( ge3 = triH.adjoint().template triangularView<Lower>() * ge2, upTri.adjoint() * ge2);
  VERIFY_IS_APPROX( ge3 = ge2 * triH.adjoint().template triangularView<Lower>(), ge2 * upTri.adjoint());
  VERIFY_IS_APPROX( ge3 = triV.template triangularView<Lower>() * ge1.adjoint(), loTri * ge1.adjoint());
  VERIFY_IS_APPROX( ge3 = ge1.adjoint() * triV.template triangularView<Lower>(), ge1.adjoint() * loTri);
  VERIFY_IS_APPROX( ge3 = triH.template triangularView<Upper>() * ge2.adjoint(), upTri * ge2.adjoint());
  VERIFY_IS_APPROX(rge3.noalias() = triH.template triangularView<Upper>() * ge2.adjoint(), upTri * ge2.adjoint());
  VERIFY_IS_APPROX( ge3 = (s1*triV).adjoint().template triangularView<Upper>() * ge2.adjoint(), internal::conj(s1) * loTri.adjoint() * ge2.adjoint());
  VERIFY_IS_APPROX(rge3.noalias() = triV.adjoint().template triangularView<Upper>() * ge2.adjoint(), loTri.adjoint() * ge2.adjoint());
  VERIFY_IS_APPROX( ge3 = triH.adjoint().template triangularView<Lower>() * ge1.adjoint(), upTri.adjoint() * ge1.adjoint());
  VERIFY_IS_APPROX(rge3.noalias() = triH.adjoint().template triangularView<Lower>() * ge1.adjoint(), upTri.adjoint() * ge1.adjoint());

  VERIFY_IS_APPROX( ge3 = triV.template triangularView<UnitLower>() * ge2, unitLoTri * ge2);
  VERIFY_IS_APPROX( rge3.noalias() = ge2 * triV.template triangularView<UnitLower>(), ge2 * unitLoTri);
  VERIFY_IS_APPROX( ge3 = ge2 * triV.template triangularView<UnitLower>(), ge2 * unitLoTri);
  VERIFY_IS_APPROX( ge3 = (s1*triV).adjoint().template triangularView<UnitUpper>() * ge2.adjoint(), internal::conj(s1) * unitLoTri.adjoint() * ge2.adjoint());

  VERIFY_IS_APPROX( ge3 = triV.template triangularView<StrictlyLower>() * ge2, strictlyLoTri * ge2);
  VERIFY_IS_APPROX( rge3.noalias() = ge2 * triV.template triangularView<StrictlyLower>(), ge2 * strictlyLoTri);
  VERIFY_IS_APPROX( ge3 = ge2 * triV.template triangularView<StrictlyLower>(), ge2 * strictlyLoTri);
  VERIFY_IS_APPROX( ge3 = (s1*triV).adjoint().template triangularView<StrictlyUpper>() * ge2.adjoint(), internal::conj(s1) * strictlyLoTri.adjoint() * ge2.adjoint());
}

void test_product_trmm()
{
  for(int i = 0; i < g_repeat ; i++)
  {
    CALL_SUBTEST_1((trmm<float>(internal::random<int>(1,320),internal::random<int>(1,320))));
    CALL_SUBTEST_2((trmm<double>(internal::random<int>(1,320),internal::random<int>(1,320))));
    CALL_SUBTEST_3((trmm<std::complex<float> >(internal::random<int>(1,200),internal::random<int>(1,200))));
    CALL_SUBTEST_4((trmm<std::complex<double> >(internal::random<int>(1,200),internal::random<int>(1,200))));
  }
}
