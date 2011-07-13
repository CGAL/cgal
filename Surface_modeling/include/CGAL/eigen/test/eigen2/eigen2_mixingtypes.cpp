// This file is part of Eigen, a lightweight C++ template library
// for linear algebra. Eigen itself is part of the KDE project.
//
// Copyright (C) 2008 Gael Guennebaud <g.gael@free.fr>
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

#ifndef EIGEN_NO_STATIC_ASSERT
#define EIGEN_NO_STATIC_ASSERT // turn static asserts into runtime asserts in order to check them
#endif

#ifndef EIGEN_DONT_VECTORIZE
#define EIGEN_DONT_VECTORIZE // SSE intrinsics aren't designed to allow mixing types
#endif

#include "main.h"


template<int SizeAtCompileType> void mixingtypes(int size = SizeAtCompileType)
{
  typedef Matrix<float, SizeAtCompileType, SizeAtCompileType> Mat_f;
  typedef Matrix<double, SizeAtCompileType, SizeAtCompileType> Mat_d;
  typedef Matrix<std::complex<float>, SizeAtCompileType, SizeAtCompileType> Mat_cf;
  typedef Matrix<std::complex<double>, SizeAtCompileType, SizeAtCompileType> Mat_cd;
  typedef Matrix<float, SizeAtCompileType, 1> Vec_f;
  typedef Matrix<double, SizeAtCompileType, 1> Vec_d;
  typedef Matrix<std::complex<float>, SizeAtCompileType, 1> Vec_cf;
  typedef Matrix<std::complex<double>, SizeAtCompileType, 1> Vec_cd;

  Mat_f mf(size,size);
  Mat_d md(size,size);
  Mat_cf mcf(size,size);
  Mat_cd mcd(size,size);
  Vec_f vf(size,1);
  Vec_d vd(size,1);
  Vec_cf vcf(size,1);
  Vec_cd vcd(size,1);

  mf+mf;
  VERIFY_RAISES_ASSERT(mf+md);
  VERIFY_RAISES_ASSERT(mf+mcf);
  VERIFY_RAISES_ASSERT(vf=vd);
  VERIFY_RAISES_ASSERT(vf+=vd);
  VERIFY_RAISES_ASSERT(mcd=md);

  mf*mf;
  md*mcd;
  mcd*md;
  mf*vcf;
  mcf*vf;
  mcf *= mf;
  vcd = md*vcd;
  vcf = mcf*vf;
#if 0
  // these are know generating hard build errors in eigen3
  VERIFY_RAISES_ASSERT(mf*md);
  VERIFY_RAISES_ASSERT(mcf*mcd);
  VERIFY_RAISES_ASSERT(mcf*vcd);
  VERIFY_RAISES_ASSERT(vcf = mf*vf);

  vf.eigen2_dot(vf);
  VERIFY_RAISES_ASSERT(vd.eigen2_dot(vf));
  VERIFY_RAISES_ASSERT(vcf.eigen2_dot(vf)); // yeah eventually we should allow this but i'm too lazy to make that change now in Dot.h
  // especially as that might be rewritten as cwise product .sum() which would make that automatic.
#endif
}

void test_eigen2_mixingtypes()
{
  // check that our operator new is indeed called:
  CALL_SUBTEST_1(mixingtypes<3>());
  CALL_SUBTEST_2(mixingtypes<4>());
  CALL_SUBTEST_3(mixingtypes<Dynamic>(20));
}
