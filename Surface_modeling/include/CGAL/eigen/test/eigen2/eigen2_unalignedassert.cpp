// This file is part of Eigen, a lightweight C++ template library
// for linear algebra. Eigen itself is part of the KDE project.
//
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

#include "main.h"

struct Good1
{
  MatrixXd m; // good: m will allocate its own array, taking care of alignment.
  Good1() : m(20,20) {}
};

struct Good2
{
  Matrix3d m; // good: m's size isn't a multiple of 16 bytes, so m doesn't have to be aligned
};

struct Good3
{
  Vector2f m; // good: same reason
};

struct Bad4
{
  Vector2d m; // bad: sizeof(m)%16==0 so alignment is required
};

struct Bad5
{
  Matrix<float, 2, 6> m; // bad: same reason
};

struct Bad6
{
  Matrix<double, 3, 4> m; // bad: same reason
};

struct Good7
{
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  Vector2d m;
  float f; // make the struct have sizeof%16!=0 to make it a little more tricky when we allow an array of 2 such objects
};

struct Good8
{
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  float f; // try the f at first -- the EIGEN_ALIGN_128 attribute of m should make that still work
  Matrix4f m;
};

struct Good9
{
  Matrix<float,2,2,DontAlign> m; // good: no alignment requested
  float f;
};

template<bool Align> struct Depends
{
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF(Align)
  Vector2d m;
  float f;
};

template<typename T>
void check_unalignedassert_good()
{
  T *x, *y;
  x = new T;
  delete x;
  y = new T[2];
  delete[] y;
}

#if EIGEN_ARCH_WANTS_ALIGNMENT
template<typename T>
void check_unalignedassert_bad()
{
  float buf[sizeof(T)+16];
  float *unaligned = buf;
  while((reinterpret_cast<std::size_t>(unaligned)&0xf)==0) ++unaligned; // make sure unaligned is really unaligned
  T *x = ::new(static_cast<void*>(unaligned)) T;
  x->~T();
}
#endif

void unalignedassert()
{
  check_unalignedassert_good<Good1>();
  check_unalignedassert_good<Good2>();
  check_unalignedassert_good<Good3>();
#if EIGEN_ARCH_WANTS_ALIGNMENT
  VERIFY_RAISES_ASSERT(check_unalignedassert_bad<Bad4>());
  VERIFY_RAISES_ASSERT(check_unalignedassert_bad<Bad5>());
  VERIFY_RAISES_ASSERT(check_unalignedassert_bad<Bad6>());
#endif

  check_unalignedassert_good<Good7>();
  check_unalignedassert_good<Good8>();
  check_unalignedassert_good<Good9>();
  check_unalignedassert_good<Depends<true> >();

#if EIGEN_ARCH_WANTS_ALIGNMENT
  VERIFY_RAISES_ASSERT(check_unalignedassert_bad<Depends<false> >());
#endif
}

void test_eigen2_unalignedassert()
{
  CALL_SUBTEST(unalignedassert());
}
