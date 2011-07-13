// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2008 Gael Guennebaud <gael.guennebaud@inria.fr>
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

#if EIGEN_ALIGN
#define ALIGNMENT 16
#else
#define ALIGNMENT 1
#endif

void check_handmade_aligned_malloc()
{
  for(int i = 1; i < 1000; i++)
  {
    char *p = (char*)internal::handmade_aligned_malloc(i);
    VERIFY(size_t(p)%ALIGNMENT==0);
    // if the buffer is wrongly allocated this will give a bad write --> check with valgrind
    for(int j = 0; j < i; j++) p[j]=0;
    internal::handmade_aligned_free(p);
  }
}

void check_aligned_malloc()
{
  for(int i = 1; i < 1000; i++)
  {
    char *p = (char*)internal::aligned_malloc(i);
    VERIFY(size_t(p)%ALIGNMENT==0);
    // if the buffer is wrongly allocated this will give a bad write --> check with valgrind
    for(int j = 0; j < i; j++) p[j]=0;
    internal::aligned_free(p);
  }
}

void check_aligned_new()
{
  for(int i = 1; i < 1000; i++)
  {
    float *p = internal::aligned_new<float>(i);
    VERIFY(size_t(p)%ALIGNMENT==0);
    // if the buffer is wrongly allocated this will give a bad write --> check with valgrind
    for(int j = 0; j < i; j++) p[j]=0;
    internal::aligned_delete(p,i);
  }
}

void check_aligned_stack_alloc()
{
  for(int i = 1; i < 1000; i++)
  {
    float *p = ei_aligned_stack_new(float,i);
    VERIFY(size_t(p)%ALIGNMENT==0);
    // if the buffer is wrongly allocated this will give a bad write --> check with valgrind
    for(int j = 0; j < i; j++) p[j]=0;
    ei_aligned_stack_delete(float,p,i);
  }
}


// test compilation with both a struct and a class...
struct MyStruct
{
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  char dummychar;
  Vector4f avec;
};

class MyClassA
{
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    char dummychar;
    Vector4f avec;
};

template<typename T> void check_dynaligned()
{
  T* obj = new T;
  VERIFY(size_t(obj)%ALIGNMENT==0);
  delete obj;
}

void test_dynalloc()
{
  // low level dynamic memory allocation
  CALL_SUBTEST(check_handmade_aligned_malloc());
  CALL_SUBTEST(check_aligned_malloc());
  CALL_SUBTEST(check_aligned_new());
  CALL_SUBTEST(check_aligned_stack_alloc());

  for (int i=0; i<g_repeat*100; ++i)
  {
    CALL_SUBTEST(check_dynaligned<Vector4f>() );
    CALL_SUBTEST(check_dynaligned<Vector2d>() );
    CALL_SUBTEST(check_dynaligned<Matrix4f>() );
    CALL_SUBTEST(check_dynaligned<Vector4d>() );
    CALL_SUBTEST(check_dynaligned<Vector4i>() );
  }
  
  // check static allocation, who knows ?
  #if EIGEN_ALIGN_STATICALLY
  {
    MyStruct foo0;  VERIFY(size_t(foo0.avec.data())%ALIGNMENT==0);
    MyClassA fooA;  VERIFY(size_t(fooA.avec.data())%ALIGNMENT==0);
  }
  
  // dynamic allocation, single object
  for (int i=0; i<g_repeat*100; ++i)
  {
    MyStruct *foo0 = new MyStruct();  VERIFY(size_t(foo0->avec.data())%ALIGNMENT==0);
    MyClassA *fooA = new MyClassA();  VERIFY(size_t(fooA->avec.data())%ALIGNMENT==0);
    delete foo0;
    delete fooA;
  }

  // dynamic allocation, array
  const int N = 10;
  for (int i=0; i<g_repeat*100; ++i)
  {
    MyStruct *foo0 = new MyStruct[N];  VERIFY(size_t(foo0->avec.data())%ALIGNMENT==0);
    MyClassA *fooA = new MyClassA[N];  VERIFY(size_t(fooA->avec.data())%ALIGNMENT==0);
    delete[] foo0;
    delete[] fooA;
  }
  #endif
  
}
