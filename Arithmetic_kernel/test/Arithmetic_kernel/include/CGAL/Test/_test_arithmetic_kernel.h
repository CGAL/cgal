// Copyright (c) 2006-2009 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL:$
// $Id:$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Michael Hemmer <hemmer@mpi-inf.mpg.de>
//
// ============================================================================
//
//    \brief provide test for Arithmetic_kernel
//

#ifndef CGAL_TEST_ARITHMETIC_KERNEL_H
#define CGAL_TEST_ARITHMETIC_KERNEL_H


#include <CGAL/Test/_test_coercion_traits.h>

namespace CGAL {

  typedef CGAL::Interval_nt<false> Interval;

  void test_coercion_from_to(CGAL::Null_tag, CGAL::Null_tag){}
  template<class A> void test_coercion_from_to(A, CGAL::Null_tag){}
  template<class B> void test_coercion_from_to(CGAL::Null_tag, B){}
  template<class A, class B> void test_coercion_from_to(A, B){
    CGAL::test_explicit_interoperable_from_to<A,B>();
  }

template <class ARK>
void test_coercion_arithmetic_kernel(){


  typedef typename ARK::Integer Integer;
  typedef typename ARK::Rational Rational;
  typedef typename ARK::Field_with_sqrt Field_with_sqrt;
  typedef typename ARK::Field_with_kth_root Field_with_kth_root;
  typedef typename ARK::Field_with_root_of  Field_with_root_of;
  typedef typename ARK::Bigfloat Bigfloat;
  typedef typename ARK::Bigfloat_interval Bigfloat_interval;


  test_coercion_from_to(int(),Integer());
  test_coercion_from_to(short(),Integer());
  test_coercion_from_to(Integer(),Integer());

  test_coercion_from_to(int(),Rational());
  test_coercion_from_to(short(),Rational());
  test_coercion_from_to(float(),Rational());
  test_coercion_from_to(double(),Rational());
  test_coercion_from_to(Integer(),Rational());
  // This is currently not consistent with LEDA and GMP, for CORE it is not even defined.
  // test_coercion_from_to(Bigfloat(),Rational());
  test_coercion_from_to(Rational(),Rational());

  test_coercion_from_to(int(),Field_with_sqrt());
  test_coercion_from_to(short(),Field_with_sqrt());
  test_coercion_from_to(float(),Field_with_sqrt());
  test_coercion_from_to(double(),Field_with_sqrt());
  test_coercion_from_to(Integer(),Field_with_sqrt());
  test_coercion_from_to(Bigfloat(),Field_with_sqrt());
  test_coercion_from_to(Rational(),Field_with_sqrt());
  test_coercion_from_to(Field_with_sqrt(),Field_with_sqrt());

  test_coercion_from_to(int(),Field_with_kth_root());
  test_coercion_from_to(short(),Field_with_kth_root());
  test_coercion_from_to(float(),Field_with_kth_root());
  test_coercion_from_to(double(),Field_with_kth_root());
  test_coercion_from_to(Integer(),Field_with_kth_root());
  test_coercion_from_to(Bigfloat(),Field_with_kth_root());
  test_coercion_from_to(Rational(),Field_with_kth_root());
  test_coercion_from_to(Field_with_sqrt(),Field_with_kth_root());
  test_coercion_from_to(Field_with_kth_root(),Field_with_kth_root());


  test_coercion_from_to(int(),Field_with_root_of());
  test_coercion_from_to(short(),Field_with_root_of());
  test_coercion_from_to(float(),Field_with_root_of());
  test_coercion_from_to(double(),Field_with_root_of());
  test_coercion_from_to(Integer(),Field_with_root_of());
  test_coercion_from_to(Bigfloat(),Field_with_root_of());
  test_coercion_from_to(Rational(),Field_with_root_of());
  test_coercion_from_to(Field_with_sqrt(),Field_with_root_of());
  test_coercion_from_to(Field_with_kth_root(),Field_with_root_of());
  test_coercion_from_to(Field_with_root_of(),Field_with_root_of());


  test_coercion_from_to(int(),Bigfloat());
  test_coercion_from_to(short(),Bigfloat());
  test_coercion_from_to(float(),Bigfloat());
  test_coercion_from_to(double(),Bigfloat());
  test_coercion_from_to(Integer(),Bigfloat());
  test_coercion_from_to(Bigfloat(),Bigfloat());


  test_coercion_from_to(int(),Bigfloat_interval());
  test_coercion_from_to(short(),Bigfloat_interval());
  test_coercion_from_to(float(),Bigfloat_interval());
  test_coercion_from_to(double(),Bigfloat_interval());
  test_coercion_from_to(Integer(),Bigfloat_interval());
  test_coercion_from_to(Bigfloat(),Bigfloat_interval());
  test_coercion_from_to(Rational(),Bigfloat_interval());
  test_coercion_from_to(Field_with_sqrt(),Bigfloat_interval());
  test_coercion_from_to(Field_with_kth_root(),Bigfloat_interval());
  test_coercion_from_to(Field_with_root_of(),Bigfloat_interval());
  test_coercion_from_to(Bigfloat_interval(),Bigfloat_interval());
}

template <typename ARK>
void test_arithmetic_kernel(){
  test_coercion_arithmetic_kernel<ARK>();
}
} //namespace CGAL
#endif // CGAL_TEST_ARITHMETIC_KERNEL_H
