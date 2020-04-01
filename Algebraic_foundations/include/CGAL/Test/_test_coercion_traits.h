// Copyright (c) 2006-2008 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael Hemmer    <hemmer@mpi-inf.mpg.de>
//
// =============================================================================


/*! \file test_coercion_traits.h
    \brief test functions for class NiX::Coercion_traits
*/

#include <cassert>
#include <CGAL/tags.h>
#include <CGAL/use.h>
#include <CGAL/number_type_config.h>
#include <CGAL/number_utils.h>
#include <CGAL/Algebraic_structure_traits.h>

// These are test functions for the Coercion_traits
namespace CGAL {

// this test implicit interoperable
template< class A, class B, class Type > void test_implicit_interoperable();
template< class FROM, class TO > void test_implicit_interoperable_from_to();

// this is testing explicit interoperability only
template <class A, class B, class Type> void test_explicit_interoperable();
template  <class FROM, class TO> void test_explicit_interoperable_from_to();

namespace INTERN_COERCION_TRAITS {

template <typename A, typename B>
void test_implicit_interoperable_for_real_embeddable (CGAL::Tag_false){}

template <typename A, typename B>
void test_implicit_interoperable_for_real_embeddable (CGAL::Tag_true){
  // two sided test for interoperability with int
  A a;
  B b;

  volatile int value_a = -5;
  volatile int value_b = -2;
  // MSVC optimizer (at least VC9 and VC10) has problems with the following
  // code with /O2 and /fp:strict (it does constant propagation but
  // produces erroneous assembler code). Using volatile variables prevents
  // the constant propagation.
  /*
    int main(){
      int i = 3;
      float f = 3.f;
      bool b = (f>= i);
      return b ? 0 : 1;
    }
  */

  a = A(value_a);
  b = B(value_b);
  // a < b
  assert (!(a == b));
  assert ( (a != b));
  assert ( (a <  b));
  assert ( (a <= b));
  assert (!(a >  b));
  assert (!(a >= b));

  assert (!(b == a));
  assert ( (b != a));
  assert (!(b <  a));
  assert (!(b <= a));
  assert ( (b >  a));
  assert ( (b >= a));

  value_a = 5;
  value_b = 2;
  // a > b
  a = A(value_a);
  b = B(value_b);
  assert (!(a == b));
  assert ( (a != b));
  assert (!(a <  b));
  assert (!(a <= b));
  assert ( (a >  b));
  assert ( (a >= b));

  assert (!(b == a));
  assert ( (b != a));
  assert ( (b <  a));
  assert ( (b <= a));
  assert (!(b >  a));
  assert (!(b >= a));

  // a == b
  value_a = 3;
  value_b = 3;
  a = A(value_a);
  b = B(value_b);
  assert ( (a == b));
  assert (!(a != b));
  assert (!(a <  b));
  assert ( (a <= b));
  assert (!(a >  b));
  assert ( (a >= b));

  assert ( (b == a));
  assert (!(b != a));
  assert (!(b <  a));
  assert ( (b <= a));
  assert (!(b >  a));
  assert ( (b >= a));
}


template <typename A, typename B>
void test_implicit_interoperable_for_algebraic_structure
(CGAL::Null_tag){}

template <typename A, typename B>
void test_implicit_interoperable_for_algebraic_structure
(CGAL::Integral_domain_without_division_tag){
  typedef CGAL::Coercion_traits<A,B> CT;
  typedef typename CT::Type C;
  A a(6); B b(2);
  assert(a + b == C(8));
  assert(a - b == C(4));
  assert(a * b == C(12));

  assert(b + a == C(8));
  assert(b - a == C(-4));
  assert(b * a == C(12));

  C c;
  c = C(4); assert((c+= A(3)) == C(7));
  c = C(4); assert((c-= A(3)) == C(1));
  c = C(4); assert((c*= A(3)) == C(12));

  c = C(4); assert((c+= B(3)) == C(7));
  c = C(4); assert((c-= B(3)) == C(1));
  c = C(4); assert((c*= B(3)) == C(12));
}

template <typename A, typename B>
void test_implicit_interoperable_for_algebraic_structure
(CGAL::Field_tag){
  test_implicit_interoperable_for_algebraic_structure<A,B>
    (CGAL::Integral_domain_without_division_tag());

  typedef CGAL::Coercion_traits<A,B> CT;
  typedef typename CT::Type C;
  A a(6); B b(2);
  C aa = C(6);
  C bb = C(2);
  assert(a / b == C(3));
  assert(b / a == bb/aa);
  C c;
  c = C(4); assert((c /= A(2)) == C(2));
  c = C(4); assert((c /= B(2)) == C(2));
}

template< class A, class B, class Type, class Compare >
class Test_compare {
  public:
    void operator()() {
      Compare compare;
      A a(4);
      B b(2);
      typename CGAL::Coercion_traits< A, B >::Cast cast;
      Type a_ret = cast(a);
      Type b_ret = cast(b);
      assert( compare( a, b ) == CGAL_NTS compare( a_ret, b_ret ) );
    }
};

template< class A, class B, class Type, class Integral_division >
class Test_integral_division {
  public:
    void operator()() {
      Integral_division integral_division;
      A a(4);
      B b(2);
      typename CGAL::Coercion_traits< A, B >::Cast cast;
      Type a_ret = cast(a);
      Type b_ret = cast(b);
      assert( integral_division( a, b ) ==
                        CGAL_NTS integral_division( a_ret, b_ret ) );
    }
};

template< class A, class B, class Type >
class Test_integral_division< A, B, Type, CGAL::Null_functor > {
  public:
  // Nothing to test
    void operator()(){}
};

template< class A, class B, class Type, class Gcd >
class Test_gcd {
  public:
    void operator()() {
      Gcd gcd;
      A a(4);
      B b(2);
      typename CGAL::Coercion_traits< A, B >::Cast cast;
      Type a_ret = cast(a);
      Type b_ret = cast(b);
      assert( gcd( a, b ) ==
                        CGAL_NTS gcd( a_ret, b_ret ) );
    }
};

template< class A, class B, class Type >
class Test_gcd< A, B, Type, CGAL::Null_functor > {
  public:
  // Nothing to test
    void operator()(){}
};

template< class A, class B, class Type, class Div_mod >
class Test_div_mod {
  public:
    void operator()() {
      Div_mod div_mod;
      A a(4);
      B b(2);
      typename CGAL::Coercion_traits< A, B >::Cast cast;
      Type a_ret = cast(a);
      Type b_ret = cast(b);
      Type q;
      Type r;
      Type q_to_compare;
      Type r_to_compare;
      div_mod( a, b, q, r );
      CGAL_NTS div_mod( a_ret, b_ret, q_to_compare, r_to_compare );
      assert( q == q_to_compare );
      assert( r == r_to_compare );
    }
};

template< class A, class B, class Type >
class Test_div_mod< A, B, Type, CGAL::Null_functor > {
  public:
  // Nothing to test
    void operator()(){}
};

template< class A, class B, class Type, class Div >
class Test_div {
  public:
    void operator()() {
      Div div;
      A a(4);
      B b(2);
      typename CGAL::Coercion_traits< A, B >::Cast cast;
      Type a_ret = cast(a);
      Type b_ret = cast(b);
      assert( div( a, b ) ==
                        CGAL_NTS div( a_ret, b_ret ) );
    }
};

template< class A, class B, class Type >
class Test_div< A, B, Type, CGAL::Null_functor > {
  public:
  // Nothing to test
    void operator()(){}
};

template< class A, class B, class Type, class Mod >
class Test_mod {
  public:
    void operator()() {
      Mod mod;
      A a(4);
      B b(2);
      typename CGAL::Coercion_traits< A, B >::Cast cast;
      Type a_ret = cast(a);
      Type b_ret = cast(b);
      assert( mod( a, b ) ==
                        CGAL_NTS mod( a_ret, b_ret ) );
    }
};

template< class A, class B, class Type >
class Test_mod< A, B, Type, CGAL::Null_functor > {
  public:
  // Nothing to test
    void operator()(){}
};


template< class Type > void test_implicit_construction(Type) {}

template< class A, class B, class Type, class Are_implicit_interoperable >
class Implicit_test_implicit_interoperable {
  public:
    void operator()() {
      // enforce implicit construction from A/B to Type
      // Results in 'no matching function for call to...' compile error, if type Type
      // is not implicit constructable from A and B (which is part of the concept)
      test_implicit_construction<Type>(A(1));
      test_implicit_construction<Type>(B(2));

      // test explicit construction
      Type test_var = Type(A(1));
      test_var = Type(B(2));
    }
};

template< class A, class B, class Type >
class Implicit_test_implicit_interoperable<A,B,Type, CGAL::Tag_false > {
  public:
    void operator()() {}
};

template< class A, class B, class Type >
void test_implicit_interoperable_one_way() {

  typedef CGAL::Coercion_traits<A,B> CT;
  typedef typename CT::Type C;
  typedef typename CT::Are_implicit_interoperable Are_implicit_interoperable;

  CGAL_static_assertion(
      (::boost::is_same<Are_implicit_interoperable, CGAL::Tag_true>::value));
  assert((::boost::is_same<Are_implicit_interoperable, CGAL::Tag_true>::value));

  typename CGAL::Real_embeddable_traits<C>::Is_real_embeddable is_real_embeddable;
  test_implicit_interoperable_for_real_embeddable<A,B>(is_real_embeddable);
  typename CGAL::Algebraic_structure_traits<C>::Algebraic_category category;
  test_implicit_interoperable_for_algebraic_structure<A,B>(category);
}



// test for explicit interoperable types
template <class A, class B, class RT>
void test_explicit_interoperable_one_way(){
  typedef CGAL::Coercion_traits<A,B> CT;
  typedef typename CT::Type Type;
  typedef typename CT::Cast Cast;
  typedef typename Cast::result_type result_type;
  CGAL_USE_TYPE(result_type);
  CGAL_static_assertion((::boost::is_same<result_type,Type>::value));
  CGAL_static_assertion((::boost::is_same< typename CT::Are_explicit_interoperable,CGAL::Tag_true>::value));
  CGAL_static_assertion((::boost::is_same<Type,RT>::value));
  typename CT::Cast cast;

  A a(3);
  B b(3);
  RT  rt(3);
  assert(rt==cast(a));
  assert(rt==cast(b));

  // Binary Functors should support explicit interoperable types
  typedef typename CGAL::Algebraic_structure_traits<Type>::Integral_division Idiv;
  Test_integral_division< A, B, Type, Idiv>()();
  typedef typename CGAL::Algebraic_structure_traits<Type>::Gcd Gcd;
  Test_gcd< A, B, Type, Gcd >()();
  typedef typename CGAL::Algebraic_structure_traits<Type>::Div_mod Div_mod;
  Test_div_mod< A, B, Type, Div_mod>()();
  typedef  typename CGAL::Algebraic_structure_traits<Type>::Div Div;
  Test_div< A, B, Type, Div >()();
  typedef typename CGAL::Algebraic_structure_traits<Type>::Mod Mod;
  Test_mod< A, B, Type, Mod >()();
  typedef typename CGAL::Real_embeddable_traits<Type>::Compare Compare;
  Test_compare< A, B, Type, Compare >()();
}

}// namespace INTERN_COERCION_TRAITS

// this test implicit interoperable
template< class A, class B, class Type >
void test_implicit_interoperable() {
  INTERN_COERCION_TRAITS::test_implicit_interoperable_one_way< A, B, Type >();
  INTERN_COERCION_TRAITS::test_implicit_interoperable_one_way< B, A, Type >();
  INTERN_COERCION_TRAITS::test_explicit_interoperable_one_way<A,B,Type>();
  INTERN_COERCION_TRAITS::test_explicit_interoperable_one_way<B,A,Type>();
}

template< class FROM, class TO >
void test_implicit_interoperable_from_to() {
  test_implicit_interoperable< FROM, TO, TO >();
}

// this is testing explicit interoperability only
template <class A, class B, class Type>
void test_explicit_interoperable(){
    INTERN_COERCION_TRAITS::test_explicit_interoperable_one_way<A,B,Type>();
    INTERN_COERCION_TRAITS::test_explicit_interoperable_one_way<B,A,Type>();
}

template  <class FROM, class TO>
void test_explicit_interoperable_from_to(){
    test_explicit_interoperable<FROM,TO,TO>();
}
} //namespace CGAL
