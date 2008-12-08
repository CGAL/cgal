// Copyright (c) 2006-2008 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Michael Hemmer    <hemmer@mpi-inf.mpg.de>
//
// =============================================================================


/*! \file test_coercion_traits.h 
    \brief test functions for class NiX::Coercion_traits
*/

#include <cassert>

// These are test functions for the Coercion_traits
CGAL_BEGIN_NAMESPACE
    namespace INTERN_COERCION_TRAITS {


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


template< class Type >
void test_implicit_construction( Type a ) {
  (void)a;
}

template< class A, class B, class Type, class Are_implicit_interoperable >
class Implicit_interoperability_test {
  public:
    void operator()() {
      // test implicit construction
      // Results in 'no matching function for call to...' compile error, if type Type
      //  is not implicit constructable from A and B (which is part of the concept)
      test_implicit_construction<Type>(A(1));
      test_implicit_construction<Type>(B(2));
      
      // test explicit construction
      Type test_var = Type(A(1));
      test_var = Type(B(2));
    }
};

template< class A, class B, class Type >
class Implicit_interoperability_test<A,B,Type, CGAL::Tag_false > {
  public:
    void operator()() {}
};

template< class A, class B, class Type >
void interoperability_test_one_way() {
  typedef CGAL::Coercion_traits< A, B > CT;

  assert((::boost::is_same< typename CT::Are_implicit_interoperable,
                                      CGAL::Tag_true 
                                    >::value));                                    
  assert((::boost::is_same< typename CT::Type, Type >::value));
  
  // Implicit_interoperability_test
  Implicit_interoperability_test< A, B, Type, 
        typename CT::Are_implicit_interoperable >()();
  
  Test_integral_division< A, B, Type, 
        typename CGAL::Algebraic_structure_traits<Type>::Integral_division >()();
  
  Test_gcd< A, B, Type, 
        typename CGAL::Algebraic_structure_traits<Type>::Gcd >()();

  Test_div_mod< A, B, Type, 
          typename CGAL::Algebraic_structure_traits<Type>::Div_mod >()();

  Test_div< A, B, Type, 
          typename CGAL::Algebraic_structure_traits<Type>::Div >()();
  Test_mod< A, B, Type, 
          typename CGAL::Algebraic_structure_traits<Type>::Mod >()();
  Test_compare< A, B, Type,
          typename CGAL::Real_embeddable_traits<Type>::Compare >()();
}

template< class A, class B, class Type > 
void interoperability_test() {
  interoperability_test_one_way< A, B, Type >();
  interoperability_test_one_way< B, A, Type >();
}

template< class FROM, class TO >
void direct_interoperability_from_to_test() {
  interoperability_test< FROM, TO, TO >();
}

template <class A, class B, class RT>
void coercion_traits_test_one_way(){
    typedef CGAL::Coercion_traits<A,B> CT;
    {
        typedef typename CT::Type Type;
        typename CT::Cast cast;
        assert((::boost::is_same< 
                                      typename CT::Are_explicit_interoperable,
                                      CGAL::Tag_true 
                                          >::value));
        assert((::boost::is_same<Type,RT>::value));
        A a(3);
        B b(3);
        RT  rt(3);
        assert(rt==cast(a));
        assert(rt==cast(b)); 
    }
}

template <class A, class B, class Type>
void coercion_traits_test(){
    coercion_traits_test_one_way<A,B,Type>();
    coercion_traits_test_one_way<B,A,Type>();   
    
    // TODO: Is this here OK?
  if((::boost::is_same< typename CGAL::Coercion_traits< A, B >::Are_implicit_interoperable,
                                      CGAL::Tag_true 
                                    >::value)) {
    interoperability_test< A, B, Type >();
  }
  
}
template  <class FROM, class TO>
void direct_coercion_from_to_test(){
    coercion_traits_test<FROM,TO,TO>();
}

}// namespace INTERN_COERCION_TRAITS
CGAL_END_NAMESPACE
