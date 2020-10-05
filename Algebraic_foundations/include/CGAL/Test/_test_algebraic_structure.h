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
// Author(s)     : Lutz Kettner  <kettner@mpi-inf.mpg.de>
//                 Michael Hemmer <hemmer@mpi-inf.mpg.de>
//
// =============================================================================

/*   \brief provides test functions for the various ring concepts of number
 *          types.
*/

// within this file AS ^= Type

#include <CGAL/Algebraic_structure_traits.h>
//#include <CGAL/Real_embeddable_traits.h>

#include <CGAL/number_utils.h>
#include <CGAL/assertions.h>
#include <CGAL/use.h>
#include <boost/type_traits.hpp>
#include <CGAL/tags.h>
#include <cassert>
#include <functional>

#include <vector>

#include <CGAL/Testsuite/use.h>


#ifndef CGAL_TEST_ALGEBRAIC_STRUCTURE_H
#define CGAL_TEST_ALGEBRAIC_STRUCTURE_H

#include <CGAL/disable_warnings.h>

// checks the result type of a functor
template <typename AdaptableFunctor, typename ResultType>
void check_result_type(AdaptableFunctor, ResultType){
  typedef typename AdaptableFunctor::result_type result_type;
  CGAL_static_assertion((::boost::is_same<result_type,ResultType>::value));
  CGAL_USE_TYPE(result_type);
}
// check nothing for CGAL::Null_functor
template <typename ResultType>
void check_result_type(CGAL::Null_functor, ResultType){}

#define CGAL_SNAP_AST_FUNCTORS(Traits)                            \
  typedef typename Traits::Simplify Simplify ;                    \
  typedef typename Traits::Unit_part Unit_part;                   \
  typedef typename Traits::Integral_division Integral_division;   \
  typedef typename Traits::Divides Divides;                       \
  typedef typename Traits::Is_square Is_square;                   \
  typedef typename Traits::Gcd Gcd;                               \
  typedef typename Traits::Div_mod Div_mod;                       \
  typedef typename Traits::Div Div;                               \
  typedef typename Traits::Mod Mod;                               \
  typedef typename Traits::Square Square;                         \
  typedef typename Traits::Is_zero Is_zero;                       \
  typedef typename Traits::Is_one Is_one;                         \
  typedef typename Traits::Sqrt Sqrt;                             \
  typedef typename Traits::Kth_root Kth_root;                     \
  typedef typename Traits::Root_of Root_of;                       \
  CGAL_USE_TYPE(Simplify) ;                                       \
  CGAL_USE_TYPE(Unit_part);                                       \
  CGAL_USE_TYPE(Integral_division);                               \
  CGAL_USE_TYPE(Divides);                                         \
  CGAL_USE_TYPE(Is_square);                                       \
  CGAL_USE_TYPE(Gcd);                                             \
  CGAL_USE_TYPE(Div_mod);                                         \
  CGAL_USE_TYPE(Div);                                             \
  CGAL_USE_TYPE(Mod);                                             \
  CGAL_USE_TYPE(Square);                                          \
  CGAL_USE_TYPE(Is_zero);                                         \
  CGAL_USE_TYPE(Is_one);                                          \
  CGAL_USE_TYPE(Sqrt);                                            \
  CGAL_USE_TYPE(Kth_root);                                        \
  CGAL_USE_TYPE(Root_of);

namespace CGAL {

template< class  AS  >
bool test_equality_epsilon( const  AS & a,
                            const  AS & b,
                            const  AS & epsilon ) {
    typedef Algebraic_structure_traits<AS> AST;
    typedef typename AST::Is_exact Is_exact;
    if( Is_exact::value )
        return a == b;
    else {
        return ( a < (b + epsilon) ) &&
            ( a > (b - epsilon) );
    }
}

template< class AS >
AS unit_normal( const AS& x ) {
    typedef Algebraic_structure_traits<AS> AST;
    CGAL_SNAP_AST_FUNCTORS(AST);
    const Unit_part unit_part = Unit_part();
    const Integral_division integral_division = Integral_division();

    return integral_division( x, unit_part(x) );
}

//Syntax tests
template< class  AS  >
void test_algebraic_structure_intern(
        const CGAL::Integral_domain_without_division_tag& ) {}

template< class  AS  >
void test_algebraic_structure_intern( const CGAL::Integral_domain_tag& ) {

    // test of functors
    typedef CGAL::Algebraic_structure_traits< AS > AST;
    CGAL_SNAP_AST_FUNCTORS(AST);

    using CGAL::Null_functor;
    CGAL_static_assertion(
        (!::boost::is_same< Integral_division, Null_functor >::value));
    CGAL_static_assertion((!::boost::is_same< Divides, Null_functor >::value));
    CGAL_static_assertion((!::boost::is_same< Is_zero, Null_functor >::value));
    CGAL_static_assertion((!::boost::is_same< Is_one, Null_functor >::value));
    CGAL_static_assertion((!::boost::is_same< Square, Null_functor >::value));

    // functor
    const Is_zero is_zero = Is_zero();
    assert(  is_zero( AS (0)) );
    assert(! is_zero( AS (23)) );
    assert(  is_zero( AS (0) + AS(0) ) );
    assert(  CGAL_NTS is_zero( AS (0)) );
    assert(! CGAL_NTS is_zero( AS (23)) );
    assert(  CGAL_NTS is_zero( AS (0) + AS(0) ) );

    const Is_one is_one = Is_one();
    assert(  is_one( AS (1)) );
    assert(! is_one( AS (23)) );
    assert(  is_one( AS (1) + AS(0) ) );
    assert(  CGAL_NTS is_one( AS (1)) );
    assert(! CGAL_NTS is_one( AS (23)) );
    assert(  CGAL_NTS is_one( AS (1) + AS(0) ) );

    const Square square = Square();
    assert( square( AS (23)) ==  AS (23*23) );
    assert( CGAL_NTS square( AS (23)) ==  AS (23*23) );

    const Integral_division integral_division = Integral_division();
    AS a(6);
    AS b(2);
    AS c(3);
    assert( integral_division( a,b) == c);
    assert( integral_division( a,c) == b);
    assert( integral_division( a+a-a,c*b-c) == b);
    assert( CGAL_NTS integral_division( a,b) == c);
    assert( CGAL_NTS integral_division( a,c) == b);
    assert( CGAL_NTS integral_division( a+a-a,c*b-c) == b);


    const Divides divides = Divides();
    assert( divides(b,AS(0)));
    assert( divides(b,a));
    assert( divides(c,a));
    //assert( divides(c*b-c,a+a-a));
    assert( CGAL_NTS divides(b,AS(0)));
    assert( CGAL_NTS divides(b,a));
    assert( CGAL_NTS divides(c,a));
    //assert( CGAL_NTS divides(c*b-c,a+a-a));

    typedef typename AST::Is_exact Is_exact;
    // VC7 produced an ICE for
    // assert( ! Is_exact::value || ... );
    bool ie = Is_exact::value; (void) ie;
    AS tmp;
    assert( divides(b,AS(0),tmp));
    assert( !ie || tmp == integral_division(AS(0),b));
    assert( divides(b,a,tmp));
    assert( !ie || tmp == integral_division(a,b));
    assert( divides(c,a,tmp));
    assert( !ie || tmp == integral_division(a,c));
    assert( divides(c*b-c,a+a-a,tmp));
    assert( !ie || tmp == integral_division(a+a-a,c*b-c));



    assert( CGAL_NTS divides(b,AS(0),tmp));
    assert( !ie || tmp == integral_division(AS(0),b));
    assert( CGAL_NTS divides(b,a,tmp));
    assert( !ie || tmp == integral_division(a,b));
    assert( CGAL_NTS divides(c,a,tmp));
    assert( !ie || tmp == integral_division(a,c));
    assert( CGAL_NTS divides(AS(c*b-c),AS(a+a-a),tmp));
    assert( !ie || tmp == integral_division(a+a-a,c*b-c));
}

template< class  AS  >
void test_algebraic_structure_intern(
        const CGAL::Unique_factorization_domain_tag& ) {
    test_algebraic_structure_intern<  AS  >(CGAL::Integral_domain_tag());
    typedef CGAL::Algebraic_structure_traits<  AS  > AST;
    CGAL_SNAP_AST_FUNCTORS(AST);

    using CGAL::Null_functor;
    CGAL_static_assertion((!::boost::is_same< Gcd, Null_functor>::value));

    const Gcd gcd = Gcd();
    assert( gcd(  AS ( 0), AS ( 0))  ==  unit_normal( AS (0) ) );
    assert( gcd(  AS ( 7), AS ( 0))  ==  unit_normal( AS (7) ) );
    assert( gcd(  AS (-7), AS ( 0))  ==  unit_normal( AS (7) ) );
    assert( gcd(  AS ( 0), AS ( 7))  ==  unit_normal( AS (7) ) );
    assert( gcd(  AS ( 0), AS (-7))  ==  unit_normal( AS (7) ) );

    assert( gcd(  AS (-7), AS ( 1))  ==  unit_normal( AS (1) ) );
    assert( gcd(  AS ( 1), AS (-7))  ==  unit_normal( AS (1) ) );

    assert( gcd(  AS ( 15), AS ( 12))  ==  unit_normal( AS (3) ) );
    assert( gcd(  AS (-15), AS ( 12))  ==  unit_normal( AS (3) ) );
    assert( gcd(  AS ( 15), AS (-12))  ==  unit_normal( AS (3) ) );
    assert( gcd(  AS (-15), AS (-12))  ==  unit_normal( AS (3) ) );

    // special test for expression template, e.g. mpz_class
    assert( gcd(AS(-10)+AS(-5),AS(-4)*AS(-3))
            ==  unit_normal( AS (3) ) );


    assert( CGAL_NTS gcd(  AS ( 0), AS ( 0))
            ==  unit_normal( AS (0) ) );
    assert( CGAL_NTS gcd(  AS ( 7), AS ( 0))
            ==  unit_normal( AS (7) ) );
    assert( CGAL_NTS gcd(  AS (-7), AS ( 0))
            ==  unit_normal( AS (7) ) );
    assert( CGAL_NTS gcd(  AS ( 0), AS ( 7))
            ==  unit_normal( AS (7) ) );
    assert( CGAL_NTS gcd(  AS ( 0), AS (-7))
            ==  unit_normal( AS (7) ) );

    assert( CGAL_NTS gcd(  AS (-7), AS ( 1))
            ==  unit_normal( AS (1) ) );
    assert( CGAL_NTS gcd(  AS ( 1), AS (-7))
            ==  unit_normal( AS (1) ) );

    assert( CGAL_NTS gcd(  AS ( 15), AS ( 12))
            ==  unit_normal( AS (3) ) );
    assert( CGAL_NTS gcd(  AS (-15), AS ( 12))
            ==  unit_normal( AS (3) ) );
    assert( CGAL_NTS gcd(  AS ( 15), AS (-12))
            ==  unit_normal( AS (3) ) );
    assert( CGAL_NTS gcd(  AS (-15), AS (-12))
            ==  unit_normal( AS (3) ) );


    // special test for expression template, e.g. mpz_class
    assert( CGAL_NTS gcd(AS(-10)+AS(-5),AS(-4)*AS(-3))
            ==  unit_normal( AS (3) ) );
}

template< class  AS  >
void test_algebraic_structure_intern( const CGAL::Euclidean_ring_tag&) {
    test_algebraic_structure_intern<  AS  >
        ( CGAL::Unique_factorization_domain_tag() );

    typedef CGAL::Algebraic_structure_traits<  AS  > AST;
    CGAL_SNAP_AST_FUNCTORS(AST);

    using CGAL::Null_functor;
    CGAL_static_assertion((!::boost::is_same< Div,     Null_functor>::value));
    CGAL_static_assertion((!::boost::is_same< Mod,     Null_functor>::value));
    CGAL_static_assertion((!::boost::is_same< Div_mod, Null_functor>::value));

    const Div     div=Div();
    const Mod     mod=Mod();
    const Div_mod div_mod=Div_mod();

    // Rounding mode for div is: to zero
    assert( div(  AS ( 3),  AS (3)) ==  AS (1));
    assert( div(  AS ( 2),  AS (3)) ==  AS (0));
    assert( div(  AS ( 1),  AS (3)) ==  AS (0));
    assert( div(  AS ( 0),  AS (3)) ==  AS (0));
    assert( div(  AS (-1),  AS (3)) ==  AS (0));
    assert( div(  AS (-2),  AS (3)) ==  AS (0));
    assert( div(  AS (-3),  AS (3)) ==  AS (-1));

    assert( mod(  AS ( 3),  AS (3)) ==  AS (0));
    assert( mod(  AS ( 2),  AS (3)) ==  AS (2));
    assert( mod(  AS ( 1),  AS (3)) ==  AS (1));
    assert( mod(  AS ( 0),  AS (3)) ==  AS (0));
    assert( mod(  AS (-1),  AS (3)) ==  AS (-1));
    assert( mod(  AS (-2),  AS (3)) ==  AS (-2));
    assert( mod(  AS (-3),  AS (3)) ==  AS (0));

    assert( div(  AS ( 3),  AS(-3)) ==  AS (-1));
    assert( div(  AS ( 2),  AS(-3)) ==  AS (0));
    assert( div(  AS ( 1),  AS(-3)) ==  AS (0));
    assert( div(  AS ( 0),  AS(-3)) ==  AS (0));
    assert( div(  AS (-1),  AS(-3)) ==  AS (0));
    assert( div(  AS (-2),  AS(-3)) ==  AS (0));
    assert( div(  AS (-3),  AS(-3)) ==  AS (1));

    assert( mod(  AS ( 3),  AS(-3)) ==  AS (0));
    assert( mod(  AS ( 2),  AS(-3)) ==  AS (2));
    assert( mod(  AS ( 1),  AS(-3)) ==  AS (1));
    assert( mod(  AS ( 0),  AS(-3)) ==  AS (0));
    assert( mod(  AS (-1),  AS(-3)) ==  AS (-1));
    assert( mod(  AS (-2),  AS(-3)) ==  AS (-2));
    assert( mod(  AS (-3),  AS(-3)) ==  AS (0));

    for (int i = -12; i <= 12; i++){
        for (int j = 1; j < 10; j++){
            AS q,r;
            {
                AS a(i),b(j);
                div_mod(a,b,q,r);
                assert(q == div(a,b));
                assert(r == mod(a,b));
            }{
                AS a(i),b(-j);
                div_mod(a,b,q,r);
                assert(q == div(a,b));
                assert(r == mod(a,b));
            }
        }
    }

    // special syntax test for expression template, e.g. mpz_class
    assert( div(AS(-4)+AS(-4),AS(3)*AS(1)) ==  AS (-2));
    assert( mod(AS(-4)+AS(-4),AS(3)*AS(1)) ==  AS (-2));
    AS q,r;
    div_mod(AS(-4)+AS(-4),AS(3)*AS(1),q,r);
}

template< class  AS  >
void test_algebraic_structure_intern( const CGAL::Field_tag& ) {
    test_algebraic_structure_intern< AS >( CGAL::Integral_domain_tag());
     AS  a(1);
     AS  b(3);
     AS  c = a / b;
    (void)c;  // avoid warnings for unused variables

    typedef Algebraic_structure_traits<AS> AST;
    typedef typename AST::Is_exact Is_exact;
   // VC7 produced an ICE for
   // assert( ! Is_exact::value || ... );
    bool ie = Is_exact::value;
    assert( ! ie || c * b == a );
    a =  AS (1);
    a /=  AS (2);
    a /=  AS (2); // that must work correctly also for float types
    assert( a *  AS (4) ==  AS (1));

    typename AST::Divides divides;
    assert(divides(AS(2),AS(0)));
    assert(divides(AS(2),AS(5)));
    assert(divides(AS(5),AS(2)));
    AS tmp;
    assert(divides(AS(2),AS(0),tmp));
    assert(!ie || tmp == AS(0));
    assert(divides(AS(2),AS(5),tmp));
    assert(!ie || tmp == AS(5)/AS(2));
    assert(divides(AS(5),AS(2),tmp));
    assert(!ie || tmp == AS(2)/AS(5));

    assert(CGAL_NTS divides(AS(2),AS(0),tmp));
    assert(!ie || tmp == AS(0));
    assert(CGAL_NTS divides(AS(2),AS(5),tmp));
    assert(!ie || tmp == AS(5)/AS(2));
    assert(CGAL_NTS divides(AS(5),AS(2),tmp));
    assert(!ie || tmp == AS(2)/AS(5));

    typename AST::Inverse inverse;
    assert(AS(1)/AS(2) == inverse(AS(2)));
    assert(AS(1)/AS(2) == CGAL::inverse(AS(2)));


}

template <class  AS >
void test_algebraic_structure_intern( const CGAL::Field_with_sqrt_tag& ) {
    test_algebraic_structure_intern< AS >( CGAL::Field_tag());
    typedef CGAL::Algebraic_structure_traits<  AS  > AST;
    typedef typename AST::Is_exact Is_exact;

    CGAL_SNAP_AST_FUNCTORS(AST);

    CGAL_static_assertion((!::boost::is_same< Sqrt, Null_functor>::value));
    const Sqrt sqrt =Sqrt();
    AS  a(4);

    AS  c = sqrt( a);
    bool ie =  Is_exact::value;
    assert( ! ie ||   (c ==  AS (2))  );
    c =  AS (5);
    assert( !ie || sqrt(c) * sqrt(c) == c );
    (void)c;  // avoid warnings for unused variables
    // #### more involved square root and root tests
}


//semantic test
template <class  AS >
void test_algebraic_structure_intern(
        const  AS & a ,const  AS & b, const  AS & c,
        const CGAL::Integral_domain_without_division_tag&) {
    assert( a !=  AS (0));
    assert( b !=  AS (0));
    assert( c !=  AS (0));
    //  AS (0) == nullptr of IntegralDomain
    assert(a* AS (0)== AS (0));
    assert(a+ AS (0)==a);
    assert(b* AS (0)== AS (0));
    assert(b+ AS (0)==b);
    assert(c* AS (0)== AS (0));
    assert(c+ AS (0)==c);
    //  AS (1) == ONE of IntegralDomain
    assert(a* AS (1)==a);
    assert(b* AS (1)==b);
    assert(c* AS (1)==c);
    assert( AS (-1)* AS (-1)== AS (1));
    //associative
    assert((a+b)+c==a+(b+c));
    assert((a*b)*c==a*(b*c));
    //commutative
    assert(a+b+c==c+b+a);
    assert(a*b*c==c*b*a);
    //distributiv
    assert((a-b)*c==a*c-b*c);
    assert((a+b)*c==a*c+b*c);
    //binom
    assert((a+b)*(a+b)==a*a+ AS (2)*a*b+b*b);
    assert((a-b)*(a-b)==a*a- AS (2)*a*b+b*b);
    assert((a-b)*(a+b)==a*a-b*b);
    // unary operators
    assert(a==+a);
    assert(b==+b);
    assert(c==+c);
    assert(-a* AS (-1)==a);
    assert(-b* AS (-1)==b);
    assert(-c* AS (-1)==c);
}

 template <class  AS >
 void test_algebraic_structure_intern(
         const  AS & a ,const  AS & b, const  AS & c,
         const CGAL::Integral_domain_tag&) {
     assert( a !=  AS (0));
     assert( b !=  AS (0));
     assert( c !=  AS (0));
     test_algebraic_structure_intern(a,b,c,
             CGAL::Integral_domain_without_division_tag());

     typedef CGAL::Algebraic_structure_traits<  AS  > AST;
     CGAL_SNAP_AST_FUNCTORS(AST);

     //Integral_div
     const Integral_division integral_division = Integral_division();
     assert(integral_division(a*b,a)==b);
     assert(integral_division(a*c,a)==c);
     assert(integral_division(b*a,b)==a);
     assert(integral_division(b*c,b)==c);
     assert(integral_division(c*a,c)==a);
     assert(integral_division(c*b,c)==b);
     assert(CGAL_NTS integral_division(c*b,c)==b);

     const Divides divides = Divides();
     assert(divides(a,a*b));
     assert(divides(a,a*c));
     assert(divides(b,b*a));
     assert(divides(b,b*c));
     assert(divides(c,c*a));
     assert(divides(c,c*b));
     assert(CGAL_NTS divides(c,c*b));
 }

template <class  AS >
void test_algebraic_structure_intern(
        const  AS & a ,
        const  AS & b,
        const  AS & c,
        const CGAL::Unique_factorization_domain_tag&) {

    assert( a !=  AS (0));
    assert( b !=  AS (0));
    assert( c !=  AS (0));
    test_algebraic_structure_intern(a,b,c,CGAL::Integral_domain_tag());

    typedef CGAL::Algebraic_structure_traits<  AS  > AST;
    CGAL_SNAP_AST_FUNCTORS(AST);
    const Gcd gcd = Gcd();
    const Unit_part unit_part = Unit_part();
    assert(unit_part(a)*gcd(a*b,a*c)==a*gcd(b,c));
    assert(unit_part(b)*gcd(a*b,b*c)==b*gcd(a,c));
    assert(unit_part(c)*gcd(a*c,b*c)==c*gcd(b,a));
    assert(gcd( AS (0),a)*unit_part(a)==a);
    assert(gcd( AS (0),b)*unit_part(b)==b);
    assert(gcd( AS (0),c)*unit_part(c)==c);
    assert(gcd(-a,a)*unit_part(a)==a);
    assert(gcd(-b,b)*unit_part(b)==b);
    assert(gcd(-c,c)*unit_part(c)==c);
}

template <class  AS >
void test_algebraic_structure_intern(
        const  AS & a ,
        const  AS & b,
        const  AS & c,
        const CGAL::Euclidean_ring_tag&) {

    assert( a !=  AS (0));
    assert( b !=  AS (0));
    assert( c !=  AS (0));

    test_algebraic_structure_intern(a,b,c,
            CGAL::Unique_factorization_domain_tag());

    typedef CGAL::Algebraic_structure_traits<  AS  > AST;
    CGAL_SNAP_AST_FUNCTORS(AST);
    const Div div = Div();
    const Mod mod = Mod();
    const Div_mod div_mod = Div_mod();

    // do we have any
    AS  tmp_mod,tmp_div;
    div_mod(a,b,tmp_div,tmp_mod);
    assert(tmp_div==div(a,b));
    assert(tmp_mod==mod(a,b));
    assert(tmp_div*b+tmp_mod==a);

    div_mod(a,c,tmp_div,tmp_mod);
    assert(tmp_div==div(a,c));
    assert(tmp_mod==mod(a,c));
    assert(tmp_div*c+tmp_mod==a);

    div_mod(c,b,tmp_div,tmp_mod);
    assert(tmp_div==div(c,b));
    assert(tmp_mod==mod(c,b));
    assert(tmp_div*b+tmp_mod==c);
}

template <class  AS >
void test_algebraic_structure_intern(
        const  AS & a,
        const  AS & b,
        const  AS & c,
        const CGAL::Field_tag&) {

    assert( a !=  AS (0));
    assert( b !=  AS (0));
    assert( c !=  AS (0));

    test_algebraic_structure_intern(a,b,c,CGAL::Integral_domain_tag());

    AS  epsilon =  AS (1)/ AS (128);

    assert( test_equality_epsilon(  AS ((a/b)*b),
                                              AS ( a ), epsilon ) );
    assert( test_equality_epsilon(  AS ( (a/c)*c ),
                                              AS ( a ), epsilon ) );
    assert( test_equality_epsilon(  AS ( (b/a)*a ),
                                              AS ( b ), epsilon ) );
    assert( test_equality_epsilon(  AS ( (b/c)*c ),
                                              AS ( b ), epsilon ) );
    assert( test_equality_epsilon(  AS ( (c/b)*b ),
                                              AS ( c ), epsilon ) );
    assert( test_equality_epsilon(  AS ( (c/a)*a ),
                                              AS ( c ), epsilon ) );
}

template <class  AS >
void test_algebraic_structure_intern(
        const  AS & a ,
        const  AS & b,
        const  AS & c,
        const CGAL::Field_with_sqrt_tag&) {

    assert( a !=  AS (0));
    assert( b !=  AS (0));
    assert( c !=  AS (0));

    test_algebraic_structure_intern(a,b,c,CGAL::Field_tag());

    typedef CGAL::Algebraic_structure_traits<  AS  > AST;
    CGAL_SNAP_AST_FUNCTORS(AST);
    const Sqrt sqrt = Sqrt();

    AS  tmp;
    AS  epsilon =  AS (1);

    tmp=CGAL_NTS unit_part(a)*a;
    assert( test_equality_epsilon( sqrt(tmp)*sqrt(tmp),
                                             tmp, epsilon ) );

    tmp=CGAL_NTS unit_part(b)*b;
    assert( test_equality_epsilon( sqrt(tmp)*sqrt(tmp),
                                             tmp, epsilon ) );

    tmp=CGAL_NTS unit_part(c)*c;
    assert( test_equality_epsilon( sqrt(tmp)*sqrt(tmp),
                                             tmp, epsilon ) );
}

template< class AS, class Is_square >
class Test_is_square {
  public:
    void operator()( const Is_square& is_square ) {
        typedef typename Is_square::first_argument_type First_argument_type;
        typedef typename Is_square::second_argument_type Second_argument_type;
        typedef typename Is_square::result_type   Result_type;
        CGAL_USE_TYPE(First_argument_type);
        CGAL_USE_TYPE(Second_argument_type);

        CGAL_static_assertion(
                ( ::boost::is_same< AS , First_argument_type>::value));
        CGAL_static_assertion(
                ( ::boost::is_same< AS& , Second_argument_type>::value));
        //CGAL_static_assertion(( ::boost::is_same< bool , Result_type>::value));
        bool b = Result_type(true); CGAL_USE(b);

        AS test_number = AS(3)*AS(3);
        AS result;
        assert( is_square( test_number));
        assert( is_square( test_number, result ));
        assert( test_equality_epsilon( result , AS(3), AS(1) ) );

        assert( CGAL_NTS is_square( test_number));
        assert( CGAL_NTS is_square( test_number, result ));
        assert( test_equality_epsilon( result , AS(3), AS(1) ) );
    }
};

template<class  AS >
class Test_is_square< AS , CGAL::Null_functor> {
public:
    void operator() (CGAL::Null_functor) {
        // ok, nothing to test
    }
};

template<class  AS , class Sqrt>
class Test_sqrt {
public:
    void operator() (const Sqrt& sqrt) {
        CGAL_USE(sqrt);
        typedef typename Sqrt::argument_type Argument_type;
        typedef typename Sqrt::result_type   Result_type;
        CGAL_USE_TYPE(Argument_type);
        CGAL_USE_TYPE(Result_type);
        CGAL_static_assertion(( ::boost::is_same< AS , Argument_type>::value));
        CGAL_static_assertion(( ::boost::is_same< AS , Result_type>::value));
        typedef Algebraic_structure_traits<AS> AST;
        typedef typename AST::Is_exact Is_exact;
        assert( !Is_exact::value ||  AS (3) == sqrt( AS (9)));
    }
};

template<class  AS >
class Test_sqrt< AS , CGAL::Null_functor> {
public:
    void operator() (CGAL::Null_functor) {
        // ok, nothing to test
    }
};

template<class  AS , class Root>
class Test_root {
public:
    void operator() (const Root& root) {
        typedef typename Root::first_argument_type  First_argument_type;
        typedef typename Root::second_argument_type Second_argument_type;
        typedef typename Root::result_type          Result_type;
        CGAL_USE_TYPE(First_argument_type);
        CGAL_USE_TYPE(Second_argument_type);
        CGAL_USE_TYPE(Result_type);
        CGAL_static_assertion(
                ( ::boost::is_same<int, First_argument_type>::value));
        CGAL_static_assertion(
                ( ::boost::is_same< AS , Second_argument_type>::value));
        CGAL_static_assertion(
                ( ::boost::is_same< AS , Result_type>::value));
         AS  epsilon(1);
        assert( test_equality_epsilon(  AS (2),
                                root( 4,  AS (16) ), epsilon ) );
        assert( test_equality_epsilon(  AS (3),
                                root( 3,  AS (27) ), epsilon ) );
    }
};

template<class  AS >
class Test_root< AS , CGAL::Null_functor> {
public:
    void operator() (CGAL::Null_functor) {
        // ok, nothing to test
    }
};

// Type_functions -----------------------------------------------
template <class  AS >
void test_Type_functions(
        const CGAL::Integral_domain_without_division_tag&) {
     AS  x(-15);
    CGAL_NTS simplify(x);
    assert(x== AS (-15));
    CGAL_NTS unit_part(x);
    assert( CGAL_NTS is_zero(  AS (0)) );
    assert( CGAL_NTS is_one(  AS (1)) );
    assert( CGAL_NTS square(  AS (23)) ==  AS (23*23) );
}

template <class  AS >
void test_Type_functions( const CGAL::Integral_domain_tag&) {
    test_Type_functions< AS >
        (CGAL::Integral_domain_without_division_tag());
    assert(CGAL_NTS integral_division( AS (10), AS (2))== AS (5));
}

template <class  AS >
void test_Type_functions(
        const CGAL::Unique_factorization_domain_tag&) {
    test_Type_functions< AS >(CGAL::Integral_domain_tag());

    assert(CGAL_NTS gcd( AS (21), AS (15)) == unit_normal(AS (3)));

}

template <class  AS >
void test_Type_functions( const CGAL::Euclidean_ring_tag&) {
    test_Type_functions< AS >(
            CGAL::Unique_factorization_domain_tag());
    //std::cerr << CGAL_NTS mod( AS(14), AS(5) ) << std::endl;
    AS  q,r,a,b;
    a = AS(14);
    b = AS(5);
    r = CGAL_NTS mod(a,b);
    q = CGAL_NTS div(a,b);
    assert( a == b*q+r);
    CGAL_NTS div_mod(a,b,q,r);
    assert( a == b*q+r);
}

template <class  AS >
void test_Type_functions( const CGAL::Field_tag&) {
    test_Type_functions< AS >(CGAL::Integral_domain_tag());
    assert(CGAL_NTS unit_part( AS (-15))== AS (-15));
    assert(CGAL_NTS unit_part( AS (1  ))== AS (  1));
    assert(CGAL_NTS unit_part( AS (0  ))== AS (  1));
}

template <class  AS >
void test_Type_functions( const CGAL::Field_with_sqrt_tag&) {
   test_Type_functions< AS >(CGAL::Field_tag());
   typedef Algebraic_structure_traits<AS> AST;
   typedef typename AST::Is_exact Is_exact;
   AS  c = CGAL_NTS sqrt( AS (4));
   bool ie = Is_exact::value;
   assert( !ie || c ==  AS (2) );
}
template <class  AS >
void test_Type_functions( const CGAL::Field_with_root_of_tag&) {
   test_Type_functions< AS >(CGAL::Field_with_sqrt_tag());
   std::vector< AS > coeffs(4);
   coeffs[0]= AS (-27);
   coeffs[1]= AS (0);
   coeffs[2]= AS (0);
   coeffs[3]= AS (1);

   typedef CGAL::Algebraic_structure_traits<  AS  > AST;
   CGAL_SNAP_AST_FUNCTORS(AST);

//   typedef typename Root_of::Boundary Boundary;

   const Root_of root_of = Root_of();
   AS  real = root_of(1,coeffs.begin(),coeffs.end());
   assert( real ==  AS (3));
   assert(real*real ==  AS (9));
   assert(real-real ==  AS (0));
   assert(CGAL_NTS sqrt(real) == CGAL_NTS sqrt( AS (3)));

  // Test the function
  assert( CGAL_NTS root_of(1, coeffs.begin(), coeffs.end() ) ==
                    root_of( 1, coeffs.begin(), coeffs.end() ) );

/*
     AS  real2 = root_of(  AS (2),
     AS (4),
     coeffs.begin(), coeffs.end() );

     assert( real2 ==  AS (3));
     assert(real2*real2 ==  AS (9));
     assert(real2-real2 ==  AS (0));
     assert(CGAL::sqrt(real2) == CGAL::sqrt( AS (3)));
*/
}

template <class  AS , class Algebraic_category, class Is_exact>
void test_algebraic_structure(){

    test_Type_functions< AS >(Algebraic_category());

    typedef CGAL::Algebraic_structure_traits<  AS  > AST;
    CGAL_SNAP_AST_FUNCTORS(AST);

    CGAL_static_assertion((::boost::is_same<AS,typename AST::Type>::value));

    typedef typename AST::Boolean Boolean;
    assert(!Boolean());
    check_result_type(Is_zero(), Boolean());
    check_result_type(Is_one(), Boolean());
    check_result_type(Divides(), Boolean());
    check_result_type(Is_square(), Boolean());

    typedef typename AST::Algebraic_category  Tag;
    using CGAL::Integral_domain_without_division_tag;
    using CGAL::Null_functor;
    // Test for desired exactness
    CGAL_static_assertion(
            ( ::boost::is_same< typename AST::Is_exact, Is_exact >::value));

    CGAL_static_assertion(( ::boost::is_convertible< Tag,
                    Integral_domain_without_division_tag >::value ));
    CGAL_static_assertion(( ::boost::is_same< Tag, Algebraic_category>::value));
    CGAL_static_assertion((!::boost::is_same< Simplify, Null_functor>::value));
    CGAL_static_assertion((!::boost::is_same< Unit_part, Null_functor>::value));
    const Simplify   simplify=Simplify();;
    const Unit_part  unit_part= Unit_part();

    // the other functors must exist as well, but they might be Null_functor's
    Integral_division  integral_div;
    Gcd                gcd;
    Div                div;
    Mod                mod;
    Div_mod            div_mod;
    Sqrt               sqrt;
    Kth_root           root;
    Is_square          is_square;
//    typename Traits::Find_only_zero_element find_only_zero_element;
//    typename Traits::Find_only_equal_pair   find_only_equal_pair;

    (void)integral_div; // avoid warnings for unused variables
    (void)gcd;
    (void)div;
    (void)mod;
    (void)div_mod;
    (void)sqrt;
    (void)root;
    (void)is_square;
//    (void)find_only_zero_element;
//    (void)find_only_equal_pair;

     AS  a; // DefaultConstructible
     AS  b(127); // construction from small integers
     AS  c(-127);
    a = b; // Assignable
     AS  d(a);
    assert( a == b); // EqualityComparable
    assert( a == d); // EqualityComparable
    assert( c != b);
    assert( c != d);
    a = +b;
    assert( a == b); // == 127
    a = -c;
    assert( a == b); // == 127
    a =  AS (1);
    b =  AS (5);
    c = a + b;
    assert( c ==  AS (6));
    c = a - b;
    assert( c ==  AS (-4));
    c = a * b;
    assert( c ==  AS (5));
    c = (a +  AS (1)) * b;
    assert( c ==  AS (10));
    c = b *  AS (-3);
    assert( c ==  AS (-15));
    c = a;
    c += b;
    assert( c ==  AS (6));
    c = a;
    c -= b;
    assert( c ==  AS (-4));
    c = a;
    c *= b;
    assert( c ==  AS (5));
    c = a +  AS (1);
    c *= b;
    assert( c ==  AS (10));
    simplify(c);
    assert( c ==  AS (10));
    unit_part( c);
    test_algebraic_structure_intern< AS >( Tag());
    Test_sqrt< AS , typename AST::Sqrt> tsr;
    tsr(sqrt);
    Test_root< AS , typename AST::Kth_root> trt;
    trt(root);
    Test_is_square< AS, typename AST::Is_square > tis;
    tis(is_square);
    {
        std::vector< AS > v(5);
        v[0]= AS (30);
        v[1]= AS (-2);
        v[2]= AS (0);
        v[3]= AS (5);
        v[4]= AS (5);
        // functor
//        assert(2 == find_only_zero_element(v.begin(),v.end()));
        // function
        // CGAL::find_only_zero_element is not yet implemented
//        assert(2 == CGAL::find_only_zero_element(v.begin(),v.end()));

        std::vector< AS > w(3);
        w[0]= AS (10);
        w[1]= AS (-8);
        w[2]= AS (-2);

//      functor
//      find_only_equal_pair is not yet implemented
//      std::pair<int,int> equal_pair(0,0);
//
//      equal_pair = find_only_equal_pair(v.begin(),v.end(),
//                                          w.begin(),w.end());
//      assert(1 == equal_pair.first);
//      assert(2 == equal_pair.second);
        //function
//      equal_pair = NiX::find_only_equal_pair(v.begin(),v.end(),
//                                               w.begin(),w.end());
//      assert(1 == equal_pair.first);
//      assert(2 == equal_pair.second);

    }
}

template <class  AS , class Algebraic_category, class Is_exact >
void test_algebraic_structure( const  AS & a, const  AS & b, const  AS & c) {

    assert( a !=  AS (0));
    assert( b !=  AS (0));
    assert( c !=  AS (0));
    test_algebraic_structure< AS ,Algebraic_category, Is_exact>();
    test_algebraic_structure_intern(a,b,c,Algebraic_category());

    typedef CGAL::Algebraic_structure_traits<AS> AST;
    typedef typename AST::Is_numerical_sensitive Is_numerical_sensitive;
    CGAL_static_assertion(
            !(::boost::is_same<Is_numerical_sensitive, CGAL::Null_tag>::value));
    CGAL_USE_TYPE(Is_numerical_sensitive);
}

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_TEST_ALGEBRAIC_STRUCTURE_H
