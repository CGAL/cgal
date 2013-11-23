// Copyright (c) 2006-2008 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
// Author(s)     : Lutz Kettner  <kettner@mpi-inf.mpg.de>
//                 Michael Hemmer <hemmer@mpi-inf.mpg.de>
//
// =============================================================================

/*   \brief provides test functions for the various ring concepts of number
 *          types.
*/

// within this file AS ^= Type

#include <CGAL/basic.h>
#include <CGAL/Algebraic_structure_traits.h>
//#include <CGAL/Real_embeddable_traits.h>

#include <CGAL/assertions.h>
#include <CGAL/use.h>
#include <boost/type_traits.hpp>
#include <CGAL/tags.h>
#include <functional>

#include <vector>

#include <CGAL/Testsuite/use.h>

#ifndef CGAL_TEST_ALGEBRAIC_STRUCTURE_H
#define CGAL_TEST_ALGEBRAIC_STRUCTURE_H

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
    CGAL_assertion(  is_zero( AS (0)) );
    CGAL_assertion(! is_zero( AS (23)) );
    CGAL_assertion(  is_zero( AS (0) + AS(0) ) );  
    CGAL_assertion(  CGAL_NTS is_zero( AS (0)) );
    CGAL_assertion(! CGAL_NTS is_zero( AS (23)) );
    CGAL_assertion(  CGAL_NTS is_zero( AS (0) + AS(0) ) );      
    
    const Is_one is_one = Is_one();
    CGAL_assertion(  is_one( AS (1)) );
    CGAL_assertion(! is_one( AS (23)) );
    CGAL_assertion(  is_one( AS (1) + AS(0) ) );
    CGAL_assertion(  CGAL_NTS is_one( AS (1)) );
    CGAL_assertion(! CGAL_NTS is_one( AS (23)) );
    CGAL_assertion(  CGAL_NTS is_one( AS (1) + AS(0) ) );
    
    const Square square = Square();
    CGAL_assertion( square( AS (23)) ==  AS (23*23) );
    CGAL_assertion( CGAL_NTS square( AS (23)) ==  AS (23*23) );
    
    const Integral_division integral_division = Integral_division();
    AS a(6);
    AS b(2);
    AS c(3);
    CGAL_assertion( integral_division( a,b) == c);
    CGAL_assertion( integral_division( a,c) == b);
    CGAL_assertion( integral_division( a+a-a,c*b-c) == b);
    CGAL_assertion( CGAL_NTS integral_division( a,b) == c);
    CGAL_assertion( CGAL_NTS integral_division( a,c) == b);
    CGAL_assertion( CGAL_NTS integral_division( a+a-a,c*b-c) == b);


    const Divides divides = Divides();
    CGAL_assertion( divides(b,AS(0)));
    CGAL_assertion( divides(b,a));
    CGAL_assertion( divides(c,a));
    //CGAL_assertion( divides(c*b-c,a+a-a));
    CGAL_assertion( CGAL_NTS divides(b,AS(0))); 
    CGAL_assertion( CGAL_NTS divides(b,a)); 
    CGAL_assertion( CGAL_NTS divides(c,a));
    //CGAL_assertion( CGAL_NTS divides(c*b-c,a+a-a));
    
    typedef typename AST::Is_exact Is_exact;
    // VC7 produced an ICE for
    // CGAL_assertion( ! Is_exact::value || ... );
    bool ie = Is_exact::value; (void) ie;
    AS tmp; 
    CGAL_assertion( divides(b,AS(0),tmp));
    CGAL_assertion( !ie || tmp == integral_division(AS(0),b));
    CGAL_assertion( divides(b,a,tmp)); 
    CGAL_assertion( !ie || tmp == integral_division(a,b));
    CGAL_assertion( divides(c,a,tmp));
    CGAL_assertion( !ie || tmp == integral_division(a,c));
    CGAL_assertion( divides(c*b-c,a+a-a,tmp));
    CGAL_assertion( !ie || tmp == integral_division(a+a-a,c*b-c));
    

    
    CGAL_assertion( CGAL_NTS divides(b,AS(0),tmp)); 
    CGAL_assertion( !ie || tmp == integral_division(AS(0),b));
    CGAL_assertion( CGAL_NTS divides(b,a,tmp)); 
    CGAL_assertion( !ie || tmp == integral_division(a,b));
    CGAL_assertion( CGAL_NTS divides(c,a,tmp));
    CGAL_assertion( !ie || tmp == integral_division(a,c));
    CGAL_assertion( CGAL_NTS divides(AS(c*b-c),AS(a+a-a),tmp));
    CGAL_assertion( !ie || tmp == integral_division(a+a-a,c*b-c));   
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
    CGAL_assertion( gcd(  AS ( 0), AS ( 0))  ==  unit_normal( AS (0) ) );
    CGAL_assertion( gcd(  AS ( 7), AS ( 0))  ==  unit_normal( AS (7) ) );
    CGAL_assertion( gcd(  AS (-7), AS ( 0))  ==  unit_normal( AS (7) ) );
    CGAL_assertion( gcd(  AS ( 0), AS ( 7))  ==  unit_normal( AS (7) ) );
    CGAL_assertion( gcd(  AS ( 0), AS (-7))  ==  unit_normal( AS (7) ) );
    
    CGAL_assertion( gcd(  AS (-7), AS ( 1))  ==  unit_normal( AS (1) ) );
    CGAL_assertion( gcd(  AS ( 1), AS (-7))  ==  unit_normal( AS (1) ) );
    
    CGAL_assertion( gcd(  AS ( 15), AS ( 12))  ==  unit_normal( AS (3) ) );
    CGAL_assertion( gcd(  AS (-15), AS ( 12))  ==  unit_normal( AS (3) ) );
    CGAL_assertion( gcd(  AS ( 15), AS (-12))  ==  unit_normal( AS (3) ) );
    CGAL_assertion( gcd(  AS (-15), AS (-12))  ==  unit_normal( AS (3) ) );

    // special test for expression template, e.g. mpz_class 
    CGAL_assertion( gcd(AS(-10)+AS(-5),AS(-4)*AS(-3))  
            ==  unit_normal( AS (3) ) );
 

    CGAL_assertion( CGAL_NTS gcd(  AS ( 0), AS ( 0))  
            ==  unit_normal( AS (0) ) );
    CGAL_assertion( CGAL_NTS gcd(  AS ( 7), AS ( 0))  
            ==  unit_normal( AS (7) ) );
    CGAL_assertion( CGAL_NTS gcd(  AS (-7), AS ( 0))  
            ==  unit_normal( AS (7) ) );
    CGAL_assertion( CGAL_NTS gcd(  AS ( 0), AS ( 7))  
            ==  unit_normal( AS (7) ) );
    CGAL_assertion( CGAL_NTS gcd(  AS ( 0), AS (-7))  
            ==  unit_normal( AS (7) ) );
    
    CGAL_assertion( CGAL_NTS gcd(  AS (-7), AS ( 1))  
            ==  unit_normal( AS (1) ) );
    CGAL_assertion( CGAL_NTS gcd(  AS ( 1), AS (-7))  
            ==  unit_normal( AS (1) ) );
    
    CGAL_assertion( CGAL_NTS gcd(  AS ( 15), AS ( 12))  
            ==  unit_normal( AS (3) ) );
    CGAL_assertion( CGAL_NTS gcd(  AS (-15), AS ( 12)) 
            ==  unit_normal( AS (3) ) );
    CGAL_assertion( CGAL_NTS gcd(  AS ( 15), AS (-12))  
            ==  unit_normal( AS (3) ) );
    CGAL_assertion( CGAL_NTS gcd(  AS (-15), AS (-12))  
            ==  unit_normal( AS (3) ) );


    // special test for expression template, e.g. mpz_class 
    CGAL_assertion( CGAL_NTS gcd(AS(-10)+AS(-5),AS(-4)*AS(-3))  
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
    CGAL_assertion( div(  AS ( 3),  AS (3)) ==  AS (1));
    CGAL_assertion( div(  AS ( 2),  AS (3)) ==  AS (0));
    CGAL_assertion( div(  AS ( 1),  AS (3)) ==  AS (0));
    CGAL_assertion( div(  AS ( 0),  AS (3)) ==  AS (0));
    CGAL_assertion( div(  AS (-1),  AS (3)) ==  AS (0));
    CGAL_assertion( div(  AS (-2),  AS (3)) ==  AS (0));
    CGAL_assertion( div(  AS (-3),  AS (3)) ==  AS (-1));

    CGAL_assertion( mod(  AS ( 3),  AS (3)) ==  AS (0));
    CGAL_assertion( mod(  AS ( 2),  AS (3)) ==  AS (2));
    CGAL_assertion( mod(  AS ( 1),  AS (3)) ==  AS (1));
    CGAL_assertion( mod(  AS ( 0),  AS (3)) ==  AS (0));
    CGAL_assertion( mod(  AS (-1),  AS (3)) ==  AS (-1));
    CGAL_assertion( mod(  AS (-2),  AS (3)) ==  AS (-2));
    CGAL_assertion( mod(  AS (-3),  AS (3)) ==  AS (0));

    CGAL_assertion( div(  AS ( 3),  AS(-3)) ==  AS (-1));
    CGAL_assertion( div(  AS ( 2),  AS(-3)) ==  AS (0));
    CGAL_assertion( div(  AS ( 1),  AS(-3)) ==  AS (0));
    CGAL_assertion( div(  AS ( 0),  AS(-3)) ==  AS (0));
    CGAL_assertion( div(  AS (-1),  AS(-3)) ==  AS (0));
    CGAL_assertion( div(  AS (-2),  AS(-3)) ==  AS (0));
    CGAL_assertion( div(  AS (-3),  AS(-3)) ==  AS (1));

    CGAL_assertion( mod(  AS ( 3),  AS(-3)) ==  AS (0));
    CGAL_assertion( mod(  AS ( 2),  AS(-3)) ==  AS (2));
    CGAL_assertion( mod(  AS ( 1),  AS(-3)) ==  AS (1));
    CGAL_assertion( mod(  AS ( 0),  AS(-3)) ==  AS (0));
    CGAL_assertion( mod(  AS (-1),  AS(-3)) ==  AS (-1));
    CGAL_assertion( mod(  AS (-2),  AS(-3)) ==  AS (-2));
    CGAL_assertion( mod(  AS (-3),  AS(-3)) ==  AS (0));

    for (int i = -12; i <= 12; i++){
        for (int j = 1; j < 10; j++){
            AS q,r;
            {               
                AS a(i),b(j);
                div_mod(a,b,q,r);
                CGAL_assertion(q == div(a,b));
                CGAL_assertion(r == mod(a,b));
            }{
                AS a(i),b(-j);
                div_mod(a,b,q,r);
                CGAL_assertion(q == div(a,b));
                CGAL_assertion(r == mod(a,b));
            }
        }
    }
    
    // special syntax test for expression template, e.g. mpz_class
    CGAL_assertion( div(AS(-4)+AS(-4),AS(3)*AS(1)) ==  AS (-2));
    CGAL_assertion( mod(AS(-4)+AS(-4),AS(3)*AS(1)) ==  AS (-2));
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
   // CGAL_assertion( ! Is_exact::value || ... );
    CGAL_assertion_code(bool ie = Is_exact::value;)
    CGAL_assertion( ! ie || c * b == a );
    a =  AS (1);
    a /=  AS (2);
    a /=  AS (2); // that must work correctly also for float types
    CGAL_assertion( a *  AS (4) ==  AS (1));
    
    CGAL_assertion_code(typename AST::Divides divides;)
    CGAL_assertion(divides(AS(2),AS(0)));
    CGAL_assertion(divides(AS(2),AS(5)));
    CGAL_assertion(divides(AS(5),AS(2)));
    AS tmp;
    CGAL_assertion(divides(AS(2),AS(0),tmp));
    CGAL_assertion(!ie || tmp == AS(0));
    CGAL_assertion(divides(AS(2),AS(5),tmp));
    CGAL_assertion(!ie || tmp == AS(5)/AS(2));
    CGAL_assertion(divides(AS(5),AS(2),tmp));
    CGAL_assertion(!ie || tmp == AS(2)/AS(5));

    CGAL_assertion(CGAL_NTS divides(AS(2),AS(0),tmp));
    CGAL_assertion(!ie || tmp == AS(0));
    CGAL_assertion(CGAL_NTS divides(AS(2),AS(5),tmp));
    CGAL_assertion(!ie || tmp == AS(5)/AS(2));
    CGAL_assertion(CGAL_NTS divides(AS(5),AS(2),tmp));
    CGAL_assertion(!ie || tmp == AS(2)/AS(5));

    CGAL_assertion_code(typename AST::Inverse inverse;)
    CGAL_assertion(AS(1)/AS(2) == inverse(AS(2)));
    CGAL_assertion(AS(1)/AS(2) == CGAL::inverse(AS(2)));
    
    
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
    CGAL_assertion_code(bool ie =  Is_exact::value;)
    CGAL_assertion( ! ie ||   (c ==  AS (2))  );   
    c =  AS (5);
    CGAL_assertion( !ie || sqrt(c) * sqrt(c) == c );
    (void)c;  // avoid warnings for unused variables
    // #### more involved square root and root tests
}
	
	
//semantic test
template <class  AS >
void test_algebraic_structure_intern( 
        const  AS & a ,const  AS & b, const  AS & c,
        const CGAL::Integral_domain_without_division_tag&) {
    CGAL_assertion( a !=  AS (0));
    CGAL_assertion( b !=  AS (0));
    CGAL_assertion( c !=  AS (0));
    //  AS (0) == NULL of IntegralDomain
    CGAL_assertion(a* AS (0)== AS (0));
    CGAL_assertion(a+ AS (0)==a);
    CGAL_assertion(b* AS (0)== AS (0));
    CGAL_assertion(b+ AS (0)==b);
    CGAL_assertion(c* AS (0)== AS (0));
    CGAL_assertion(c+ AS (0)==c);
    //  AS (1) == ONE of IntegralDomain
    CGAL_assertion(a* AS (1)==a);
    CGAL_assertion(b* AS (1)==b);
    CGAL_assertion(c* AS (1)==c);
    CGAL_assertion( AS (-1)* AS (-1)== AS (1));
    //associative
    CGAL_assertion((a+b)+c==a+(b+c));
    CGAL_assertion((a*b)*c==a*(b*c));
    //commutative
    CGAL_assertion(a+b+c==c+b+a);
    CGAL_assertion(a*b*c==c*b*a);
    //distributiv
    CGAL_assertion((a-b)*c==a*c-b*c);
    CGAL_assertion((a+b)*c==a*c+b*c);
    //binom
    CGAL_assertion((a+b)*(a+b)==a*a+ AS (2)*a*b+b*b);
    CGAL_assertion((a-b)*(a-b)==a*a- AS (2)*a*b+b*b);
    CGAL_assertion((a-b)*(a+b)==a*a-b*b);
    // unary operators
    CGAL_assertion(a==+a);  
    CGAL_assertion(b==+b);  
    CGAL_assertion(c==+c);
    CGAL_assertion(-a* AS (-1)==a);  
    CGAL_assertion(-b* AS (-1)==b);  
    CGAL_assertion(-c* AS (-1)==c);
}
 
 template <class  AS >
 void test_algebraic_structure_intern( 
         const  AS & a ,const  AS & b, const  AS & c,
         const CGAL::Integral_domain_tag&) {
     CGAL_assertion( a !=  AS (0));
     CGAL_assertion( b !=  AS (0));
     CGAL_assertion( c !=  AS (0));
     test_algebraic_structure_intern(a,b,c,
             CGAL::Integral_domain_without_division_tag());

     typedef CGAL::Algebraic_structure_traits<  AS  > AST;
     CGAL_SNAP_AST_FUNCTORS(AST);
     
     //Integral_div
     const Integral_division integral_division = Integral_division();
     CGAL_assertion(integral_division(a*b,a)==b);
     CGAL_assertion(integral_division(a*c,a)==c);
     CGAL_assertion(integral_division(b*a,b)==a);
     CGAL_assertion(integral_division(b*c,b)==c);
     CGAL_assertion(integral_division(c*a,c)==a);
     CGAL_assertion(integral_division(c*b,c)==b);
     CGAL_assertion(CGAL_NTS integral_division(c*b,c)==b);
     
     const Divides divides = Divides();
     CGAL_assertion(divides(a,a*b));
     CGAL_assertion(divides(a,a*c));
     CGAL_assertion(divides(b,b*a));
     CGAL_assertion(divides(b,b*c));
     CGAL_assertion(divides(c,c*a));
     CGAL_assertion(divides(c,c*b));
     CGAL_assertion(CGAL_NTS divides(c,c*b));
 }
    
template <class  AS >
void test_algebraic_structure_intern( 
        const  AS & a ,
        const  AS & b, 
        const  AS & c,
        const CGAL::Unique_factorization_domain_tag&) {

    CGAL_assertion( a !=  AS (0));
    CGAL_assertion( b !=  AS (0));
    CGAL_assertion( c !=  AS (0));
    test_algebraic_structure_intern(a,b,c,CGAL::Integral_domain_tag());
    
    typedef CGAL::Algebraic_structure_traits<  AS  > AST;
    CGAL_SNAP_AST_FUNCTORS(AST);
    const Gcd gcd = Gcd();
    const Unit_part unit_part = Unit_part();
    CGAL_assertion(unit_part(a)*gcd(a*b,a*c)==a*gcd(b,c));
    CGAL_assertion(unit_part(b)*gcd(a*b,b*c)==b*gcd(a,c));
    CGAL_assertion(unit_part(c)*gcd(a*c,b*c)==c*gcd(b,a));
    CGAL_assertion(gcd( AS (0),a)*unit_part(a)==a);
    CGAL_assertion(gcd( AS (0),b)*unit_part(b)==b);
    CGAL_assertion(gcd( AS (0),c)*unit_part(c)==c);
    CGAL_assertion(gcd(-a,a)*unit_part(a)==a);
    CGAL_assertion(gcd(-b,b)*unit_part(b)==b);
    CGAL_assertion(gcd(-c,c)*unit_part(c)==c);
}

template <class  AS >
void test_algebraic_structure_intern( 
        const  AS & a ,
        const  AS & b, 
        const  AS & c,
        const CGAL::Euclidean_ring_tag&) {
    
    CGAL_assertion( a !=  AS (0));
    CGAL_assertion( b !=  AS (0));
    CGAL_assertion( c !=  AS (0));
    
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
    CGAL_assertion(tmp_div==div(a,b));
    CGAL_assertion(tmp_mod==mod(a,b));
    CGAL_assertion(tmp_div*b+tmp_mod==a);

    div_mod(a,c,tmp_div,tmp_mod);
    CGAL_assertion(tmp_div==div(a,c));
    CGAL_assertion(tmp_mod==mod(a,c));
    CGAL_assertion(tmp_div*c+tmp_mod==a);

    div_mod(c,b,tmp_div,tmp_mod);
    CGAL_assertion(tmp_div==div(c,b));
    CGAL_assertion(tmp_mod==mod(c,b));
    CGAL_assertion(tmp_div*b+tmp_mod==c);
}

template <class  AS >
void test_algebraic_structure_intern( 
        const  AS & a,
        const  AS & b, 
        const  AS & c,
        const CGAL::Field_tag&) {
    
    CGAL_assertion( a !=  AS (0));
    CGAL_assertion( b !=  AS (0));
    CGAL_assertion( c !=  AS (0));
    
    test_algebraic_structure_intern(a,b,c,CGAL::Integral_domain_tag());
    
    AS  epsilon =  AS (1)/ AS (128);

    CGAL_assertion( test_equality_epsilon(  AS ((a/b)*b), 
                                              AS ( a ), epsilon ) );
    CGAL_assertion( test_equality_epsilon(  AS ( (a/c)*c ),
                                              AS ( a ), epsilon ) );
    CGAL_assertion( test_equality_epsilon(  AS ( (b/a)*a ),
                                              AS ( b ), epsilon ) );
    CGAL_assertion( test_equality_epsilon(  AS ( (b/c)*c ),
                                              AS ( b ), epsilon ) );
    CGAL_assertion( test_equality_epsilon(  AS ( (c/b)*b ),
                                              AS ( c ), epsilon ) );
    CGAL_assertion( test_equality_epsilon(  AS ( (c/a)*a ),
                                              AS ( c ), epsilon ) );
}
  
template <class  AS >
void test_algebraic_structure_intern( 
        const  AS & a ,
        const  AS & b, 
        const  AS & c,
        const CGAL::Field_with_sqrt_tag&) {
    
    CGAL_assertion( a !=  AS (0));
    CGAL_assertion( b !=  AS (0));
    CGAL_assertion( c !=  AS (0));
    
    test_algebraic_structure_intern(a,b,c,CGAL::Field_tag());
    
    typedef CGAL::Algebraic_structure_traits<  AS  > AST;
    CGAL_SNAP_AST_FUNCTORS(AST);
    const Sqrt sqrt = Sqrt();
   
    AS  tmp;
    AS  epsilon =  AS (1);
    
    tmp=CGAL_NTS unit_part(a)*a;
    CGAL_assertion( test_equality_epsilon( sqrt(tmp)*sqrt(tmp),
                                             tmp, epsilon ) );
   
    tmp=CGAL_NTS unit_part(b)*b;
    CGAL_assertion( test_equality_epsilon( sqrt(tmp)*sqrt(tmp),
                                             tmp, epsilon ) );
   
    tmp=CGAL_NTS unit_part(c)*c;
    CGAL_assertion( test_equality_epsilon( sqrt(tmp)*sqrt(tmp),
                                             tmp, epsilon ) );
}
   
template< class AS, class Is_square > 
class Test_is_square {
  public:
    void operator()( const Is_square& is_square ) {
        typedef typename Is_square::first_argument_type First_argument_type;
        typedef typename Is_square::second_argument_type Second_argument_type;
        typedef typename Is_square::result_type   Result_type;
        CGAL_static_assertion(
                ( ::boost::is_same< AS , First_argument_type>::value));
        CGAL_static_assertion(
                ( ::boost::is_same< AS& , Second_argument_type>::value));
        //CGAL_static_assertion(( ::boost::is_same< bool , Result_type>::value));
        bool b = Result_type(true); CGAL_USE(b);
 
        AS test_number = AS(3)*AS(3);
        AS result;
        CGAL_assertion( is_square( test_number));
        CGAL_assertion( is_square( test_number, result ));
        CGAL_assertion( test_equality_epsilon( result , AS(3), AS(1) ) );

        CGAL_assertion( CGAL_NTS is_square( test_number));
        CGAL_assertion( CGAL_NTS is_square( test_number, result ));
        CGAL_assertion( test_equality_epsilon( result , AS(3), AS(1) ) );
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
        typedef typename Sqrt::argument_type Argument_type;
        typedef typename Sqrt::result_type   Result_type;
        CGAL_static_assertion(( ::boost::is_same< AS , Argument_type>::value));
        CGAL_static_assertion(( ::boost::is_same< AS , Result_type>::value));
        typedef Algebraic_structure_traits<AS> AST;
        typedef typename AST::Is_exact Is_exact; 
        CGAL_assertion( !Is_exact::value ||  AS (3) == sqrt( AS (9)));
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
        CGAL_static_assertion(
                ( ::boost::is_same<int, First_argument_type>::value));
        CGAL_static_assertion(
                ( ::boost::is_same< AS , Second_argument_type>::value));
        CGAL_static_assertion(
                ( ::boost::is_same< AS , Result_type>::value));
         AS  epsilon(1);
        CGAL_assertion( test_equality_epsilon(  AS (2),
                                root( 4,  AS (16) ), epsilon ) );
        CGAL_assertion( test_equality_epsilon(  AS (3),
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
    CGAL_assertion(x== AS (-15));
    CGAL_NTS unit_part(x);
    CGAL_assertion( CGAL_NTS is_zero(  AS (0)) );
    CGAL_assertion( CGAL_NTS is_one(  AS (1)) );
    CGAL_assertion( CGAL_NTS square(  AS (23)) ==  AS (23*23) );
}
 
template <class  AS >
void test_Type_functions( const CGAL::Integral_domain_tag&) {
    test_Type_functions< AS >
        (CGAL::Integral_domain_without_division_tag());
    CGAL_assertion(CGAL_NTS integral_division( AS (10), AS (2))== AS (5));
}
    
template <class  AS >
void test_Type_functions( 
        const CGAL::Unique_factorization_domain_tag&) {
    test_Type_functions< AS >(CGAL::Integral_domain_tag());
        
    CGAL_assertion(CGAL_NTS gcd( AS (21), AS (15)) == unit_normal(AS (3)));
    
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
    CGAL_assertion( a == b*q+r);
    CGAL_NTS div_mod(a,b,q,r);
    CGAL_assertion( a == b*q+r);
}
 
template <class  AS >
void test_Type_functions( const CGAL::Field_tag&) {
    test_Type_functions< AS >(CGAL::Integral_domain_tag());
    CGAL_assertion(CGAL_NTS unit_part( AS (-15))== AS (-15));
    CGAL_assertion(CGAL_NTS unit_part( AS (1  ))== AS (  1));
    CGAL_assertion(CGAL_NTS unit_part( AS (0  ))== AS (  1));
}
  
template <class  AS >
void test_Type_functions( const CGAL::Field_with_sqrt_tag&) {   
   test_Type_functions< AS >(CGAL::Field_tag());  
   typedef Algebraic_structure_traits<AS> AST;
   typedef typename AST::Is_exact Is_exact; 
   AS  c = CGAL_NTS sqrt( AS (4));
   CGAL_assertion_code(bool ie = Is_exact::value;)
   CGAL_assertion( !ie || c ==  AS (2) );
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
   CGAL_assertion( real ==  AS (3));
   CGAL_assertion(real*real ==  AS (9));
   CGAL_assertion(real-real ==  AS (0));
   CGAL_assertion(CGAL_NTS sqrt(real) == CGAL_NTS sqrt( AS (3)));
      
  // Test the function
  CGAL_assertion( CGAL_NTS root_of(1, coeffs.begin(), coeffs.end() ) ==
                    root_of( 1, coeffs.begin(), coeffs.end() ) );      
      
/*   
     AS  real2 = root_of(  AS (2), 
     AS (4),
     coeffs.begin(), coeffs.end() );
     
     CGAL_assertion( real2 ==  AS (3));
     CGAL_assertion(real2*real2 ==  AS (9));
     CGAL_assertion(real2-real2 ==  AS (0));
     CGAL_assertion(CGAL::sqrt(real2) == CGAL::sqrt( AS (3)));
*/
}

template <class  AS , class Algebraic_category, class Is_exact>
void test_algebraic_structure(){
   
    test_Type_functions< AS >(Algebraic_category());
  
    typedef CGAL::Algebraic_structure_traits<  AS  > AST;
    CGAL_SNAP_AST_FUNCTORS(AST);
    
    CGAL_static_assertion((::boost::is_same<AS,typename AST::Type>::value));
    
    typedef typename AST::Boolean Boolean;
    CGAL_assertion(!Boolean());
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
    CGAL_assertion( a == b); // EqualityComparable
    CGAL_assertion( a == d); // EqualityComparable
    CGAL_assertion( c != b);
    CGAL_assertion( c != d);
    a = +b;
    CGAL_assertion( a == b); // == 127
    a = -c;
    CGAL_assertion( a == b); // == 127
    a =  AS (1);
    b =  AS (5);
    c = a + b;
    CGAL_assertion( c ==  AS (6));
    c = a - b;
    CGAL_assertion( c ==  AS (-4));
    c = a * b;
    CGAL_assertion( c ==  AS (5));
    c = (a +  AS (1)) * b;
    CGAL_assertion( c ==  AS (10));
    c = b *  AS (-3);
    CGAL_assertion( c ==  AS (-15));
    c = a;
    c += b;
    CGAL_assertion( c ==  AS (6));
    c = a;
    c -= b;
    CGAL_assertion( c ==  AS (-4));
    c = a;
    c *= b;
    CGAL_assertion( c ==  AS (5));
    c = a +  AS (1);
    c *= b;
    CGAL_assertion( c ==  AS (10));
    simplify(c);
    CGAL_assertion( c ==  AS (10));
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
//        CGAL_assertion(2 == find_only_zero_element(v.begin(),v.end()));
        // function
        // CGAL::find_only_zero_element is not yet implemented
//        CGAL_assertion(2 == CGAL::find_only_zero_element(v.begin(),v.end()));

        std::vector< AS > w(3);
        w[0]= AS (10);
        w[1]= AS (-8);
        w[2]= AS (-2);
        std::pair<int,int> equal_pair(0,0);
        // functor
        // find_only_equal_pair is not yet implemented
//        equal_pair = find_only_equal_pair(v.begin(),v.end(),
//                                          w.begin(),w.end());
//        CGAL_assertion(1 == equal_pair.first);
//        CGAL_assertion(2 == equal_pair.second);
        //function
//        equal_pair = NiX::find_only_equal_pair(v.begin(),v.end(),
//                                               w.begin(),w.end());
//        CGAL_assertion(1 == equal_pair.first);
//        CGAL_assertion(2 == equal_pair.second);
        
    } 
}

template <class  AS , class Algebraic_category, class Is_exact >
void test_algebraic_structure( const  AS & a, const  AS & b, const  AS & c) {
 
    CGAL_assertion( a !=  AS (0));
    CGAL_assertion( b !=  AS (0));
    CGAL_assertion( c !=  AS (0));
    test_algebraic_structure< AS ,Algebraic_category, Is_exact>();
    test_algebraic_structure_intern(a,b,c,Algebraic_category());
    
    typedef CGAL::Algebraic_structure_traits<AS> AST;
    typedef typename AST::Is_numerical_sensitive Is_numerical_sensitive;
    CGAL_static_assertion(
            !(::boost::is_same<Is_numerical_sensitive, CGAL::Null_tag>::value));
    CGAL_USE_TYPE(Is_numerical_sensitive);
}
  
} //namespace CGAL

#endif // CGAL_TEST_ALGEBRAIC_STRUCTURE_H
