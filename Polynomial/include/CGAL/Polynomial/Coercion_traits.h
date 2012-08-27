// Copyright (c) 2008 Max-Planck-Institute Saarbruecken (Germany)
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
// Author(s)     :  Michael Hemmer <hemmer@mpi-inf.mpg.de>
//
// ============================================================================

#ifndef CGAL_POLYNOMIAL_COERCION_TRAITS_H
#define CGAL_POLYNOMIAL_COERCION_TRAITS_H

// The coercion type of two polynomials is a polynomial in d=max(d1,d2) 
// variables, where d1 and d2 are the number of variables the two
// polynomials. (This also includes the case of d1 = 0 or d2 = 0.) 
// Moreover, the new Innermost_coefficient_type is the coercion type of the 
// two Innermost_coefficient_types of the two involved polynomials. 
// (Again, this is generalized if one of the involved types is just a scalar 
// type)
// Though the coercion type is clear, the problem is how to match the 
// variables. The recursive definition of Polynomial<Coeff> suggest that 
// the coercion type of two polynomial types Polynomial<A> and Polynomial<B>
// is defined as Polynomial<C>, where C is the coercion type. 
// However, this is not in line with the fact that a Polynomial<A> 
// is interoperable with its coefficient type A, that is, if A is a polynomial 
// the variables of A should not be moved outward while casting A to 
// Polynomial<A>. 

#include <CGAL/Polynomial/misc.h>

namespace CGAL {

namespace internal{

// A has less variables than B
template <typename A, typename B, bool less >
struct Coercion_traits_for_polynomial_comp_d
  :public Coercion_traits_for_polynomial_comp_d< B, A , false >{};

// Polynomial<A> has more variables than B 
template <typename A, typename B >
struct Coercion_traits_for_polynomial_comp_d< Polynomial<A>, B , false>{
  typedef Coercion_traits<A,B> CT;

    typedef CGAL::Tag_true  Are_explicit_interoperable;
    typedef CGAL::Tag_false Are_implicit_interoperable;

    typedef Polynomial<typename CT::Type> Type;
    struct Cast{                                      
        typedef Type result_type;                               
        Type operator()(const Polynomial<A>& poly) const {
            typename CT::Cast cast;
            return Type(::boost::make_transform_iterator(poly.begin(),cast),
                       ::boost::make_transform_iterator(poly.end()  ,cast));
        } 
        Type operator()(const B& x) const {
            typename CT::Cast cast;
            return Type(cast(x));
        } 
    };                          
};

// number of variables is different 
template <typename A, typename B, int a, int b>
struct Coercion_traits_for_polynomial_equal_d
  :public Coercion_traits_for_polynomial_comp_d <A,B, a < b >{};

// number of variables is equal and at least one.
template <class A,class B, int d>
struct Coercion_traits_for_polynomial_equal_d<Polynomial<A>, Polynomial<B>, d, d >{
    typedef Coercion_traits<A,B> CT;            

    typedef CGAL::Tag_true  Are_explicit_interoperable;
    typedef CGAL::Tag_false Are_implicit_interoperable;
    typedef Polynomial<typename CT::Type> Type;
    struct Cast{                                      
        typedef Type result_type;                               
        Type operator()(const Polynomial<A>& poly) const { 
            typename CT::Cast cast; 
            return Type(::boost::make_transform_iterator(poly.begin(),cast),
                    ::boost::make_transform_iterator(poly.end()  ,cast));
        } 
        Type operator()(const Polynomial<B>& poly) const {  
            typename CT::Cast cast;  
            return Type(::boost::make_transform_iterator(poly.begin(),cast),
                    ::boost::make_transform_iterator(poly.end()  ,cast));
        } 
    }; 
};

// determine number of variables in each polynomial
template <typename A, typename B>
struct Coercion_traits_for_polynomial
  : public Coercion_traits_for_polynomial_equal_d
   < A , B , Dimension<A>::value, Dimension<B>::value >{};

}// namespace internal 

template <class A,class B>
struct Coercion_traits_for_level< Polynomial<A> , Polynomial<B>, CTL_POLYNOMIAL >
  :public internal::Coercion_traits_for_polynomial< Polynomial<A>, Polynomial<B> >
{};
template <class A,class B>
struct Coercion_traits_for_level< Polynomial<A> , B , CTL_POLYNOMIAL >
  :public internal::Coercion_traits_for_polynomial< Polynomial<A>, B >
{};
template <class A,class B>
struct Coercion_traits_for_level< A , Polynomial<B> , CTL_POLYNOMIAL >
  :public internal::Coercion_traits_for_polynomial< A , Polynomial<B> >
{};




#if 0 
// COERCION_TRAITS BEGIN 

//Coercion_traits_polynomial-----------------------------------
// If there is a Polynomial_traits, valid for more than one Polynomial
// class this part should be adapted, using a Polynomial_traits 
// and the nesting_depth 
template <class A,class B>
struct Coercion_traits_for_level<Polynomial<A>, Polynomial<B>, CTL_POLYNOMIAL >{
    typedef Coercion_traits<A,B> CT;            

    typedef CGAL::Tag_true  Are_explicit_interoperable;
    typedef CGAL::Tag_false Are_implicit_interoperable;
    typedef Polynomial<typename CT::Type> Type;
    struct Cast{                                      
        typedef Type result_type;                               
        Type operator()(const Polynomial<A>& poly) const { 
            typename CT::Cast cast; 
            return Type(::boost::make_transform_iterator(poly.begin(),cast),
                    ::boost::make_transform_iterator(poly.end()  ,cast));
        } 
        Type operator()(const Polynomial<B>& poly) const {  
            typename CT::Cast cast;  
            return Type(::boost::make_transform_iterator(poly.begin(),cast),
                    ::boost::make_transform_iterator(poly.end()  ,cast));
        } 
    }; 
};
        
template <class A,class B>
struct Coercion_traits_for_level<Polynomial<A>,B ,CTL_POLYNOMIAL >{
    typedef Coercion_traits<A,B> CT;

    typedef CGAL::Tag_true  Are_explicit_interoperable;
    typedef CGAL::Tag_false Are_implicit_interoperable;

    typedef Polynomial<typename CT::Type> Type;
    struct Cast{                                      
        typedef Type result_type;                               
        Type operator()(const Polynomial<A>& poly) const {
            typename CT::Cast cast;
            return Type(::boost::make_transform_iterator(poly.begin(),cast),
                       ::boost::make_transform_iterator(poly.end()  ,cast));
        } 
        Type operator()(const B& x) const {
            typename CT::Cast cast;
            return Type(cast(x));
        } 
    };                                                        
}; 
template <class A,class B> 
struct Coercion_traits_for_level<B,Polynomial<A>,CTL_POLYNOMIAL  >
    :public Coercion_traits_for_level<Polynomial<A>,B,CTL_POLYNOMIAL >
{};

#endif // 0 

// COERCION_TRAITS END

} //namespace CGAL

#endif // CGAL_POLYNOMIAL_COERCION_TRAITS_H
