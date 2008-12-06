// Copyright (c) 2008 Max-Planck-Institute Saarbruecken (Germany)
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
// Author(s)     : Michael Hemmer <mhemmer@uni-mainz.de>
                   Dominik Huelse <dominik.huelse@gmx.de>
//
// ============================================================================

/*! \file CGAL/Algebraic_number_traits.h
 *  \brief Defines traits class CGAL::Algebraic_structure_traits_extended. 
 */

#ifndef ALGEBRAIC_STRUCTURE_TRAITS_EXTENDED_H
#define ALGEBRAIC_STRUCTURE_TRAITS_EXTENDED_H

#include <CGAL/basic.h>
#include <CGAL/type_traits.h>
#include <CGAL/Algebraic_structure_traits.h>
#include <CGAL/Polynomial.h>        

namespace CGAL {


// The algebraic structure traits extended base class
// =========================================================================
template< class Type, class Algebraic_category >
class Algebraic_structure_traits_extended_base;

template< class Type, class Algebraic_category >
class Sqrt_algebraic_structure_traits_extended;

template< class Type, class Algebraic_category >
class Polynomial_algebraic_structure_traits_extended;

template< class Type_ > 
class Algebraic_structure_traits_extended;

//! The template specialization that can be used for types that are not any
//! of the number type concepts. The divides functor is set to \c Null_functor.
template< class Type_ >
class Algebraic_structure_traits_extended_base< Type_, CGAL::Null_tag > 
    : public CGAL::Algebraic_structure_traits< Type_> {
public:
    typedef Type_       Type;
    typedef CGAL::Null_functor Divides;
};

//! The template specialization that is used if the number type is
//! a model of the \c IntegralDomainWithoutDiv concept.
template< class Type_ > 
class Algebraic_structure_traits_extended_base< Type_, 
                                       Integral_domain_without_division_tag > 
    : public Algebraic_structure_traits_extended_base< Type_, 
                                              CGAL::Null_tag > {
public:
    typedef Type_                   Type;

};

//! The template specialization that is used if the number type is
//! a model of the \c IntegralDomain concept.
template< class Type_ >
class Algebraic_structure_traits_extended_base< Type_, 
                                       Integral_domain_tag >
    : public Algebraic_structure_traits_extended_base< Type_, 
                                       Integral_domain_without_division_tag > { 

};

//! The template specialization that is used if the number type is
//! a model of the \c UFDomain concept. It is equivalent to the specialization
//! for the \c IntegralDomain concept.
template< class Type_ >
class Algebraic_structure_traits_extended_base< Type_, 
                                       Unique_factorization_domain_tag >
    : public Algebraic_structure_traits_extended_base< Type_, 
                                              Integral_domain_tag > {

public:
    typedef Type_  Type;
  
    class Divides { 
    public:
        typedef Type    first_argument_type;
        typedef Type    second_argument_type;
        typedef Type&   third_argument_type;
        typedef bool  result_type;
        bool operator()( const Type& x, 
                const Type& y, 
                Type& q) const {            
            typedef CGAL::Algebraic_structure_traits<Type> AST;
            typename AST::Gcd gcd; 
//            std::cout<<" UFD gcd computation"<<std::endl;
            if(gcd(y,x) == x/unit_part(x)){
                q = y/x;
                return true;
            }
            else
                return false;
        }
    };

};

//! The template specialization that is used if the number type is
//! a model of the \c EuclideanRing concept.
template< class Type_ >
class Algebraic_structure_traits_extended_base< Type_, 
                                       Euclidean_ring_tag >
    : public Algebraic_structure_traits_extended_base< Type_, 
                                              Unique_factorization_domain_tag > {
public:
    typedef Type_        Type;
    
    class Divides {
    public:
        typedef Type    first_argument_type;
        typedef Type    second_argument_type;
        typedef Type&   third_argument_type;
        typedef bool  result_type;
        bool operator()( const Type& x, 
                const Type& y, 
                Type& q) const {
            typedef CGAL::Algebraic_structure_traits<Type> AST; 
            typename AST::Div div;
//            std::cout<<"euclidean ring "<<std::endl;
            q = div(y,x);
            if(q*x==y)
                return true;
            else
                return false;
        }
    };

};

//! The template specialization that is used if the number type is
//! a model of the \c Field concept.
template< class Type_ >
class Algebraic_structure_traits_extended_base< Type_, Field_tag >
    : public Algebraic_structure_traits_extended_base< Type_, 
                                              Integral_domain_tag > {
public:
    typedef Type_        Type;

    class Divides { 
    public:
        typedef Type    first_argument_type;
        typedef Type    second_argument_type;
        typedef Type&   third_argument_type;
        typedef bool  result_type;
        bool operator()( const Type& x, 
                const Type& y, 
                Type& q) const {
//            std::cout<<" divides field"<<std::endl;
            q = y/x;
            return true;
        }
    };
};

//! The template specialization that is used if the number type is a model
//! of the \c FieldWithSqrt concept. It is equivalent to the 
//! specialization for the \c Field concept.
template< class Type_ >
class Algebraic_structure_traits_extended_base< Type_, 
                                       Field_with_sqrt_tag>
    : public Algebraic_structure_traits_extended_base< Type_, 
                                              Field_tag> {
public:
    typedef Type_        Type;
};

//! The template specialization that is used if the number type is a model
//! of the \c FieldWithKthRoot concept. It is equivalent to the 
//! specialization for the \c Field concept.
template< class Type_ >
class Algebraic_structure_traits_extended_base< Type_, 
                                       Field_with_kth_root_tag>
    : public Algebraic_structure_traits_extended_base< Type_, 
                                              Field_with_sqrt_tag> {    
    
public:
    typedef Type_        Type;
};

//! The template specialization that is used if the number type is a model
//! of the \c FieldWithRootOf concept. It is equivalent to the 
//! specialization for the \c FieldWithKthRoot concept.
template< class Type_ >
class Algebraic_structure_traits_extended_base< Type_, 
                                       Field_with_root_of_tag >
    : public Algebraic_structure_traits_extended_base< Type_, 
                                              Field_with_kth_root_tag > {
public:
    typedef Type_           Type;
};


//! The template specialization that is used for Sqrt_extensions if the COEFF type is a model
//! of the \c IntegralDomain concept.
template< class Type_, class ROOT >
class Sqrt_algebraic_structure_traits_extended< Sqrt_extension< Type_, ROOT >, 
                                                        CGAL::Integral_domain_tag > 
    : public Algebraic_structure_traits_extended_base< Sqrt_extension< Type_, ROOT >, 
                                                              CGAL::Integral_domain_tag >{
public:
    typedef Type_ COEFF;
    typedef Sqrt_extension< COEFF, ROOT >  Type;   

private:
    typedef typename CGAL::Coercion_traits< ROOT, COEFF >::Cast Root_nt_cast;

public:
    class Divides { 
    public:
        typedef Type    first_argument_type;
        typedef Type    second_argument_type;
        typedef Type&   third_argument_type;
        typedef bool  result_type;
        bool operator()( const Type& x, 
                const Type& y, 
                Type& q) const {       
            typedef CGAL::Algebraic_structure_traits_extended<COEFF> ASTE;
            typename ASTE::Divides divides;   
//            std::cout<<"integral domain for sqrt"<<std::endl;
            bool result;
            COEFF q1, q2;
            if(x.is_extended()){
//                std::cout<<" y is extended "<<std::endl;
                COEFF denom = x.a0()*x.a0() - x.a1()*x.a1() * Root_nt_cast()(x.root());
                if ( denom == COEFF(0) ) {   
                    // this is for the rare case in which root is a square
                    // and the (pseudo) algebraic conjugate of p becomes zero 
                    result = divides(COEFF(2)*x.a0(),y.a0(),q1);
                    if(!result) return false;
                    result = divides(COEFF(2)*x.a1(),y.a1(),q2);
                    if(!result) return false;
                    q = Type(q1 + q2); 
                }else{
                    q = y;
                    q *= Type(x.a0(),-x.a1(),x.root());
                    result = divides(denom,q.a0(),q1);
                    if(!result) return false;
                    result = divides(denom,q.a1(),q2);
                    if(!result) return false;
                    q =  Type( q1, q2, y.root());
                }
            }else{
//                std::cout<<" x is not extended "<<std::endl;
                if(y.is_extended()){
                    result = divides(x.a0(),y.a0(),q1);
                    if(!result) return false;
                    result = divides(x.a0(),y.a1(),q2);  
                    if(!result) return false;
                    q = Type(q1, q2, y.root());
                }else{
                    result = divides(x.a0(),y.a0(),q1);
                    if(!result) return false;
                    q = q1;
                }
            }
            if(q*x==y)
                return true; 
            else
                return false;
        }
    }; 
    
};

//! The template specialization that is used for Sqrt_extensions if the COEFF type is a model
//! of the \c UFDomain concept.
template< class Type_, class ROOT >
class Sqrt_algebraic_structure_traits_extended< Sqrt_extension< Type_, ROOT >, 
                                                   CGAL::Unique_factorization_domain_tag >
    : public Sqrt_algebraic_structure_traits_extended< Sqrt_extension< Type_, ROOT >, 
                                                              CGAL::Integral_domain_tag >
{ };


//! The template specialization that is used for Sqrt_extensions if the COEFF type is a model
//! of the \c EuclideanRing concept.
template< class Type_, class ROOT >
class Sqrt_algebraic_structure_traits_extended< Sqrt_extension< Type_, ROOT >, 
                                                          CGAL::Euclidean_ring_tag >
  : public Sqrt_algebraic_structure_traits_extended< Sqrt_extension< Type_, ROOT >, 
                                                             CGAL::Integral_domain_tag >
{
public:
    typedef Type_ COEFF;
    typedef Sqrt_extension< COEFF, ROOT >  Type;   
};


//! The template specialization that is used for polynomials if the number type is a model
//! of the \c IntegralDomain concept.
template< class Type_ >
class Polynomial_algebraic_structure_traits_extended< CGAL::Polynomial< Type_ >, 
                                                            Integral_domain_tag > 
    : public Algebraic_structure_traits_extended_base< CGAL::Polynomial< Type_ >, 
                                                              CGAL::Integral_domain_tag >{
  
public:
    typedef Type_ NT;
    typedef CGAL::Polynomial< NT >  Type;
//    typedef Unique_factorization_domain_tag Algebraic_category;

    class Divides { 
    public:
        typedef Type    first_argument_type;
        typedef Type    second_argument_type;
        typedef Type&   third_argument_type;
        typedef bool  result_type;
        bool operator()( const Type& p1, 
                const Type& p2, 
                Type& q) const {
            q=Type(0);
//            std::cout<<"new Sqrt POLY"<<std::endl;
            
            typedef CGAL::Algebraic_structure_traits_extended<NT> ASTE;
            typename ASTE::Divides divides;   
            NT q1;
            bool result;
            int d1 = p1.degree();
            int d2 = p2.degree();
            if (p2.is_zero()) {
                q=Type(0); 
                return true;
            }

            if ( d2 < d1 ) {
                q = Type(0);
                return false;
            }

            typedef std::vector<NT> Vector;
            Vector V_R, V_Q;    
            V_Q.reserve(d2);
            if(d1==0){
                for(int i=d2;i>=0;--i){   
                    result=divides(p1[0],p2[i],q1);
                    if(!result) return false;
                    V_Q.push_back(q1);
                }
                V_R.push_back(NT(0));
            }
            else{        
                V_R.reserve(d2);
                V_R=Vector(p2.begin(),p2.end());
                Vector tmp1;
                tmp1.reserve(d1);
                for(int k=0; k<=d2-d1; ++k){
                    result=divides(p1[d1],V_R[d2-k],q1);
                    if(!result) return false;
                    V_Q.push_back(q1);  
                    for(int j=0;j<d1;++j){   
                        tmp1.push_back(p1[j]*V_Q[k]);
                    }   
                    V_R[d2-k]=0;            
                    for(int i=d2-d1-k;i<=d2-k-1;++i){
                        V_R[i]=V_R[i]-tmp1[i-(d2-d1-k)];
                    }   
                    tmp1.clear();
                }
        

            }
            q = Type(V_Q.rbegin(),V_Q.rend());
            Type r = Type(V_R.begin(),V_R.end());
            if(r == Type(0))
                return true;
            else
                return false;
        }
    };
};

//! The template specialization that is used for polynomials if the COEFF type is a model
//! of the \c UFDomain concept.
template< class Type_ >
class Polynomial_algebraic_structure_traits_extended< CGAL::Polynomial< Type_ >, 
                                                   CGAL::Unique_factorization_domain_tag >
    : public Polynomial_algebraic_structure_traits_extended< CGAL::Polynomial< Type_ >, 
                                                              CGAL::Integral_domain_tag >
{ };


//! The template specialization that is used for polynomials if the number type is a model
//! of the \c EuclideanRing concept.
template< class Type_ >
class Polynomial_algebraic_structure_traits_extended< CGAL::Polynomial< Type_ >, 
                                                               Euclidean_ring_tag > 
    : public Polynomial_algebraic_structure_traits_extended< CGAL::Polynomial< Type_ >, 
                                                              Integral_domain_tag >{
  
public:
    typedef Type_ NT;
    typedef CGAL::Polynomial< Type_ >  Type;
//    typedef Unique_factorization_domain_tag Algebraic_category;
};

//! The template specialization that is used for polynomials if the number type is a model
//! of the \c Field concept.
template< class Type_ >
class Polynomial_algebraic_structure_traits_extended< CGAL::Polynomial< Type_ >, Field_tag > 
    : public Polynomial_algebraic_structure_traits_extended< CGAL::Polynomial< Type_ >, 
                                                              Integral_domain_tag >{
  
public:
    typedef CGAL::Polynomial< Type_ >  Type;
//    typedef Euclidean_ring_tag Algebraic_category;

    class Divides { 
    public:
        typedef Type    first_argument_type;
        typedef Type    second_argument_type;
        typedef Type&   third_argument_type;
        typedef bool  result_type;
        bool operator()( const Type& p1, 
                const Type& p2, 
                Type& q) const {

//            std::cout<<"Field POLY"<<std::endl;
            
            if (p2.is_zero()) {
                q = Type(0);
                return true;
            }

            Type r;
            Type::euclidean_division(p2, p1, q, r);
            if(r == Type(0))
                return true;
            else
                return false;
        }
    };

};

template< class Type_ > 
class Algebraic_structure_traits_extended
    : public Algebraic_structure_traits_extended_base< Type_, 
             typename CGAL::Algebraic_structure_traits<Type_>::Algebraic_category > {

};

// The actual algebraic structure extended traits class
template< class Type_, class ROOT > 
class Algebraic_structure_traits_extended< Sqrt_extension< Type_, ROOT > >
    : public Sqrt_algebraic_structure_traits_extended< Sqrt_extension< Type_, ROOT >, 
 typename CGAL::Algebraic_structure_traits<Type_>::Algebraic_category > {

};

// The actual algebraic structure extended traits class
template< class Type_ > 
class Algebraic_structure_traits_extended< CGAL::Polynomial< Type_ > >
    : public Polynomial_algebraic_structure_traits_extended< CGAL::Polynomial< Type_ >, 
 typename CGAL::Algebraic_structure_traits<Type_>::Algebraic_category > {

};

} //namespace CGAL

#endif  // CGAL_ALGEBRAIC_STRUCTURE_TRAITS_EXTENDED_H
// EOF
