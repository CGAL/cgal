// Copyright (c) 2005,2006  INRIA Sophia-Antipolis (France)
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
// SPDX-License-Identifier: LGPL-3.0+
//
//
// Author(s)     : Sylvain Pion, Monique Teillaud, Athanasios Kakargias, Michael Hemmer

#ifndef CGAL_ROOT_OF_TRAITS_H
#define CGAL_ROOT_OF_TRAITS_H

#include <CGAL/number_type_basic.h>
#include <CGAL/Get_arithmetic_kernel.h>
#include <CGAL/Sqrt_extension.h>
#include <CGAL/Quotient.h>
#include <boost/mpl/has_xxx.hpp>

namespace CGAL {

template < typename NT >
struct Root_of_traits;

template < class NT >
inline
typename Root_of_traits< NT >::Root_of_2
make_root_of_2(const NT &a, const NT &b, const NT &c)
{
    typename Root_of_traits<NT>::Make_root_of_2 make_root_of_2;
    return make_root_of_2(a,b,c);
}

template < class NT >
inline
typename Root_of_traits< NT >::Root_of_2
make_root_of_2(const NT &a, const NT &b, const NT &c,const bool smaller)
{
    typename Root_of_traits<NT>::Make_root_of_2 make_root_of_2;
    return make_root_of_2(a,b,c,smaller);
}

template < class NT >
inline 
typename Root_of_traits< NT >::Root_of_2
make_sqrt(const NT &a)
{
  typename Root_of_traits<NT>::Make_sqrt make_sqrt;
  return make_sqrt(a);
}

template < class NT , class OutputIterator>
inline 
OutputIterator
compute_roots_of_2(const NT &a_, const NT &b_, const NT &c_, OutputIterator oit)
{
  typedef typename Root_of_traits<NT>::Root_of_1 Root_of_1;
  typedef typename Root_of_traits<NT>::Root_of_2 Root_of_2;
  typename CGAL::Coercion_traits<Root_of_1,NT>::Cast cast; 
  Root_of_1 a(cast(a_)), b(cast(b_)), c(cast(c_));
    
  if ( a != 0 ) {
    Root_of_1 a0_  (-b/(2*a));
    Root_of_1 root_(CGAL_NTS square(a0_) - c/a);
    switch(CGAL::sign(root_)){
    case CGAL::NEGATIVE: return oit; 
    case CGAL::ZERO: *oit++ = Root_of_2(a0_);  return oit;
    default:
      // two roots 
      *oit++ = make_root_of_2(a0_,Root_of_1(-1),root_);
      *oit++ = make_root_of_2(a0_,Root_of_1( 1),root_);
      return oit; 
    }
  }
  else { 
    *oit++ = -c/b; return oit;   
  }
}


namespace internal {

BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Has_typedef_Arithmetic_kernel,Arithmetic_kernel,false)  

template <class NT,bool has_AK=Has_typedef_Arithmetic_kernel<Get_arithmetic_kernel<NT> >::value>
struct Get_rational_type{
  typedef Quotient<NT> type;
};

template <class NT>
struct Get_rational_type<NT,true>{
  typedef typename Get_arithmetic_kernel<NT>::Arithmetic_kernel::Rational type;
};

  
//Default or not a field.
//If no specialization of Get_arithmetic_kernel is available, a field type compatible with NT 
//is made using CGAL::Quotient
template < typename NT, class Algebraic_category>
struct Root_of_traits_helper{
//    typedef Quotient<NT> Root_of_1;
    typedef typename Get_rational_type<NT>::type Root_of_1;
    typedef CGAL::Sqrt_extension<Root_of_1,Root_of_1,::CGAL::Tag_true,::CGAL::Tag_true>             Root_of_2;
//    typedef CGAL::Root_of_2<NT> Root_of_2;
    struct Make_root_of_2{
        typedef Root_of_2 result_type;
        Root_of_2 operator()(const NT& a, const NT& b, const NT& c) const {
            return Root_of_2(a,b,c);
        }
        Root_of_2 operator()(const NT& a, const NT& b, const NT& c, bool s) const {
            return Root_of_2(a,b,c,s);
        }
        Root_of_2 operator()(const Root_of_1& a,
                             const Root_of_1& b,
                             const Root_of_1& c) const {
            return Root_of_2(a,b,c);
        }
        Root_of_2 operator()(const Root_of_1& a,
                             const Root_of_1& b,
                             const Root_of_1& c,
                             bool s) const {
            return Root_of_2(a,b,c,s);
        }
    };
  
private:
  typedef CGAL::Algebraic_structure_traits<Root_of_2> AST;
public:
  typedef typename AST::Square  Square; 
  typedef typename AST::Inverse Inverse;
  
  struct Make_sqrt {
    typedef Root_of_2 result_type;
    Root_of_2 operator()(const NT& x) const {
      return Root_of_2(x,true);
    }
  };
};

template < typename FT>
struct Root_of_traits_helper < FT, Field_tag >
{
    typedef FT               Root_of_1;
private:
    typedef Fraction_traits<FT> FrT;
    // Field must be a Type (Decomposable)
    // We have the typedef as VC10 fails with 
    // static_assert(FrT::Is_fraction::value)
    typedef typename FrT::Is_fraction ISF;
    CGAL_static_assertion((ISF::value));


    typedef typename FrT::Numerator_type      RT;
    typedef typename FrT::Decompose Decompose;
public:
    typedef CGAL::Sqrt_extension<Root_of_1,Root_of_1,::CGAL::Tag_true,::CGAL::Tag_true>             Root_of_2;

    struct Make_root_of_2{
        typedef Root_of_2 result_type;
        Root_of_2
        operator()(const FT& a, const FT& b, const FT& c) const {
            return Root_of_2(a,b,c);
        }
        Root_of_2
        operator()(const FT& a, const FT& b, const FT& c, bool smaller) const {
            Decompose decompose;
            RT a_num,b_num,c_num;
            RT a_den,b_den,c_den; // Denomiantor same as RT

            decompose(a,a_num,a_den);
            decompose(b,b_num,b_den);
            decompose(c,c_num,c_den);

            RT a_ = a_num * b_den * c_den;
            RT b_ = b_num * a_den * c_den;
            RT c_ = c_num * a_den * b_den;

            return make_root_of_2<RT>(a_,b_,c_,smaller);
        } 
    };

private:
  typedef CGAL::Algebraic_structure_traits<Root_of_2> AST;
public:
  typedef typename AST::Square  Square; 
  typedef typename AST::Inverse Inverse;
  
  struct Make_sqrt{
    typedef Root_of_2 result_type;
    Root_of_2 operator()(const FT& x) const {
      return Root_of_2( FT(0),FT(1),x);
    }
  };
};

template < typename NT >
struct Root_of_traits_helper < NT, Field_with_sqrt_tag >
{
    typedef NT  Root_of_1;
    typedef NT  Root_of_2;

    struct Make_root_of_2{
        typedef NT result_type;
        // just a copy, not sure about the semantic of smaller
        NT operator()(const NT& a, const NT& b, const NT& c, bool smaller) const {
            // former make_root_of_2_sqrt()
            CGAL_assertion( a != 0 );
            NT discriminant = CGAL_NTS square(b) - a*c*4;
            CGAL_assertion( discriminant >= 0 );
            NT d = CGAL_NTS sqrt(discriminant);
            if ((smaller && a>0) || (!smaller && a<0))
                d = -d;
            return (d-b)/(a*2);
        }
        // it's so easy
        NT operator()(const NT& a, const NT& b, const NT& c) const {
            return a + b * CGAL_NTS sqrt(c) ;
        }
    };

private:
  typedef CGAL::Algebraic_structure_traits<Root_of_2> AST;
public:
  typedef typename AST::Square  Square; 
  typedef typename AST::Inverse Inverse;
  
  struct Make_sqrt{
    typedef Root_of_2 result_type;
    Root_of_2 operator()(const NT& x) const {
      return CGAL::sqrt(x);
    }
  };
};

template < typename NT >
struct Root_of_traits_helper < NT, Field_with_kth_root_tag >
    :public Root_of_traits_helper < NT, Field_with_sqrt_tag>{};

template < typename NT >
struct Root_of_traits_helper < NT, Field_with_root_of_tag >
    :public Root_of_traits_helper < NT, Field_with_sqrt_tag>{};


} // namespace internal



// Default Traits class for NT types
template < typename NT >
struct Root_of_traits
    : public internal::Root_of_traits_helper<NT,
      typename Algebraic_structure_traits<NT>::Algebraic_category> {
    typedef internal::Root_of_traits_helper<NT,
      typename Algebraic_structure_traits<NT>::Algebraic_category> Base;
    typedef typename Base::Root_of_1 RootOf_1;
    typedef typename Base::Root_of_2 RootOf_2;
};

template <bool B>
struct Root_of_traits<Interval_nt<B> >{
  typedef Interval_nt<B> Root_of_1;
  typedef Interval_nt<B> Root_of_2;
  typedef Root_of_1 RootOf_1;
  typedef Root_of_2 RootOf_2;
  struct Make_root_of_2{
    typedef Interval_nt<B> result_type;
    // just a copy, not sure about the semantic of smaller
    Interval_nt<B> operator()(const Interval_nt<B>& a, const Interval_nt<B>& b, const Interval_nt<B>& c, bool smaller) const {
        // former make_root_of_2_sqrt()
        if (CGAL::possibly(a==0))
          return Interval_nt<B>::largest();
        Interval_nt<B> discriminant = CGAL_NTS square(b) - a*c*4;
        CGAL_assertion(discriminant >= 0);
        Interval_nt<B> d = CGAL_NTS sqrt(discriminant);
        if ((smaller && a>0) || (!smaller && a<0))
            d = -d;
        return (d-b)/(a*2);
    }
    // it's so easy
    Interval_nt<B> operator()(const Interval_nt<B>& a, const Interval_nt<B>& b, const Interval_nt<B>& c) const {
        return a + b * CGAL_NTS sqrt(c) ;
    }
  };
  
private:
  typedef CGAL::Algebraic_structure_traits<Interval_nt<B> > AST;
public:
  typedef typename AST::Square  Square; 
  typedef typename AST::Inverse Inverse;
  typedef typename AST::Sqrt    Make_sqrt;
  
};


} //namespace CGAL

#include <CGAL/Root_of_traits_specializations.h>

#endif // CGAL_ROOT_OF_TRAITS_H
