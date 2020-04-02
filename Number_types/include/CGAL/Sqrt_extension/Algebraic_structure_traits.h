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
// Author(s)     : Michael Hemmer   <hemmer@mpi-inf.mpg.de>


#ifndef CGAL_SQRT_EXTENSION_ALGEBRAIC_STRUCTURE_TRAITS_H
#define CGAL_SQRT_EXTENSION_ALGEBRAIC_STRUCTURE_TRAITS_H

#include <CGAL/basic.h>

namespace CGAL {

namespace internal {

// Algebraic structure traits
template< class Type, class Algebraic_type >
class Sqrt_extension_algebraic_structure_traits_base;

template< class Type >
class Sqrt_extension_algebraic_structure_traits_base< Type,
                                        CGAL::Integral_domain_without_division_tag >
  : public Algebraic_structure_traits_base< Type,
                                      CGAL::Integral_domain_without_division_tag > {
  public:
    typedef CGAL::Integral_domain_without_division_tag Algebraic_category;

    class Simplify
      : public CGAL::cpp98::unary_function< Type&, void > {
      public:
        typedef void result_type;
        typedef Type& argument_type;

        void operator()( Type& x ) const {
          x.simplify();
        }
    };
};

template< class Type >
class Sqrt_extension_algebraic_structure_traits_base< Type,
                                                    CGAL::Integral_domain_tag >
    : public Sqrt_extension_algebraic_structure_traits_base< Type,
                                      CGAL::Integral_domain_without_division_tag > {
public:
    typedef CGAL::Integral_domain_tag Algebraic_category;

    class Integral_division
      : public CGAL::cpp98::binary_function< Type, Type, Type > {
    public:
      Type operator()( const Type& x,const Type& y ) const {
        return x/y;
      }
      CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR( Type )
    };

private:
  typedef typename Type::NT COEFF;
  typedef typename Type::ROOT ROOT;
  typedef typename CGAL::Coercion_traits< ROOT, COEFF >::Cast Root_nt_cast;
  typedef CGAL::Algebraic_structure_traits<COEFF> AST_COEFF;
  typedef typename AST_COEFF::Divides Divides_coeff;
public:
  class Divides
    : public CGAL::cpp98::binary_function<Type,Type,typename Divides_coeff::result_type>{
    typedef typename Divides_coeff::result_type BOOL;
  public:
    BOOL operator()( const Type& x, const Type& y) const {
      Type q;
      return (*this)(x,y,q);
    }

    BOOL operator()( const Type& x, const Type& y, Type& q) const {
      Divides_coeff divides;

//            std::cout<<"integral domain for sqrt"<<std::endl;
      BOOL result;
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
    CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR_WITH_RT(Type,BOOL)
  };
};

template< class Type >
class Sqrt_extension_algebraic_structure_traits_base< Type,
                                                    CGAL::Unique_factorization_domain_tag >
  : public Sqrt_extension_algebraic_structure_traits_base< Type,
                                      CGAL::Integral_domain_tag > {
  // Nothing new
};

template< class Type >
class Sqrt_extension_algebraic_structure_traits_base< Type,
                                                    CGAL::Euclidean_ring_tag >
  : public Sqrt_extension_algebraic_structure_traits_base< Type,
                                      CGAL::Integral_domain_tag > {
  // Nothing new
};

template< class Type >
class Sqrt_extension_algebraic_structure_traits_base< Type,
                                                    CGAL::Field_tag >
  : public Sqrt_extension_algebraic_structure_traits_base< Type,
                                      CGAL::Integral_domain_tag > {
  public:
    typedef Field_tag Algebraic_category;

    class Unit_part
      : public CGAL::cpp98::unary_function< Type, Type > {
      public:
        Type operator()( const Type& x ) const {
          return( x == Type(0) ? Type(1) : x );
        }
    };
  class Inverse
    : public CGAL::cpp98::unary_function< Type, Type > {
  public:
    Type operator()( const Type& x ) const {
      return Type(1)/x ;
    }
  };
};

template< class Type >
class Sqrt_extension_algebraic_structure_traits_base< Type,
                                                    CGAL::Field_with_sqrt_tag >
  : public Sqrt_extension_algebraic_structure_traits_base< Type,
                                      CGAL::Field_tag > {
  // Nothing new
};

template< class Type >
class Sqrt_extension_algebraic_structure_traits_base< Type,
                                                CGAL::Field_with_kth_root_tag >
  : public Sqrt_extension_algebraic_structure_traits_base< Type,
                                      // TODO: Why not Fiel_tag?
                                      CGAL::Field_with_sqrt_tag > {
  // Nothing new
};

template< class Type >
class Sqrt_extension_algebraic_structure_traits_base< Type,
                                                CGAL::Field_with_root_of_tag >
  : public Sqrt_extension_algebraic_structure_traits_base< Type,
                                      // TODO: Why not Fiel_tag?
                                      CGAL::Field_with_sqrt_tag > {
  // Nothing new
};

} // namespace internal


template< class COEFF_, class ROOT_, class ACDE_TAG,class FP_TAG>
class Algebraic_structure_traits< Sqrt_extension< COEFF_, ROOT_, ACDE_TAG,FP_TAG > >
    : public internal::Sqrt_extension_algebraic_structure_traits_base<
      Sqrt_extension< COEFF_, ROOT_, ACDE_TAG,FP_TAG >,
      typename Algebraic_structure_traits< COEFF_ >::Algebraic_category > {
public:
  typedef Sqrt_extension< COEFF_, ROOT_, ACDE_TAG,FP_TAG > Type;

    // Tag_true if COEFF and ROOT are exact
    typedef typename ::boost::mpl::if_c<
       bool( ::boost::is_same<typename CGAL::Algebraic_structure_traits<ROOT_ >::Is_exact,::CGAL::Tag_true>::value )&&
       bool( ::boost::is_same<typename CGAL::Algebraic_structure_traits<COEFF_>::Is_exact,::CGAL::Tag_true>::value )
           ,::CGAL::Tag_true,::CGAL::Tag_false>::type Is_exact;

    typedef typename Algebraic_structure_traits<COEFF_>::Is_numerical_sensitive
    Is_numerical_sensitive;
};

} //namespace CGAL

#endif
