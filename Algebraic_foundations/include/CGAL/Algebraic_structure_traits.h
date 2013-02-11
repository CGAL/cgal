// Copyright (c) 2006-2007 Max-Planck-Institute Saarbruecken (Germany).
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
// Author(s)     : Michael Hemmer    <hemmer@mpi-inf.mpg.de>
//
// =============================================================================


#ifndef CGAL_ALGEBRAIC_STRUCTURE_TRAITS_H
#define CGAL_ALGEBRAIC_STRUCTURE_TRAITS_H

#include <functional>
#include <CGAL/tags.h>
#include <CGAL/type_traits.h>
#include <CGAL/Coercion_traits.h>
#include <CGAL/assertions.h>
#include <CGAL/use.h>

namespace CGAL {

// REMARK: Some of the following comments and references are just copy & pasted
//         from EXACUS and have to be adapted/removed in the future.

// The tags for Algebra_type corresponding to the number type concepts
// ===================================================================

//! corresponds to the \c IntegralDomainWithoutDiv concept.
struct Integral_domain_without_division_tag {};

//! corresponds to the \c IntegralDomain concept.
struct Integral_domain_tag : public Integral_domain_without_division_tag {};

//! corresponds to the \c UFDomain concept.
struct Unique_factorization_domain_tag : public Integral_domain_tag {};

//! corresponds to the \c EuclideanRing concept.
struct Euclidean_ring_tag : public Unique_factorization_domain_tag {};

//! corresponds to the \c Field concept.
struct Field_tag : public Integral_domain_tag {};

//! corresponds to the \c FieldWithSqrt concept.
struct Field_with_sqrt_tag : public Field_tag {};

//! corresponds to the \c FieldWithKthRoot concept
struct Field_with_kth_root_tag : public Field_with_sqrt_tag {};

//! corresponds to the \c FieldWithRootOF concept.
struct Field_with_root_of_tag : public Field_with_kth_root_tag {};


// The algebraic structure traits template
// =========================================================================
template< class Type_ > 
class Algebraic_structure_traits  {
  public:
    typedef Type_  Type;
    typedef Null_tag       Algebraic_category;
    typedef Null_tag       Is_exact;
    typedef Null_tag       Is_numerical_sensitive;

    typedef Null_functor Simplify;
    typedef Null_functor Unit_part;
    typedef Null_functor Integral_division;
    typedef Null_functor Is_square;    
    typedef Null_functor Gcd;
    typedef Null_functor Div_mod;
    typedef Null_functor Div;
    typedef Null_functor Mod;
    typedef Null_functor Square;
    typedef Null_functor Is_zero;
    typedef Null_functor Is_one;    
    typedef Null_functor Sqrt;
    typedef Null_functor Kth_root;
    typedef Null_functor Root_of; 
    typedef Null_functor Divides; 
    typedef Null_functor Inverse; 
};

// The algebraic structure traits base class
// =========================================================================
template< class Type, class Algebra_type >
class Algebraic_structure_traits_base;

//! The template specialization that can be used for types that are not any
//! of the number type concepts. All functors are set to \c Null_functor
//! or suitable defaults. The \c Simplify functor does nothing by default.
template< class Type_ >
class Algebraic_structure_traits_base< Type_, Null_tag > {
  public:
    typedef Type_  Type;
    typedef Null_tag       Algebraic_category;
    typedef Tag_false      Is_exact;
    typedef Null_tag       Is_numerical_sensitive;
    typedef Null_tag       Boolean; 

    // does nothing by default
    class Simplify 
      : public std::unary_function< Type&, void > {
      public:
        void operator()( Type& ) const {}
    };

    typedef Null_functor Unit_part;
    typedef Null_functor Integral_division;
    typedef Null_functor Is_square;    
    typedef Null_functor Gcd;
    typedef Null_functor Div_mod;
    typedef Null_functor Div;
    typedef Null_functor Mod;
    typedef Null_functor Square;
    typedef Null_functor Is_zero;
    typedef Null_functor Is_one;    
    typedef Null_functor Sqrt;
    typedef Null_functor Kth_root;
    typedef Null_functor Root_of; 
    typedef Null_functor Divides;
    typedef Null_functor Inverse;
};

//! The template specialization that is used if the number type is
//! a model of the \c IntegralDomainWithoutDiv concept. The \c Simplify
//! does nothing by default and the \c Unit_part is equal to
//! \c Type(-1) for negative numbers and 
//! \c Type(1) otherwise
template< class Type_ >
class Algebraic_structure_traits_base< Type_, 
                                       Integral_domain_without_division_tag > 
    : public Algebraic_structure_traits_base< Type_, 
                                              Null_tag > {
  public:
    typedef Type_                                 Type;
    typedef Integral_domain_without_division_tag  Algebraic_category;
    typedef bool                                  Boolean;

    // returns Type(1) by default
    class Unit_part 
      : public std::unary_function< Type, Type > { 
      public:
        Type operator()( const Type& x ) const {
          return( x < Type(0)) ? 
                  Type(-1) : Type(1); 
        }
    };
    
    class Square 
      : public std::unary_function< Type, Type > {
      public:        
        Type operator()( const Type& x ) const {
          return x*x;
        }
    };
    
    class Is_zero 
      : public std::unary_function< Type, bool > {
      public:        
        bool operator()( const Type& x ) const {
          return x == Type(0);
        }
    };

    class Is_one 
      : public std::unary_function< Type, bool > {
      public:        
        bool operator()( const Type& x ) const {
          return x == Type(1);
        }
    };

};


//! The template specialization that is used if the number type is
//! a model of the \c IntegralDomain concept. It is equivalent to the 
//! specialization
//! for the \c IntegralDomainWithoutDiv concept. The additionally required 
//! \c Integral_division functor needs to be implemented in the 
//! \c Algebraic_structure_traits itself.
template< class Type_ >
class Algebraic_structure_traits_base< Type_, 
                                       Integral_domain_tag >
    : public Algebraic_structure_traits_base< Type_, 
                                       Integral_domain_without_division_tag > {
  public:
    typedef Type_       Type;
    typedef Integral_domain_tag  Algebraic_category;
};


//! The template specialization that is used if the number type is
//! a model of the \c UFDomain concept. It is equivalent to the specialization
//! for the \c IntegralDomain concept. The additionally required 
//! \c Integral_div functor
//! and \c Gcd functor need to be implemented in the 
//! \c Algebraic_structure_traits itself.
template< class Type_ >
class Algebraic_structure_traits_base< Type_, 
                                       Unique_factorization_domain_tag >
    : public Algebraic_structure_traits_base< Type_, 
                                              Integral_domain_tag > {
  public:
    typedef Type_  Type;
    typedef Unique_factorization_domain_tag    Algebraic_category;

  // Default implementation of Divides functor for unique factorization domains
  // x divides y if gcd(y,x) equals x up to inverses 
  class Divides 
    : public std::binary_function<Type,Type,bool>{ 
  public:
    bool operator()( const Type& x,  const Type& y) const {  
      typedef CGAL::Algebraic_structure_traits<Type> AST;
      typename AST::Gcd gcd;
      typename AST::Unit_part unit_part;
      typename AST::Integral_division idiv;
      return gcd(y,x) == idiv(x,unit_part(x));
    }
    // second operator computing q = x/y 
    bool operator()( const Type& x,  const Type& y, Type& q) const {    
      typedef CGAL::Algebraic_structure_traits<Type> AST;
      typename AST::Integral_division idiv;
      bool result = (*this)(x,y);
      if( result == true ) 
        q = idiv(x,y);
      return result; 
    }
    CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR_WITH_RT(Type,bool)
  };
};


//! The template specialization that is used if the number type is
//! a model of the \c EuclideanRing concept.
template< class Type_ >
class Algebraic_structure_traits_base< Type_, 
                                       Euclidean_ring_tag >
    : public Algebraic_structure_traits_base< Type_, 
                                              Unique_factorization_domain_tag > {
  public:
    typedef Type_        Type;
    typedef Euclidean_ring_tag    Algebraic_category;

    // maps to \c Div by default.
    class Integral_division 
      : public std::binary_function< Type, Type,
                                Type > { 
      public:
        Type operator()( 
                const Type& x, 
                const Type& y) const { 
            typedef Algebraic_structure_traits<Type> AST; 
            typedef typename AST::Is_exact Is_exact;
            CGAL_USE_TYPE(Is_exact);
            typename AST::Div actual_div;
            
            CGAL_precondition_msg( 
                    !Is_exact::value || actual_div( x, y) * y == x,
                    "'x' must be divisible by 'y' in "
                    "Algebraic_structure_traits<...>::Integral_div()(x,y)" );
            return actual_div( x, y);          
        }
        CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR( Type )
    };

    // Algorithm from NiX/euclids_algorithm.h
    class Gcd 
      : public std::binary_function< Type, Type,
                                Type > { 
      public:
        Type operator()( 
                const Type& x, 
                const Type& y) const {
            typedef Algebraic_structure_traits<Type> AST;
            typename AST::Mod mod;
            typename AST::Unit_part unit_part;
            typename AST::Integral_division integral_div;
            // First: the extreme cases and negative sign corrections.
            if (x == Type(0)) {
                if (y == Type(0))  
                    return Type(0);
                return integral_div( y, unit_part(y) );
            }
            if (y == Type(0))
                return integral_div(x, unit_part(x) );
            Type u = integral_div( x, unit_part(x) );
            Type v = integral_div( y, unit_part(y) );
            // Second: assuming mod is the most expensive op here, 
            // we don't compute it unnecessarily if u < v
            if (u < v) {
                v = mod(v,u);
                // maintain invariant of v > 0 for the loop below
                if ( v == Type(0) )
                    return u;
            }
            // Third: generic case of two positive integer values and u >= v.
            // The standard loop would be:
            //      while ( v != 0) {
            //          int tmp = mod(u,v);
            //          u = v;
            //          v = tmp;
            //      }
            //      return u;
            //
            // But we want to save us all the variable assignments and unroll
            // the loop. Before that, we transform it into a do {...} while()
            // loop to reduce branching statements.
            Type w;
            do {
                w = mod(u,v);
                if ( w == Type(0))
                    return v;
                u = mod(v,w);
                if ( u == Type(0))
                    return w;
                v = mod(w,u);
            } while (v != Type(0));
            return u;
        }  
        CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR( Type )
    };

    // based on \c Div and \c Mod.
    class Div_mod { 
    public:
        typedef Type    first_argument_type;
        typedef Type    second_argument_type;
        typedef Type&   third_argument_type;
        typedef Type&   fourth_argument_type;
        typedef void  result_type;
        void operator()( const Type& x, 
                const Type& y, 
                Type& q, Type& r) const {
            typedef Algebraic_structure_traits<Type> Traits;
            typename Traits::Div  actual_div;
            typename Traits::Mod  actual_mod;
            q = actual_div( x, y );
            r = actual_mod( x, y );          
            return;
        }
        
        template < class NT1, class NT2 >
        void operator()( 
                const NT1& x, 
                const NT2& y,
                Type& q, 
                Type& r ) const {
            typedef Coercion_traits< NT1, NT2 > CT;
            typedef typename CT::Type Type; 
            CGAL_USE_TYPE(Type);
            CGAL_static_assertion(( 
              ::boost::is_same<Type , Type >::value));
            
            typename Coercion_traits< NT1, NT2 >::Cast cast;
            operator()( cast(x), cast(y), q, r );          
        }
    };
    
    // based on \c Div_mod.
    class Div 
      : public std::binary_function< Type, Type,
                                Type > {
      public:
        Type operator()( const Type& x, 
                                        const Type& y) const {
          typename Algebraic_structure_traits<Type>
                                                    ::Div_mod actual_div_mod;
          Type q;     
          Type r;
          actual_div_mod( x, y, q, r );
          return q;
        };
        
        CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR( Type )
    };

    // based on \c Div_mod.
    class Mod 
      : public std::binary_function< Type, Type,
                                Type > { 
      public:
        Type operator()( const Type& x, 
                                        const Type& y) const {
          typename Algebraic_structure_traits<Type>
                                                    ::Div_mod actual_div_mod;
          Type q;     
          Type r;
          actual_div_mod( x, y, q, r );
          return r;
        };
        
        CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR( Type )
    };

  // Divides for Euclidean Ring 
  class Divides 
    : public std::binary_function<Type, Type, bool>{
  public:
    bool operator()( const Type& x, const Type& y) const {
      typedef Algebraic_structure_traits<Type> AST;
      typename AST::Mod mod;
      CGAL_precondition(typename AST::Is_zero()(x) == false );
      return typename AST::Is_zero()(mod(y,x));
    }
    // second operator computing q 
    bool operator()( const Type& x, const Type& y, Type& q) const {
      typedef Algebraic_structure_traits<Type> AST;
      typename AST::Div_mod div_mod;
      CGAL_precondition(typename AST::Is_zero()(x) == false );
      Type r;
      div_mod(y,x,q,r);
      return (typename AST::Is_zero()(r));
    }
    CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR_WITH_RT(Type,bool)
  };
};


//! The template specialization that is used if the number type is
//! a model of the \c Field concept. \c Unit_part ()(x)
//! returns \c NT(1) if the value \c x is equal to \c NT(0) and
//! otherwise the value \c x itself. The \c Integral_div
//! maps to the \c operator/.
//! See also \link NiX_NT_traits_functors concept NT_traits \endlink .
//! \ingroup NiX_NT_traits_bases
//
template< class Type_ >
class Algebraic_structure_traits_base< Type_, Field_tag >
    : public Algebraic_structure_traits_base< Type_, 
                                              Integral_domain_tag > {
  public:
    typedef Type_        Type;
    typedef Field_tag             Algebraic_category;

    // returns the argument \a a by default
    class Unit_part 
      : public std::unary_function< Type, Type > { 
      public:
        Type operator()( const Type& x ) const {
            return( x == Type(0)) ? Type(1) : x;
        }
    };
    // maps to \c operator/ by default.
    class Integral_division 
      : public std::binary_function< Type, Type,
                                Type > { 
      public:
        Type operator()( const Type& x, 
                                        const Type& y) const { 
            typedef Algebraic_structure_traits<Type> AST; 
            typedef typename AST::Is_exact Is_exact;
            CGAL_USE_TYPE(Is_exact);
	    CGAL_precondition_code( bool ie = Is_exact::value; )
            CGAL_precondition_msg( !ie || (x / y) * y  == x,
                    "'x' must be divisible by 'y' in "
                    "Algebraic_structure_traits<...>::Integral_div()(x,y)" );
            return x / y;
        }
      CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR( Type )
    };
  
  // maps to \c 1/x by default.
  class Inverse 
    : public std::unary_function< Type, Type > { 
  public:
    Type operator()( const Type& x ) const { 
      return Type(1)/x;
    }
  };
  

  // Default implementation of Divides functor for Field: 
  // returns always true
  // \pre: x != 0
  class Divides
    : public std::binary_function< Type, Type, bool > { 
  public:
    bool operator()( const Type& CGAL_precondition_code(x), const Type& /* y */) const {
      CGAL_precondition_code( typedef Algebraic_structure_traits<Type> AST);
      CGAL_precondition( typename AST::Is_zero()(x) == false );
      return true;
    } 
    // second operator computing q
    bool operator()( const Type& x, const Type& y, Type& q) const {
      CGAL_precondition_code(typedef Algebraic_structure_traits<Type> AST);
      CGAL_precondition( typename AST::Is_zero()(x) == false );
      q = y/x;
      return true;
    }
    CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR_WITH_RT(Type,bool)
  };
};


//! The template specialization that is used if the number type is a model
//! of the \c FieldWithSqrt concept. It is equivalent to the 
//! specialization for the \c Field concept. The additionally required 
//! \c NiX::NT_traits::Sqrt functor need to be
//! implemented in the \c NT_traits itself.
//! \ingroup NiX_NT_traits_bases
//
template< class Type_ >
class Algebraic_structure_traits_base< Type_, 
                                       Field_with_sqrt_tag>
    : public Algebraic_structure_traits_base< Type_, 
                                              Field_tag> {
  public:
    typedef Type_        Type;
    typedef Field_with_sqrt_tag   Algebraic_category;

    struct Is_square
        :public std::binary_function<Type,Type&,bool>
    {
        bool operator()(const Type& ) const {return true;}
        bool operator()(
                const Type& x,
                Type      & result) const {
            typename Algebraic_structure_traits<Type>::Sqrt sqrt;
            result = sqrt(x);
            return true;
        }
    };
};

//! The template specialization that is used if the number type is a model
//! of the \c FieldWithKthRoot concept. It is equivalent to the 
//! specialization for the \c Field concept. The additionally required 
//! \c NiX::NT_traits::Kth_root functor need to be
//! implemented in the \c Algebraic_structure_traits itself.
//! \ingroup NiX_NT_traits_bases
//
template< class Type_ >
class Algebraic_structure_traits_base< Type_, 
                                       Field_with_kth_root_tag>
    : public Algebraic_structure_traits_base< Type_, 
                                              Field_with_sqrt_tag> {
    
    
    
  public:
    typedef Type_        Type;
    typedef Field_with_kth_root_tag   Algebraic_category;
};




//! The template specialization that is used if the number type is a model
//! of the \c FieldWithRootOf concept. It is equivalent to the 
//! specialization for the \c FieldWithKthRoot concept. The additionally 
//! required \c NiX::NT_traits::Root_of functor need to be
//! implemented in the \c NT_traits itself.
//! \ingroup NiX_NT_traits_bases
//
template< class Type_ >
class Algebraic_structure_traits_base< Type_, 
                                       Field_with_root_of_tag >
    : public Algebraic_structure_traits_base< Type_, 
                                              Field_with_kth_root_tag > {
  public:
    typedef Type_           Type;
    typedef Field_with_root_of_tag   Algebraic_category;
};

// Some common functors to be used by AST specializations
namespace INTERN_AST {
  template< class Type >
  class Div_per_operator 
    : public std::binary_function< Type, Type, 
                              Type > {
    public:      
      Type operator()( const Type& x, 
                                      const Type& y ) const {
        return x / y;
      }
      
      CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR( Type )
  };
  
  template< class Type >
  class Mod_per_operator 
    : public std::binary_function< Type, Type,
                              Type > {
    public:
      Type operator()( const Type& x, 
                                      const Type& y ) const {
        return x % y;
      }
      
      CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR( Type )
  };
  
  template< class Type >
  class Is_square_per_sqrt
    : public std::binary_function< Type, Type&,
                              bool > {
    public:      
      bool operator()( const Type& x, 
                       Type& y ) const {
          typename Algebraic_structure_traits< Type >::Sqrt
              actual_sqrt;
          y = actual_sqrt( x );
          return y * y == x;
      }
      bool operator()( const Type& x) const {
          Type dummy;
          return operator()(x,dummy);
      }
  };
} // INTERN_AST
} //namespace CGAL

#endif  // CGAL_ALGEBRAIC_STRUCTURE_TRAITS_H
