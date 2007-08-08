// Copyright (c) 2006-2007 Max-Planck-Institute Saarbruecken (Germany).
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


#ifndef CGAL_REAL_EMBEDDABLE_TRAITS_H
#define CGAL_REAL_EMBEDDABLE_TRAITS_H

#include <CGAL/number_type_basic.h>

CGAL_BEGIN_NAMESPACE


template< class Type_ , 
          typename Is_real_embeddable_ = Tag_false > 
class Real_embeddable_traits {
  public:
    typedef Type_  Type;
    typedef Tag_false   Is_real_embeddable;

    typedef Null_functor Abs;
    typedef Null_functor Sign;
    typedef Null_functor Is_finite;
    typedef Null_functor Is_positive;
    typedef Null_functor Is_negative;
    typedef Null_functor Is_zero;
    typedef Null_functor Compare;
    typedef Null_functor To_double;
    typedef Null_functor To_interval;
};

namespace INTERN_RET {
    template< class Type, class AST_is_zero >
    class Is_zero_selector
      : public Unary_function< Type, bool > {
      public:        
        //! the function call.
        bool operator()( const Type& x ) const {
            return AST_is_zero()(x);
        }
    };
    
    template< class Type >
    class Is_zero_selector< Type, Null_functor >
      : public Unary_function< Type, bool > {
      public:        
        //! the function call.
        bool operator()( const Type& x ) const {
            return x == Type(0);
        }
    };
    
} // INTERN_RET

template< class Type_ >
class Real_embeddable_traits_base {
  public:
    typedef Type_  Type;
    typedef Tag_true    Is_real_embeddable; 
        
    //! The generic \c Is_zero functor implementation uses one comparison
    typedef INTERN_RET::Is_zero_selector< Type, 
                typename Algebraic_structure_traits< Type >::Is_zero 
                                        > Is_zero;
    
    //! The generic \c Is_finite functor returns true
    class Is_finite
      : public Unary_function< Type, bool > {
      public:
        bool operator()( const Type& ) const {
          return true;
        }
    };
    //! The generic \c Abs functor implementation
    //! uses one comparisons and the unary minus if necessary.
    class Abs
      : public Unary_function< Type, Type > {
      public:
        //! the function call.
        Type  operator()( const Type& x ) const {
          return( x < Type(0) ) ? -x : x;
        }
    };
    
    //! The generic \c Sign functor implementation uses two comparisons.
    class Sign 
      : public Unary_function< Type, ::CGAL::Sign > {
      public:
        //! the function call.
        ::CGAL::Sign operator()( const Type& x ) const {
          if ( x < Type(0))
            return NEGATIVE;
          if ( x > Type(0))
            return POSITIVE;
          return ZERO;
        }
    };
    
    //! The generic \c Is_positive functor implementation uses one comparison.
    class Is_positive 
      : public Unary_function< Type, bool > {
      public:        
        //! the function call.
        bool operator()( const Type& x ) const {
          return x > Type(0);
        }
    };
    
    //! The generic \c Is_negative functor implementation uses one comparison.
    class Is_negative 
      : public Unary_function< Type, bool > {
      public:        
        //! the function call.
        bool operator()( const Type& x ) const {
          return x < Type(0);
        }
    };
        
    //! The generic \c Compare functor implementation uses two comparisons.
    class Compare 
      : public Binary_function< Type, Type, 
                                Comparison_result > {
      public:
        //! the function call.
        Comparison_result operator()( const Type& x, 
                                            const Type& y) const {
          if( x < y )
            return SMALLER;
          if( x > y )
            return LARGER;
         return EQUAL;
        }
        
        CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR_WITH_RT( Type,
                                                      Comparison_result )
    };

    typedef Null_functor To_double;
    typedef Null_functor To_interval;
};

// Some common functors to be used by RET specializations
namespace INTERN_RET {
  template< class Type, class Is_real_embeddable >
  class Real_embeddable_traits_base_selector;
  
  template< class Type >
  class Real_embeddable_traits_base_selector< Type, Tag_false >
    : public Real_embeddable_traits< Null_tag > {};
    
  template< class Type >
  class Real_embeddable_traits_base_selector< Type, Tag_true >
    : public Real_embeddable_traits_base< Type > {};

  template< class Type >
  class To_double_by_conversion 
    : public Unary_function< Type, double > {
    public:      
      //! the function call.
      double operator()( const Type& x ) const {
        return static_cast<double>(x);
      }
  };
  
  template< class Type >
  class To_interval_by_conversion 
    : public Unary_function< Type, std::pair< double, double > > {
    public:      
        //! the function call.
      std::pair<double, double> operator()( const Type& x ) const {

        return std::pair<double,double>( x, x );
      }
  };  
} // INTERN_RET

CGAL_END_NAMESPACE

#endif  // CGAL_REAL_EMBEDDABLE_TRAITS_H
