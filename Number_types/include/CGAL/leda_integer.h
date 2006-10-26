// Copyright (c) 1999  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
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
// Author(s)     : Andreas Fabri, Michael Hemmer
 
#ifndef CGAL_LEDA_INTEGER_H
#define CGAL_LEDA_INTEGER_H

#include <CGAL/basic.h>
#include <CGAL/leda_coercion_traits.h>
#include <CGAL/Number_type_traits.h>
#include <CGAL/Interval_nt.h>
#include <CGAL/Algebraic_structure_traits.h>
#include <CGAL/Real_embeddable_traits.h>
#include <CGAL/utils.h>
#include <CGAL/Needs_parens_as_product.h>

#include <CGAL/functional_base.h> // Unary_function, Binary_function

#include <utility>

#include <CGAL/LEDA_basic.h>
#include <CGAL/LEDA_basic.h>
#if CGAL_LEDA_VERSION < 500
#include <LEDA/integer.h>
#include <LEDA/bigfloat.h>// for To_interval
#else
#include <LEDA/numbers/integer.h>
#include <LEDA/numbers/bigfloat.h>// for To_interval
#endif

CGAL_BEGIN_NAMESPACE

template <> struct Number_type_traits<leda_integer> {
  typedef Tag_true  Has_gcd;
  typedef Tag_true  Has_division;
  typedef Tag_true  Has_sqrt;

  typedef Tag_true  Has_exact_ring_operations;
  typedef Tag_false Has_exact_division;
  typedef Tag_false Has_exact_sqrt;
};

template <> class Algebraic_structure_traits< leda_integer >
  : public Algebraic_structure_traits_base< leda_integer, 
                                            CGAL::Euclidean_ring_tag >  {
  public:
    typedef CGAL::Tag_true            Is_exact;
                
    typedef CGAL::INTERN_AST::Is_square_per_sqrt< Algebraic_structure >
                                                                 Is_square;
                                                                 
    class Gcd 
      : public Binary_function< Algebraic_structure, Algebraic_structure,
                                Algebraic_structure > {
      public:
        Algebraic_structure operator()( const Algebraic_structure& x, 
                                        const Algebraic_structure& y ) const {
          // By definition gcd(0,0) == 0
          if( x == Algebraic_structure(0) && y == Algebraic_structure(0) )
            return Algebraic_structure(0);
            
          return CGAL_LEDA_SCOPE::gcd( x, y );
        }
        
        CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR( Algebraic_structure )
    };
    
    typedef CGAL::INTERN_AST::Div_per_operator< Algebraic_structure > Div;
    
    class Mod 
      : public Binary_function< Algebraic_structure, Algebraic_structure,
                                Algebraic_structure > {
      public:
        Algebraic_structure operator()( const Algebraic_structure& x, 
                                        const Algebraic_structure& y ) const {
          Algebraic_structure m = x % y;
          
          // Fix wrong lede result if first operand is negative
          if( x < 0 && m != 0 )
            m -= y;

          return m;
        }
        
        CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR( Algebraic_structure )
    };
    
    class Sqrt 
      : public Unary_function< Algebraic_structure, Algebraic_structure > {
      public:
        Algebraic_structure operator()( const Algebraic_structure& x ) const {
          return CGAL_LEDA_SCOPE::sqrt( x );
        }
    };        
};

template <> class Real_embeddable_traits< leda_integer > 
  : public Real_embeddable_traits_base< leda_integer > {
  public:
      
    class Abs 
      : public Unary_function< Real_embeddable, Real_embeddable > {
      public:
        Real_embeddable operator()( const Real_embeddable& x ) const {
            return CGAL_LEDA_SCOPE::abs( x );
        }
    };
    
    class Sign 
      : public Unary_function< Real_embeddable, CGAL::Sign > {
      public:
        CGAL::Sign operator()( const Real_embeddable& x ) const {
          return (CGAL::Sign) CGAL_LEDA_SCOPE::sign( x );
        }        
    };
    
    class Compare 
      : public Binary_function< Real_embeddable, Real_embeddable,
                                CGAL::Comparison_result > {
      public:
        CGAL::Comparison_result operator()( const Real_embeddable& x, 
                                            const Real_embeddable& y ) const {
          return (CGAL::Comparison_result) CGAL_LEDA_SCOPE::compare( x, y );
        }
        
    };
    
    class To_double 
      : public Unary_function< Real_embeddable, double > {
      public:
        double operator()( const Real_embeddable& x ) const {
          return x.to_double();
        }
    };
    
    class To_interval 
      : public Unary_function< Real_embeddable, std::pair< double, double > > {
      public:
        std::pair<double, double> operator()( const Real_embeddable& x ) const {

          Protect_FPU_rounding<true> P (CGAL_FE_TONEAREST);
          double cn = CGAL_NTS to_double(x);
          leda_integer pn = ( x>0 ? x : -x);
          if ( pn.iszero() || log(pn) < 53 )
              return CGAL_NTS to_interval(cn);
          else {
            FPU_set_cw(CGAL_FE_UPWARD);
            Interval_nt_advanced ina(cn);
            ina += Interval_nt_advanced::smallest();
            return ina.pair();
          }
          
/*        CGAL_LEDA_SCOPE::bigfloat h(x);
          CGAL_LEDA_SCOPE::bigfloat low = 
                        CGAL_LEDA_SCOPE::round(h,53,CGAL_LEDA_SCOPE::TO_N_INF);
          CGAL_LEDA_SCOPE::bigfloat high = 
                        CGAL_LEDA_SCOPE::round(h,53,CGAL_LEDA_SCOPE::TO_P_INF);
          return Double_interval(low.to_double(), high.to_double());                    
        }*/
        }
    };
};

//
// Needs_parens_as_product
//
template <> 
struct Needs_parens_as_product<leda_integer> {
  bool operator()(const leda_integer& x) {
    return CGAL_NTS is_negative(x);
  } 
};


inline
io_Operator
io_tag(const leda_integer &)
{ return io_Operator(); }

// missing mixed operators
inline
bool
operator==(int a, const leda_integer& b)
{ return b == a; }

inline
bool
operator!=(int a, const leda_integer& b)
{ return b != a; }

CGAL_END_NAMESPACE

// Unary + is missing for leda::integer
namespace leda {
    inline integer operator+( const integer& i) { return i; }
} // namespace leda


#endif // CGAL_LEDA_INTEGER_H
