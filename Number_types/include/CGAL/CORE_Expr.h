// Copyright (c) 2002-2004  Utrecht University (The Netherlands),
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
// Author(s)     : Sylvain Pion, Michael Hemmer
 
#ifndef CGAL_CORE_EXPR_H
#define CGAL_CORE_EXPR_H

#include <CGAL/basic.h>
#include <CGAL/Number_type_traits.h>
#include <CGAL/Interval_nt.h>
#include <CGAL/CORE_Expr_fwd.h>
#include <CGAL/functional_base.h> // Unary_function, Binary_function

#include <CGAL/Algebraic_structure_traits.h>
#include <CGAL/Real_embeddable_traits.h>
#include <CGAL/CORE_coercion_traits.h>

#include <utility>

#define CORE_LEVEL 4
#include <CORE/CORE.h>

CGAL_BEGIN_NAMESPACE

template <>
struct Number_type_traits<CORE::Expr> {
  typedef Tag_false Has_gcd;
  typedef Tag_true  Has_division;
  typedef Tag_true  Has_sqrt;

  typedef Tag_true  Has_exact_ring_operations;
  typedef Tag_true  Has_exact_division;
  typedef Tag_true  Has_exact_sqrt;
};

template <> class Algebraic_structure_traits< CORE::Expr >
  : public Algebraic_structure_traits_base< CORE::Expr, 
                                            CGAL::Field_with_root_of_tag >  {
  public:
    typedef CGAL::Tag_true            Is_exact;

    class Sqrt 
      : public Unary_function< Algebraic_structure, Algebraic_structure > {
      public:
        Algebraic_structure operator()( const Algebraic_structure& x ) const {
          return CORE::sqrt( x );
        }
    };
    
    class Kth_root 
      : public Binary_function<int, Algebraic_structure, Algebraic_structure> {
      public:
        Algebraic_structure operator()( int k, 
                                        const Algebraic_structure& x) const {
          CGAL_precondition_msg( k > 0, "'k' must be positive for k-th roots");
          // CORE::radical isn't implemented for negative values of x, so we
          //  have to handle this case separately
          if( x < 0 && k%2 != 0)
            return -CORE::radical( -x, k );
  
          return CORE::radical( x, k );
        }
    };
    
    class Root_of {
      public:
//        typedef CORE::BigRat Boundary;
        typedef Algebraic_structure   result_type;

        typedef Arity_tag< 3 >         Arity;
        
      public:
        // constructs the kth roots of the polynomial 
        // given by the iterator range, starting from 0. 
        template< class ForwardIterator >
        Algebraic_structure operator()( int k, 
                                        ForwardIterator begin, 
                                        ForwardIterator end) const {
            std::vector<Algebraic_structure> coeffs;
            for(ForwardIterator it = begin; it != end; it++){
                coeffs.push_back(*it);
            }
            CORE::Polynomial<Algebraic_structure> polynomial(coeffs);
            return Algebraic_structure(polynomial,k);
        };
        
// TODO: Need to be fixed: polynomial<CORE::Expr>.eval() cannot return 
//       CORE::BigFloat, so this does not compile.

/*        template <class ForwardIterator>
        Algebraic_structure operator()( CORE::BigRat lower,
                                        CORE::BigRat upper,
                                        ForwardIterator begin, 
                                        ForwardIterator end) const {
            std::vector<Algebraic_structure> coeffs;
            for(ForwardIterator it = begin; it != end; it++){
                 coeffs.push_back(*it);
            }
            CORE::Polynomial<Algebraic_structure> polynomial(coeffs);
            CORE::BigFloat lower_bf, upper_bf;
            CORE::BigFloat eval_at_lower(0), eval_at_upper(0);
            
            CORE::extLong r(16),a(16);
            while((eval_at_lower.isZeroIn() || 
                   eval_at_upper.isZeroIn())){
                //std::cout << "while"<<std::endl;
                r*=2; 
                a*=2;
                lower_bf.approx(lower,r,a);
                upper_bf.approx(upper,r,a);
                // The most expensive precond I've ever seen :)),
                // since the coefficients of the polynomial are CORE::Expr 
                // TODO: be sure that lower_bf, upper_bf contain exactly one root
                //NiX_expensive_precond(
                //     CORE::Sturm(polynomial).numberOfRoots(lower_bf,upper_bf)==1);
                eval_at_lower = polynomial.eval(lower_bf);
                eval_at_upper = polynomial.eval(upper_bf);
            }  
            CORE::BFInterval interval(lower_bf,upper_bf);
   
            return Algebraic_structure(polynomial,interval);
        };  */
    };
    
};

template <> class Real_embeddable_traits< CORE::Expr > 
  : public Real_embeddable_traits_base< CORE::Expr > {
  public:
      
    class Abs 
      : public Unary_function< Real_embeddable, Real_embeddable > {
      public:
        Real_embeddable operator()( const Real_embeddable& x ) const {
            return CORE::abs( x );
        }
    };
    
    class Sign 
      : public Unary_function< Real_embeddable, CGAL::Sign > {
      public:
        CGAL::Sign operator()( const Real_embeddable& x ) const {
          return (CGAL::Sign) CORE::sign( x );
        }        
    };
    
    class Compare 
      : public Binary_function< Real_embeddable, Real_embeddable,
                                CGAL::Comparison_result > {
      public:
        CGAL::Comparison_result operator()( const Real_embeddable& x, 
                                            const Real_embeddable& y ) const {
          return (CGAL::Comparison_result) CORE::cmp( x, y );
        }
        
        CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR_WITH_RT( Real_embeddable,
                                                      CGAL::Comparison_result );
        
    };
    
    class To_double 
      : public Unary_function< Real_embeddable, double > {
      public:
        double operator()( const Real_embeddable& x ) const {
          // this call is required to get reasonable values for the double
          // approximation 
          x.approx( 53, 1 );
          return x.doubleValue();
        }
    };
    
    class To_interval 
      : public Unary_function< Real_embeddable, std::pair< double, double > > {
      public:
        std::pair<double, double> operator()( const Real_embeddable& x ) const {

            std::pair<double,double> result;
            x.doubleInterval(result.first, result.second);
            CGAL_expensive_assertion(result.first <= x);
            CGAL_expensive_assertion(result.second >= x);
            return result;
        }          
    };
};

inline
io_Operator
io_tag(const CORE::Expr &)
{ return io_Operator(); }

CGAL_END_NAMESPACE

#endif // CGAL_CORE_EXPR_H
