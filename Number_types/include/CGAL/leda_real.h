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
// Author(s)     : Stefan Schirra, Michael Hemmer
 
#ifndef CGAL_LEDA_REAL_H
#define CGAL_LEDA_REAL_H

#include <CGAL/number_type_basic.h>
#include <CGAL/leda_coercion_traits.h>

#include <CGAL/utils.h>
#include <CGAL/Interval_nt.h>

#include <utility>

#include <CGAL/LEDA_basic.h>
#if CGAL_LEDA_VERSION < 500
#include <LEDA/real.h>
#else
#include <LEDA/numbers/real.h>
#endif


CGAL_BEGIN_NAMESPACE

template <> struct Number_type_traits<leda_real> {
  typedef Tag_false Has_gcd;
  typedef Tag_true  Has_division;
  typedef Tag_true  Has_sqrt;

  typedef Tag_true  Has_exact_ring_operations;
  typedef Tag_true  Has_exact_division;
  typedef Tag_true  Has_exact_sqrt;
};

template <> class Algebraic_structure_traits< leda_real >

#if CGAL_LEDA_VERSION >= 500 
  : public Algebraic_structure_traits_base< leda_real, 
                                            Field_with_root_of_tag >  {
#else
  : public Algebraic_structure_traits_base< leda_real, 
                                            Field_with_kth_root_tag >  {
#endif

  public:
    typedef Tag_true           Is_exact;
                                                                             
    class Sqrt 
      : public Unary_function< Algebraic_structure, Algebraic_structure > {
      public:
        Algebraic_structure operator()( const Algebraic_structure& x ) const {
          return CGAL_LEDA_SCOPE::sqrt( x );
        }
    };
    
    class Kth_root 
      : public Binary_function<int, Algebraic_structure, Algebraic_structure> {
      public:
        Algebraic_structure operator()( int k, 
                                        const Algebraic_structure& x) const {
            CGAL_precondition_msg(k > 0, "'k' must be positive for k-th roots");
            return CGAL_LEDA_SCOPE::root( x, k);
        }
    };

// Root_of is only available for LEDA versions >= 5.0
#if CGAL_LEDA_VERSION >= 500 
    class Root_of {
      public:
        typedef Algebraic_structure result_type;
        typedef Arity_tag< 3 >         Arity;

//        typedef leda_rational Boundary;
      private:
        template< class ForwardIterator >
        inline 
        CGAL_LEDA_SCOPE::polynomial<Algebraic_structure>
        make_polynomial(ForwardIterator begin, 
                        ForwardIterator end) const {
          CGAL_LEDA_SCOPE::growing_array<Algebraic_structure> coeffs;
          for(ForwardIterator it = begin; it < end; it++) 
              coeffs.push_back(*it);
          return CGAL_LEDA_SCOPE::polynomial<Algebraic_structure>(coeffs);
        }
      public:
        template <class ForwardIterator>
        Algebraic_structure operator()( int k, 
                       ForwardIterator begin, 
                       ForwardIterator end) const {
            return CGAL_LEDA_SCOPE::diamond(k,make_polynomial(begin,end));
        };
/*        template <class ForwardIterator>
        Algebraic_structure operator()( leda_rational lower,
                                        leda_rational upper,
                                        ForwardIterator begin, 
                                        ForwardIterator end) const {
            return CGAL_LEDA_SCOPE::diamond(lower,upper,
                                             make_polynomial(begin,end));
        };*/
    };
    
#endif

    
};

template <> class Real_embeddable_traits< leda_real > 
  : public Real_embeddable_traits_base< leda_real > {
  public:
      
    class Abs 
      : public Unary_function< Real_embeddable, Real_embeddable > {
      public:
        Real_embeddable operator()( const Real_embeddable& x ) const {
            return CGAL_LEDA_SCOPE::abs( x );
        }
    };
    
    class Sign 
      : public Unary_function< Real_embeddable, ::CGAL::Sign > {
      public:
        ::CGAL::Sign operator()( const Real_embeddable& x ) const {
          return (::CGAL::Sign) CGAL_LEDA_SCOPE::sign( x );
        }        
    };
    
    class Compare 
      : public Binary_function< Real_embeddable, Real_embeddable,
                                Comparison_result > {
      public:
        Comparison_result operator()( const Real_embeddable& x, 
                                            const Real_embeddable& y ) const {
          return (Comparison_result) CGAL_LEDA_SCOPE::compare( x, y );
        }
        
        CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR_WITH_RT( Real_embeddable,
                                                      Comparison_result );
        
    };
    
    class To_double 
      : public Unary_function< Real_embeddable, double > {
      public:
        double operator()( const Real_embeddable& x ) const {
          // this call is required to get reasonable values for the double
          // approximation (as of LEDA-4.3.1)
          x.improve_approximation_to(53);
          return x.to_double();
        }
    };
    
    class To_interval 
      : public Unary_function< Real_embeddable, std::pair< double, double > > {
      public:
        std::pair<double, double> operator()( const Real_embeddable& x ) const {

#if CGAL_LEDA_VERSION >= 501
            leda_bigfloat bnum = x.to_bigfloat();  
            leda_bigfloat berr = x.get_bigfloat_error();
            
            double dummy;
            double low = (bnum - berr).to_double(dummy, 
                                                     CGAL_LEDA_SCOPE::TO_N_INF);
            double upp = (bnum + berr).to_double(dummy, 
                                                     CGAL_LEDA_SCOPE::TO_P_INF);
            
            std::pair<double, double> result(low, upp);
            CGAL_postcondition(Real_embeddable(result.first)<=x);
            CGAL_postcondition(Real_embeddable(result.second)>=x);
            return result;
#else
            CGAL_LEDA_SCOPE::interval temp(x); //bug in leda
            std::pair<double, double> result(temp.lower_bound(),temp.upper_bound());
            CGAL_postcondition_msg(Real_embeddable(result.first)<=x, 
                                                    "Known bug in LEDA <=5.0");
            CGAL_postcondition_msg(Real_embeddable(result.first)>=x, 
                                                    "Known bug in LEDA <=5.0");
            return result;
            // If x is very small and we look closer at x 
            // (i.e. comparison or to_double() or to_bigfloat())
            // then x gets 0, which is really bad.
            // Therefore we do not touch x. 
            // The LEDA interval above returns (-inf, inf) for
            // very small x, which is also bad and leads to 
            // problems lateron. The postcondition fails in this
            // situation.
#endif
              // Original CGAL to_interval:
            //  Protect_FPU_rounding<true> P (CGAL_FE_TONEAREST);
            //  double approx = z.to_double();
            //  double rel_error = z.get_double_error();
            //  FPU_set_cw(CGAL_FE_UPWARD);
            //  Interval_nt_advanced ina(-rel_error,rel_error);
            //  ina += 1;
            //  ina *= approx;
            //  return ina.pair();
        }
    };
};

inline
io_Operator
io_tag(const leda_real &)
{ return io_Operator(); }

CGAL_END_NAMESPACE

// Unary + is missing for leda::real

namespace leda {
    inline real operator+( const real& i) { return i; }
} // namespace leda


#endif // CGAL_LEDA_REAL_H
