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
 

#ifndef CGAL_LEDA_BIGFLOAT_H
#define CGAL_LEDA_BIGFLOAT_H

#include <CGAL/basic.h>
#include <CGAL/leda_coercion_traits.h>
#include <CGAL/Number_type_traits.h>
#include <CGAL/Interval_nt.h>
#include <CGAL/Algebraic_structure_traits.h>
#include <CGAL/Real_embeddable_traits.h>
#include <CGAL/utils.h>
#include <CGAL/functional_base.h> // Unary_function, Binary_function

#include <utility>


#include <CGAL/LEDA_basic.h>
#if CGAL_LEDA_VERSION < 500
#include <LEDA/bigfloat.h>
#else
#include <LEDA/numbers/bigfloat.h>
#endif


CGAL_BEGIN_NAMESPACE

template <> struct Number_type_traits<leda_bigfloat> {
  typedef Tag_false Has_gcd;
  typedef Tag_true  Has_division;
  typedef Tag_true  Has_sqrt;

  typedef Tag_true  Has_exact_ring_operations;
  typedef Tag_false Has_exact_division;
  typedef Tag_false Has_exact_sqrt;
};

template <> class Algebraic_structure_traits< leda_bigfloat >
  : public Algebraic_structure_traits_base< leda_bigfloat, 
                                            Field_with_kth_root_tag >  {
  public:
    typedef Tag_false           Is_exact;
                                         
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
            // heuristic: we ask for as many precision as the argument has
            long d = x.get_significant_length();
            if ( d < 53) // O.K. we want at least double precision
                d = 53;
            return CGAL_LEDA_SCOPE::sqrt_d( x, d, k);
        }
    };
    
};

template <> class Real_embeddable_traits< leda_bigfloat > 
  : public Real_embeddable_traits_base< leda_bigfloat > {
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

          // assuming leda_bigfloat guarantee 1 bit error max
          Protect_FPU_rounding<true> P (CGAL_FE_TONEAREST);
          Interval_nt_advanced approx (CGAL_LEDA_SCOPE::to_double(x));
          FPU_set_cw(CGAL_FE_UPWARD);
          approx += Interval_nt<false>::smallest();
          return approx.pair();          
        }
    };
    
    class Is_finite
      : public Unary_function< Real_embeddable, bool > {
      public:
        bool operator()( const Real_embeddable& x )  const {
          return !( CGAL_LEDA_SCOPE::isInf(x) || CGAL_LEDA_SCOPE::isNaN(x) ); 
        }
    };
};

template<>
class Is_valid< leda_bigfloat > 
  : public Unary_function< leda_bigfloat, bool > {
  public :
    bool operator()( const leda_bigfloat& x ) const {
      return !( CGAL_LEDA_SCOPE::isNaN(x) );
    }  
};

inline
io_Operator
io_tag(const leda_bigfloat &)
{ return io_Operator(); }

CGAL_END_NAMESPACE

// Unary + is missing for leda::bigfloat
namespace leda {
    inline bigfloat operator+( const bigfloat& i) { return i; }
} // namespace leda

#endif // CGAL_LEDA_BIGFLOAT_H
