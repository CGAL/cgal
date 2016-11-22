// Copyright (c) 1999,2007  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
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
// Author(s)     : Stefan Schirra, Michael Hemmer

#ifndef CGAL_LEDA_BIGFLOAT_H
#define CGAL_LEDA_BIGFLOAT_H

#include <CGAL/basic.h>

#include <utility>
#include <CGAL/leda_coercion_traits.h>
#include <CGAL/Interval_nt.h>

#include <CGAL/LEDA_basic.h>
#include <LEDA/numbers/bigfloat.h>

namespace CGAL {

template <> class Algebraic_structure_traits< leda_bigfloat >
  : public Algebraic_structure_traits_base< leda_bigfloat,
                                            Field_with_kth_root_tag >  {
  public:
    typedef Tag_false           Is_exact;
    typedef Tag_true            Is_numerical_sensitive;

    class Sqrt
      : public std::unary_function< Type, Type > {
      public:
        Type operator()( const Type& x ) const {
          return CGAL_LEDA_SCOPE::sqrt( x );
        }
    };

    class Kth_root
      : public std::binary_function<int, Type, Type> {
      public:
        Type operator()( int k,
                                        const Type& x) const {
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
  : public INTERN_RET::Real_embeddable_traits_base< leda_bigfloat , CGAL::Tag_true > {
public:
  
    class Abs
      : public std::unary_function< Type, Type > {
      public:
        Type operator()( const Type& x ) const {
            return CGAL_LEDA_SCOPE::abs( x );
        }
    };

    class Sgn
      : public std::unary_function< Type, ::CGAL::Sign > {
      public:
        ::CGAL::Sign operator()( const Type& x ) const {
          return (::CGAL::Sign) CGAL_LEDA_SCOPE::sign( x );
        }
    };

    class Compare
      : public std::binary_function< Type, Type,
                                Comparison_result > {
      public:
        Comparison_result operator()( const Type& x,
                                            const Type& y ) const {
          return (Comparison_result) CGAL_LEDA_SCOPE::compare( x, y );
        }
        
        CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR_WITH_RT( Type, 
                Comparison_result )
    };

    class To_double
      : public std::unary_function< Type, double > {
      public:
        double operator()( const Type& x ) const {
          return x.to_double();
        }
    };

    class To_interval
      : public std::unary_function< Type, std::pair< double, double > > {
      public:
        std::pair<double, double> operator()( const Type& x ) const {

          // assuming leda_bigfloat guarantee 1 bit error max
          Protect_FPU_rounding<true> P (CGAL_FE_TONEAREST);
          Interval_nt_advanced approx (CGAL_LEDA_SCOPE::to_double(x));
          FPU_set_cw(CGAL_FE_UPWARD);
          approx += Interval_nt<false>::smallest();
          return approx.pair();
        }
    };

    class Is_finite
      : public std::unary_function< Type, bool > {
      public:
        bool operator()( const Type& x )  const {
          return !( CGAL_LEDA_SCOPE::isInf(x) || CGAL_LEDA_SCOPE::isNaN(x) );
        }
    };
};

template<>
class Is_valid< leda_bigfloat >
  : public std::unary_function< leda_bigfloat, bool > {
  public :
    bool operator()( const leda_bigfloat& x ) const {
      return !( CGAL_LEDA_SCOPE::isNaN(x) );
    }
};


} //namespace CGAL

// Unary + is missing for leda::bigfloat
namespace leda {
    inline bigfloat operator+( const bigfloat& i) { return i; }
} // namespace leda

//since types are included by LEDA_coercion_traits.h:
#include <CGAL/leda_integer.h>
#include <CGAL/leda_rational.h>
#include <CGAL/leda_real.h>
#include <CGAL/LEDA_arithmetic_kernel.h>

#endif // CGAL_LEDA_BIGFLOAT_H
