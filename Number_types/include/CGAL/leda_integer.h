// Copyright (c) 1999,2007  Utrecht University (The Netherlands),
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

#include <CGAL/number_type_basic.h>

#include <utility>

#include <CGAL/leda_coercion_traits.h>
#include <CGAL/Interval_nt.h>

#include <CGAL/LEDA_basic.h>
#if CGAL_LEDA_VERSION < 500
#include <LEDA/integer.h>
#include <LEDA/bigfloat.h>// for To_interval
#else
#include <LEDA/numbers/integer.h>
#include <LEDA/numbers/bigfloat.h>// for To_interval
#endif

#include <CGAL/Modular.h>
#include <CGAL/Modular_traits.h>

CGAL_BEGIN_NAMESPACE

template <> class Algebraic_structure_traits< leda_integer >
  : public Algebraic_structure_traits_base< leda_integer,
                                            Euclidean_ring_tag >  {
  public:
    typedef Tag_true            Is_exact;
    typedef Tag_false           Is_numerical_sensitive;

    typedef INTERN_AST::Is_square_per_sqrt< Type >
                                                                 Is_square;

    class Gcd
      : public Binary_function< Type, Type,
                                Type > {
      public:
        Type operator()( const Type& x,
                                        const Type& y ) const {
          // By definition gcd(0,0) == 0
          if( x == Type(0) && y == Type(0) )
            return Type(0);

          return CGAL_LEDA_SCOPE::gcd( x, y );
        }

        CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR( Type )
    };

    typedef INTERN_AST::Div_per_operator< Type > Div;

    class Mod
      : public Binary_function< Type, Type,
                                Type > {
      public:
        Type operator()( const Type& x,
                                        const Type& y ) const {
          Type m = x % y;

#if CGAL_LEDA_VERSION < 520
          // Fix wrong leda result
          if( x < 0 && m != 0 )
            m -= y;
#else
          // Fix another wrong leda result
          // TODO: be careful for future improvements of LEDA
          if( x < 0 && y > 0 && m != 0 )
            m -= y;
#endif
          return m;
        }

        CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR( Type )
    };

    class Sqrt
      : public Unary_function< Type, Type > {
      public:
        Type operator()( const Type& x ) const {
          return CGAL_LEDA_SCOPE::sqrt( x );
        }
    };
};

template <> class Real_embeddable_traits< leda_integer >
  : public Real_embeddable_traits_base< leda_integer > {
  public:

    class Abs
      : public Unary_function< Type, Type > {
      public:
        Type operator()( const Type& x ) const {
            return CGAL_LEDA_SCOPE::abs( x );
        }
    };

    class Sign
      : public Unary_function< Type, ::CGAL::Sign > {
      public:
        ::CGAL::Sign operator()( const Type& x ) const {
            return (::CGAL::Sign) CGAL_LEDA_SCOPE::sign( x );
        }
    };

    class Compare
      : public Binary_function< Type, Type,
                                Comparison_result > {
      public:
        Comparison_result operator()( const Type& x,
                                            const Type& y ) const {
          return (Comparison_result) CGAL_LEDA_SCOPE::compare( x, y );
        }

    };

    class To_double
      : public Unary_function< Type, double > {
      public:
        double operator()( const Type& x ) const {
          return x.to_double();
        }
    };

    class To_interval
      : public Unary_function< Type, std::pair< double, double > > {
      public:
        std::pair<double, double> operator()( const Type& x ) const {

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

template<>
class Modular_traits< ::leda::integer > {
    typedef Modular MOD;
 public:
    typedef ::leda::integer NT;
    typedef ::CGAL::Tag_true Is_modularizable;
    typedef MOD Modular_NT;

    struct Modular_image{
        Modular_NT operator()(const NT& a){
            return Modular_NT ((a%NT(MOD::get_current_prime())).to_long());
        }
    };
    struct Modular_image_inv{
        NT operator()(const Modular& x){
            return NT(x.get_value());
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

// missing mixed operators
inline
bool
operator==(int a, const leda_integer& b)
{ return b == a; }

inline
bool
operator!=(int a, const leda_integer& b)
{ return b != a; }


template <>
struct Split_double<leda_integer>
{
  void operator()(double d, leda_integer &num, leda_integer &den) const
  {
    std::pair<double, double> p = split_numerator_denominator(d);
    num = leda_integer(p.first);
    den = leda_integer(p.second);
  }
};

CGAL_END_NAMESPACE

// Unary + is missing for leda::integer
namespace leda {
    inline integer operator+( const integer& i) { return i; }
} // namespace leda

//since types are included by leda_coercion_traits.h:
#include <CGAL/leda_integer.h>
#include <CGAL/leda_rational.h>
#include <CGAL/leda_bigfloat.h>
#include <CGAL/leda_real.h>

#endif // CGAL_LEDA_INTEGER_H
