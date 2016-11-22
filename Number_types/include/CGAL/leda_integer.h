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
// Author(s)     : Andreas Fabri, Michael Hemmer

#ifndef CGAL_LEDA_INTEGER_H
#define CGAL_LEDA_INTEGER_H

#include <CGAL/number_type_basic.h>

#include <utility>

#include <CGAL/leda_coercion_traits.h>
#include <CGAL/Interval_nt.h>

#include <CGAL/LEDA_basic.h>
#include <LEDA/numbers/integer.h>
#include <LEDA/numbers/bigfloat.h>// for To_interval

#include <CGAL/Residue.h>
#include <CGAL/Modular_traits.h>

namespace CGAL {


template <> class Algebraic_structure_traits< leda_integer >
  : public Algebraic_structure_traits_base< leda_integer,
                                            Euclidean_ring_tag >  {
  public:
    typedef Tag_true            Is_exact;
    typedef Tag_false           Is_numerical_sensitive;

    typedef INTERN_AST::Is_square_per_sqrt< Type >
                                                                 Is_square;

    class Gcd
      : public std::binary_function< Type, Type,
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

    // Unfortunately the behaviour of leda has changed here several times
    // The following Div_mod is invariant under these changes
    // However, the Div and Mod defined below might be more efficient 
    // TODO: recover Div Mod implementation for all leda versions
    class Div_mod {
    public: 
        typedef Type first_argument_type;
        typedef Type second_argument_type; 
        typedef Type& third_argument_type; 
        typedef Type& fourth_argument_type; 
        typedef void result_type;
        
        void operator()(const Type& x, const Type& y, Type& q, Type& r) const {
            
            q = x / y;             
            r = x - q*y;
            CGAL_postcondition(x == y*q + r);  
            
            if (r == 0) return;   
             
            // round q towards zero 
            if ( r.sign() != x.sign() ){
                q -= x.sign();
                r -= x.sign()*y;
            }

            CGAL_postcondition(x == y*q + r);            
            CGAL_postcondition(r.sign() == x.sign());
        }  
    };
    // Div defined via base using Div_mod
    // Mod defined via base using Div_mod

    // This code results in an inconsisten div/mod for some leda versions 
    // TODO: reactivate this code 

//     typedef INTERN_AST::Div_per_operator< Type > Div;
//     class Mod
//       : public std::binary_function< Type, Type,
//                                 Type > {
//       public:
//         Type operator()( const Type& x, const Type& y ) const {
//           Type m = x % y;
//           return m;
//         }
//         CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR( Type )
//     };

    class Sqrt
      : public std::unary_function< Type, Type > {
      public:
        Type operator()( const Type& x ) const {
          return CGAL_LEDA_SCOPE::sqrt( x );
        }
    };
};

template <> class Real_embeddable_traits< leda_integer >
  : public INTERN_RET::Real_embeddable_traits_base< leda_integer , CGAL::Tag_true > {
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
        leda::bigfloat h(x);
	double abs_err = 0;
	double  low =h.to_double(abs_err, leda::TO_N_INF);
	double high =h.to_double(abs_err, leda::TO_P_INF);
	return std::make_pair(low,high);
      }
    };
};

template<>
class Modular_traits< ::leda::integer > {
    typedef Residue MOD;
 public:
    typedef ::leda::integer NT;
    typedef ::CGAL::Tag_true Is_modularizable;
    typedef MOD Residue_type;

    struct Modular_image{
        Residue_type operator()(const NT& a){
            return Residue_type ((a%NT(MOD::get_current_prime())).to_long());
        }
    };
    struct Modular_image_representative{
        NT operator()(const Residue_type& x){
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

// Benchmark_rep specialization 
template<>
class Benchmark_rep< leda_integer > {
    const leda_integer& t;
public:
    //! initialize with a const reference to \a t.
    Benchmark_rep( const leda_integer& tt) : t(tt) {}
    //! perform the output, calls \c operator\<\< by default.
    std::ostream& operator()( std::ostream& out) const { 
            out << t;
            return out;
    }
    
    static std::string get_benchmark_name() {
        return "Integer";
    }
};


} //namespace CGAL

// Unary + is missing for leda::integer
namespace leda {
    inline integer operator+( const integer& i) { return i; }
} // namespace leda

//since types are included by LEDA_coercion_traits.h:
#include <CGAL/leda_rational.h>
#include <CGAL/leda_bigfloat.h>
#include <CGAL/leda_real.h>
#include <CGAL/LEDA_arithmetic_kernel.h>

#endif // CGAL_LEDA_INTEGER_H
