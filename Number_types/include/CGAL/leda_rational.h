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


#ifndef CGAL_LEDA_RATIONAL_H
#define CGAL_LEDA_RATIONAL_H

#include <CGAL/basic.h>
#include <CGAL/leda_coercion_traits.h>
#include <CGAL/Number_type_traits.h>
#include <CGAL/Interval_nt.h>
#include <CGAL/Algebraic_structure_traits.h>
#include <CGAL/Real_embeddable_traits.h>
#include <CGAL/utils.h>

#include <CGAL/functional_base.h> // Unary_function, Binary_function
#include <CGAL/Needs_parens_as_product.h>

#include <utility>

#include <CGAL/LEDA_basic.h>
#if CGAL_LEDA_VERSION < 500
#include <LEDA/rational.h>
#include <LEDA/interval.h>
#else
#include <LEDA/numbers/rational.h>
#include <LEDA/numbers/interval.h>
#endif


#include <LEDA/rational.h>
#include <LEDA/interval.h>

CGAL_BEGIN_NAMESPACE

template <>
struct Number_type_traits<leda_rational> {
  typedef Tag_false Has_gcd;
  typedef Tag_true  Has_division;
  typedef Tag_false Has_sqrt;

  typedef Tag_true  Has_exact_ring_operations;
  typedef Tag_true  Has_exact_division;
  typedef Tag_false Has_exact_sqrt;
};

template <>
struct Rational_traits<leda_rational> {
  typedef leda_integer RT;
  RT numerator   (const leda_rational & r) const { return r.numerator(); }
  RT denominator (const leda_rational & r) const { return r.denominator(); }
  leda_rational make_rational(const RT & n, const RT & d) const
  { return leda_rational(n, d); }
  leda_rational make_rational(const leda_rational & n,
                              const leda_rational & d) const
  { return n / d; }
};

template <> class Algebraic_structure_traits< leda_rational >
  : public Algebraic_structure_traits_base< leda_rational, 
                                            CGAL::Field_tag >  {
  public:
    typedef CGAL::Tag_true            Is_exact;
                
//    TODO: How to implement this without having sqrt?
//    typedef CGAL::INTERN_AST::Is_square_per_sqrt< Algebraic_structure >
//                                                                 Is_square;

    class Simplify 
      : public Unary_function< Algebraic_structure&, void > {
      public:        
        void operator()( Algebraic_structure& x) const {
            x.normalize();
        }
    };
                                                                 
};

template <> class Real_embeddable_traits< leda_rational > 
  : public Real_embeddable_traits_base< leda_rational > {
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

#if CGAL_LEDA_VERSION >= 501
          CGAL_LEDA_SCOPE::interval temp(x);
          std::pair<double, double> result(temp.lower_bound(),temp.upper_bound());
          CGAL_postcondition(Real_embeddable(result.first)<=x);
          CGAL_postcondition(Real_embeddable(result.second)>=x);
          return result; 
#else 
          CGAL_LEDA_SCOPE::bigfloat xnum = x.numerator();
          CGAL_LEDA_SCOPE::bigfloat xden = x.denominator();
          CGAL_LEDA_SCOPE::bigfloat xupp = 
                                    div(xnum,xden,53,CGAL_LEDA_SCOPE::TO_P_INF);
          CGAL_LEDA_SCOPE::bigfloat xlow = 
                                    div(xnum,xden,53,CGAL_LEDA_SCOPE::TO_N_INF);
    
          // really smallest positive double
          double MinDbl = CGAL_LEDA_SCOPE::fp::compose_parts(0,0,0,1); 
   
          double low = xlow.to_double();
          while(Real_embeddable(low) > x) low = low - MinDbl;
       
          double upp = xupp.to_double();
          while(Real_embeddable(upp) < x) upp = upp + MinDbl;
         
          std::pair<double, double> result(low,upp);
          CGAL_postcondition(Real_embeddable(result.first)<=x);
          CGAL_postcondition(Real_embeddable(result.second)>=x);
          return result; 
#endif
          // Original CGAL to_interval (seemed to be inferior)
          //  // There's no guarantee about the error of to_double(), so I add 
          //  //  3 ulps...
          //  Protect_FPU_rounding<true> P (CGAL_FE_TONEAREST);
          //  Interval_nt_advanced approx (z.to_double());
          //  FPU_set_cw(CGAL_FE_UPWARD);
          //
          //  approx += Interval_nt<false>::smallest();
          //  approx += Interval_nt<false>::smallest();
          //  approx += Interval_nt<false>::smallest();
          //  return approx.pair();
        
        }
    };
};

inline
io_Operator
io_tag(const leda_rational &)
{ return io_Operator(); }

template <class F>
class Output_rep< leda_rational, F> {
    const leda_rational& t;
public:
    //! initialize with a const reference to \a t.
    Output_rep( const leda_rational& tt) : t(tt) {}
    //! perform the output, calls \c operator\<\< by default.
    std::ostream& operator()( std::ostream& out) const {
        switch (CGAL::get_mode(out)) {
        case CGAL::IO::BENCHMARK:
            return out << "Rational(" 
                       << t.numerator()<< "," 
                       << t.denominator() << ")";
            break;
        case CGAL::IO::PRETTY:{
            if(t.denominator() == leda_integer(1))
                return out <<t.numerator();
            else
                return out << t.numerator()
                           << "/" 
                           << t.denominator();
            break;
        }
            
        default:
            return out << t.numerator()
                       << "/" 
                       << t.denominator();
        }
    }
};

template <>
struct Needs_parens_as_product< leda_rational >{
    bool operator()( leda_rational t){
        if (t.denominator() != 1 ) 
            return true;
        else
            return needs_parens_as_product(t.numerator()) ;
    }
};

template <>
class Output_rep< leda_rational, CGAL::Parens_as_product_tag > {
    const leda_rational& t;
public:
    // Constructor 
    Output_rep( const leda_rational& tt) : t(tt) {}
    // operator 
    std::ostream& operator()( std::ostream& out) const { 
        Needs_parens_as_product< leda_rational > needs_parens_as_product;
        if (needs_parens_as_product(t))
            return out <<"("<< oformat(t) <<")";
        else
            return out << oformat(t);
    }
};
////////////////////////////////////////////////////////////////////////////////
// Additional Fraction traits (not part of official concept)
template <>
class Fraction_traits< leda_rational > {
  public:
    typedef leda_rational Fraction;
    typedef CGAL::Tag_true Is_fraction;
    typedef leda_integer Numerator;
    typedef leda_integer Denominator;
    typedef Algebraic_structure_traits< leda_integer >::Gcd Common_factor;
    class Decompose {
      public:
        typedef leda_rational first_argument_type;
        typedef leda_integer& second_argument_type;
        typedef leda_integer& third_argument_type;
        void operator () (const leda_rational& rat,
                          leda_integer& num,
                          leda_integer& den
        ) {
            num = rat.numerator();
            den = rat.denominator();
        }
    };
    class Compose {
    public:
        typedef leda_integer first_argument_type;
        typedef leda_integer second_argument_type;
        typedef leda_rational result_type;
        leda_rational operator () (const leda_integer& num,
                                   const leda_integer& den
        ) {
            return leda_rational(num, den);
        }
    };
};
////////////////////////////////////////////////////////////////////////////////



CGAL_END_NAMESPACE

// Unary + is missing for leda::rational
namespace leda{
inline rational operator+( const rational& i) { return i; }
}


#endif  // CGAL_LEDA_RATIONAL_H
