// Copyright (c) 2002,2003  Utrecht University (The Netherlands),
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
// Author(s)     : Andreas Fabri, Sylvain Pion


#ifndef CGAL_GMPQ_H
#define CGAL_GMPQ_H

#include <CGAL/basic.h>
#include <CGAL/GMP/Gmpq_type.h>
#include <CGAL/Gmpz.h>
#include <CGAL/Gmp_coercion_traits.h>

CGAL_BEGIN_NAMESPACE

inline
io_Operator
io_tag(const Gmpq&)
{ return io_Operator(); }


// AST for Gmpq-class
template <> class Algebraic_structure_traits< Gmpq >
  : public Algebraic_structure_traits_base< Gmpq, Field_tag >  {
  public:
    typedef Tag_true            Is_exact;
    
    class Is_square
      : public Binary_function< Algebraic_structure, Algebraic_structure&, 
                                bool > {
      public:
        bool operator()( const Algebraic_structure& x_, Algebraic_structure& y ) const {
          Gmpq x( x_ );
          mpq_canonicalize( x.mpq() );
          Algebraic_structure_traits< Gmpz >::Sqrt sqrt;
          y = Gmpq( sqrt( x.numerator() ), sqrt( x.denominator() ) );
          return y*y == x;
        }
        bool operator()( const Algebraic_structure& x) const {
            Algebraic_structure y;
            return operator()(x,y);
        }
        
    };

    class Simplify 
      : public Unary_function< Algebraic_structure&, void > {
      public:
        void operator()( Algebraic_structure& x) const {
          mpq_canonicalize( x.mpq() );
        }
    };
                                                                 
};

// RET for Gmpq-class

template <> class Real_embeddable_traits< Gmpq > 
  : public Real_embeddable_traits_base< Gmpq > {
  public:

    class Sign 
      : public Unary_function< Real_embeddable, ::CGAL::Sign > {
      public:
        ::CGAL::Sign operator()( const Real_embeddable& x ) const {
          return x.sign();
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
          mpfr_t y;
          mpfr_init2 (y, 53); /* Assume IEEE-754 */
          mpfr_set_q (y, x.mpq(), GMP_RNDD);
          double i = mpfr_get_d (y, GMP_RNDD); /* EXACT but can overflow */
          mpfr_set_q (y, x.mpq(), GMP_RNDU);
          double s = mpfr_get_d (y, GMP_RNDU); /* EXACT but can overflow */
          mpfr_clear (y);
          return std::pair<double, double>(i, s);
        }
    };
};

/*! \ingroup NiX_Fraction_traits_spec
 *  \brief Specialization of Fraction_traits for Gmpq
 */
template <>
class Fraction_traits< Gmpq > {
public:
    typedef Gmpq Fraction;
    typedef ::CGAL::Tag_true Is_fraction;
    typedef Gmpz Numerator;
    typedef Gmpz Denominator;
    typedef Algebraic_structure_traits< Gmpz >::Gcd Common_factor;
    class Decompose {
    public:
        typedef Gmpq first_argument_type;
        typedef Gmpz& second_argument_type;
        typedef Gmpz& third_argument_type;
        void operator () (const Gmpq& rat, Gmpz& num,Gmpz& den) {
            num = rat.numerator();
            den = rat.denominator();
        }
    };
    class Compose {
    public:
        typedef Gmpz first_argument_type;
        typedef Gmpz second_argument_type;
        typedef Gmpq result_type;
        Gmpq operator () (const Gmpz& num,const Gmpz& den) {
            return Gmpq(num, den);
        }
    };
};





/* FIX ME: this not compile
inline
Root_of_2< Gmpz >
make_root_of_2(const Gmpq &a, const CGAL::Gmpq &b,
               const Gmpq &c, bool d)
{
  return CGALi::make_root_of_2_rational< Gmpz, CGAL::Gmpq >(a,b,c,d);
}
*/

#include <CGAL/make_root_of_2.h>
#include <CGAL/Root_of_traits.h>
#include <CGAL/Root_of_2.h>
// Gmpq is the same as Root_of_traits< CGAL::Gmpz >::RootOf_1
template <>
struct Root_of_traits< Gmpq >
{
  typedef Gmpq               RootOf_1;
  typedef Root_of_2< Gmpz >  RootOf_2;
};

CGAL_END_NAMESPACE

#endif // CGAL_GMPQ_H
