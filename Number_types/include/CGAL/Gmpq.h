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
#include <CGAL/Gmpq_type.h>
#include <CGAL/Gmp_coercion_traits.h>

CGAL_BEGIN_NAMESPACE

inline
io_Operator
io_tag(const Gmpq&)
{ return io_Operator(); }


// AST for Gmpq-class
template <> class Algebraic_structure_traits< Gmpq >
  : public Algebraic_structure_traits_base< Gmpq, CGAL::Field_tag >  {
  public:
    typedef CGAL::Tag_true            Is_exact;
    
    class Is_square
      : public Binary_function< Algebraic_structure, Algebraic_structure&, 
                                bool > {
      public:
        bool operator()( const Algebraic_structure& x_, Algebraic_structure& y ) {
          Gmpq x( x_ );
          mpq_canonicalize( x.mpq() );
          Algebraic_structure_traits< Gmpz >::Sqrt sqrt;
          y = Gmpq( sqrt( x.numerator() ), sqrt( x.denominator() ) );
          return y*y == x;
        }
        bool operator()( const Algebraic_structure& x) {
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
      : public Unary_function< Real_embeddable, CGAL::Sign > {
      public:
        CGAL::Sign operator()( const Real_embeddable& x ) const {
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


template <>
struct Rational_traits<Gmpq> {
  typedef Gmpz RT;
  RT   numerator     (const Gmpq & r) const { return r.numerator(); }
  RT   denominator   (const Gmpq & r) const { return r.denominator(); }
  Gmpq make_rational (const RT & n, const RT & d) const
  { return Gmpq(n, d); }
  Gmpq make_rational (const Gmpq & n, const Gmpq & d) const
  { return n / d; }
};

/* FIX ME: this not compile
inline
Root_of_2< CGAL::Gmpz >
make_root_of_2(const CGAL::Gmpq &a, const CGAL::Gmpq &b,
               const CGAL::Gmpq &c, bool d)
{
  return CGALi::make_root_of_2_rational< CGAL::Gmpz, CGAL::Gmpq >(a,b,c,d);
}
*/

#include <CGAL/make_root_of_2.h>
#include <CGAL/Root_of_traits.h>
#include <CGAL/Root_of_2.h>
// CGAL::Gmpq is the same as Root_of_traits< CGAL::Gmpz >::RootOf_1
template <>
struct Root_of_traits< CGAL::Gmpq >
{
  typedef CGAL::Gmpq               RootOf_1;
  typedef Root_of_2< CGAL::Gmpz >  RootOf_2;
};

CGAL_END_NAMESPACE

#endif // CGAL_GMPQ_H
