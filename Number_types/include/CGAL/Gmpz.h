// Copyright (c) 1999,2003,2004  Utrecht University (The Netherlands),
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
// Author(s)     : Andreas Fabri, Stefan Schirra, Sylvain Pion


#ifndef CGAL_GMPZ_H
#define CGAL_GMPZ_H

#include <CGAL/basic.h>
#include <CGAL/Gmpz_type.h>
#include <CGAL/Gmp_coercion_traits.h>

#include <CGAL/Quotient.h> // spec of AST for Quotient<Gmpz>

#include <string>
#ifndef CGAL_CFG_NO_LOCALE
#  include <locale>
#else
#  include <cctype>
#endif

//#include <CGAL/Root_of_traits.h>
//#include <CGAL/Root_of_2_fwd.h>

CGAL_BEGIN_NAMESPACE

inline
io_Operator
io_tag(const Gmpz&)
{ return io_Operator(); }


// Algebraic structure traits
template <> class Algebraic_structure_traits< Gmpz >
    : public Algebraic_structure_traits_base< Gmpz, 
                                            CGAL::Euclidean_ring_tag >  {
public:
    typedef CGAL::Tag_true            Is_exact;
                
    typedef CGAL::INTERN_AST::Is_square_per_sqrt< Algebraic_structure >
    Is_square;
    class Integral_division 
        : public Binary_function< Algebraic_structure, Algebraic_structure,
                                Algebraic_structure > {
    public:
        Algebraic_structure operator()( const Algebraic_structure& x,
                const Algebraic_structure& y ) const {
            Gmpz result;
            mpz_divexact(result.mpz(), x.mpz(), y.mpz());
            CGAL_postcondition_msg(result * y == x, "exact_division failed\n");
            return result;          
        }
    };
                                                                 
    class Gcd 
        : public Binary_function< Algebraic_structure, Algebraic_structure, 
                                Algebraic_structure > {
    public:
        Algebraic_structure operator()( const Algebraic_structure& x, 
                const Algebraic_structure& y ) const {
            Gmpz result;
            mpz_gcd(result.mpz(), x.mpz(), y.mpz());
            return result;

        }
        
        Algebraic_structure operator()( const Algebraic_structure& x,
                                        const int& y ) const {
          if (y > 0)
          {
              Gmpz Res;
              mpz_gcd_ui(Res.mpz(), x.mpz(), y);
              return Res;
          }
          return CGAL_NTS gcd(x, Gmpz(y));                                          
        }
        
        Algebraic_structure operator()( const int& x,
                                        const Algebraic_structure& y ) const {
          return CGAL_NTS gcd(Gmpz(x), y );
        }
    };
    
    typedef CGAL::INTERN_AST::Div_per_operator< Algebraic_structure > Div;
    typedef CGAL::INTERN_AST::Mod_per_operator< Algebraic_structure > Mod;
    
    class Sqrt 
        : public Unary_function< Algebraic_structure, Algebraic_structure > {
    public:
        Algebraic_structure operator()( const Algebraic_structure& x ) const {
            Gmpz result;
            mpz_sqrt(result.mpz(), x.mpz());
            return result;
        }
    };        
};

template <> class Real_embeddable_traits< Gmpz > 
    : public Real_embeddable_traits_base< Gmpz > {
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
            mpfr_set_z (y, x.mpz(), GMP_RNDD);
            double i = mpfr_get_d (y, GMP_RNDD); /* EXACT but can overflow */
            mpfr_set_z (y, x.mpz(), GMP_RNDU);
            double s = mpfr_get_d (y, GMP_RNDU); /* EXACT but can overflow */
            mpfr_clear (y);
            return std::pair<double, double>(i, s);
        }
    };
};

template<> class Algebraic_structure_traits< Quotient<Gmpz> >
    : public INTERN_QUOTIENT::Algebraic_structure_traits_quotient_base<Quotient<Gmpz> >{
    // specialization of to double functor
    typedef Quotient<Gmpz> Algebraic_structure;
    
    struct To_double: public Unary_function<Quotient<Gmpz>, double>{
        double operator()(const Quotient<Gmpz>& quot){
            mpq_t  mpQ;
            mpq_init(mpQ);
            const Gmpz& n = quot.numerator();
            const Gmpz& d = quot.denominator();
            mpz_set(mpq_numref(mpQ), n.mpz());
            mpz_set(mpq_denref(mpQ), d.mpz());
            
            mpq_canonicalize(mpQ);
            
            double ret = mpq_get_d(mpQ);
            mpq_clear(mpQ);
            return ret;
        }
    };
};

CGAL_END_NAMESPACE



#include <CGAL/Root_of_2.h>
  
CGAL_BEGIN_NAMESPACE
  
class Gmpq;

template <>
struct Root_of_traits< CGAL::Gmpz >
{
    typedef CGAL::Gmpq               RootOf_1;
    typedef Root_of_2< CGAL::Gmpz >  RootOf_2;
};

// FIX ME: This not compile
// inline
// Root_of_2<Gmpz>
// make_root_of_2(const Gmpz &a, const Gmpz &b, const Gmpz &c, bool smaller)
// {
//   CGAL_assertion( a != 0 );
//   return Root_of_2<Gmpz>(a, b, c, smaller);
// }
CGAL_END_NAMESPACE


#endif // CGAL_GMPZ_H
