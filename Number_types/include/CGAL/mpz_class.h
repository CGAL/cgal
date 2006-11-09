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
// $URL: svn+ssh://hemmer@scm.gforge.inria.fr/svn/cgal/branches/CGAL_with_EXACUS/Algebraic_foundations/include/CGAL/gmpxx.h $
// $Id: gmpxx.h 34439 2006-09-21 08:56:06Z slimbach $
// 
//
// Author(s)     : Sylvain Pion
 
#ifndef CGAL_MPZ_CLASS_H
#define CGAL_MPZ_CLASS_H

#include <CGAL/number_type_basic.h>
//#include <gmpxx.h>
#include <CGAL/gmpxx_coercion_traits.h>

// This file gathers the necessary adaptors so that the following
// C++ number types that come with GMP can be used by CGAL :
// - mpz_class

// Note that GMP++ use the expression template mechanism, which makes things
// a little bit complicated in order to make square(x+y) work for example.
// Reading gmpxx.h shows that ::__gmp_expr<T, T> is the mp[zqf]_class proper,
// while ::__gmp_expr<T, U> is the others "expressions".

CGAL_BEGIN_NAMESPACE


template <>
struct Number_type_traits<mpz_class> {
    typedef Tag_false Has_gcd;
    typedef Tag_true  Has_division;
    typedef Tag_true  Has_sqrt;

    typedef Tag_true  Has_exact_ring_operations;
    typedef Tag_false Has_exact_division;
    typedef Tag_false Has_exact_sqrt;
};

// AST for mpz_class

template <class U> 
class Algebraic_structure_traits< ::__gmp_expr< ::__gmpz_value,U>  >
    :public Algebraic_structure_traits_base<  ::__gmp_expr< ::__gmpz_value,U>  , Null_tag > {
public:
    typedef Euclidean_ring_tag  Algebraic_structure_tag;
    typedef Tag_true            Is_exact;
    typedef mpz_class                 Algebraic_structure;
    
    struct Is_zero: public Unary_function< mpz_class , bool > {
        template <class U2> 
        bool operator()( const ::__gmp_expr< ::__gmpz_value,U2>& x) const {
            return ::sgn(x) == 0;
        }
    };  

    struct Is_one: public Unary_function< mpz_class , bool > {
        template <class U2> 
        bool operator()( const ::__gmp_expr< ::__gmpz_value,U2>& x) const {
            return x == 1;
        }
    }; 

    struct Simplify: public Unary_function< mpz_class , void > {
        template <class U2> 
        void operator()( const ::__gmp_expr< ::__gmpz_value,U2>& x) const {}
    }; 
    
    struct Square: public Unary_function< mpz_class , mpz_class > {
        template <class U2> 
        mpz_class operator()( const ::__gmp_expr< ::__gmpz_value,U2>& x) const {
            return x*x;
        }
    }; 

    struct Unit_part: public Unary_function< mpz_class , mpz_class > {
        template <class U2> 
        mpz_class operator()( const ::__gmp_expr< ::__gmpz_value,U2>& x) const {
            return( x < mpz_class(0)) ?  mpz_class(-1) : mpz_class(1); 
        }
    }; 

   

    struct Integral_division: public Binary_function< mpz_class , mpz_class, mpz_class > {
        template <class U2, class U3> 
        mpz_class operator()( 
                const ::__gmp_expr< ::__gmpz_value,U2>& x,
                const ::__gmp_expr< ::__gmpz_value,U3>& y) const {
            mpz_class result = x / y;
            CGAL_precondition_msg( result * y == x,
            "'x' must be divisible by 'y' in "
            "Algebraic_structure_traits<mpz_class>::Integral_div()(x,y)" );
            return result;         
        } 
        CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR( Algebraic_structure )
    }; 
    
    struct Gcd : public Binary_function< mpz_class, mpz_class, mpz_class > {
        template <class U2, class U3> 
        mpz_class operator()( 
                const ::__gmp_expr< ::__gmpz_value,U2>& x,
                const ::__gmp_expr< ::__gmpz_value,U3>& y) const {
            mpz_class c;
            mpz_gcd( c.get_mpz_t(), mpz_class(x).get_mpz_t(), mpz_class(y).get_mpz_t() );
            return c;
        } 
        CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR( Algebraic_structure )
    };
    
    struct Div : public Binary_function< mpz_class, mpz_class, mpz_class > {
        template <class U2, class U3> 
        mpz_class operator()( 
                const ::__gmp_expr< ::__gmpz_value,U2>& x,
                const ::__gmp_expr< ::__gmpz_value,U3>& y) const {
            return x / y; 
        } 
        CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR( Algebraic_structure )
    };
    
    struct Mod : public Binary_function< mpz_class, mpz_class, mpz_class > {
        template <class U2, class U3> 
        mpz_class operator()( 
                const ::__gmp_expr< ::__gmpz_value,U2>& x,
                const ::__gmp_expr< ::__gmpz_value,U3>& y) const {
            return x % y; 
        }
        CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR( Algebraic_structure )
    };
    struct Div_mod {
        typedef mpz_class    first_argument_type;
        typedef mpz_class    second_argument_type;
        typedef mpz_class&   third_argument_type;
        typedef mpz_class&   fourth_argument_type;
        typedef Arity_tag< 4 >         Arity;
        typedef void         result_type;
        template <class U2, class U3> 
        void operator()( 
                const ::__gmp_expr< ::__gmpz_value,U2>& x,
                const ::__gmp_expr< ::__gmpz_value,U3>& y,
                mpz_class& q, 
                mpz_class& r
        ) const {
            typedef Algebraic_structure_traits<mpz_class> Traits;
                typename Traits::Div  actual_div;
                typename Traits::Mod  actual_mod;
                q = actual_div( x, y );
                r = actual_mod( x, y );          
                return;
            };  
    };
    
    
    struct Sqrt: public Unary_function< mpz_class , mpz_class > {
        template <class U2> 
        mpz_class operator()( const ::__gmp_expr< ::__gmpz_value,U2>& x) const {
            return ::sqrt(x);
        }
    }; 
   

    /*struct Is_square: public Binary_function< mpz_class , mpz_class& , bool > {
        template <class U2> 
        bool operator()(
                const ::__gmp_expr< ::__gmpz_value,U2>& x,
                mpz_class&                              r){
            r = ::sqrt(x);
            return (r*r==x) ? true : false; 
        }
        template <class U2> 
        bool operator()(const ::__gmp_expr< ::__gmpz_value,U2>& x){
            mpz_class r = ::sqrt(x);
            return (r*r==x) ? true : false; 
        }
        };*/ 
};

// RET for mpz_class

template <class U> 
class Real_embeddable_traits< ::__gmp_expr< ::__gmpz_value,U>  > 
    : public Real_embeddable_traits_base< ::__gmp_expr< ::__gmpz_value,U> > {
public:
    typedef ::__gmp_expr< ::__gmpz_value, ::__gmpz_value>  mpz_class;
    typedef mpz_class  Real_embeddable;

    struct Is_zero: public Unary_function< mpz_class , bool > {
        template <class U2> 
        bool operator()( const ::__gmp_expr< ::__gmpz_value,U2>& x) const {
            return ::sgn(x) == 0;
        }
    }; 
    struct Is_finite: public Unary_function<mpz_class,bool> {
        template <class U2> 
        bool operator()( const ::__gmp_expr< ::__gmpz_value,U2>& x) const {
            return true;
        }
    };

    struct Is_positive: public Unary_function< mpz_class , bool > {
        template <class U2> 
        bool operator()( const ::__gmp_expr< ::__gmpz_value,U2>& x) const {
            return ::sgn(x) > 0;
        }
    }; 
  
    struct Is_negative: public Unary_function< mpz_class , bool > {
        template <class U2> 
        bool operator()( const ::__gmp_expr< ::__gmpz_value,U2>& x) const {
            return ::sgn(x) < 0;
        }
    };

    struct Abs: public Unary_function< mpz_class , mpz_class > {
        template <class U2> 
        mpz_class operator()( const ::__gmp_expr< ::__gmpz_value,U2>& x) const {
            return ::abs(x);
        }
    };
    
    struct Sign 
        : public Unary_function< mpz_class, ::CGAL::Sign > {
    public:
        template <class U2> 
        ::CGAL::Sign operator()( const ::__gmp_expr< ::__gmpz_value,U2>& x ) const {
            return (::CGAL::Sign) ::sgn( x );
        }        
    };
    
    struct Compare 
        : public Binary_function< mpz_class, mpz_class, Comparison_result > {
        template <class U2, class U3>
        Comparison_result operator()( const ::__gmp_expr< ::__gmpz_value,U2>& x, 
                const ::__gmp_expr< ::__gmpz_value,U3>& y ) const {
            // cmp returns any int value, not just -1/0/1...
            return (Comparison_result) CGAL_NTS sign( ::cmp(x, y) );
        }
        CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR_WITH_RT
        ( Real_embeddable, Comparison_result);
    };
    
    struct To_double 
        : public Unary_function< mpz_class, double > {
        template <class U2>
        double operator()( const ::__gmp_expr< ::__gmpz_value,U2>& x ) const {
            return mpz_class(x).get_d();
        }
    };
    
    struct To_interval 
    
        : public Unary_function< mpz_class, std::pair< double, double > > {
        template <class U2>
        std::pair<double, double> 
        operator()( const ::__gmp_expr< ::__gmpz_value,U2>& x_ ) const {
            mpz_class x = mpz_class(x_);
            mpfr_t y;
            mpfr_init2 (y, 53); /* Assume IEEE-754 */
            mpfr_set_z (y, x.get_mpz_t(), GMP_RNDD);
            double i = mpfr_get_d (y, GMP_RNDD); /* EXACT but can overflow */
            mpfr_set_z (y, x.get_mpz_t(), GMP_RNDU);
            double s = mpfr_get_d (y, GMP_RNDU); /* EXACT but can overflow */
            mpfr_clear (y);
            return std::pair<double, double>(i, s);
        }
    };
};

template <class T, typename U >
inline
io_Operator
io_tag(const ::__gmp_expr< T , U> &)
{ return io_Operator(); }


/* FIX ME: THERE IS NO CONSTRUCTOR FT(x,y)
   AVALIABLE FOR THIS TYPE

   namespace CGALi {

   inline
   Root_of_2< ::mpz_class >
   make_root_of_2_gmpxx(const ::mpz_class & a,
   const ::mpz_class & b,
   const ::mpz_class & c,
   bool d)
   {
   return Root_of_2< ::mpz_class >(a, b, c, d);
   }

   inline
   Root_of_2< ::mpz_class >
   make_root_of_2_gmpxx(const ::mpq_class & a,
   const ::mpq_class & b,
   const ::mpq_class & c,
   bool d)
   {
   typedef Rational_traits< ::mpq_class > Rational;

   Rational r;
   CGAL_assertion( r.denominator(a) > 0 );
   CGAL_assertion( r.denominator(b) > 0 );
   CGAL_assertion( r.denominator(c) > 0 );
*/
/*   const RT lcm = ( r.denominator(a) * r.denominator(b) * r.denominator(c)          )/
     ( gcd( r.denominator(a), gcd(r.denominator(b), r.denominator(c)) ) );

     RT a_ = r.numerator(a) * ( lcm / r.denominator(a) );
     RT b_ = r.numerator(b) * ( lcm / r.denominator(b) );
     RT c_ = r.numerator(c) * ( lcm / r.denominator(c) );
*/
/*    ::mpz_class a_ = r.numerator(a) * r.denominator(b) * r.denominator(c);
      ::mpz_class b_ = r.numerator(b) * r.denominator(a) * r.denominator(c);
      ::mpz_class c_ = r.numerator(c) * r.denominator(a) * r.denominator(b);

      return Root_of_2< ::mpz_class >(a, b, c, d);
      }

      } // CGALi

      template < typename T, typename U1, typename U2, typename U3 >
      inline
      typename Root_of_traits< ::__gmp_expr<T, T> >::RootOf_2
      make_root_of_2(const ::__gmp_expr< T, U1> & a,
      const ::__gmp_expr< T, U2> & b,
      const ::__gmp_expr< T, U3> & c,
      bool d)
      {
      return CGALi::make_root_of_2_gmpxx(a, b, c, d);
      }

      template < typename T, typename U >
      struct Root_of_traits< ::__gmp_expr<T, U> >
      {
      typedef ::mpq_class               RootOf_1;
      typedef Root_of_2< ::mpz_class >  RootOf_2;
      };*/

CGAL_END_NAMESPACE

// XXX : These seem necessary.
// I don't know why doing them in namespace CGAL is not enough.
// using to_double;
// using is_valid;


#endif // CGAL_MPZ_CLASS_H
