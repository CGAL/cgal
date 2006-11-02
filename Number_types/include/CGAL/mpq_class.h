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
 
#ifndef CGAL_MPQ_CLASS_H
#define CGAL_MPQ_CLASS_H

#include <CGAL/basic.h>
#include <CGAL/gmpxx.h>
#include <CGAL/gmpxx_coercion_traits.h>


// This file gathers the necessary adaptors so that the following
// C++ number types that come with GMP can be used by CGAL :
// - mpq_class

// Note that GMP++ use the expression template mechanism, which makes things
// a little bit complicated in order to make square(x+y) work for example.
// Reading gmpxx.h shows that ::__gmp_expr<T, T> is the mp[zqf]_class proper,
// while ::__gmp_expr<T, U> is the others "expressions".



CGAL_BEGIN_NAMESPACE

template <>
struct Number_type_traits<mpq_class> {
  typedef Tag_false Has_gcd;
  typedef Tag_true  Has_division;
  typedef Tag_false Has_sqrt;

  typedef Tag_true  Has_exact_ring_operations;
  typedef Tag_true  Has_exact_division;
  typedef Tag_false Has_exact_sqrt;
};

template <>
struct Rational_traits<mpq_class> {
  typedef mpz_class RT;
  RT numerator   (const mpq_class & r) const { return r.get_num(); }
  RT denominator (const mpq_class & r) const { return r.get_den(); }

  mpq_class make_rational(const RT & n, const RT & d) const
  { return mpq_class(n, d); } 
  mpq_class make_rational(const mpq_class & n, const mpq_class & d) const
  { return n / d; } 
};

// AST for mpq_class
template <class U> 
class Algebraic_structure_traits< ::__gmp_expr< ::__gmpq_value,U> >
  : public Algebraic_structure_traits_base< ::__gmp_expr< ::__gmpq_value,U>, 
                                            Null_tag >  {
  public:
    typedef Field_tag           Algebraic_structure_tag;
    typedef Tag_true            Is_exact;
    typedef mpq_class                 Algebraic_structure;

    struct Is_zero: public Unary_function< mpq_class , bool > {
        template <class U2> 
        bool operator()( const ::__gmp_expr< ::__gmpq_value,U2>& x) const {
            return ::sgn(x) == 0;
        }
    };  

    struct Is_one: public Unary_function< mpq_class , bool > {
        template <class U2> 
        bool operator()( const ::__gmp_expr< ::__gmpq_value,U2>& x) const {
            return x == 1;
        }
    }; 

    struct Simplify: public Unary_function< mpq_class , void > {
        template <class U2> 
        void operator()( ::__gmp_expr< ::__gmpq_value,U2>& x) const {
          //TODO: cast x to (mpq_class)??
          x.canonicalize();
        }
    }; 
    
    struct Square: public Unary_function< mpq_class , mpq_class > {
        template <class U2> 
        mpq_class operator()( const ::__gmp_expr< ::__gmpq_value,U2>& x) const {
            return x*x;
        }
    }; 

    struct Unit_part: public Unary_function< mpq_class , mpq_class > {
        template <class U2> 
        mpq_class operator()( const ::__gmp_expr< ::__gmpq_value,U2>& x) const {
            return( x == mpq_class(0)) ? mpq_class(1) : x;
        }
    }; 

    struct Integral_division
        : public Binary_function< mpq_class , mpq_class, mpq_class > {
        template <class U2, class U3> 
        mpq_class operator()( 
                const ::__gmp_expr< ::__gmpq_value,U2>& x,
                const ::__gmp_expr< ::__gmpq_value,U3>& y) const {
            mpq_class result = x / y;
            CGAL_precondition_msg( result * y == x,
            "'x' must be divisible by 'y' in "
            "Algebraic_structure_traits<mpq_class>::Integral_div()(x,y)" );
            return result;         
        }
    }; 
    
    class Is_square
        : public Binary_function< mpq_class, mpq_class&, bool > {
    public:
        template <class U2>
        bool operator()( 
                const ::__gmp_expr< ::__gmpq_value,U2>& x_, 
                mpq_class& y ) const {
            mpq_class x( x_ );
            y = mpq_class (::sqrt( x.get_num() ), ::sqrt( x.get_den() )) ;
            
            return y*y == x;
        }
        
        template <class U2>
        bool operator()( const ::__gmp_expr< ::__gmpq_value,U2>& x ) const {
            mpq_class y;
            return operator()(x,y);
        }
    };   
};

// RET for mpq_class

template < class U > 
class Real_embeddable_traits< ::__gmp_expr< ::__gmpq_value,U> > 
  : public Real_embeddable_traits_base< ::__gmp_expr< ::__gmpq_value,U> > {
  public:
    typedef mpq_class  Real_embeddable;
      
    struct Is_zero: public Unary_function< mpq_class , bool > {
        template <class U2> 
        bool operator()( const ::__gmp_expr< ::__gmpq_value,U2>& x) const {
            return ::sgn(x) == 0;
        }
    }; 
    struct Is_finite: public Unary_function<mpq_class,bool> {
        template <class U2> 
        bool operator()( const ::__gmp_expr< ::__gmpq_value,U2>& x) const {
            return true;
        }
    };

    struct Is_positive: public Unary_function< mpq_class , bool > {
        template <class U2> 
        bool operator()( const ::__gmp_expr< ::__gmpq_value,U2>& x) const {
            return ::sgn(x) > 0;
        }
    }; 
  
    struct Is_negative: public Unary_function< mpq_class , bool > {
        template <class U2> 
        bool operator()( const ::__gmp_expr< ::__gmpq_value,U2>& x) const {
            return ::sgn(x) < 0;
        }
    };

    struct Abs: public Unary_function< mpq_class , mpq_class > {
        template <class U2> 
        mpq_class operator()( const ::__gmp_expr< ::__gmpq_value,U2>& x) const {
            return ::abs(x);
        }
    };
    
    struct Sign 
        : public Unary_function< mpq_class, ::CGAL::Sign > {
    public:
        template <class U2> 
        ::CGAL::Sign 
        operator()( const ::__gmp_expr< ::__gmpq_value,U2>& x ) const {
            return (::CGAL::Sign) ::sgn( x );
        }        
    };
    
    struct Compare 
        : public Binary_function< mpq_class, mpq_class, Comparison_result>
    {
        template <class U2, class U3>
        Comparison_result operator()( 
                const ::__gmp_expr< ::__gmpq_value,U2>& x, 
                const ::__gmp_expr< ::__gmpq_value,U3>& y ) const {
            // cmp returns any int value, not just -1/0/1...
            return (Comparison_result) CGAL_NTS sign( ::cmp(x, y) );
        }
        CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR_WITH_RT
        ( Real_embeddable, Comparison_result);
    };
    
    struct To_double 
        : public Unary_function< mpq_class, double > {
        template <class U2>
        double operator()( const ::__gmp_expr< ::__gmpq_value,U2>& x ) const {
            return mpq_class(x).get_d();
        }
    };
    
    struct To_interval 
    
        : public Unary_function< mpq_class, std::pair< double, double > > {
        template <class U2>
        std::pair<double, double> 
        operator()( const ::__gmp_expr< ::__gmpq_value,U2>& x_ ) const {
            mpq_class x = mpq_class(x_);
            mpfr_t y;
            mpfr_init2 (y, 53); /* Assume IEEE-754 */
            mpfr_set_q (y, x.get_mpq_t(), GMP_RNDD);
            double i = mpfr_get_d (y, GMP_RNDD); /* EXACT but can overflow */
            mpfr_set_q (y, x.get_mpq_t(), GMP_RNDU);
            double s = mpfr_get_d (y, GMP_RNDU); /* EXACT but can overflow */
            mpfr_clear (y);
            return std::pair<double, double>(i, s);
        }
    };
};

CGAL_END_NAMESPACE

#endif // CGAL_MPQ_CLASS_H
