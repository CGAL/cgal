// Copyright (c) 2002,2003,2007
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sylvain Pion, Michael Hemmer

#ifndef CGAL_MPZ_CLASS_H
#define CGAL_MPZ_CLASS_H

// This file gathers the necessary adaptors so that the following
// C++ number types that come with GMP can be used by CGAL :
// - mpz_class

// Note that GMP++ use the expression template mechanism, which makes things
// a little bit complicated in order to make square(x+y) work for example.
// Reading gmpxx.h shows that ::__gmp_expr<T, T> is the mp[zqf]_class proper,
// while ::__gmp_expr<T, U> is the others "expressions".

#include <CGAL/number_type_config.h>
#include <CGAL/functional.h>
#include <CGAL/Algebraic_structure_traits.h>
#include <CGAL/Real_embeddable_traits.h>
#include <CGAL/number_utils.h>
#include <CGAL/double.h>
#include <boost/type_traits/is_same.hpp>
#include <mpfr.h>
#include <gmpxx.h>

#define CGAL_CHECK_GMP_EXPR                                             \
    CGAL_static_assertion(                                                \
            (::boost::is_same< ::__gmp_expr< T , T >,Type>::value ));

namespace CGAL {

// AST for mpz_class

template<>
class Algebraic_structure_traits< mpz_class >
  :public Algebraic_structure_traits_base<  mpz_class , Euclidean_ring_tag > {

public:
    typedef mpz_class           Type;
    typedef Euclidean_ring_tag  Algebraic_category;
    typedef Tag_true            Is_exact;
    typedef Tag_false           Is_numerical_sensitive;

    struct Is_zero: public CGAL::cpp98::unary_function< mpz_class , bool > {
        template <typename T, typename U>
        bool operator()( const ::__gmp_expr< T , U >& x) const {
            CGAL_CHECK_GMP_EXPR;
            return ::sgn(x) == 0;
        }
    };

    struct Is_one: public CGAL::cpp98::unary_function< mpz_class , bool > {
        template <typename T, typename U>
        bool operator()( const ::__gmp_expr< T , U >& x) const {
            CGAL_CHECK_GMP_EXPR;
            return x == 1;
        }
    };

    struct Simplify: public CGAL::cpp98::unary_function< mpz_class , void > {
        template <class T, class U>
        void operator()( const ::__gmp_expr< T ,U >&) const {
            CGAL_CHECK_GMP_EXPR;
        }
    };

    struct Square: public CGAL::cpp98::unary_function< mpz_class , mpz_class > {
        mpz_class operator()( const mpz_class& x) const {
            return x*x;
        }
    };

    struct Unit_part: public CGAL::cpp98::unary_function< mpz_class , mpz_class > {
        mpz_class operator()( const mpz_class& x) const {
            return (x < 0) ? -1 : 1;
        }
    };



    struct Integral_division:
        public CGAL::cpp98::binary_function< mpz_class , mpz_class, mpz_class > {
        template <typename T,  typename U1, typename U2>
        mpz_class operator()(
                const ::__gmp_expr< T , U1 >& x,
                const ::__gmp_expr< T , U2 >& y) const {
            CGAL_CHECK_GMP_EXPR;
            mpz_class result = x / y;
            CGAL_precondition_msg( result * y == x,
            "'x' must be divisible by 'y' in "
            "Algebraic_structure_traits<mpz_class>::Integral_div()(x,y)" );
            return result;
        }
        CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR( Type )
    };

    struct Gcd : public CGAL::cpp98::binary_function< mpz_class, mpz_class, mpz_class > {
        mpz_class operator()(
                const mpz_class& x,
                const mpz_class& y) const {
            mpz_class c;
            mpz_gcd(c.get_mpz_t(),x.get_mpz_t(), y.get_mpz_t() );
            return c;
        }
        CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR( Type )
    };

    struct Div : public CGAL::cpp98::binary_function< mpz_class, mpz_class, mpz_class > {
        template <typename T,  typename U1, typename U2>
        mpz_class operator()(
                const ::__gmp_expr< T , U1 >& x,
                const ::__gmp_expr< T , U2 >& y) const {
            CGAL_CHECK_GMP_EXPR;
            return x / y;
        }
        CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR( Type )
    };

    struct Mod : public CGAL::cpp98::binary_function< mpz_class, mpz_class, mpz_class > {
        template <typename T,  typename U1, typename U2>
        mpz_class operator()(
                const ::__gmp_expr< T , U1 >& x,
                const ::__gmp_expr< T , U2 >& y) const {
            CGAL_CHECK_GMP_EXPR;
            return x % y;
        }
        CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR( Type )
    };
    struct Div_mod {
        typedef mpz_class    first_argument_type;
        typedef mpz_class    second_argument_type;
        typedef mpz_class&   third_argument_type;
        typedef mpz_class&   fourth_argument_type;
        typedef void         result_type;
        void operator()(
                const mpz_class& x,
                const mpz_class& y,
                mpz_class& q,
                mpz_class& r
        ) const {
            typedef Algebraic_structure_traits<mpz_class> Traits;
                Traits::Div  actual_div;
                Traits::Mod  actual_mod;
                q = actual_div( x, y );
                r = actual_mod( x, y );
                // use mpz_tdiv_qr to do both at once
                return;
            }
    };


    struct Sqrt: public CGAL::cpp98::unary_function< mpz_class , mpz_class > {
        template <typename T, typename U>
        mpz_class operator()( const ::__gmp_expr< T , U >& x) const {
            CGAL_CHECK_GMP_EXPR;
            return ::sqrt(x);
        }
    };


    /*struct Is_square: public CGAL::cpp98::binary_function< mpz_class , mpz_class& , bool > {
        template <typename T, typename U>
        bool operator()(
                const ::__gmp_expr< T , U >& x,
                mpz_class&                              r){
            r = ::sqrt(x);
            return (r*r==x) ? true : false;
        }
        template <typename T, typename U>
        bool operator()(const ::__gmp_expr< T , U >& x){
            mpz_class r = ::sqrt(x);
            return (r*r==x) ? true : false;
        }
        };*/
};

} //namespace CGAL

#include <CGAL/int.h> // for `sign( ::cmp(x, y) )`, below

namespace CGAL {

// RET for mpz_class
template<>
class Real_embeddable_traits< mpz_class  >
    : public INTERN_RET::Real_embeddable_traits_base< mpz_class , CGAL::Tag_true > {
public:

    struct Is_zero: public CGAL::cpp98::unary_function< mpz_class , bool > {
        template <typename T, typename U>
        bool operator()( const ::__gmp_expr< T , U >& x) const {
            CGAL_CHECK_GMP_EXPR;
            return ::sgn(x) == 0;
        }
    };
    struct Is_finite: public CGAL::cpp98::unary_function<mpz_class,bool> {
        template <typename T, typename U>
        bool operator()( const ::__gmp_expr< T , U >& ) const {
            CGAL_CHECK_GMP_EXPR;
            return true;
        }
    };

    struct Is_positive: public CGAL::cpp98::unary_function< mpz_class , bool > {
        template <typename T, typename U>
        bool operator()( const ::__gmp_expr< T , U >& x) const {
            CGAL_CHECK_GMP_EXPR;
            return ::sgn(x) > 0;
        }
    };

    struct Is_negative: public CGAL::cpp98::unary_function< mpz_class , bool > {
        template <typename T, typename U>
        bool operator()( const ::__gmp_expr< T , U >& x) const {
            CGAL_CHECK_GMP_EXPR;
            return ::sgn(x) < 0;
        }
    };

    struct Abs: public CGAL::cpp98::unary_function< mpz_class , mpz_class > {
        template <typename T, typename U>
        mpz_class operator()( const ::__gmp_expr< T , U >& x) const {
            CGAL_CHECK_GMP_EXPR;
            return ::abs(x);
        }
    };

    struct Sgn
        : public CGAL::cpp98::unary_function< mpz_class, ::CGAL::Sign > {
    public:
        template <typename T, typename U>
        ::CGAL::Sign operator()( const ::__gmp_expr< T , U >& x ) const {
            CGAL_CHECK_GMP_EXPR;
            return (::CGAL::Sign) ::sgn( x );
        }
    };

    struct Compare
        : public CGAL::cpp98::binary_function< mpz_class, mpz_class, Comparison_result > {
        template <typename T,  typename U1, typename U2>
        Comparison_result operator()(
                const ::__gmp_expr< T , U1 >& x,
                const ::__gmp_expr< T , U2 >& y ) const {
            CGAL_CHECK_GMP_EXPR;
            // cmp returns any int value, not just -1/0/1...
            return (Comparison_result) CGAL_NTS sign( ::cmp(x, y) );
        }
        CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR_WITH_RT
        ( Type, Comparison_result)
    };

    struct To_double
        : public CGAL::cpp98::unary_function< mpz_class, double > {
        double operator()( const mpz_class& x ) const {
            return x.get_d();
        }
    };

    struct To_interval

        : public CGAL::cpp98::unary_function< mpz_class, std::pair< double, double > > {
        std::pair<double, double>
        operator()( const mpz_class& x ) const {
#if MPFR_VERSION_MAJOR >= 3
          MPFR_DECL_INIT (y, 53); /* Assume IEEE-754 */
          int r = mpfr_set_z (y, x.get_mpz_t(), MPFR_RNDA);
          double i = mpfr_get_d (y, MPFR_RNDA); /* EXACT but can overflow */
          if (r == 0 && is_finite (i))
            return std::pair<double, double>(i, i);
          else
            {
              double s = nextafter (i, 0);
              if (i < 0)
                return std::pair<double, double>(i, s);
              else
                return std::pair<double, double>(s, i);
            }
#else
          mpfr_t y;
          mpfr_init2 (y, 53); /* Assume IEEE-754 */
          mpfr_set_z (y, x.get_mpz_t(), GMP_RNDD);
          double i = mpfr_get_d (y, GMP_RNDD); /* EXACT but can overflow */
          mpfr_set_z (y, x.get_mpz_t(), GMP_RNDU);
          double s = mpfr_get_d (y, GMP_RNDU); /* EXACT but can overflow */
          mpfr_clear (y);
          return std::pair<double, double>(i, s);
#endif
        }
    };
};

} //namespace CGAL

#include <CGAL/gmpxx.h>
#include <CGAL/gmpxx_coercion_traits.h>
#include <CGAL/Residue.h>
#include <CGAL/Modular_traits.h>

namespace CGAL {

/*! \ingroup NiX_Modular_traits_spec
 *  \brief a model of concept ModularTraits,
 *  specialization of NiX::Modular_traits.
 */

template<>
class Modular_traits< mpz_class > {
public:
  typedef mpz_class NT;
  typedef CGAL::Tag_true Is_modularizable;
  typedef Residue Residue_type;

  struct Modular_image{
    Residue_type operator()(const mpz_class& a){
      NT tmp(CGAL::mod(a,NT(Residue::get_current_prime())));
      return CGAL::Residue(int(mpz_get_si(tmp.get_mpz_t())));
    }
  };
  struct Modular_image_representative{
    NT operator()(const Residue_type& x){
      return NT(x.get_value());
    }
  };
};
} //namespace CGAL

#include <CGAL/Quotient.h>

namespace CGAL {
template <>
struct Split_double<mpz_class>
{
  void operator()(double d, mpz_class &num, mpz_class &den) const
  {
    std::pair<double, double> p = split_numerator_denominator(d);
    num = mpz_class(p.first);
    den = mpz_class(p.second);
  }
};

} //namespace CGAL

#undef CGAL_CHECK_GMP_EXPR

#endif // CGAL_MPZ_CLASS_H
