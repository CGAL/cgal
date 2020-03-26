// Copyright (c) 2002,2003
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

#ifndef CGAL_MPQ_CLASS_H
#define CGAL_MPQ_CLASS_H

#include <CGAL/Algebraic_structure_traits.h>
#include <CGAL/Real_embeddable_traits.h>
#include <CGAL/number_utils.h>
#include <CGAL/double.h>
#include <CGAL/IO/io.h>
#include <mpfr.h>
#include <gmpxx.h>

// This file gathers the necessary adaptors so that the following
// C++ number types that come with GMP can be used by CGAL :
// - mpq_class

// Note that GMP++ use the expression template mechanism, which makes things
// a little bit complicated in order to make square(x+y) work for example.
// Reading gmpxx.h shows that ::__gmp_expr<T, T> is the mp[zqf]_class proper,
// while ::__gmp_expr<T, U> is the others "expressions".

#define CGAL_CHECK_GMP_EXPR_MPQ_CLASS                                     \
    CGAL_static_assertion(                                                \
            (::boost::is_same< ::__gmp_expr< T , T >,Type>::value ));

namespace CGAL {


// AST for mpq_class
template<>
class Algebraic_structure_traits< mpq_class  >
  : public Algebraic_structure_traits_base< mpq_class , Field_tag >  {
  public:
    typedef mpq_class           Type;
    typedef Field_tag           Algebraic_category;
    typedef Tag_true            Is_exact;
    typedef Tag_false           Is_numerical_sensitive;

    struct Is_zero: public CGAL::cpp98::unary_function< mpq_class , bool > {
        template <class T, class U>
        bool operator()( const ::__gmp_expr< T , U >& x) const {
            CGAL_CHECK_GMP_EXPR_MPQ_CLASS;
            return ::sgn(x) == 0;
        }
    };

    struct Is_one: public CGAL::cpp98::unary_function< mpq_class , bool > {
        template <typename T, typename U>
        bool operator()( const ::__gmp_expr< T , U >& x) const {
            CGAL_CHECK_GMP_EXPR_MPQ_CLASS;
            return x == 1;
        }
    };

    struct Simplify: public CGAL::cpp98::unary_function< mpq_class , void > {
        void operator()( mpq_class& x) const {
            // do nothing because x is already canonical?
            x.canonicalize();
        }
    };

    struct Square: public CGAL::cpp98::unary_function< mpq_class , mpq_class > {
        mpq_class operator()( const mpq_class& x) const {
            return x*x;
        }
    };

    struct Unit_part: public CGAL::cpp98::unary_function< mpq_class , mpq_class > {
        mpq_class operator()( const mpq_class& x) const {
            return( x == 0) ? mpq_class(1) : x;
        }
    };

    struct Integral_division
        : public CGAL::cpp98::binary_function< mpq_class , mpq_class, mpq_class > {
        template <typename T,  typename U1, typename U2>
        mpq_class operator()(
                const ::__gmp_expr< T , U1 >& x,
                const ::__gmp_expr< T , U2 > & y) const {
            CGAL_CHECK_GMP_EXPR_MPQ_CLASS;
            mpq_class result = x / y;
            CGAL_precondition_msg( result * y == x,
            "'x' must be divisible by 'y' in "
            "Algebraic_structure_traits<mpq_class>::Integral_div()(x,y)" );
            return result;
        }
        CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR( Type )
    };

    class Is_square
        : public CGAL::cpp98::binary_function< mpq_class, mpq_class&, bool > {
    public:
        bool operator()( const mpq_class& x, mpq_class& y ) const {
            y = mpq_class (::sqrt( x.get_num() ), ::sqrt( x.get_den() )) ;
            return y*y == x;
            // for efficiency, only handle den if num is a square
        }

        bool operator()( const mpq_class& x ) const {
            mpq_class y;
            return operator()(x,y);
        }
    };
};

// RET for mpq_class

template < >
class Real_embeddable_traits< mpq_class >
  : public INTERN_RET::Real_embeddable_traits_base< mpq_class , CGAL::Tag_true > {
  public:

    struct Is_zero: public CGAL::cpp98::unary_function< mpq_class , bool > {
        template <typename T, typename U>
        bool operator()( const ::__gmp_expr< T , U >& x) const {
            CGAL_CHECK_GMP_EXPR_MPQ_CLASS;
            return ::sgn(x) == 0;
        }
    };
    struct Is_finite: public CGAL::cpp98::unary_function<mpq_class,bool> {
        template <typename T, typename U>
        bool operator()( const ::__gmp_expr< T , U >&) const {
            CGAL_CHECK_GMP_EXPR_MPQ_CLASS;
            return true;
        }
    };

    struct Is_positive: public CGAL::cpp98::unary_function< mpq_class , bool > {
        template <typename T, typename U>
        bool operator()( const ::__gmp_expr< T , U >& x) const {
            CGAL_CHECK_GMP_EXPR_MPQ_CLASS;
            return ::sgn(x) > 0;
        }
    };

    struct Is_negative: public CGAL::cpp98::unary_function< mpq_class , bool > {
        template <typename T, typename U>
        bool operator()( const ::__gmp_expr< T , U >& x) const {
            CGAL_CHECK_GMP_EXPR_MPQ_CLASS;
            return ::sgn(x) < 0;
        }
    };

    struct Abs: public CGAL::cpp98::unary_function< mpq_class , mpq_class > {
        template <typename T, typename U>
        mpq_class operator()( const ::__gmp_expr< T , U >& x) const {
            CGAL_CHECK_GMP_EXPR_MPQ_CLASS;
            return ::abs(x);
        }
    };

    struct Sgn
        : public CGAL::cpp98::unary_function< mpq_class, ::CGAL::Sign > {
    public:
        template <typename T, typename U>
        ::CGAL::Sign
        operator()( const ::__gmp_expr< T , U >& x ) const {
            CGAL_CHECK_GMP_EXPR_MPQ_CLASS;
            return (::CGAL::Sign) ::sgn( x );
        }
    };

    struct Compare
        : public CGAL::cpp98::binary_function< mpq_class, mpq_class, Comparison_result>
    {
        template <typename T, typename U1, typename U2>
        Comparison_result operator()(
                const ::__gmp_expr< T , U1 >& x,
                const ::__gmp_expr< T , U2 >& y ) const {
            CGAL_CHECK_GMP_EXPR_MPQ_CLASS;
            // cmp returns any int value, not just -1/0/1...
            return (Comparison_result) CGAL_NTS sign( ::cmp(x, y) );
        }
        CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR_WITH_RT
        ( Type, Comparison_result)
    };

    struct To_double
        : public CGAL::cpp98::unary_function< mpq_class, double > {
        double operator()( const mpq_class& x ) const {
            return x.get_d();
        }
    };

    struct To_interval

        : public CGAL::cpp98::unary_function< mpq_class, std::pair< double, double > > {
        std::pair<double, double>
        operator()( const mpq_class& x ) const {
#if MPFR_VERSION_MAJOR >= 3
          mpfr_exp_t emin = mpfr_get_emin();
          mpfr_set_emin(-1073);
          MPFR_DECL_INIT (y, 53); /* Assume IEEE-754 */
          int r = mpfr_set_q (y, x.get_mpq_t(), MPFR_RNDA);
          r = mpfr_subnormalize (y, r, MPFR_RNDA); /* Round subnormals */
          double i = mpfr_get_d (y, MPFR_RNDA); /* EXACT but can overflow */
          mpfr_set_emin(emin); /* Restore old value, users may care */
          // With mpfr_set_emax(1024) we could drop the is_finite test
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
          mpfr_set_q (y, x.get_mpq_t(), GMP_RNDD);
          double i = mpfr_get_d (y, GMP_RNDD); /* EXACT but can overflow */
          mpfr_set_q (y, x.get_mpq_t(), GMP_RNDU);
          double s = mpfr_get_d (y, GMP_RNDU); /* EXACT but can overflow */
          mpfr_clear (y);
          return std::pair<double, double>(i, s);
#endif
        }
    };
};

} // namespace CGAL

#include <CGAL/gmpxx_coercion_traits.h>
#include <CGAL/mpz_class.h> // for GCD in Type traits
#include <CGAL/Fraction_traits.h>

namespace CGAL {

/*! \ingroup NiX_Fraction_traits_spec
 *  \brief Specialization of Fraction_traits for mpq_class
 */
template <>
class Fraction_traits< mpq_class > {
public:
    typedef mpq_class Type;

    typedef ::CGAL::Tag_true Is_fraction;
    typedef mpz_class Numerator_type;
    typedef mpz_class Denominator_type;

    typedef Algebraic_structure_traits< mpz_class >::Gcd Common_factor;

    class Decompose {
    public:
        typedef mpq_class first_argument_type;
        typedef mpz_class& second_argument_type;
        typedef mpz_class& third_argument_type;
        void operator () (
                const mpq_class& rat,
                mpz_class& num,
                mpz_class& den) {
            num = rat.get_num();
            den = rat.get_den();
        }
    };

    class Compose {
    public:
        typedef mpz_class first_argument_type;
        typedef mpz_class second_argument_type;
        typedef mpq_class result_type;
        mpq_class operator ()(
                const mpz_class& num ,
                const mpz_class& den ) {
            mpq_class result(num, den);
            result.canonicalize();
            return result;
        }
    };
};

template <>
class Input_rep<mpq_class> : public IO_rep_is_specialized {
    mpq_class& q;
public:
    Input_rep( mpq_class& qq) : q(qq) {}
    std::istream& operator()( std::istream& in) const {
      internal::read_float_or_quotient<mpz_class,mpq_class>(in, q);
      return in;
    }
};

// Copied from leda_rational.h
namespace internal {
  // See: Stream_support/include/CGAL/IO/io.h
  template <typename ET>
  void read_float_or_quotient(std::istream & is, ET& et);

  template <>
  inline void read_float_or_quotient(std::istream & is, mpq_class& et)
  {
    internal::read_float_or_quotient<mpz_class,mpq_class>(is, et);
  }
} // namespace internal

} //namespace CGAL

#undef CGAL_CHECK_GMP_EXPR_MPQ_CLASS

#endif // CGAL_MPQ_CLASS_H
