// Copyright (c) 2017
// INRIA Saclay-Ile de France (France),
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Marc Glisse

#ifndef CGAL_BOOST_MP_TYPE_H
#define CGAL_BOOST_MP_TYPE_H

#include <CGAL/config.h>
#include <CGAL/number_utils.h>

// CGAL_USE_BOOST_MP is defined in
// CGAL/Installation/internal/enable_third_party_libraries.h
#if CGAL_USE_BOOST_MP

#include <CGAL/Quotient.h>
#include <CGAL/functional.h> // *ary_function
#include <CGAL/number_type_basic.h>
#include <CGAL/Modular_traits.h>
// We can't just include all Boost.Multiprecision here...
#include <boost/multiprecision/number.hpp>
#include <boost/type_traits/common_type.hpp>
// ... but we kind of have to :-(
#include <boost/multiprecision/cpp_int.hpp>
#ifdef CGAL_USE_GMP
// Same dance as in CGAL/gmp.h
# include <CGAL/disable_warnings.h>
# if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable: 4127 4244 4146 4267) // conversion with loss of data
                                     // warning on - applied on unsigned number
# endif

# include <boost/multiprecision/gmp.hpp>

# if defined(BOOST_MSVC)
#  pragma warning(pop)
# endif

# include <CGAL/enable_warnings.h>
#endif
#ifdef CGAL_USE_MPFR
# include <mpfr.h>
#endif

// TODO: work on the coercions (end of the file)

namespace CGAL {
template<>
struct Needs_parens_as_product<typename boost::multiprecision::cpp_int>{
  bool operator()(const typename boost::multiprecision::cpp_int& x){
    return x < 0;
  }
};

template<>
struct Needs_parens_as_product<typename boost::multiprecision::cpp_rational>{
  bool operator()(const typename boost::multiprecision::cpp_rational& x){
    if (denominator(x) != 1 )
      return true;
    else
      return needs_parens_as_product(numerator(x)) ;
  }

};

#ifdef CGAL_USE_GMP
template<>
struct Needs_parens_as_product<typename boost::multiprecision::mpz_int>{
  bool operator()(const typename boost::multiprecision::mpz_int& x){
    return x < 0;
  }
};

template<>
struct Needs_parens_as_product<typename boost::multiprecision::mpq_rational>{
  bool operator()(const typename boost::multiprecision::mpq_rational& x){
    if (denominator(x) != 1 )
      return true;
    else
      return needs_parens_as_product(numerator(x)) ;
  }
};
#endif



// Algebraic_structure_traits

template <class T, class = boost::mpl::int_<boost::multiprecision::number_category<T>::value> >
struct AST_boost_mp;

template <class NT>
struct AST_boost_mp <NT, boost::mpl::int_<boost::multiprecision::number_kind_integer> >
  : Algebraic_structure_traits_base< NT, Euclidean_ring_tag > {
    typedef NT Type;
    typedef Euclidean_ring_tag Algebraic_category;
    typedef Boolean_tag<std::numeric_limits<Type>::is_exact> Is_exact;
    typedef Tag_false Is_numerical_sensitive;

    typedef INTERN_AST::Is_square_per_sqrt< Type > Is_square;

    struct Is_zero: public CGAL::cpp98::unary_function<Type ,bool> {
        bool operator()( const Type& x) const {
            return x.is_zero();
        }
    };

    struct Div:
        public CGAL::cpp98::binary_function<Type ,Type, Type> {
        template <typename T, typename U>
        Type operator()(const T& x, const U& y) const {
            return x / y;
        }
    };

    struct Mod:
        public CGAL::cpp98::binary_function<Type ,Type, Type> {
        template <typename T, typename U>
        Type operator()(const T& x, const U& y) const {
            return x % y;
        }
    };

    struct Gcd : public CGAL::cpp98::binary_function<Type, Type, Type> {
        template <typename T, typename U>
        Type operator()(const T& x, const U& y) const {
            return boost::multiprecision::gcd(x, y);
        }
    };

    struct Sqrt : public CGAL::cpp98::unary_function<Type, Type> {
        template <typename T>
        Type operator()(const T& x) const {
            return boost::multiprecision::sqrt(x);
        }
    };
};

template <class NT>
struct AST_boost_mp <NT, boost::mpl::int_<boost::multiprecision::number_kind_rational> >
  : public Algebraic_structure_traits_base< NT , Field_tag >  {
  public:
    typedef NT                  Type;
    typedef Field_tag           Algebraic_category;
    typedef Tag_true            Is_exact;
    typedef Tag_false           Is_numerical_sensitive;

    struct Is_zero: public CGAL::cpp98::unary_function<Type ,bool> {
        bool operator()( const Type& x) const {
            return x.is_zero();
        }
    };

    struct Div:
        public CGAL::cpp98::binary_function<Type ,Type, Type> {
        template <typename T, typename U>
        Type operator()(const T& x, const U& y) const {
            return x / y;
        }
    };
};

template <class Backend, boost::multiprecision::expression_template_option Eto>
struct Algebraic_structure_traits<boost::multiprecision::number<Backend, Eto> >
: AST_boost_mp <boost::multiprecision::number<Backend, Eto> > {};
template <class T1,class T2,class T3,class T4,class T5>
struct Algebraic_structure_traits<boost::multiprecision::detail::expression<T1,T2,T3,T4,T5> >
: Algebraic_structure_traits<typename boost::multiprecision::detail::expression<T1,T2,T3,T4,T5>::result_type > {};

// Real_embeddable_traits

namespace Boost_MP_internal {

  // here we know that `intv` contains int64 numbers such that their msb is std::numeric_limits<double>::digits-1
  // TODO: possibly return denormals sometimes...
  inline
  std::pair<double,double> shift_positive_interval( const std::pair<double,double>& intv, const int e ) {
    CGAL_assertion(intv.first > 0.0);
    CGAL_assertion(intv.second > 0.0);

#ifdef CGAL_LITTLE_ENDIAN
    CGAL_assertion_code(
    union {
      struct { uint64_t man:52; uint64_t exp:11; uint64_t sig:1; } s;
      double d;
    } conv;

    conv.d = intv.first;
    )
#else
    //WARNING: untested!
    CGAL_assertion_code(
    union {

      struct { uint64_t sig:1; uint64_t exp:11; uint64_t man:52; } s;
      double d;
    } conv;

    conv.d = intv.first;
    )
#endif
    // Check that the exponent of intv.inf is 52, which corresponds to a 53 bit integer
    CGAL_assertion(conv.s.exp - ((1 << (11 - 1)) - 1) == std::numeric_limits<double>::digits - 1);

    typedef std::numeric_limits<double> limits;

    // warning: min_exponent and max_exponent are 1 more than what the name suggests
    if (e < limits::min_exponent - limits::digits)
      return {0, (limits::min)()};
    if (e > limits::max_exponent - limits::digits)
      return {(limits::max)(), limits::infinity()}; // intv is positive

    const double scale = std::ldexp(1.0, e); // ldexp call is exact
    return { scale * intv.first, scale * intv.second }; // cases that would require a rounding mode have been handled above
  }

  // This function checks if the computed interval is correct and if it is tight.
  template<typename Type>
  bool are_bounds_correct( const double l, const double u, const Type& x ) {
    typedef std::numeric_limits<double> limits;

    const double inf = std::numeric_limits<double>::infinity();
    if ( u!=l && (l==-inf || u==inf
                          || (u==0 && l >=  -(limits::min)())
                          || (l==0 && u <= (limits::min)())) )
    {
      return x >  Type((limits::max)()) ||
             x < Type(-(limits::max)()) ||
             (x > Type(-(limits::min)()) && x < Type((limits::min)()));
    }

    if (!(u == l || u == std::nextafter(l, +inf))) return false;
    //TODO: Type(nextafter(l,inf))>x && Type(nextafter(u,-inf))<x

    const Type lb(l), ub(u);
    const bool are_bounds_respected = (lb <= x && x <= ub);
    return are_bounds_respected;
  }

  // This one returns zero length interval that is inf = sup.
  inline
  std::pair<double, double> get_0ulp_interval( const int shift, const uint64_t p ) {

    const double pp_dbl = static_cast<double>(p);
    const std::pair<double,double> intv(pp_dbl, pp_dbl);

    return shift_positive_interval(intv, -shift);
  }

  // This one returns 1 unit length interval.
  inline
  std::pair<double, double> get_1ulp_interval( const int shift, const uint64_t p ) {

    const double pp_dbl = static_cast<double>(p);
    const double qq_dbl = pp_dbl+1;
    const std::pair<double,double> intv(pp_dbl, qq_dbl);
    return shift_positive_interval(intv, -shift);
  }

  template<typename ET>
  std::pair<double, double> to_interval( ET x, int extra_shift = 0 );

  // This is a version of to_interval that converts a rational type into a
  // double tight interval.
  template<typename Type, typename ET>
  std::pair<double, double> to_interval( ET xnum, ET xden ) {

    CGAL_assertion(!CGAL::is_zero(xden));
    CGAL_assertion_code(const Type input(xnum, xden));
    double l = 0.0, u = 0.0;
    if (CGAL::is_zero(xnum)) { // return [0.0, 0.0]
      CGAL_assertion(are_bounds_correct(l, u, input));
      return std::make_pair(l, u);
    }
    CGAL_assertion(!CGAL::is_zero(xnum));

    // Handle signs.
    bool change_sign = false;
    const bool is_num_pos = CGAL::is_positive(xnum);
    const bool is_den_pos = CGAL::is_positive(xden);
    if (!is_num_pos && !is_den_pos) {
      xnum = -xnum;
      xden = -xden;
    } else if (!is_num_pos && is_den_pos) {
      change_sign = true;
      xnum = -xnum;
    } else if (is_num_pos && !is_den_pos) {
      change_sign = true;
      xden = -xden;
    }
    CGAL_assertion(CGAL::is_positive(xnum) && CGAL::is_positive(xden));

    const int64_t num_dbl_digits = std::numeric_limits<double>::digits - 1;
    const int64_t msb_num = static_cast<int64_t>(boost::multiprecision::msb(xnum));
    const int64_t msb_den = static_cast<int64_t>(boost::multiprecision::msb(xden));

#if 0 // Optimization for the case of input that are double
    // An alternative strategy would be to convert numerator and denominator to
    // intervals, then divide. However, this would require setting the rounding
    // mode (and dividing intervals is not completely free). An important
    // special case is when the rational is exactly equal to a double
    // (fit_in_double). Then the denominator is a power of 2, so we can skip
    // the division and it becomes unnecessary to set the rounding mode, we
    // just need to modify the exponent correction for the denominator.
    if(msb_den == static_cast<int64_t>(lsb(xden))) {
      std::tie(l,u)=to_interval(xnum, msb_den);
      if (change_sign) {
        CGAL_assertion(are_bounds_correct(-u, -l, input));
        return {-u, -l};
      }
      CGAL_assertion(are_bounds_correct(l, u, input));
      return {u, l};
    }
#endif

    const int64_t msb_diff = msb_num - msb_den;
    // Shift so the division result has at least 53 (and at most 54) bits
    int shift = static_cast<int>(num_dbl_digits - msb_diff + 1);
    CGAL_assertion(shift == num_dbl_digits - msb_diff + 1);

    if (shift > 0) {
      xnum <<= +shift;
    } else if (shift < 0) {
      xden <<= -shift;
    }
    CGAL_assertion(num_dbl_digits + 1 ==
      static_cast<int64_t>(boost::multiprecision::msb(xnum)) -
      static_cast<int64_t>(boost::multiprecision::msb(xden)));

    ET p, r;
    boost::multiprecision::divide_qr(xnum, xden, p, r);
    uint64_t uip = static_cast<uint64_t>(p);
    const int64_t p_bits = static_cast<int64_t>(boost::multiprecision::msb(p));
    bool exact = r.is_zero();

    if (p_bits > num_dbl_digits) { // case 54 bits
      exact &= ((uip & 1) == 0);
      uip>>=1;
      --shift;
    }
    std::tie(l, u) = exact ? get_0ulp_interval(shift, uip) : get_1ulp_interval(shift, uip);

    if (change_sign) {
      const double t = l;
      l = -u;
      u = -t;
    }

    CGAL_assertion(are_bounds_correct(l, u, input));
    return std::make_pair(l, u);
  }

  // This is a version of to_interval that converts an integer type into a
  // double tight interval.
  template<typename ET>
  std::pair<double, double> to_interval( ET x, int extra_shift) {

    CGAL_assertion_code(const ET input = x);
    double l = 0.0, u = 0.0;
    if (CGAL::is_zero(x)) { // return [0.0, 0.0]
      CGAL_assertion(are_bounds_correct(l, u, input));
      return std::make_pair(l, u);
    }
    CGAL_assertion(!CGAL::is_zero(x));

    bool change_sign = false;
    const bool is_pos = CGAL::is_positive(x);
    if (!is_pos) {
      change_sign = true;
      x = -x;
    }
    CGAL_assertion(CGAL::is_positive(x));

    const int64_t n = static_cast<int64_t>(boost::multiprecision::msb(x)) + 1;
    const int64_t num_dbl_digits = std::numeric_limits<double>::digits;

    if (n > num_dbl_digits) {
      const int64_t mindig = static_cast<int64_t>(boost::multiprecision::lsb(x));
      int e = static_cast<int>(n - num_dbl_digits);
      x >>= e;
      if (n - mindig > num_dbl_digits)
        std::tie(l, u) = get_1ulp_interval(-e+extra_shift, static_cast<uint64_t>(x));
      else
        std::tie(l, u) = get_0ulp_interval(-e+extra_shift, static_cast<uint64_t>(x));
    } else {
      l = u = extra_shift==0 ? static_cast<double>(static_cast<uint64_t>(x))
                             : std::ldexp(static_cast<double>(static_cast<uint64_t>(x)),-extra_shift);
    }

    if (change_sign) {
      const double t = l;
      l = -u;
      u = -t;
    }

    CGAL_assertion(extra_shift != 0 || are_bounds_correct(l, u, input));
    return std::make_pair(l, u);
  }

} // Boost_MP_internal

template <class NT>
struct RET_boost_mp_base
    : public INTERN_RET::Real_embeddable_traits_base< NT , CGAL::Tag_true > {

    typedef NT Type;

    struct Is_zero: public CGAL::cpp98::unary_function<Type ,bool> {
        bool operator()( const Type& x) const {
            return x.is_zero();
        }
    };

    struct Is_positive: public CGAL::cpp98::unary_function<Type ,bool> {
        bool operator()( const Type& x) const {
            return x.sign() > 0;
        }
    };

    struct Is_negative: public CGAL::cpp98::unary_function<Type ,bool> {
        bool operator()( const Type& x) const {
            return x.sign() < 0;
        }
    };

    struct Abs : public CGAL::cpp98::unary_function<Type, Type> {
        template <typename T>
        Type operator()(const T& x) const {
            return boost::multiprecision::abs(x);
        }
    };

    struct Sgn : public CGAL::cpp98::unary_function<Type, ::CGAL::Sign> {
        ::CGAL::Sign operator()(Type const& x) const {
            return CGAL::sign(x.sign());
        }
    };

    struct Compare
        : public CGAL::cpp98::binary_function<Type, Type, Comparison_result> {
        Comparison_result operator()(const Type& x, const Type& y) const {
            return CGAL::sign(x.compare(y));
        }
    };

    struct To_double
        : public CGAL::cpp98::unary_function<Type, double> {
        double operator()(const Type& x) const {
            return x.template convert_to<double>();
        }
    };

    struct To_interval
        : public CGAL::cpp98::unary_function< Type, std::pair< double, double > > {

        std::pair<double, double>
        operator()(const Type& x) const {

          // See if https://github.com/boostorg/multiprecision/issues/108 suggests anything better
          // assume the conversion is within 1 ulp
          // adding IA::smallest() doesn't work because inf-e=inf, even rounded down.

          // We must use to_nearest here.
          double i;
          const double inf = std::numeric_limits<double>::infinity();
          {
            Protect_FPU_rounding<true> P(CGAL_FE_TONEAREST);
            i = static_cast<double>(x);
            if (i == +inf) {
              return std::make_pair((std::numeric_limits<double>::max)(), i);
            } else if (i == -inf) {
              return std::make_pair(i, std::numeric_limits<double>::lowest());
            }
          }
          double s = i;
          CGAL_assertion(CGAL::abs(i) != inf && CGAL::abs(s) != inf);

          // Throws uncaught exception: Cannot convert a non-finite number to an integer.
          // We can catch it earlier by using the CGAL_assertion() one line above.
          const int cmp = x.compare(i);
          if (cmp > 0) {
            s = nextafter(s, +inf);
            CGAL_assertion(x.compare(s) < 0);
          }
          else if (cmp < 0) {
            i = nextafter(i, -inf);
            CGAL_assertion(x.compare(i) > 0);
          }
          return std::pair<double, double>(i, s);
        }
    };
};

template <class T, class = boost::mpl::int_<boost::multiprecision::number_category<T>::value> >
struct RET_boost_mp;

template <class NT>
struct RET_boost_mp <NT, boost::mpl::int_<boost::multiprecision::number_kind_integer> >
    : RET_boost_mp_base <NT> {
    typedef NT Type;
    struct To_interval
        : public CGAL::cpp98::unary_function< Type, std::pair< double, double > > {

        std::pair<double, double> operator()( const Type& x ) const {
          return Boost_MP_internal::to_interval(x);
        }
    };
};

template <class NT>
struct RET_boost_mp <NT, boost::mpl::int_<boost::multiprecision::number_kind_rational> >
    : RET_boost_mp_base <NT> {
    typedef NT Type;
    struct To_interval
        : public CGAL::cpp98::unary_function< Type, std::pair< double, double > > {

        std::pair<double, double> operator()( const Type& x ) const {
          return Boost_MP_internal::to_interval<Type>(
            boost::multiprecision::numerator(x), boost::multiprecision::denominator(x));
        }
    };
};

#ifdef CGAL_USE_MPFR
// Because of these full specializations, things get instantiated more eagerly. Make them artificially partial if necessary.
template <>
struct RET_boost_mp <boost::multiprecision::mpz_int>
    : RET_boost_mp_base <boost::multiprecision::mpz_int> {
    typedef boost::multiprecision::mpz_int Type;
    struct To_interval
        : public CGAL::cpp98::unary_function< Type, std::pair< double, double > > {
        std::pair<double, double>
        operator()(const Type& x) const {
#if MPFR_VERSION_MAJOR >= 3
          MPFR_DECL_INIT (y, 53); /* Assume IEEE-754 */
          int r = mpfr_set_z (y, x.backend().data(), MPFR_RNDA);
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
          mpfr_set_z (y, x.backend().data(), GMP_RNDD);
          double i = mpfr_get_d (y, GMP_RNDD); /* EXACT but can overflow */
          mpfr_set_z (y, x.backend().data(), GMP_RNDU);
          double s = mpfr_get_d (y, GMP_RNDU); /* EXACT but can overflow */
          mpfr_clear (y);
          return std::pair<double, double>(i, s);
#endif
        }
    };
};
template <>
struct RET_boost_mp <boost::multiprecision::mpq_rational>
    : RET_boost_mp_base <boost::multiprecision::mpq_rational> {
    typedef boost::multiprecision::mpq_rational Type;
    struct To_interval
        : public CGAL::cpp98::unary_function< Type, std::pair< double, double > > {
        std::pair<double, double>
        operator()(const Type& x) const {
# if MPFR_VERSION_MAJOR >= 3
            mpfr_exp_t emin = mpfr_get_emin();
            mpfr_set_emin(-1073);
            MPFR_DECL_INIT (y, 53); /* Assume IEEE-754 */
            int r = mpfr_set_q (y, x.backend().data(), MPFR_RNDA);
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
# else
            mpfr_t y;
            mpfr_init2 (y, 53); /* Assume IEEE-754 */
            mpfr_set_q (y, x.backend().data(), GMP_RNDD);
            double i = mpfr_get_d (y, GMP_RNDD); /* EXACT but can overflow */
            mpfr_set_q (y, x.backend().data(), GMP_RNDU);
            double s = mpfr_get_d (y, GMP_RNDU); /* EXACT but can overflow */
            mpfr_clear (y);
            return std::pair<double, double>(i, s);
# endif
        }
    };
};
#endif

template <class Backend, boost::multiprecision::expression_template_option Eto>
struct Real_embeddable_traits<boost::multiprecision::number<Backend, Eto> >
: RET_boost_mp <boost::multiprecision::number<Backend, Eto> > {};
template <class T1,class T2,class T3,class T4,class T5>
struct Real_embeddable_traits<boost::multiprecision::detail::expression<T1,T2,T3,T4,T5> >
: Real_embeddable_traits<typename boost::multiprecision::detail::expression<T1,T2,T3,T4,T5>::result_type > {};

// Modular_traits

template <class T, class = boost::mpl::int_<boost::multiprecision::number_category<T>::value> >
struct MT_boost_mp {
  typedef T NT;
  typedef ::CGAL::Tag_false Is_modularizable;
  typedef ::CGAL::Null_functor Residue_type;
  typedef ::CGAL::Null_functor Modular_image;
  typedef ::CGAL::Null_functor Modular_image_representative;
};

template <class T>
struct MT_boost_mp <T, boost::mpl::int_<boost::multiprecision::number_kind_integer> > {
  typedef T NT;
  typedef CGAL::Tag_true Is_modularizable;
  typedef Residue Residue_type;

  struct Modular_image{
    Residue_type operator()(const NT& a){
      NT tmp(CGAL::mod(a,NT(Residue::get_current_prime())));
      return CGAL::Residue(tmp.template convert_to<int>());
    }
  };
  struct Modular_image_representative{
    NT operator()(const Residue_type& x){
      return NT(x.get_value());
    }
  };
};

template <class Backend, boost::multiprecision::expression_template_option Eto>
struct Modular_traits<boost::multiprecision::number<Backend, Eto> >
: MT_boost_mp <boost::multiprecision::number<Backend, Eto> > {};
template <class T1,class T2,class T3,class T4,class T5>
struct Modular_traits<boost::multiprecision::detail::expression<T1,T2,T3,T4,T5> >
: Modular_traits<typename boost::multiprecision::detail::expression<T1,T2,T3,T4,T5>::result_type > {};

// Split_double

template <class NT, class = boost::mpl::int_<boost::multiprecision::number_category<NT>::value> >
struct SD_boost_mp {
  void operator()(double d, NT &num, NT &den) const
  {
    num = d;
    den = 1;
  }
};

template <class NT>
struct SD_boost_mp <NT, boost::mpl::int_<boost::multiprecision::number_kind_integer> >
{
  void operator()(double d, NT &num, NT &den) const
  {
    std::pair<double, double> p = split_numerator_denominator(d);
    num = NT(p.first);
    den = NT(p.second);
  }
};

template <class Backend, boost::multiprecision::expression_template_option Eto>
struct Split_double<boost::multiprecision::number<Backend, Eto> >
: SD_boost_mp <boost::multiprecision::number<Backend, Eto> > {};
template <class T1,class T2,class T3,class T4,class T5>
struct Split_double<boost::multiprecision::detail::expression<T1,T2,T3,T4,T5> >
: Split_double<typename boost::multiprecision::detail::expression<T1,T2,T3,T4,T5>::result_type > {};


// Fraction_traits

template <class T, class = boost::mpl::int_<boost::multiprecision::number_category<T>::value> >
struct FT_boost_mp {
  typedef T Type;
  typedef Tag_false Is_fraction;
  typedef Null_tag Numerator_type;
  typedef Null_tag Denominator_type;
  typedef Null_functor Common_factor;
  typedef Null_functor Decompose;
  typedef Null_functor Compose;
};

template <class NT>
struct FT_boost_mp <NT, boost::mpl::int_<boost::multiprecision::number_kind_rational> > {
    typedef NT Type;

    typedef ::CGAL::Tag_true Is_fraction;
    typedef typename boost::multiprecision::component_type<NT>::type Numerator_type;
    typedef Numerator_type Denominator_type;

    typedef typename Algebraic_structure_traits< Numerator_type >::Gcd Common_factor;

    class Decompose {
    public:
        typedef Type first_argument_type;
        typedef Numerator_type& second_argument_type;
        typedef Denominator_type& third_argument_type;
        void operator () (
                const Type& rat,
                Numerator_type& num,
                Denominator_type& den) {
            num = numerator(rat);
            den = denominator(rat);
        }
    };

    class Compose {
    public:
        typedef Numerator_type first_argument_type;
        typedef Denominator_type second_argument_type;
        typedef Type result_type;
        Type operator ()(
                const Numerator_type& num ,
                const Denominator_type& den ) {
            return Type(num, den);
        }
    };
};

template <class Backend, boost::multiprecision::expression_template_option Eto>
struct Fraction_traits<boost::multiprecision::number<Backend, Eto> >
: FT_boost_mp <boost::multiprecision::number<Backend, Eto> > {};
template <class T1,class T2,class T3,class T4,class T5>
struct Fraction_traits<boost::multiprecision::detail::expression<T1,T2,T3,T4,T5> >
: Fraction_traits<typename boost::multiprecision::detail::expression<T1,T2,T3,T4,T5>::result_type > {};


// Coercions

namespace internal { namespace boost_mp { BOOST_MPL_HAS_XXX_TRAIT_DEF(type) } }

template <class B1, boost::multiprecision::expression_template_option E1, class B2, boost::multiprecision::expression_template_option E2>
struct Coercion_traits<boost::multiprecision::number<B1, E1>, boost::multiprecision::number<B2, E2> >
{
  typedef boost::common_type<boost::multiprecision::number<B1, E1>, boost::multiprecision::number<B2, E2> > CT;
  typedef Boolean_tag<internal::boost_mp::has_type<CT>::value> Are_implicit_interoperable;
  // FIXME: the implicit/explicit answers shouldn't be the same...
  typedef Are_implicit_interoperable Are_explicit_interoperable;
  // FIXME: won't compile when they are not interoperable.
  typedef typename CT::type Type;
  struct Cast{
    typedef Type result_type;
    template <class U>
      Type operator()(const U& x) const {
        return Type(x);
      }
  };
};
// Avoid ambiguity with the specialization for <A,A> ...
template <class B1, boost::multiprecision::expression_template_option E1>
struct Coercion_traits<boost::multiprecision::number<B1, E1>, boost::multiprecision::number<B1, E1> >
{
  typedef boost::multiprecision::number<B1, E1> Type;
  typedef Tag_true Are_implicit_interoperable;
  typedef Tag_true Are_explicit_interoperable;
  struct Cast{
    typedef Type result_type;
    template <class U>
      Type operator()(const U& x) const {
        return Type(x);
      }
  };
};

template <class T1, class T2, class T3, class T4, class T5, class U1, class U2, class U3, class U4, class U5>
struct Coercion_traits <
boost::multiprecision::detail::expression<T1,T2,T3,T4,T5>,
boost::multiprecision::detail::expression<U1,U2,U3,U4,U5> >
: Coercion_traits <
typename boost::multiprecision::detail::expression<T1,T2,T3,T4,T5>::result_type,
typename boost::multiprecision::detail::expression<U1,U2,U3,U4,U5>::result_type>
{ };
// Avoid ambiguity with the specialization for <A,A> ...
template <class T1, class T2, class T3, class T4, class T5>
struct Coercion_traits <
boost::multiprecision::detail::expression<T1,T2,T3,T4,T5>,
boost::multiprecision::detail::expression<T1,T2,T3,T4,T5> >
: Coercion_traits <
typename boost::multiprecision::detail::expression<T1,T2,T3,T4,T5>::result_type,
typename boost::multiprecision::detail::expression<T1,T2,T3,T4,T5>::result_type>
{ };

template <class B, boost::multiprecision::expression_template_option E, class T1, class T2, class T3, class T4, class T5>
struct Coercion_traits<boost::multiprecision::number<B, E>, boost::multiprecision::detail::expression<T1,T2,T3,T4,T5> >
: Coercion_traits <
boost::multiprecision::number<B, E>,
typename boost::multiprecision::detail::expression<T1,T2,T3,T4,T5>::result_type>
{ };

template <class B, boost::multiprecision::expression_template_option E, class T1, class T2, class T3, class T4, class T5>
struct Coercion_traits<boost::multiprecision::detail::expression<T1,T2,T3,T4,T5>, boost::multiprecision::number<B, E> >
: Coercion_traits <
typename boost::multiprecision::detail::expression<T1,T2,T3,T4,T5>::result_type,
boost::multiprecision::number<B, E> >
{ };

// TODO: fix existing coercions
// (double -> rational is implicit only for 1.56+, see ticket #10082)
// The real solution would be to avoid specializing Coercion_traits for all pairs of number types and let it auto-detect what works, so only broken types need an explicit specialization.

// Ignore types smaller than long
#define CGAL_COERCE_INT(int) \
template <class B1, boost::multiprecision::expression_template_option E1> \
struct Coercion_traits<boost::multiprecision::number<B1, E1>, int> { \
  typedef boost::multiprecision::number<B1, E1> Type; \
  typedef Tag_true Are_implicit_interoperable; \
  typedef Tag_true Are_explicit_interoperable; \
  struct Cast{ \
    typedef Type result_type; \
    template <class U> Type operator()(const U& x) const { return Type(x); } \
  }; \
}; \
template <class B1, boost::multiprecision::expression_template_option E1> \
struct Coercion_traits<int, boost::multiprecision::number<B1, E1> > \
: Coercion_traits<boost::multiprecision::number<B1, E1>, int> {}; \
template <class T1, class T2, class T3, class T4, class T5> \
struct Coercion_traits<boost::multiprecision::detail::expression<T1,T2,T3,T4,T5>, int> \
: Coercion_traits<typename boost::multiprecision::detail::expression<T1,T2,T3,T4,T5>::result_type, int>{}; \
template <class T1, class T2, class T3, class T4, class T5> \
struct Coercion_traits<int, boost::multiprecision::detail::expression<T1,T2,T3,T4,T5> > \
: Coercion_traits<typename boost::multiprecision::detail::expression<T1,T2,T3,T4,T5>::result_type, int>{}

CGAL_COERCE_INT(short);
CGAL_COERCE_INT(int);
CGAL_COERCE_INT(long);
#undef CGAL_COERCE_INT

// Ignore bounded-precision rationals
#define CGAL_COERCE_FLOAT(float) \
template <class B1, boost::multiprecision::expression_template_option E1> \
struct Coercion_traits<boost::multiprecision::number<B1, E1>, float> { \
  typedef boost::multiprecision::number<B1, E1> Type; \
  typedef Boolean_tag<boost::multiprecision::number_category<Type>::value != boost::multiprecision::number_kind_integer> Are_implicit_interoperable; \
  typedef Are_implicit_interoperable Are_explicit_interoperable; \
  struct Cast{ \
    typedef Type result_type; \
    template <class U> Type operator()(const U& x) const { return Type(x); } \
  }; \
}; \
template <class B1, boost::multiprecision::expression_template_option E1> \
struct Coercion_traits<float, boost::multiprecision::number<B1, E1> > \
: Coercion_traits<boost::multiprecision::number<B1, E1>, float> {}; \
template <class T1, class T2, class T3, class T4, class T5> \
struct Coercion_traits<boost::multiprecision::detail::expression<T1,T2,T3,T4,T5>, float> \
: Coercion_traits<typename boost::multiprecision::detail::expression<T1,T2,T3,T4,T5>::result_type, float>{}; \
template <class T1, class T2, class T3, class T4, class T5> \
struct Coercion_traits<float, boost::multiprecision::detail::expression<T1,T2,T3,T4,T5> > \
: Coercion_traits<typename boost::multiprecision::detail::expression<T1,T2,T3,T4,T5>::result_type, float>{}

CGAL_COERCE_FLOAT(float);
CGAL_COERCE_FLOAT(double);
#undef CGAL_COERCE_FLOAT

// Because of https://github.com/boostorg/multiprecision/issues/29 , this is not perfect and fails to read some KDS files.

template <>
class Input_rep<boost::multiprecision::cpp_rational> : public IO_rep_is_specialized {
    boost::multiprecision::cpp_rational& q;
public:
    Input_rep(boost::multiprecision::cpp_rational& qq) : q(qq) {}
    std::istream& operator()(std::istream& in) const {
      internal::read_float_or_quotient<boost::multiprecision::cpp_int,boost::multiprecision::cpp_rational>(in, q);
      return in;
    }
};
#ifdef CGAL_USE_GMP
template <>
class Input_rep<boost::multiprecision::mpq_rational> : public IO_rep_is_specialized {
    boost::multiprecision::mpq_rational& q;
public:
    Input_rep(boost::multiprecision::mpq_rational& qq) : q(qq) {}
    std::istream& operator()(std::istream& in) const {
      internal::read_float_or_quotient<boost::multiprecision::mpz_int,boost::multiprecision::mpq_rational>(in, q);
      return in;
    }
};
#endif

// Copied from leda_rational.h
namespace internal {
  // See: Stream_support/include/CGAL/IO/io.h
  template <typename ET>
  void read_float_or_quotient(std::istream & is, ET& et);

  template <>
  inline void read_float_or_quotient(std::istream & is, boost::multiprecision::cpp_rational& et)
  {
    internal::read_float_or_quotient<boost::multiprecision::cpp_int,boost::multiprecision::cpp_rational>(is, et);
  }
#ifdef CGAL_USE_GMP
  template <>
  inline void read_float_or_quotient(std::istream & is, boost::multiprecision::mpq_rational& et)
  {
    internal::read_float_or_quotient<boost::multiprecision::mpz_int,boost::multiprecision::mpq_rational>(is, et);
  }
#endif
} // namespace internal

#ifdef CGAL_USE_BOOST_MP

template< > class Real_embeddable_traits< Quotient<boost::multiprecision::cpp_int> >
    : public INTERN_QUOTIENT::Real_embeddable_traits_quotient_base< Quotient<boost::multiprecision::cpp_int> > {

  public:
    typedef Quotient<boost::multiprecision::cpp_int> Type;

    class To_interval
      : public CGAL::cpp98::unary_function< Type, std::pair< double, double > > {
      public:
        std::pair<double, double> operator()( const Type& x ) const {
          return Boost_MP_internal::to_interval<Type>(x.num, x.den);
        }
    };
};

#endif // CGAL_USE_BOOST_MP

} //namespace CGAL

namespace Eigen {
  template<class> struct NumTraits;
  template<> struct NumTraits<boost::multiprecision::cpp_int>
  {
    typedef boost::multiprecision::cpp_int Real;
    typedef boost::multiprecision::cpp_rational NonInteger;
    typedef boost::multiprecision::cpp_int Nested;
    typedef boost::multiprecision::cpp_int Literal;

    static inline Real epsilon() { return 0; }
    static inline Real dummy_precision() { return 0; }

    enum {
      IsInteger = 1,
      IsSigned = 1,
      IsComplex = 0,
      RequireInitialization = 1,
      ReadCost = 6,
      AddCost = 30,
      MulCost = 50
    };
  };

  template<> struct NumTraits<boost::multiprecision::cpp_rational>
  {
    typedef boost::multiprecision::cpp_rational Real;
    typedef boost::multiprecision::cpp_rational NonInteger;
    typedef boost::multiprecision::cpp_rational Nested;
    typedef boost::multiprecision::cpp_rational Literal;

    static inline Real epsilon() { return 0; }
    static inline Real dummy_precision() { return 0; }

    enum {
      IsInteger = 0,
      IsSigned = 1,
      IsComplex = 0,
      RequireInitialization = 1,
      ReadCost = 6,
      AddCost = 150,
      MulCost = 100
    };
  };

#ifdef CGAL_USE_GMP

  template<> struct NumTraits<boost::multiprecision::mpz_int>
  {
    typedef boost::multiprecision::mpz_int Real;
    typedef boost::multiprecision::mpq_rational NonInteger;
    typedef boost::multiprecision::mpz_int Nested;
    typedef boost::multiprecision::mpz_int Literal;

    static inline Real epsilon() { return 0; }
    static inline Real dummy_precision() { return 0; }

    enum {
      IsInteger = 1,
      IsSigned = 1,
      IsComplex = 0,
      RequireInitialization = 1,
      ReadCost = 6,
      AddCost = 30,
      MulCost = 50
    };
  };

  template<> struct NumTraits<boost::multiprecision::mpq_rational>
  {
    typedef boost::multiprecision::mpq_rational Real;
    typedef boost::multiprecision::mpq_rational NonInteger;
    typedef boost::multiprecision::mpq_rational Nested;
    typedef boost::multiprecision::mpq_rational Literal;

    static inline Real epsilon() { return 0; }
    static inline Real dummy_precision() { return 0; }

    enum {
      IsInteger = 0,
      IsSigned = 1,
      IsComplex = 0,
      RequireInitialization = 1,
      ReadCost = 6,
      AddCost = 150,
      MulCost = 100
    };
  };
#endif // CGAL_USE_GMP


} // namespace Eigen

#endif // BOOST_VERSION
#endif
