//
//  Copyright (c) 2000-2002
//  Joerg Walter, Mathias Koch
//
//  Permission to use, copy, modify, distribute and sell this software
//  and its documentation for any purpose is hereby granted without fee,
//  provided that the above copyright notice appear in all copies and
//  that both that copyright notice and this permission notice appear
//  in supporting documentation.  The authors make no representations
//  about the suitability of this software for any purpose.
//  It is provided "as is" without express or implied warranty.
//
//  The authors gratefully acknowledge the support of
//  GeNeSys mbH & Co. KG in producing this work.
//

#ifndef BOOST_UBLAS_TRAITS_H
#define BOOST_UBLAS_TRAITS_H

#include <algorithm>
#include <iterator>
#include <complex>
#include <cmath>

#include <boost/numeric/ublas/iterator.hpp>
#if !defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION) && !defined(BOOST_NO_SFINAE)
#include <boost/numeric/ublas/returntype_deduction.hpp>
#endif

namespace boost { namespace numeric { namespace ublas {

    template<class T>
    struct type_traits {
        typedef type_traits<T> self_type;
        typedef T value_type;
        typedef const T &const_reference;
        typedef T &reference;

        /*
         * Don't define unknown properties
         * 
        typedef T real_type;
        typedef T precision_type;

        BOOST_STATIC_CONSTANT (unsigned, plus_complexity = 0);
        BOOST_STATIC_CONSTANT (unsigned, multiplies_complexity = 0);

        static
        BOOST_UBLAS_INLINE
        real_type real (const_reference) {
        }
        static
        BOOST_UBLAS_INLINE
        real_type imag (const_reference) {
        }
        static
        BOOST_UBLAS_INLINE
        value_type conj (const_reference) {
        }

        static
        BOOST_UBLAS_INLINE
        real_type abs (const_reference) {
        }
        static
        BOOST_UBLAS_INLINE
        value_type sqrt (const_reference) {
        }

        static
        BOOST_UBLAS_INLINE
        real_type norm_1 (const_reference t) {
        }
        static
        BOOST_UBLAS_INLINE
        real_type norm_2 (const_reference t) {
        }
        static
        BOOST_UBLAS_INLINE
        real_type norm_inf (const_reference t) {
        }

        static
        BOOST_UBLAS_INLINE
        bool equals (const_reference t1, const_reference t2) {
        }
        */
        // Dummy definition for compilers that error if undefined even though it is never used
#ifdef BOOST_NO_SFINAE
        typedef void real_type;
        typedef void precision_type;
#endif
    };

    template<>
    struct type_traits<float> {
        typedef type_traits<float> self_type;
        typedef float value_type;
        typedef const value_type &const_reference;
        typedef value_type &reference;
        typedef value_type real_type;
        typedef double precision_type;

        BOOST_STATIC_CONSTANT (unsigned, plus_complexity = 1);
        BOOST_STATIC_CONSTANT (unsigned, multiplies_complexity = 1);

        static
        BOOST_UBLAS_INLINE
        real_type real (const_reference t) {
                return t;
        }
        static
        BOOST_UBLAS_INLINE
        real_type imag (const_reference /*t*/) {
                return 0;
        }
        static
        BOOST_UBLAS_INLINE
        value_type conj (const_reference t) {
                return t;
        }

        static
        BOOST_UBLAS_INLINE
        real_type abs (const_reference t) {
#if defined (BOOST_NO_STDC_NAMESPACE) || defined (BOOST_UBLAS_CMATH_BAD_STD)
            return ::fabsf (t);
#else
            return std::abs (t);
#endif
        }
        static
        BOOST_UBLAS_INLINE
        value_type sqrt (const_reference t) {
#if defined (BOOST_NO_STDC_NAMESPACE) || defined (BOOST_UBLAS_CMATH_BAD_STD)
            return ::sqrtf (t);
#else
            return std::sqrt (t);
#endif
        }

        static
        BOOST_UBLAS_INLINE
        real_type norm_1 (const_reference t) {
            return self_type::abs (t);
        }
        static
        BOOST_UBLAS_INLINE
        real_type norm_2 (const_reference t) {
            return self_type::abs (t);
        }
        static
        BOOST_UBLAS_INLINE
        real_type norm_inf (const_reference t) {
            return self_type::abs (t);
        }

        static
        BOOST_UBLAS_INLINE
        bool equals (const_reference t1, const_reference t2) {
            return self_type::norm_inf (t1 - t2) < BOOST_UBLAS_TYPE_CHECK_EPSILON *
                   (std::max) ((std::max) (self_type::norm_inf (t1),
                                       self_type::norm_inf (t2)),
                             BOOST_UBLAS_TYPE_CHECK_MIN);
        }
    };
    template<>
    struct type_traits<double> {
        typedef type_traits<double> self_type;
        typedef double value_type;
        typedef const value_type &const_reference;
        typedef value_type &reference;
        typedef value_type real_type;
#ifndef BOOST_UBLAS_NO_LONG_DOUBLE
        typedef long double precision_type;
#else
        typedef value_type precision_type;
#endif

        BOOST_STATIC_CONSTANT (unsigned, plus_complexity = 1);
        BOOST_STATIC_CONSTANT (unsigned, multiplies_complexity = 1);

        static
        BOOST_UBLAS_INLINE
        real_type real (const_reference t) {
                return t;
        }
        static
        BOOST_UBLAS_INLINE
        real_type imag (const_reference /*t*/) {
                return 0;
        }
        static
        BOOST_UBLAS_INLINE
        value_type conj (const_reference t) {
                return t;
        }

        static
        BOOST_UBLAS_INLINE
        real_type abs (const_reference t) {
#if defined (BOOST_NO_STDC_NAMESPACE) || defined (BOOST_UBLAS_CMATH_BAD_STD)
            return ::fabs (t);
#else
            return std::abs (t);
#endif
        }
        static
        BOOST_UBLAS_INLINE
        value_type sqrt (const_reference t) {
#if defined (BOOST_NO_STDC_NAMESPACE) || defined (BOOST_UBLAS_CMATH_BAD_STD)
            return ::sqrt (t);
#else
            return std::sqrt (t);
#endif
        }

        static
        BOOST_UBLAS_INLINE
        real_type norm_1 (const_reference t) {
            return self_type::abs (t);
        }
        static
        BOOST_UBLAS_INLINE
        real_type norm_2 (const_reference t) {
            return self_type::abs (t);
        }
        static
        BOOST_UBLAS_INLINE
        real_type norm_inf (const_reference t) {
            return self_type::abs (t);
        }

        static
        BOOST_UBLAS_INLINE
        bool equals (const_reference t1, const_reference t2) {
            return self_type::norm_inf (t1 - t2) < BOOST_UBLAS_TYPE_CHECK_EPSILON *
                   (std::max) ((std::max) (self_type::norm_inf (t1),
                                       self_type::norm_inf (t2)),
                             BOOST_UBLAS_TYPE_CHECK_MIN);
        }
    };
#ifndef BOOST_UBLAS_NO_LONG_DOUBLE
    template<>
    struct type_traits<long double> {
        typedef type_traits<long double> self_type;
        typedef long double value_type;
        typedef const value_type &const_reference;
        typedef value_type &reference;
        typedef value_type real_type;
        typedef value_type precision_type;

        BOOST_STATIC_CONSTANT (unsigned, plus_complexity = 1);
        BOOST_STATIC_CONSTANT (unsigned, multiplies_complexity = 1);

        static
        BOOST_UBLAS_INLINE
        real_type real (const_reference t) {
                return t;
        }
        static
        BOOST_UBLAS_INLINE
        real_type imag (const_reference /*t*/) {
                return 0;
        }
        static
        BOOST_UBLAS_INLINE
        value_type conj (const_reference t) {
                return t;
        }

        static
        BOOST_UBLAS_INLINE
        real_type abs (const_reference t) {
#if defined (BOOST_NO_STDC_NAMESPACE) || defined (BOOST_UBLAS_CMATH_BAD_STD)
            return ::fabsl (t);
#else
            return std::abs (t);
#endif
        }
        static
        BOOST_UBLAS_INLINE
        value_type sqrt (const_reference t) {
#if defined (BOOST_NO_STDC_NAMESPACE) || defined (BOOST_UBLAS_CMATH_BAD_STD)
            return ::sqrtl (t);
#else
            return std::sqrt (t);
#endif
        }

        static
        BOOST_UBLAS_INLINE
        real_type norm_1 (const_reference t) {
            return self_type::abs (t);
        }
        static
        BOOST_UBLAS_INLINE
        real_type norm_2 (const_reference t) {
            return self_type::abs (t);
        }
        static
        BOOST_UBLAS_INLINE
        real_type norm_inf (const_reference t) {
            return self_type::abs (t);
        }

        static
        BOOST_UBLAS_INLINE
        bool equals (const_reference t1, const_reference t2) {
            return self_type::norm_inf (t1 - t2) < BOOST_UBLAS_TYPE_CHECK_EPSILON *
                   (std::max) ((std::max) (self_type::norm_inf (t1),
                                       self_type::norm_inf (t2)),
                             BOOST_UBLAS_TYPE_CHECK_MIN);
        }
    };
#endif

    template<>
    struct type_traits<std::complex<float> > {
        typedef type_traits<std::complex<float> > self_type;
        typedef std::complex<float> value_type;
        typedef const value_type &const_reference;
        typedef value_type &reference;
        typedef float real_type;
        typedef std::complex<double> precision_type;

        BOOST_STATIC_CONSTANT (unsigned, plus_complexity = 2);
        BOOST_STATIC_CONSTANT (unsigned, multiplies_complexity = 6);

        static
        BOOST_UBLAS_INLINE
        real_type real (const_reference t) {
                return std::real (t);
        }
        static
        BOOST_UBLAS_INLINE
        real_type imag (const_reference t) {
                return std::imag (t);
        }
        static
        BOOST_UBLAS_INLINE
        value_type conj (const_reference t) {
                return std::conj (t);
        }

        static
        BOOST_UBLAS_INLINE
        real_type abs (const_reference t) {
                return std::abs (t);
        }
        static
        BOOST_UBLAS_INLINE
        value_type sqrt (const_reference t) {
                return std::sqrt (t);
        }

        static
        BOOST_UBLAS_INLINE
        real_type norm_1 (const_reference t) {
            return type_traits<real_type>::abs (self_type::real (t)) +
                   type_traits<real_type>::abs (self_type::imag (t));
        }
        static
        BOOST_UBLAS_INLINE
        real_type norm_2 (const_reference t) {
            return self_type::abs (t);
        }
        static
        BOOST_UBLAS_INLINE
        real_type norm_inf (const_reference t) {
            return (std::max) (type_traits<real_type>::abs (self_type::real (t)),
                             type_traits<real_type>::abs (self_type::imag (t)));
        }

        static
        BOOST_UBLAS_INLINE
        bool equals (const_reference t1, const_reference t2) {
            return self_type::norm_inf (t1 - t2) < BOOST_UBLAS_TYPE_CHECK_EPSILON *
                   (std::max) ((std::max) (self_type::norm_inf (t1),
                                       self_type::norm_inf (t2)),
                             BOOST_UBLAS_TYPE_CHECK_MIN);
        }
    };
    template<>
    struct type_traits<std::complex<double> > {
        typedef type_traits<std::complex<double> > self_type;
        typedef std::complex<double> value_type;
        typedef const value_type &const_reference;
        typedef value_type &reference;
        typedef double real_type;
#ifndef BOOST_UBLAS_NO_LONG_DOUBLE
        typedef std::complex<long double> precision_type;
#else
        typedef value_type precision_type;
#endif

        BOOST_STATIC_CONSTANT (unsigned, plus_complexity = 2);
        BOOST_STATIC_CONSTANT (unsigned, multiplies_complexity = 6);

        static
        BOOST_UBLAS_INLINE
        real_type real (const_reference t) {
                return std::real (t);
        }
        static
        BOOST_UBLAS_INLINE
        real_type imag (const_reference t) {
                return std::imag (t);
        }
        static
        BOOST_UBLAS_INLINE
        value_type conj (const_reference t) {
                return std::conj (t);
        }

        static
        BOOST_UBLAS_INLINE
        real_type abs (const_reference t) {
                return std::abs (t);
        }
        static
        BOOST_UBLAS_INLINE
        value_type sqrt (const_reference t) {
                return std::sqrt (t);
        }

        static
        BOOST_UBLAS_INLINE
        real_type norm_1 (const_reference t) {
            return type_traits<real_type>::abs (self_type::real (t)) +
                   type_traits<real_type>::abs (self_type::imag (t));
        }
        static
        BOOST_UBLAS_INLINE
        real_type norm_2 (const_reference t) {
            return self_type::abs (t);
        }
        static
        BOOST_UBLAS_INLINE
        real_type norm_inf (const_reference t) {
            return (std::max) (type_traits<real_type>::abs (self_type::real (t)),
                             type_traits<real_type>::abs (self_type::imag (t)));
        }

        static
        BOOST_UBLAS_INLINE
        bool equals (const_reference t1, const_reference t2) {
            return self_type::norm_inf (t1 - t2) < BOOST_UBLAS_TYPE_CHECK_EPSILON *
                   (std::max) ((std::max) (self_type::norm_inf (t1),
                                       self_type::norm_inf (t2)),
                             BOOST_UBLAS_TYPE_CHECK_MIN);
        }
    };
#ifndef BOOST_UBLAS_NO_LONG_DOUBLE
    template<>
    struct type_traits<std::complex<long double> > {
        typedef type_traits<std::complex<long double> > self_type;
        typedef std::complex<long double> value_type;
        typedef const value_type &const_reference;
        typedef value_type &reference;
        typedef long double real_type;
        typedef value_type precision_type;

        BOOST_STATIC_CONSTANT (unsigned, plus_complexity = 2);
        BOOST_STATIC_CONSTANT (unsigned, multiplies_complexity = 6);

        static
        BOOST_UBLAS_INLINE
        real_type real (const_reference t) {
                return std::real (t);
        }
        static
        BOOST_UBLAS_INLINE
        real_type imag (const_reference t) {
                return std::imag (t);
        }
        static
        BOOST_UBLAS_INLINE
        value_type conj (const_reference t) {
                return std::conj (t);
        }

        static
        BOOST_UBLAS_INLINE
        real_type abs (const_reference t) {
                return std::abs (t);
        }
        static
        BOOST_UBLAS_INLINE
        value_type sqrt (const_reference t) {
                return std::sqrt (t);
        }

        static
        BOOST_UBLAS_INLINE
        real_type norm_1 (const_reference t) {
            return type_traits<real_type>::abs (self_type::real (t)) +
                   type_traits<real_type>::abs (self_type::imag (t));
        }
        static
        BOOST_UBLAS_INLINE
        real_type norm_2 (const_reference t) {
            return self_type::abs (t);
        }
        static
        BOOST_UBLAS_INLINE
        real_type norm_inf (const_reference t) {
            return (std::max) (type_traits<real_type>::abs (self_type::real (t)),
                             type_traits<real_type>::abs (self_type::imag (t)));
        }

        static
        BOOST_UBLAS_INLINE
        bool equals (const_reference t1, const_reference t2) {
            return self_type::norm_inf (t1 - t2) < BOOST_UBLAS_TYPE_CHECK_EPSILON *
                   (std::max) ((std::max) (self_type::norm_inf (t1),
                                       self_type::norm_inf (t2)),
                             BOOST_UBLAS_TYPE_CHECK_MIN);
        }
    };
#endif

#ifdef BOOST_UBLAS_USE_INTERVAL
    template<>
    struct type_traits<boost::numeric::interval<float> > {
        typedef type_traits<boost::numeric::interval<float> > self_type;
        typedef boost::numeric::interval<float> value_type;
        typedef const value_type &const_reference;
        typedef value_type &reference;
        typedef value_type real_type;
        typedef boost::numeric::interval<double> precision_type;

        BOOST_STATIC_CONSTANT (unsigned, plus_complexity = 1);
        BOOST_STATIC_CONSTANT (unsigned, multiplies_complexity = 1);

        static
        BOOST_UBLAS_INLINE
        real_type real (const_reference t) {
                return t;
        }
        static
        BOOST_UBLAS_INLINE
        real_type imag (const_reference t) {
                return 0;
        }
        static
        BOOST_UBLAS_INLINE
        value_type conj (const_reference t) {
                return t;
        }

        static
        BOOST_UBLAS_INLINE
        real_type abs (const_reference t) {
            return boost::numeric::abs (t);
        }
        static
        BOOST_UBLAS_INLINE
        value_type sqrt (const_reference t) {
            return boost::numeric::sqrt (t);
        }

        static
        BOOST_UBLAS_INLINE
        real_type norm_1 (const_reference t) {
            return self_type::abs (t);
        }
        static
        BOOST_UBLAS_INLINE
        real_type norm_2 (const_reference t) {
            return self_type::abs (t);
        }
        static
        BOOST_UBLAS_INLINE
        real_type norm_inf (const_reference t) {
            return self_type::abs (t);
        }

        static
        BOOST_UBLAS_INLINE
        bool equals (const_reference t1, const_reference t2) {
            return self_type::norm_inf (t1 - t2) < BOOST_UBLAS_TYPE_CHECK_EPSILON *
                   (std::max) ((std::max) (self_type::norm_inf (t1),
                                       self_type::norm_inf (t2)),
                             BOOST_UBLAS_TYPE_CHECK_MIN);
        }
    };
    template<>
    struct type_traits<boost::numeric::interval<double> > {
        typedef type_traits<boost::numeric::interval<double> > self_type;
        typedef boost::numeric::interval<double> value_type;
        typedef const value_type &const_reference;
        typedef value_type &reference;
        typedef value_type real_type;
#ifndef BOOST_UBLAS_NO_LONG_DOUBLE
        typedef boost::numeric::interval<long double> precision_type;
#else
        typedef value_type precision_type;
#endif

        BOOST_STATIC_CONSTANT (unsigned, plus_complexity = 1);
        BOOST_STATIC_CONSTANT (unsigned, multiplies_complexity = 1);

        static
        BOOST_UBLAS_INLINE
        real_type real (const_reference t) {
                return t;
        }
        static
        BOOST_UBLAS_INLINE
        real_type imag (const_reference t) {
                return 0;
        }
        static
        BOOST_UBLAS_INLINE
        value_type conj (const_reference t) {
                return t;
        }

        static
        BOOST_UBLAS_INLINE
        real_type abs (const_reference t) {
            return boost::numeric::abs (t);
        }
        static
        BOOST_UBLAS_INLINE
        value_type sqrt (const_reference t) {
            return boost::numeric::sqrt (t);
        }

        static
        BOOST_UBLAS_INLINE
        real_type norm_1 (const_reference t) {
            return self_type::abs (t);
        }
        static
        BOOST_UBLAS_INLINE
        real_type norm_2 (const_reference t) {
            return self_type::abs (t);
        }
        static
        BOOST_UBLAS_INLINE
        real_type norm_inf (const_reference t) {
            return self_type::abs (t);
        }

        static
        BOOST_UBLAS_INLINE
        bool equals (const_reference t1, const_reference t2) {
            return self_type::norm_inf (t1 - t2) < BOOST_UBLAS_TYPE_CHECK_EPSILON *
                   (std::max) ((std::max) (self_type::norm_inf (t1),
                                       self_type::norm_inf (t2)),
                             BOOST_UBLAS_TYPE_CHECK_MIN);
        }
    };
#ifndef BOOST_UBLAS_NO_LONG_DOUBLE
    template<>
    struct type_traits<boost::numeric::interval<long double> > {
        typedef type_traits<boost::numeric::interval<long double> > self_type;
        typedef boost::numeric::interval<long double> value_type;
        typedef const value_type &const_reference;
        typedef value_type &reference;
        typedef value_type real_type;
        typedef value_type precision_type;

        BOOST_STATIC_CONSTANT (unsigned, plus_complexity = 1);
        BOOST_STATIC_CONSTANT (unsigned, multiplies_complexity = 1);

        static
        BOOST_UBLAS_INLINE
        real_type real (const_reference t) {
                return t;
        }
        static
        BOOST_UBLAS_INLINE
        real_type imag (const_reference t) {
                return 0;
        }
        static
        BOOST_UBLAS_INLINE
        value_type conj (const_reference t) {
                return t;
        }

        static
        BOOST_UBLAS_INLINE
        real_type abs (const_reference t) {
            return boost::numeric::abs (t);
        }
        static
        BOOST_UBLAS_INLINE
        value_type sqrt (const_reference t) {
            return boost::numeric::sqrt (t);
        }

        static
        BOOST_UBLAS_INLINE
        real_type norm_1 (const_reference t) {
            return self_type::abs (t);
        }
        static
        BOOST_UBLAS_INLINE
        real_type norm_2 (const_reference t) {
            return self_type::abs (t);
        }
        static
        BOOST_UBLAS_INLINE
        real_type norm_inf (const_reference t) {
            return self_type::abs (t);
        }

        static
        BOOST_UBLAS_INLINE
        bool equals (const_reference t1, const_reference t2) {
            return self_type::norm_inf (t1 - t2) < BOOST_UBLAS_TYPE_CHECK_EPSILON *
                   (std::max) ((std::max) (self_type::norm_inf (t1),
                                       self_type::norm_inf (t2)),
                             BOOST_UBLAS_TYPE_CHECK_MIN);
        }
    };
#endif

#ifdef BOOST_UBLAS_USE_BOOST_COMPLEX
    template<>
    struct type_traits<boost::complex<boost::numeric::interval<float> > > {
        typedef type_traits<boost::complex<boost::numeric::interval<float> > > self_type;
        typedef boost::complex<boost::numeric::interval<float> > value_type;
        typedef const value_type &const_reference;
        typedef value_type &reference;
        typedef boost::numeric::interval<float> real_type;
        typedef boost::complex<boost::numeric::interval<double> > precision_type;

        BOOST_STATIC_CONSTANT (unsigned, plus_complexity = 2);
        BOOST_STATIC_CONSTANT (unsigned, multiplies_complexity = 6);

        static
        BOOST_UBLAS_INLINE
        real_type real (const_reference t) {
                return std::real (t);
        }
        static
        BOOST_UBLAS_INLINE
        real_type imag (const_reference t) {
                return std::imag (t);
        }
        static
        BOOST_UBLAS_INLINE
        value_type conj (const_reference t) {
                return std::conj (t);
        }

        static
        BOOST_UBLAS_INLINE
        real_type abs (const_reference t) {
                return boost::numeric::abs (t);
        }
        static
        BOOST_UBLAS_INLINE
        value_type sqrt (const_reference t) {
                return boost::numeric::sqrt (t);
        }

        static
        BOOST_UBLAS_INLINE
        real_type norm_1 (const_reference t) {
            return type_traits<real_type>::abs (self_type::real (t)) +
                   type_traits<real_type>::abs (self_type::imag (t));
        }
        static
        BOOST_UBLAS_INLINE
        real_type norm_2 (const_reference t) {
            return self_type::abs (t);
        }
        static
        BOOST_UBLAS_INLINE
        real_type norm_inf (const_reference t) {
            return (std::max) (type_traits<real_type>::abs (self_type::real (t)),
                             type_traits<real_type>::abs (self_type::imag (t)));
        }

        static
        BOOST_UBLAS_INLINE
        bool equals (const_reference t1, const_reference t2) {
            return self_type::norm_inf (t1 - t2) < BOOST_UBLAS_TYPE_CHECK_EPSILON *
                   (std::max) ((std::max) (self_type::norm_inf (t1),
                                       self_type::norm_inf (t2)),
                             BOOST_UBLAS_TYPE_CHECK_MIN);
        }
    };
    template<>
    struct type_traits<boost::complex<boost::numeric::interval<double> > {
        typedef type_traits<boost::complex<boost::numeric::interval<double> > self_type;
        typedef boost::complex<boost::numeric::interval<double> > value_type;
        typedef const value_type &const_reference;
        typedef value_type &reference;
        typedef boost::numeric::interval<double> real_type;
#ifndef BOOST_UBLAS_NO_LONG_DOUBLE
        typedef boost::complex<boost::numeric::interval<long double> > precision_type;
#else
        typedef value_type precision_type;
#endif

        BOOST_STATIC_CONSTANT (unsigned, plus_complexity = 2);
        BOOST_STATIC_CONSTANT (unsigned, multiplies_complexity = 6);

        static
        BOOST_UBLAS_INLINE
        real_type real (const_reference t) {
                return std::real (t);
        }
        static
        BOOST_UBLAS_INLINE
        real_type imag (const_reference t) {
                return std::imag (t);
        }
        static
        BOOST_UBLAS_INLINE
        value_type conj (const_reference t) {
                return std::conj (t);
        }

        static
        BOOST_UBLAS_INLINE
        real_type abs (const_reference t) {
                return boost::numeric::abs (t);
        }
        static
        BOOST_UBLAS_INLINE
        value_type sqrt (const_reference t) {
                return boost::numeric::sqrt (t);
        }

        static
        BOOST_UBLAS_INLINE
        real_type norm_1 (const_reference t) {
            return type_traits<real_type>::abs (self_type::real (t)) +
                   type_traits<real_type>::abs (self_type::imag (t));
        }
        static
        BOOST_UBLAS_INLINE
        real_type norm_2 (const_reference t) {
            return self_type::abs (t);
        }
        static
        BOOST_UBLAS_INLINE
        real_type norm_inf (const_reference t) {
            return (std::max) (type_traits<real_type>::abs (self_type::real (t)),
                             type_traits<real_type>::abs (self_type::imag (t)));
        }

        static
        BOOST_UBLAS_INLINE
        bool equals (const_reference t1, const_reference t2) {
            return self_type::norm_inf (t1 - t2) < BOOST_UBLAS_TYPE_CHECK_EPSILON *
                   (std::max) ((std::max) (self_type::norm_inf (t1),
                                       self_type::norm_inf (t2)),
                             BOOST_UBLAS_TYPE_CHECK_MIN);
        }
    };
#ifndef BOOST_UBLAS_NO_LONG_DOUBLE
    template<>
    struct type_traits<boost::complex<boost::numeric::interval<long double> > > {
        typedef type_traits<boost::complex<boost::numeric::interval<long double> > > self_type;
        typedef boost::complex<boost::numeric::interval<long double> > value_type;
        typedef const value_type &const_reference;
        typedef value_type &reference;
        typedef boost::numeric::interval<long double> real_type;
        typedef value_type precision_type;

        BOOST_STATIC_CONSTANT (unsigned, plus_complexity = 2);
        BOOST_STATIC_CONSTANT (unsigned, multiplies_complexity = 6);

        static
        BOOST_UBLAS_INLINE
        real_type real (const_reference t) {
                return std::real (t);
        }
        static
        BOOST_UBLAS_INLINE
        real_type imag (const_reference t) {
                return std::imag (t);
        }
        static
        BOOST_UBLAS_INLINE
        value_type conj (const_reference t) {
                return std::conj (t);
        }

        static
        BOOST_UBLAS_INLINE
        real_type abs (const_reference t) {
                return boost::numeric::abs (t);
        }
        static
        BOOST_UBLAS_INLINE
        value_type sqrt (const_reference t) {
                return boost::numeric::sqrt (t);
        }

        static
        BOOST_UBLAS_INLINE
        real_type norm_1 (const_reference t) {
            return type_traits<real_type>::abs (self_type::real (t)) +
                   type_traits<real_type>::abs (self_type::imag (t));
        }
        static
        BOOST_UBLAS_INLINE
        real_type norm_2 (const_reference t) {
            return self_type::abs (t);
        }
        static
        BOOST_UBLAS_INLINE
        real_type norm_inf (const_reference t) {
            return (std::max) (type_traits<real_type>::abs (self_type::real (t)),
                             type_traits<real_type>::abs (self_type::imag (t)));
        }

        static
        BOOST_UBLAS_INLINE
        bool equals (const_reference t1, const_reference t2) {
            return self_type::norm_inf (t1 - t2) < BOOST_UBLAS_TYPE_CHECK_EPSILON *
                   (std::max) ((std::max) (self_type::norm_inf (t1),
                                       self_type::norm_inf (t2)),
                             BOOST_UBLAS_TYPE_CHECK_MIN);
        }
    };
#endif
#endif
#endif



#if !defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION) && !defined(BOOST_NO_SFINAE)
    // Use Joel de Guzman's return type deduction
    // uBLAS assumes a common return type for all binary arithmetic operators
    template<class X, class Y>
    struct promote_traits {
        typedef type_deduction_detail::base_result_of<X, Y> base_type;
        static typename base_type::x_type x;
        static typename base_type::y_type y;
        BOOST_STATIC_CONSTANT(int,
            size = sizeof(
                type_deduction_detail::test<
                    typename base_type::x_type
                  , typename base_type::y_type
                >(x + y)     // Use x+y to stand of all the arithmetic actions
            ));

        BOOST_STATIC_CONSTANT(int, index = (size / sizeof(char)) - 1);
        typedef typename mpl::at_c<
            typename base_type::types, index>::type id;
        typedef typename id::type promote_type;
    };
    template<class X, class Y>
    struct promote_type_multiplies {
        typedef type_deduction_detail::base_result_of<X, Y> base_type;
        static typename base_type::x_type x;
        static typename base_type::y_type y;
        BOOST_STATIC_CONSTANT(int,
            size = sizeof(
                type_deduction_detail::test<
                    typename base_type::x_type
                  , typename base_type::y_type
                >(x * y)     // Specifically the * arithmetic actions
            ));

        BOOST_STATIC_CONSTANT(int, index = (size / sizeof(char)) - 1);
        typedef typename mpl::at_c<
            typename base_type::types, index>::type id;
        typedef typename id::type promote_type;
    };


#else
    template<class T1, class T2>
    struct promote_traits {
        // Default promotion will badly fail, if the types are different.
        // Thanks to Kresimir Fresl for spotting this.
        BOOST_STATIC_ASSERT ((boost::is_same<T1, T2>::value));
        typedef T1 promote_type;
    };

    template<>
    struct promote_traits<float, double> {
        typedef double promote_type;
    };
    template<>
    struct promote_traits<double, float> {
        typedef double promote_type;
    };
#ifndef BOOST_UBLAS_NO_LONG_DOUBLE
    template<>
    struct promote_traits<float, long double> {
        typedef long double promote_type;
    };
    template<>
    struct promote_traits<long double, float> {
        typedef long double promote_type;
    };
    template<>
    struct promote_traits<double, long double> {
        typedef long double promote_type;
    };
    template<>
    struct promote_traits<long double, double> {
        typedef long double promote_type;
    };
#endif

    template<>
    struct promote_traits<float, std::complex<float> > {
        typedef std::complex<float> promote_type;
    };
    template<>
    struct promote_traits<std::complex<float>, float> {
        typedef std::complex<float> promote_type;
    };
    template<>
    struct promote_traits<float, std::complex<double> > {
        typedef std::complex<double> promote_type;
    };
    template<>
    struct promote_traits<std::complex<double>, float> {
        typedef std::complex<double> promote_type;
    };
#ifndef BOOST_UBLAS_NO_LONG_DOUBLE
    template<>
    struct promote_traits<float, std::complex<long double> > {
        typedef std::complex<long double> promote_type;
    };
    template<>
    struct promote_traits<std::complex<long double>, float> {
        typedef std::complex<long double> promote_type;
    };
#endif

    template<>
    struct promote_traits<double, std::complex<float> > {
        // Here we'd better go the conservative way.
        // typedef std::complex<float> promote_type;
        typedef std::complex<double> promote_type;
    };
    template<>
    struct promote_traits<std::complex<float>, double> {
        // Here we'd better go the conservative way.
        // typedef std::complex<float> promote_type;
        typedef std::complex<double> promote_type;
    };
    template<>
    struct promote_traits<double, std::complex<double> > {
        typedef std::complex<double> promote_type;
    };
    template<>
    struct promote_traits<std::complex<double>, double> {
        typedef std::complex<double> promote_type;
    };
#ifndef BOOST_UBLAS_NO_LONG_DOUBLE
    template<>
    struct promote_traits<double, std::complex<long double> > {
        typedef std::complex<long double> promote_type;
    };
    template<>
    struct promote_traits<std::complex<long double>, double> {
        typedef std::complex<long double> promote_type;
    };
#endif

#ifndef BOOST_UBLAS_NO_LONG_DOUBLE
    template<>
    struct promote_traits<long double, std::complex<float> > {
        // Here we'd better go the conservative way.
        // typedef std::complex<float> promote_type;
        typedef std::complex<long double> promote_type;
    };
    template<>
    struct promote_traits<std::complex<float>, long double> {
        // Here we'd better go the conservative way.
        // typedef std::complex<float> promote_type;
        typedef std::complex<long double> promote_type;
    };
    template<>
    struct promote_traits<long double, std::complex<double> > {
        // Here we'd better go the conservative way.
        // typedef std::complex<double> promote_type;
        typedef std::complex<long double> promote_type;
    };
    template<>
    struct promote_traits<std::complex<double>, long double> {
        // Here we'd better go the conservative way.
        // typedef std::complex<double> promote_type;
        typedef std::complex<long double> promote_type;
    };
    template<>
    struct promote_traits<long double, std::complex<long double> > {
        typedef std::complex<long double> promote_type;
    };
    template<>
    struct promote_traits<std::complex<long double>, long double> {
        typedef std::complex<long double> promote_type;
    };
#endif

    template<>
    struct promote_traits<std::complex<float>, std::complex<double> > {
        typedef std::complex<double> promote_type;
    };
    template<>
    struct promote_traits<std::complex<double>, std::complex<float> > {
        typedef std::complex<double> promote_type;
    };
#ifndef BOOST_UBLAS_NO_LONG_DOUBLE
    template<>
    struct promote_traits<std::complex<float>, std::complex<long double> > {
        typedef std::complex<long double> promote_type;
    };
    template<>
    struct promote_traits<std::complex<long double>, std::complex<float> > {
        typedef std::complex<long double> promote_type;
    };
    template<>
    struct promote_traits<std::complex<double>, std::complex<long double> > {
        typedef std::complex<long double> promote_type;
    };
    template<>
    struct promote_traits<std::complex<long double>, std::complex<double> > {
        typedef std::complex<long double> promote_type;
    };
#endif

#ifdef BOOST_UBLAS_USE_INTERVAL
    template<>
    struct promote_traits<boost::numeric::interval<float>, boost::numeric::interval<double> > {
        typedef boost::numeric::interval<double> promote_type;
    };
    template<>
    struct promote_traits<boost::numeric::interval<double>, boost::numeric::interval<float> > {
        typedef boost::numeric::interval<double> promote_type;
    };
#ifndef BOOST_UBLAS_NO_LONG_DOUBLE
    template<>
    struct promote_traits<boost::numeric::interval<float>, boost::numeric::interval<long double> > {
        typedef boost::numeric::interval<long double> promote_type;
    };
    template<>
    struct promote_traits<boost::numeric::interval<long double>, boost::numeric::interval<float> > {
        typedef boost::numeric::interval<long double> promote_type;
    };
    template<>
    struct promote_traits<boost::numeric::interval<double>, boost::numeric::interval<long double> > {
        typedef boost::numeric::interval<long double> promote_type;
    };
    template<>
    struct promote_traits<boost::numeric::interval<long double>, boost::numeric::interval<double> > {
        typedef boost::numeric::interval<long double> promote_type;
    };
#endif

#ifdef BOOST_UBLAS_USE_BOOST_COMPLEX
    template<>
    struct promote_traits<boost::numeric::interval<float>, boost::complex<boost::numeric::interval<float> > > {
        typedef boost::complex<boost::numeric::interval<float> > promote_type;
    };
    template<>
    struct promote_traits<boost::complex<boost::numeric::interval<float> >, boost::numeric::interval<float> > {
        typedef boost::complex<boost::numeric::interval<float> > promote_type;
    };
    template<>
    struct promote_traits<boost::numeric::interval<float>, boost::complex<boost::numeric::interval<double> > > {
        typedef boost::complex<boost::numeric::interval<double> > promote_type;
    };
    template<>
    struct promote_traits<boost::complex<boost::numeric::interval<double> >, boost::numeric::interval<float> > {
        typedef boost::complex<boost::numeric::interval<double> > promote_type;
    };
#ifndef BOOST_UBLAS_NO_LONG_DOUBLE
    template<>
    struct promote_traits<boost::numeric::interval<float>, boost::complex<boost::numeric::interval<long double> > > {
        typedef boost::complex<boost::numeric::interval<long double> > promote_type;
    };
    template<>
    struct promote_traits<boost::complex<boost::numeric::interval<long double> >, boost::numeric::interval<float> > {
        typedef boost::complex<boost::numeric::interval<long double> > promote_type;
    };
#endif

    template<>
    struct promote_traits<boost::numeric::interval<double>, boost::complex<boost::numeric::interval<float> > > {
        // Here we'd better go the conservative way.
        // typedef boost::complex<boost::numeric::interval<float> > promote_type;
        typedef boost::complex<boost::numeric::interval<double> > promote_type;
    };
    template<>
    struct promote_traits<boost::complex<boost::numeric::interval<float> >, boost::numeric::interval<double> > {
        // Here we'd better go the conservative way.
        // typedef boost::complex<boost::numeric::interval<float> > promote_type;
        typedef boost::complex<boost::numeric::interval<double> > promote_type;
    };
    template<>
    struct promote_traits<boost::numeric::interval<double>, boost::complex<boost::numeric::interval<double> > > {
        typedef boost::complex<boost::numeric::interval<double> > promote_type;
    };
    template<>
    struct promote_traits<boost::complex<boost::numeric::interval<double> >, boost::numeric::interval<double> > {
        typedef boost::complex<boost::numeric::interval<double> > promote_type;
    };
#ifndef BOOST_UBLAS_NO_LONG_DOUBLE
    template<>
    struct promote_traits<boost::numeric::interval<double>, boost::complex<boost::numeric::interval<long double> > > {
        typedef boost::complex<boost::numeric::interval<long double> > promote_type;
    };
    template<>
    struct promote_traits<boost::complex<boost::numeric::interval<long double> >, boost::numeric::interval<double> > {
        typedef boost::complex<boost::numeric::interval<long double> > promote_type;
    };
#endif

#ifndef BOOST_UBLAS_NO_LONG_DOUBLE
    template<>
    struct promote_traits<boost::numeric::interval<long double>, boost::complex<boost::numeric::interval<float> > > {
        // Here we'd better go the conservative way.
        // typedef boost::complex<boost::numeric::interval<float> > promote_type;
        typedef boost::complex<boost::numeric::interval<long double> > promote_type;
    };
    template<>
    struct promote_traits<boost::complex<boost::numeric::interval<float> >, boost::numeric::interval<long double> > {
        // Here we'd better go the conservative way.
        // typedef boost::complex<boost::numeric::interval<float> > promote_type;
        typedef boost::complex<boost::numeric::interval<long double> > promote_type;
    };
    template<>
    struct promote_traits<boost::numeric::interval<long double>, boost::complex<boost::numeric::interval<double> > > {
        // Here we'd better go the conservative way.
        // typedef boost::complex<boost::numeric::interval<double> > promote_type;
        typedef boost::complex<boost::numeric::interval<long double> > promote_type;
    };
    template<>
    struct promote_traits<boost::complex<boost::numeric::interval<double> >, boost::numeric::interval<long double> > {
        // Here we'd better go the conservative way.
        // typedef boost::complex<boost::numeric::interval<double> > promote_type;
        typedef boost::complex<boost::numeric::interval<long double> > promote_type;
    };
    template<>
    struct promote_traits<boost::numeric::interval<long double>, boost::complex<boost::numeric::interval<long double> > > {
        typedef boost::complex<boost::numeric::interval<long double> > promote_type;
    };
    template<>
    struct promote_traits<boost::complex<boost::numeric::interval<long double> >, boost::numeric::interval<long double> > {
        typedef boost::complex<boost::numeric::interval<long double> > promote_type;
    };
#endif

    template<>
    struct promote_traits<boost::complex<boost::numeric::interval<float> >, boost::complex<boost::numeric::interval<double> > > {
        typedef boost::complex<boost::numeric::interval<double> > promote_type;
    };
    template<>
    struct promote_traits<boost::complex<boost::numeric::interval<double> >, boost::complex<boost::numeric::interval<float> > > {
        typedef boost::complex<boost::numeric::interval<double> > promote_type;
    };
#ifndef BOOST_UBLAS_NO_LONG_DOUBLE
    template<>
    struct promote_traits<boost::complex<boost::numeric::interval<float> >, boost::complex<boost::numeric::interval<long double> > > {
        typedef boost::complex<boost::numeric::interval<long double> > promote_type;
    };
    template<>
    struct promote_traits<boost::complex<boost::numeric::interval<long double> >, boost::complex<boost::numeric::interval<float> > > {
        typedef boost::complex<boost::numeric::interval<long double> > promote_type;
    };
    template<>
    struct promote_traits<boost::complex<boost::numeric::interval<double> >, boost::complex<boost::numeric::interval<long double> > > {
        typedef boost::complex<boost::numeric::interval<long double> > promote_type;
    };
    template<>
    struct promote_traits<boost::complex<boost::numeric::interval<long double> >, boost::complex<boost::numeric::interval<double> > > {
        typedef boost::complex<boost::numeric::interval<long double> > promote_type;
    };
#endif
#endif
#endif
#endif

    struct unknown_storage_tag {};
    struct sparse_proxy_tag: public unknown_storage_tag {};
    struct sparse_tag: public sparse_proxy_tag {};
    struct packed_proxy_tag: public sparse_proxy_tag {};
    struct packed_tag: public packed_proxy_tag {};
    struct dense_proxy_tag: public packed_proxy_tag {};
    struct dense_tag: public dense_proxy_tag {};

    template<class S1, class S2>
    struct storage_restrict_traits {
        typedef S1 storage_category;
    };

    template<>
    struct storage_restrict_traits<sparse_tag, dense_proxy_tag> {
        typedef sparse_proxy_tag storage_category;
    };
    template<>
    struct storage_restrict_traits<sparse_tag, packed_proxy_tag> {
        typedef sparse_proxy_tag storage_category;
    };
    template<>
    struct storage_restrict_traits<sparse_tag, sparse_proxy_tag> {
        typedef sparse_proxy_tag storage_category;
    };

    template<>
    struct storage_restrict_traits<packed_tag, dense_proxy_tag> {
        typedef packed_proxy_tag storage_category;
    };
    template<>
    struct storage_restrict_traits<packed_tag, packed_proxy_tag> {
        typedef packed_proxy_tag storage_category;
    };
    template<>
    struct storage_restrict_traits<packed_tag, sparse_proxy_tag> {
        typedef sparse_proxy_tag storage_category;
    };

    template<>
    struct storage_restrict_traits<packed_proxy_tag, sparse_proxy_tag> {
        typedef sparse_proxy_tag storage_category;
    };

    template<>
    struct storage_restrict_traits<dense_tag, dense_proxy_tag> {
        typedef dense_proxy_tag storage_category;
    };
    template<>
    struct storage_restrict_traits<dense_tag, packed_proxy_tag> {
        typedef packed_proxy_tag storage_category;
    };
    template<>
    struct storage_restrict_traits<dense_tag, sparse_proxy_tag> {
        typedef sparse_proxy_tag storage_category;
    };

    template<>
    struct storage_restrict_traits<dense_proxy_tag, packed_proxy_tag> {
        typedef packed_proxy_tag storage_category;
    };
    template<>
    struct storage_restrict_traits<dense_proxy_tag, sparse_proxy_tag> {
        typedef sparse_proxy_tag storage_category;
    };

    struct sparse_bidirectional_iterator_tag : public std::bidirectional_iterator_tag {};
    struct packed_random_access_iterator_tag : public std::random_access_iterator_tag {};
    struct dense_random_access_iterator_tag : public packed_random_access_iterator_tag {};

    // Thanks to Kresimir Fresl for convincing Comeau with iterator_base_traits ;-)
    template<class IC>
    struct iterator_base_traits {};

    template<>
    struct iterator_base_traits<std::forward_iterator_tag> {
        template<class I, class T>
        struct iterator_base {
            typedef forward_iterator_base<std::forward_iterator_tag, I, T> type;
        };
    };

    template<>
    struct iterator_base_traits<std::bidirectional_iterator_tag> {
        template<class I, class T>
        struct iterator_base {
            typedef bidirectional_iterator_base<std::bidirectional_iterator_tag, I, T> type;
        };
    };
    template<>
    struct iterator_base_traits<sparse_bidirectional_iterator_tag> {
        template<class I, class T>
        struct iterator_base {
            typedef bidirectional_iterator_base<sparse_bidirectional_iterator_tag, I, T> type;
        };
    };

    template<>
    struct iterator_base_traits<std::random_access_iterator_tag> {
        template<class I, class T>
        struct iterator_base {
            typedef random_access_iterator_base<std::random_access_iterator_tag, I, T> type;
        };
    };
    template<>
    struct iterator_base_traits<packed_random_access_iterator_tag> {
        template<class I, class T>
        struct iterator_base {
            typedef random_access_iterator_base<packed_random_access_iterator_tag, I, T> type;
        };
    };
    template<>
    struct iterator_base_traits<dense_random_access_iterator_tag> {
        template<class I, class T>
        struct iterator_base {
            typedef random_access_iterator_base<dense_random_access_iterator_tag, I, T> type;
        };
    };

    template<class I1, class I2>
    struct iterator_restrict_traits {
        typedef I1 iterator_category;
    };

    template<>
    struct iterator_restrict_traits<packed_random_access_iterator_tag, sparse_bidirectional_iterator_tag> {
        typedef sparse_bidirectional_iterator_tag iterator_category;
    };
    template<>
    struct iterator_restrict_traits<sparse_bidirectional_iterator_tag, packed_random_access_iterator_tag> {
        typedef sparse_bidirectional_iterator_tag iterator_category;
    };

    template<>
    struct iterator_restrict_traits<dense_random_access_iterator_tag, sparse_bidirectional_iterator_tag> {
        typedef sparse_bidirectional_iterator_tag iterator_category;
    };
    template<>
    struct iterator_restrict_traits<sparse_bidirectional_iterator_tag, dense_random_access_iterator_tag> {
        typedef sparse_bidirectional_iterator_tag iterator_category;
    };

    template<>
    struct iterator_restrict_traits<dense_random_access_iterator_tag, packed_random_access_iterator_tag> {
        typedef packed_random_access_iterator_tag iterator_category;
    };
    template<>
    struct iterator_restrict_traits<packed_random_access_iterator_tag, dense_random_access_iterator_tag> {
        typedef packed_random_access_iterator_tag iterator_category;
    };

    template<class I>
    BOOST_UBLAS_INLINE
    void increment (I &it, const I &it_end, BOOST_UBLAS_TYPENAME I::difference_type compare, packed_random_access_iterator_tag) {
        it += (std::min) (compare, it_end - it);
    }
    template<class I>
    BOOST_UBLAS_INLINE
    void increment (I &it, const I &/* it_end */, BOOST_UBLAS_TYPENAME I::difference_type /* compare */, sparse_bidirectional_iterator_tag) {
        ++ it;
    }
    template<class I>
    BOOST_UBLAS_INLINE
    void increment (I &it, const I &it_end, BOOST_UBLAS_TYPENAME I::difference_type compare) {
        increment (it, it_end, compare, BOOST_UBLAS_TYPENAME I::iterator_category ());
    }

    template<class I>
    BOOST_UBLAS_INLINE
    void increment (I &it, const I &it_end) {
#if BOOST_UBLAS_TYPE_CHECK
        I cit (it);
        while (cit != it_end) {
            BOOST_UBLAS_CHECK (*cit == BOOST_UBLAS_TYPENAME I::value_type (0), internal_logic ());
            ++ cit;
        }
#endif
        it = it_end;
    }

}}}

#endif
