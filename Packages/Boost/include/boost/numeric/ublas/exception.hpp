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

#ifndef BOOST_UBLAS_EXCEPTION_H
#define BOOST_UBLAS_EXCEPTION_H

#if ! defined (BOOST_NO_EXCEPTIONS) && ! defined (BOOST_UBLAS_NO_EXCEPTIONS)
#include <stdexcept>
#else
#include <cstdlib>
#endif
#ifndef BOOST_UBLAS_NO_STD_CERR
#include <iostream>
#endif

#include <boost/numeric/ublas/config.hpp>

namespace boost { namespace numeric { namespace ublas {

    struct divide_by_zero
#if ! defined (BOOST_NO_EXCEPTIONS) && ! defined (BOOST_UBLAS_NO_EXCEPTIONS)
        // Inherit from standard exceptions as requested during review.
        : public std::runtime_error {
        explicit
        divide_by_zero (const char *s = "divide by zero") :
            std::runtime_error (s) {}
        void raise () {
            throw *this;
        }
#else
    {
        explicit
        divide_by_zero (const char *s = 0)
            {}
        void raise () {
#ifdef BOOST_NO_STDC_NAMESPACE
            ::abort ();
#else
            std::abort ();
#endif
        }
#endif
    };

    struct internal_logic
#if ! defined (BOOST_NO_EXCEPTIONS) && ! defined (BOOST_UBLAS_NO_EXCEPTIONS)
        // Inherit from standard exceptions as requested during review.
        : public std::logic_error {
        explicit
        internal_logic (const char *s = "internal logic") :
            std::logic_error (s) {}
        void raise () {
            throw *this;
        }
#else
    {
        explicit
        internal_logic (const char *s = 0)
            {}
        void raise () {
#ifdef BOOST_NO_STDC_NAMESPACE
            ::abort ();
#else
            std::abort ();
#endif
        }
#endif
    };

    struct external_logic
#if ! defined (BOOST_NO_EXCEPTIONS) && ! defined (BOOST_UBLAS_NO_EXCEPTIONS)
        // Inherit from standard exceptions as requested during review.
        : public std::logic_error {
        explicit
        external_logic (const char *s = "external logic") :
            std::logic_error (s) {}
        // virtual const char *what () const throw () {
        //     return "exception: external logic";
        // }
        void raise () {
            throw *this;
        }
#else
    {
        explicit
        external_logic (const char *s = 0)
            {}
        void raise () {
#ifdef BOOST_NO_STDC_NAMESPACE
            ::abort ();
#else
            std::abort ();
#endif
        }
#endif
    };

    struct bad_argument
#if ! defined (BOOST_NO_EXCEPTIONS) && ! defined (BOOST_UBLAS_NO_EXCEPTIONS)
        // Inherit from standard exceptions as requested during review.
        : public std::invalid_argument {
        explicit
        bad_argument (const char *s = "bad argument") :
            std::invalid_argument (s) {}
        void raise () {
            throw *this;
        }
#else
    {
        explicit
        bad_argument (const char *s = 0)
            {}
        void raise () {
            throw *this;
#ifdef BOOST_NO_STDC_NAMESPACE
            ::abort ();
#else
            std::abort ();
#endif
        }
#endif
    };

    struct bad_size
#if ! defined (BOOST_NO_EXCEPTIONS) && ! defined (BOOST_UBLAS_NO_EXCEPTIONS)
        // Inherit from standard exceptions as requested during review.
        : public std::domain_error {
        explicit
        bad_size (const char *s = "bad size") :
            std::domain_error (s) {}
        void raise () {
            throw *this;
        }
#else
    {
        explicit
        bad_size (const char *s = 0)
            {}
        void raise () {
#ifdef BOOST_NO_STDC_NAMESPACE
            ::abort ();
#else
            std::abort ();
#endif
        }
#endif
    };

    struct bad_index
#if ! defined (BOOST_NO_EXCEPTIONS) && ! defined (BOOST_UBLAS_NO_EXCEPTIONS)
        // Inherit from standard exceptions as requested during review.
        : public std::out_of_range {
        explicit
        bad_index (const char *s = "bad index") :
            std::out_of_range (s) {}
        void raise () {
            throw *this;
        }
#else
    {
        explicit
        bad_index (const char *s = 0)
            {}
        void raise () {
#ifdef BOOST_NO_STDC_NAMESPACE
            ::abort ();
#else
            std::abort ();
#endif
        }
#endif
    };

    struct singular
#if ! defined (BOOST_NO_EXCEPTIONS) && ! defined (BOOST_UBLAS_NO_EXCEPTIONS)
        // Inherit from standard exceptions as requested during review.
        : public std::runtime_error {
        explicit
        singular (const char *s = "singular") :
            std::runtime_error (s) {}
        void raise () {
            throw *this;
        }
#else
    {
        explicit
        singular (const char *s = 0)
            {}
        void raise () {
            throw *this;
#ifdef BOOST_NO_STDC_NAMESPACE
            ::abort ();
#else
            std::abort ();
#endif
        }
#endif
    };

    struct non_real
#if ! defined (BOOST_NO_EXCEPTIONS) && ! defined (BOOST_UBLAS_NO_EXCEPTIONS)
        // Inherit from standard exceptions as requested during review.
        : public std::domain_error {
        explicit
        non_real (const char *s = "exception: non real") :
            std::domain_error (s) {}
        void raise () {
            throw *this;
        }
#else
     {
        explicit
        non_real (const char *s = 0)
            {}
        void raise () {
#ifdef BOOST_NO_STDC_NAMESPACE
            ::abort ();
#else
            std::abort ();
#endif
        }
#endif
    };

#if BOOST_UBLAS_CHECK_ENABLE
// FIXME: for performance reasons we better use macros
//    template<class E>
//    BOOST_UBLAS_INLINE
//    void check (bool expression, const E &e) {
//        if (! expression)
//            e.raise ();
//    }
//    template<class E>
//    BOOST_UBLAS_INLINE
//    void check_ex (bool expression, const char *file, int line, const E &e) {
//        if (! expression)
//            e.raise ();
//    }
// Dan Muller reported problems with COMO in GUI applications
// So we need a new preprocessor symbol:
#ifndef BOOST_UBLAS_NO_STD_CERR
#define BOOST_UBLAS_CHECK(expression, e) \
    if (! (expression)) { \
        std::cerr << "Assertion failed in file " << __FILE__ << " at line " << __LINE__ << ":" << std::endl; \
        std::cerr << #expression << std::endl; \
        e.raise (); \
    }
#define BOOST_UBLAS_CHECK_EX(expression, file, line, e) \
    if (! (expression)) { \
        std::cerr << "Assertion failed in file " << (file) << " at line " << (line) << ":" << std::endl; \
        std::cerr << #expression << std::endl; \
        e.raise (); \
    }
#else
#define BOOST_UBLAS_CHECK(expression, e) \
    if (! (expression)) { \
        e.raise (); \
    }
#define BOOST_UBLAS_CHECK_EX(expression, file, line, e) \
    if (! (expression)) { \
        e.raise (); \
    }
#endif
#else
// FIXME: for performance reasons we better use macros
//    template<class E>
//    BOOST_UBLAS_INLINE
//    void check (bool expression, const E &e) {}
//    template<class E>
//    BOOST_UBLAS_INLINE
//    void check_ex (bool expression, const char *file, int line, const E &e) {}
#define BOOST_UBLAS_CHECK(expression, e)
#define BOOST_UBLAS_CHECK_EX(expression, file, line, e)
#endif

}}}

#endif
