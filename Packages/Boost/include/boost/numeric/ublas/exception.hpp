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

#include <cstdlib>
#include <exception>
#include <stdexcept>
#include <iostream>

#include <boost/numeric/ublas/config.hpp>

namespace boost { namespace numeric { namespace ublas {

    struct divide_by_zero:
        // Inherit from standard exceptions as requested during review.
        // public std::exception {
        public std::runtime_error {
        BOOST_UBLAS_EXPLICIT
        divide_by_zero (const std::string &s = "exception: divide by zero"):
            std::runtime_error (s) {}
        // virtual const char *what () const throw () {
        //     return "exception: divide by zero";
        // }
        virtual void raise () {
#if ! defined (BOOST_NO_EXCEPTIONS) && defined (BOOST_UBLAS_USE_EXCEPTIONS)
            throw *this;
#else
#ifdef BOOST_NO_STDC_NAMESPACE
            ::abort ();
#else
            std::abort ();
#endif
#endif
        }
    };
    struct internal_logic:
        // Inherit from standard exceptions as requested during review.
        // public std::exception {
        public std::logic_error {
        BOOST_UBLAS_EXPLICIT
        internal_logic (const std::string &s = "exception: internal logic"):
            std::logic_error (s) {}
        // virtual const char *what () const throw () {
        //     return "exception: internal logic";
        // }
        virtual void raise () {
#if ! defined (BOOST_NO_EXCEPTIONS) && defined (BOOST_UBLAS_USE_EXCEPTIONS)
            throw *this;
#else
#ifdef BOOST_NO_STDC_NAMESPACE
            ::abort ();
#else
            std::abort ();
#endif
#endif
        }
    };
    struct external_logic:
        // Inherit from standard exceptions as requested during review.
        // public std::exception {
        public std::logic_error {
        BOOST_UBLAS_EXPLICIT
        external_logic (const std::string &s = "exception: external logic"):
            std::logic_error (s) {}
        // virtual const char *what () const throw () {
        //     return "exception: external logic";
        // }
        virtual void raise () {
#if ! defined (BOOST_NO_EXCEPTIONS) && defined (BOOST_UBLAS_USE_EXCEPTIONS)
            throw *this;
#else
#ifdef BOOST_NO_STDC_NAMESPACE
            ::abort ();
#else
            std::abort ();
#endif
#endif
        }
    };
    struct bad_argument:
        // Inherit from standard exceptions as requested during review.
        // public std::exception {
        public std::invalid_argument {
        BOOST_UBLAS_EXPLICIT
        bad_argument (const std::string &s = "exception: bad argument"):
            std::invalid_argument (s) {}
        // virtual const char *what () const throw () {
        //     return "exception: bad argument";
        // }
        virtual void raise () {
#if ! defined (BOOST_NO_EXCEPTIONS) && defined (BOOST_UBLAS_USE_EXCEPTIONS)
            throw *this;
#else
#ifdef BOOST_NO_STDC_NAMESPACE
            ::abort ();
#else
            std::abort ();
#endif
#endif
        }
    };
    struct bad_size:
        // Inherit from standard exceptions as requested during review.
        // public std::exception {
        public std::domain_error {
        BOOST_UBLAS_EXPLICIT
        bad_size (const std::string &s = "exception: bad size"):
            std::domain_error (s) {}
        // virtual const char *what () const throw () {
        //     return "exception: bad size";
        // }
        virtual void raise () {
#if ! defined (BOOST_NO_EXCEPTIONS) && defined (BOOST_UBLAS_USE_EXCEPTIONS)
            throw *this;
#else
#ifdef BOOST_NO_STDC_NAMESPACE
            ::abort ();
#else
            std::abort ();
#endif
#endif
        }
    };
    struct bad_index:
        // Inherit from standard exceptions as requested during review.
        // public std::exception {
        public std::out_of_range {
        BOOST_UBLAS_EXPLICIT
        bad_index (const std::string &s = "exception: bad index"):
            std::out_of_range (s) {}
        // virtual const char *what () const throw () {
        //     return "exception: bad index";
        // }
        virtual void raise () {
#if ! defined (BOOST_NO_EXCEPTIONS) && defined (BOOST_UBLAS_USE_EXCEPTIONS)
            throw *this;
#else
#ifdef BOOST_NO_STDC_NAMESPACE
            ::abort ();
#else
            std::abort ();
#endif
#endif
        }
    };
    struct singular:
        // Inherit from standard exceptions as requested during review.
        // public std::exception {
        public std::runtime_error {
        BOOST_UBLAS_EXPLICIT
        singular (const std::string &s = "exception: singular"):
            std::runtime_error (s) {}
        // virtual const char *what () const throw () {
        //     return "exception: singular";
        // }
        virtual void raise () {
#if ! defined (BOOST_NO_EXCEPTIONS) && defined (BOOST_UBLAS_USE_EXCEPTIONS)
            throw *this;
#else
#ifdef BOOST_NO_STDC_NAMESPACE
            ::abort ();
#else
            std::abort ();
#endif
#endif
        }
    };
    struct non_real:
        // Inherit from standard exceptions as requested during review.
        // public std::exception {
        public std::domain_error {
        BOOST_UBLAS_EXPLICIT
        non_real (const std::string &s = "exception: non real"):
            std::domain_error (s) {}
        // virtual const char *what () const throw () {
        //     return "exception: non real";
        // }
        virtual void raise () {
#if ! defined (BOOST_NO_EXCEPTIONS) && defined (BOOST_UBLAS_USE_EXCEPTIONS)
            throw *this;
#else
#ifdef BOOST_NO_STDC_NAMESPACE
            ::abort ();
#else
            std::abort ();
#endif
#endif
        }
    };

#ifdef BOOST_UBLAS_BOUNDS_CHECK
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






