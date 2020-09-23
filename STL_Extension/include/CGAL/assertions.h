// Copyright (c) 1999
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
// Author(s)     : Geert-Jan Giezeman and Sven Schoenherr

#include <CGAL/config.h>

#ifndef CGAL_ASSERTIONS_H
#define CGAL_ASSERTIONS_H

#include <CGAL/export/CGAL.h>

// #include <CGAL/assertions_behaviour.h> // for backward compatibility

#ifdef NDEBUG
#  ifndef CGAL_NDEBUG
#    define CGAL_NDEBUG
#  endif
#endif

// The macro `CGAL_DEBUG` allows to force CGAL assertions, even if `NDEBUG`
// is defined,
#ifdef CGAL_DEBUG
#  ifdef CGAL_NDEBUG
#    undef CGAL_NDEBUG
#  endif
#endif

#ifdef CGAL_NDEBUG
#  define CGAL_NO_ASSERTIONS
#  define CGAL_NO_PRECONDITIONS
#  define CGAL_NO_POSTCONDITIONS
#  define CGAL_NO_WARNINGS
#endif

#ifdef CGAL_CFG_NO_CPP0X_STATIC_ASSERT
#include <boost/static_assert.hpp>
#endif

namespace CGAL {

// function declarations
// =====================
// failure functions
// -----------------
[[noreturn]] CGAL_EXPORT void assertion_fail      ( const char*, const char*, int, const char* = "") ;
[[noreturn]] CGAL_EXPORT void precondition_fail   ( const char*, const char*, int, const char* = "") ;
[[noreturn]] CGAL_EXPORT void postcondition_fail  ( const char*, const char*, int, const char* = "") ;

// warning function
// ----------------
CGAL_EXPORT
void warning_fail( const char*, const char*, int, const char* = "");


// The following declarations are needed first, before assertions are used.
// ---------------------------------------------------------------------
template < typename T > class Uncertain;
inline bool possibly(bool b);
inline bool possibly(Uncertain<bool> c);


// macro definitions
// =================
// assertions
// ----------

#if defined(CGAL_NO_ASSERTIONS)
#  define CGAL_assertion(EX) (static_cast<void>(0))
#  define CGAL_destructor_assertion(EX) (static_cast<void>(0))
#  define CGAL_assertion_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_assertion_code(CODE)
#  ifdef CGAL_ASSUME
#    define CGAL_assume(EX) CGAL_ASSUME(EX)
#    define CGAL_assume_code(CODE) CODE
#  else // not def CGAL_ASSUME
#    define CGAL_assume(EX)  CGAL_assertion(EX)
#    define CGAL_assume_code(CODE) CGAL_assertion_code(CODE)
#  endif // not def CGAL_ASSUME
#else // no CGAL_NO_ASSERTIONS
#  define CGAL_assertion(EX) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__))
#  if __cpp_lib_uncaught_exceptions || ( _MSVC_LANG >= 201703L )  // C++17
#    define CGAL_destructor_assertion(EX) \
     (CGAL::possibly(EX)||(std::uncaught_exceptions() > 0)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__))
#  else // use C++03 `std::uncaught_exception()`
#    define CGAL_destructor_assertion(EX) \
     (CGAL::possibly(EX)||std::uncaught_exception()?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__))
#  endif // use C++03 `std::uncaught_exception()`
#  define CGAL_assertion_msg(EX,MSG)                                    \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_assertion_code(CODE) CODE
#  define CGAL_assume(EX) CGAL_assertion(EX)
#  define CGAL_assume_code(CODE) CGAL_assertion_code(CODE)
#endif // no CGAL_NO_ASSERTIONS

#ifndef CGAL_CFG_NO_CPP0X_STATIC_ASSERT

#  if defined(CGAL_NO_ASSERTIONS)

#    define CGAL_static_assertion(EX) \
     static_assert(true, "")

#    define CGAL_static_assertion_msg(EX,MSG) \
     static_assert(true, "")

#  else // no CGAL_NO_ASSERTIONS

#    define CGAL_static_assertion(EX) \
     static_assert(EX, #EX)

#    define CGAL_static_assertion_msg(EX,MSG) \
     static_assert(EX, MSG)

#  endif // no CGAL_NO_ASSERTIONS

#else // if CGAL_CFG_NO_CPP0X_STATIC_ASSERT is true

#  if defined(CGAL_NO_ASSERTIONS)

#    define CGAL_static_assertion(EX) \
     BOOST_STATIC_ASSERT(true) CGAL_UNUSED

#    define CGAL_static_assertion_msg(EX,MSG) \
     BOOST_STATIC_ASSERT(true) CGAL_UNUSED

#  else // no CGAL_NO_ASSERTIONS

#    define CGAL_static_assertion(EX) \
     BOOST_STATIC_ASSERT(EX) CGAL_UNUSED

#    define CGAL_static_assertion_msg(EX,MSG) \
     BOOST_STATIC_ASSERT(EX) CGAL_UNUSED

#  endif // no CGAL_NO_ASSERTIONS

#endif // if CGAL_CFG_NO_CPP0X_STATIC_ASSERT is true

#if defined(CGAL_NO_ASSERTIONS) || !defined(CGAL_CHECK_EXACTNESS)
#  define CGAL_exactness_assertion(EX) (static_cast<void>(0))
#  define CGAL_exactness_assertion_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_exactness_assertion_code(CODE)
#else
#  define CGAL_exactness_assertion(EX) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__))
#  define CGAL_exactness_assertion_msg(EX,MSG) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_exactness_assertion_code(CODE) CODE
#endif // CGAL_NO_ASSERTIONS

#if defined(CGAL_NO_ASSERTIONS) || !defined(CGAL_CHECK_EXPENSIVE)
#  define CGAL_expensive_assertion(EX) (static_cast<void>(0))
#  define CGAL_expensive_assertion_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_expensive_assertion_code(CODE)
#else
#  define CGAL_expensive_assertion(EX) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__))
#  define CGAL_expensive_assertion_msg(EX,MSG) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_expensive_assertion_code(CODE) CODE
#endif // CGAL_NO_ASSERTIONS

#if defined(CGAL_NO_ASSERTIONS) || !defined(CGAL_CHECK_EXACTNESS) || !defined(CGAL_CHECK_EXPENSIVE)
#  define CGAL_expensive_exactness_assertion(EX) (static_cast<void>(0))
#  define CGAL_expensive_exactness_assertion_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_expensive_exactness_assertion_code(CODE)
#else
#  define CGAL_expensive_exactness_assertion(EX) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__))
#  define CGAL_expensive_exactness_assertion_msg(EX,MSG) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_expensive_exactness_assertion_code(CODE) CODE
#endif // CGAL_NO_ASSERTIONS


// preconditions
// -------------

#if defined(CGAL_NO_PRECONDITIONS)
#  define CGAL_precondition(EX) (static_cast<void>(0))
#  define CGAL_precondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_precondition_code(CODE)
#else
#  define CGAL_precondition(EX) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_precondition_msg(EX,MSG) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_precondition_code(CODE) CODE
#endif // CGAL_NO_PRECONDITIONS

#if defined(CGAL_NO_PRECONDITIONS) || !defined(CGAL_CHECK_EXACTNESS)
#  define CGAL_exactness_precondition(EX) (static_cast<void>(0))
#  define CGAL_exactness_precondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_exactness_precondition_code(CODE)
#else
#  define CGAL_exactness_precondition(EX) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_exactness_precondition_msg(EX,MSG) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_exactness_precondition_code(CODE) CODE
#endif // CGAL_NO_PRECONDITIONS

#if defined(CGAL_NO_PRECONDITIONS) || !defined(CGAL_CHECK_EXPENSIVE)
#  define CGAL_expensive_precondition(EX) (static_cast<void>(0))
#  define CGAL_expensive_precondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_expensive_precondition_code(CODE)
#else
#  define CGAL_expensive_precondition(EX) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_expensive_precondition_msg(EX,MSG) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_expensive_precondition_code(CODE) CODE
#endif // CGAL_NO_PRECONDITIONS

#if defined(CGAL_NO_PRECONDITIONS) || !defined(CGAL_CHECK_EXACTNESS) || !defined(CGAL_CHECK_EXPENSIVE)
#  define CGAL_expensive_exactness_precondition(EX) (static_cast<void>(0))
#  define CGAL_expensive_exactness_precondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_expensive_exactness_precondition_code(CODE)
#else
#  define CGAL_expensive_exactness_precondition(EX) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_expensive_exactness_precondition_msg(EX,MSG) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_expensive_exactness_precondition_code(CODE) CODE
#endif // CGAL_NO_PRECONDITIONS


// postconditions
// --------------

#if defined(CGAL_NO_POSTCONDITIONS)
#  define CGAL_postcondition(EX) (static_cast<void>(0))
#  define CGAL_postcondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_postcondition_code(CODE)
#else
#  define CGAL_postcondition(EX) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_postcondition_msg(EX,MSG) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_postcondition_code(CODE) CODE
#endif // CGAL_NO_POSTCONDITIONS

#if defined(CGAL_NO_POSTCONDITIONS) || !defined(CGAL_CHECK_EXACTNESS)
#  define CGAL_exactness_postcondition(EX) (static_cast<void>(0))
#  define CGAL_exactness_postcondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_exactness_postcondition_code(CODE)
#else
#  define CGAL_exactness_postcondition(EX) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_exactness_postcondition_msg(EX,MSG) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_exactness_postcondition_code(CODE) CODE
#endif // CGAL_NO_POSTCONDITIONS

#if defined(CGAL_NO_POSTCONDITIONS) || !defined(CGAL_CHECK_EXPENSIVE)
#  define CGAL_expensive_postcondition(EX) (static_cast<void>(0))
#  define CGAL_expensive_postcondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_expensive_postcondition_code(CODE)
#else
#  define CGAL_expensive_postcondition(EX) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_expensive_postcondition_msg(EX,MSG) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_expensive_postcondition_code(CODE) CODE
#endif // CGAL_NO_POSTCONDITIONS

#if defined(CGAL_NO_POSTCONDITIONS) || !defined(CGAL_CHECK_EXACTNESS) || !defined(CGAL_CHECK_EXPENSIVE)
#  define CGAL_expensive_exactness_postcondition(EX) (static_cast<void>(0))
#  define CGAL_expensive_exactness_postcondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_expensive_exactness_postcondition_code(CODE)
#else
#  define CGAL_expensive_exactness_postcondition(EX) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_expensive_exactness_postcondition_msg(EX,MSG) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_expensive_exactness_postcondition_code(CODE) CODE
#endif // CGAL_NO_POSTCONDITIONS


// warnings
// --------

#if defined(CGAL_NO_WARNINGS)
#  define CGAL_warning(EX) (static_cast<void>(0))
#  define CGAL_warning_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_warning_code(CODE)
#else
#  define CGAL_warning(EX) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__))
#  define CGAL_warning_msg(EX,MSG) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_warning_code(CODE) CODE
#endif // CGAL_NO_WARNINGS

#if defined(CGAL_NO_WARNINGS) || !defined(CGAL_CHECK_EXACTNESS)
#  define CGAL_exactness_warning(EX) (static_cast<void>(0))
#  define CGAL_exactness_warning_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_exactness_warning_code(CODE)
#else
#  define CGAL_exactness_warning(EX) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__))
#  define CGAL_exactness_warning_msg(EX,MSG) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_exactness_warning_code(CODE) CODE
#endif // CGAL_NO_WARNINGS

#if defined(CGAL_NO_WARNINGS) || !defined(CGAL_CHECK_EXPENSIVE)
#  define CGAL_expensive_warning(EX) (static_cast<void>(0))
#  define CGAL_expensive_warning_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_expensive_warning_code(CODE)
#else
#  define CGAL_expensive_warning(EX) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__))
#  define CGAL_expensive_warning_msg(EX,MSG) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_expensive_warning_code(CODE) CODE
#endif // CGAL_NO_WARNINGS

#if defined(CGAL_NO_WARNINGS) || !defined(CGAL_CHECK_EXACTNESS) || !defined(CGAL_CHECK_EXPENSIVE)
#  define CGAL_expensive_exactness_warning(EX) (static_cast<void>(0))
#  define CGAL_expensive_exactness_warning_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_expensive_exactness_warning_code(CODE)
#else
#  define CGAL_expensive_exactness_warning(EX) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__))
#  define CGAL_expensive_exactness_warning_msg(EX,MSG) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_expensive_exactness_warning_code(CODE) CODE
#endif // CGAL_NO_WARNINGS

// CGAL error
#define CGAL_error_msg(MSG) ::CGAL::assertion_fail( "", __FILE__, __LINE__, MSG )
#define CGAL_error()        ::CGAL::assertion_fail( "", __FILE__, __LINE__ )

} //namespace CGAL

// This comes last as it is dependant on the macros to be defined.
// But the macros need CGAL::possibly().
#include <CGAL/Uncertain.h>

#ifdef CGAL_HEADER_ONLY
#include <CGAL/assertions_impl.h>
#endif // CGAL_HEADER_ONLY

#endif // CGAL_ASSERTIONS_H
