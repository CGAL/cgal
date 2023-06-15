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

#if defined(CGAL_ENABLE_DISABLE_ASSERTIONS_AT_RUNTIME)

#include <CGAL/tss.h>
namespace CGAL{
inline bool& get_use_assertions()
{
  CGAL_STATIC_THREAD_LOCAL_VARIABLE(bool, b, true);
  return b;
}
inline void set_use_assertions(bool b)
{
  get_use_assertions() = b;
}
}

#elif !defined(CGAL_USER_DEFINED_USE_ASSERTIONS)

namespace CGAL{
inline void set_use_assertions(bool){}
inline constexpr bool get_use_assertions(){return true;}
}

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
#  define CGAL_ASSERTIONS_ENABLED false
#  define CGAL_assertion(EX) (static_cast<void>(0))
#  define CGAL_destructor_assertion(EX) (static_cast<void>(0))
#  define CGAL_destructor_assertion_catch(CODE) CODE
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
#  define CGAL_ASSERTIONS_ENABLED true
#  define CGAL_assertion(EX) \
   ((!CGAL::get_use_assertions() || CGAL::possibly(EX))?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__))
#  if __cpp_lib_uncaught_exceptions || ( _MSVC_LANG >= 201703L )  // C++17
#    define CGAL_destructor_assertion(EX) \
     ((!CGAL::get_use_assertions() || CGAL::possibly(EX)||(std::uncaught_exceptions() > 0))?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__))
#  define CGAL_destructor_assertion_catch(CODE) try{ CODE } catch(...) { if(std::uncaught_exceptions() <= 0) throw; }
#  else // use C++03 `std::uncaught_exception()`
#    define CGAL_destructor_assertion(EX) \
     ((!CGAL::get_use_assertions() || CGAL::possibly(EX)||std::uncaught_exception())?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__))
#  define CGAL_destructor_assertion_catch(CODE) try{ CODE } catch(...) { if(!std::uncaught_exception()) throw; }
#  endif // use C++03 `std::uncaught_exception()`
#  define CGAL_assertion_msg(EX,MSG)                                    \
   ((!CGAL::get_use_assertions() || CGAL::possibly(EX))?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_assertion_code(CODE) CODE
#  define CGAL_assume(EX) CGAL_assertion(EX)
#  define CGAL_assume_code(CODE) CGAL_assertion_code(CODE)
#endif // no CGAL_NO_ASSERTIONS

#  ifdef CGAL_UNREACHABLE
#    define CGAL_unreachable() CGAL_UNREACHABLE()
#  else // not def CGAL_UNREACHABLE
#    define CGAL_unreachable()  CGAL_assertion(false)
#  endif // CGAL_UNREACHABLE

#if defined(CGAL_NO_ASSERTIONS) || !defined(CGAL_CHECK_EXACTNESS)
#  define CGAL_exactness_assertion(EX) (static_cast<void>(0))
#  define CGAL_exactness_assertion_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_exactness_assertion_code(CODE)
#else
#  define CGAL_exactness_assertion(EX) \
   ((!CGAL::get_use_assertions() || CGAL::possibly(EX))?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__))
#  define CGAL_exactness_assertion_msg(EX,MSG) \
   ((!CGAL::get_use_assertions() || CGAL::possibly(EX))?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_exactness_assertion_code(CODE) CODE
#endif // CGAL_NO_ASSERTIONS

#if defined(CGAL_NO_ASSERTIONS) || !defined(CGAL_CHECK_EXPENSIVE)
#  define CGAL_expensive_assertion(EX) (static_cast<void>(0))
#  define CGAL_expensive_assertion_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_expensive_assertion_code(CODE)
#else
#  define CGAL_expensive_assertion(EX) \
   ((!CGAL::get_use_assertions() || CGAL::possibly(EX))?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__))
#  define CGAL_expensive_assertion_msg(EX,MSG) \
   ((!CGAL::get_use_assertions() || CGAL::possibly(EX))?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_expensive_assertion_code(CODE) CODE
#endif // CGAL_NO_ASSERTIONS

#if defined(CGAL_NO_ASSERTIONS) || !defined(CGAL_CHECK_EXACTNESS) || !defined(CGAL_CHECK_EXPENSIVE)
#  define CGAL_expensive_exactness_assertion(EX) (static_cast<void>(0))
#  define CGAL_expensive_exactness_assertion_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_expensive_exactness_assertion_code(CODE)
#else
#  define CGAL_expensive_exactness_assertion(EX) \
   ((!CGAL::get_use_assertions() || CGAL::possibly(EX))?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__))
#  define CGAL_expensive_exactness_assertion_msg(EX,MSG) \
   ((!CGAL::get_use_assertions() || CGAL::possibly(EX))?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_expensive_exactness_assertion_code(CODE) CODE
#endif // CGAL_NO_ASSERTIONS


// preconditions
// -------------

#if defined(CGAL_NO_PRECONDITIONS)
#  define CGAL_PRECONDITIONS_ENABLED false
#  define CGAL_precondition(EX) (static_cast<void>(0))
#  define CGAL_precondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_precondition_code(CODE)
#else
#  define CGAL_PRECONDITIONS_ENABLED true
#  define CGAL_precondition(EX) \
   ((!CGAL::get_use_assertions() || CGAL::possibly(EX))?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_precondition_msg(EX,MSG) \
   ((!CGAL::get_use_assertions() || CGAL::possibly(EX))?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_precondition_code(CODE) CODE
#endif // CGAL_NO_PRECONDITIONS

#if defined(CGAL_NO_PRECONDITIONS) || !defined(CGAL_CHECK_EXACTNESS)
#  define CGAL_exactness_precondition(EX) (static_cast<void>(0))
#  define CGAL_exactness_precondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_exactness_precondition_code(CODE)
#else
#  define CGAL_exactness_precondition(EX) \
   ((!CGAL::get_use_assertions() || CGAL::possibly(EX))?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_exactness_precondition_msg(EX,MSG) \
   ((!CGAL::get_use_assertions() || CGAL::possibly(EX))?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_exactness_precondition_code(CODE) CODE
#endif // CGAL_NO_PRECONDITIONS

#if defined(CGAL_NO_PRECONDITIONS) || !defined(CGAL_CHECK_EXPENSIVE)
#  define CGAL_expensive_precondition(EX) (static_cast<void>(0))
#  define CGAL_expensive_precondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_expensive_precondition_code(CODE)
#else
#  define CGAL_expensive_precondition(EX) \
   ((!CGAL::get_use_assertions() || CGAL::possibly(EX))?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_expensive_precondition_msg(EX,MSG) \
   ((!CGAL::get_use_assertions() || CGAL::possibly(EX))?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_expensive_precondition_code(CODE) CODE
#endif // CGAL_NO_PRECONDITIONS

#if defined(CGAL_NO_PRECONDITIONS) || !defined(CGAL_CHECK_EXACTNESS) || !defined(CGAL_CHECK_EXPENSIVE)
#  define CGAL_expensive_exactness_precondition(EX) (static_cast<void>(0))
#  define CGAL_expensive_exactness_precondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_expensive_exactness_precondition_code(CODE)
#else
#  define CGAL_expensive_exactness_precondition(EX) \
   ((!CGAL::get_use_assertions() || CGAL::possibly(EX))?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_expensive_exactness_precondition_msg(EX,MSG) \
   ((!CGAL::get_use_assertions() || CGAL::possibly(EX))?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, MSG))
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
   ((!CGAL::get_use_assertions() || CGAL::possibly(EX))?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_postcondition_msg(EX,MSG) \
   ((!CGAL::get_use_assertions() || CGAL::possibly(EX))?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_postcondition_code(CODE) CODE
#endif // CGAL_NO_POSTCONDITIONS

#if defined(CGAL_NO_POSTCONDITIONS) || !defined(CGAL_CHECK_EXACTNESS)
#  define CGAL_exactness_postcondition(EX) (static_cast<void>(0))
#  define CGAL_exactness_postcondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_exactness_postcondition_code(CODE)
#else
#  define CGAL_exactness_postcondition(EX) \
   ((!CGAL::get_use_assertions() || CGAL::possibly(EX))?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_exactness_postcondition_msg(EX,MSG) \
   ((!CGAL::get_use_assertions() || CGAL::possibly(EX))?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_exactness_postcondition_code(CODE) CODE
#endif // CGAL_NO_POSTCONDITIONS

#if defined(CGAL_NO_POSTCONDITIONS) || !defined(CGAL_CHECK_EXPENSIVE)
#  define CGAL_expensive_postcondition(EX) (static_cast<void>(0))
#  define CGAL_expensive_postcondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_expensive_postcondition_code(CODE)
#else
#  define CGAL_expensive_postcondition(EX) \
   ((!CGAL::get_use_assertions() || CGAL::possibly(EX))?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_expensive_postcondition_msg(EX,MSG) \
   ((!CGAL::get_use_assertions() || CGAL::possibly(EX))?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_expensive_postcondition_code(CODE) CODE
#endif // CGAL_NO_POSTCONDITIONS

#if defined(CGAL_NO_POSTCONDITIONS) || !defined(CGAL_CHECK_EXACTNESS) || !defined(CGAL_CHECK_EXPENSIVE)
#  define CGAL_expensive_exactness_postcondition(EX) (static_cast<void>(0))
#  define CGAL_expensive_exactness_postcondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_expensive_exactness_postcondition_code(CODE)
#else
#  define CGAL_expensive_exactness_postcondition(EX) \
   ((!CGAL::get_use_assertions() || CGAL::possibly(EX))?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_expensive_exactness_postcondition_msg(EX,MSG) \
   ((!CGAL::get_use_assertions() || CGAL::possibly(EX))?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_expensive_exactness_postcondition_code(CODE) CODE
#endif // CGAL_NO_POSTCONDITIONS


// warnings
// --------

#if defined(CGAL_NO_WARNINGS)
#  define CGAL_warning(EX) (static_cast<void>(0))
#  define CGAL_destructor_warning(EX) (static_cast<void>(0))
#  define CGAL_warning_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_warning_code(CODE)
#else
#  define CGAL_warning(EX) \
   ((!CGAL::get_use_assertions() || CGAL::possibly(EX))?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__))
#  if __cpp_lib_uncaught_exceptions || ( _MSVC_LANG >= 201703L )  // C++17
#    define CGAL_destructor_warning(EX) \
     ((!CGAL::get_use_assertions() || CGAL::possibly(EX)||(std::uncaught_exceptions() > 0))?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__))
#  else // use C++03 `std::uncaught_exception()`
#    define CGAL_destructor_warning(EX) \
     ((!CGAL::get_use_assertions() || CGAL::possibly(EX)||std::uncaught_exception())?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__))
#  endif // use C++03 `std::uncaught_exception()`
#  define CGAL_warning_msg(EX,MSG) \
   ((!CGAL::get_use_assertions() || CGAL::possibly(EX))?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_warning_code(CODE) CODE
#endif // CGAL_NO_WARNINGS

#if defined(CGAL_NO_WARNINGS) || !defined(CGAL_CHECK_EXACTNESS)
#  define CGAL_exactness_warning(EX) (static_cast<void>(0))
#  define CGAL_exactness_warning_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_exactness_warning_code(CODE)
#else
#  define CGAL_exactness_warning(EX) \
   ((!CGAL::get_use_assertions() || CGAL::possibly(EX))?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__))
#  define CGAL_exactness_warning_msg(EX,MSG) \
   ((!CGAL::get_use_assertions() || CGAL::possibly(EX))?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_exactness_warning_code(CODE) CODE
#endif // CGAL_NO_WARNINGS

#if defined(CGAL_NO_WARNINGS) || !defined(CGAL_CHECK_EXPENSIVE)
#  define CGAL_expensive_warning(EX) (static_cast<void>(0))
#  define CGAL_expensive_warning_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_expensive_warning_code(CODE)
#else
#  define CGAL_expensive_warning(EX) \
   ((!CGAL::get_use_assertions() || CGAL::possibly(EX))?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__))
#  define CGAL_expensive_warning_msg(EX,MSG) \
   ((!CGAL::get_use_assertions() || CGAL::possibly(EX))?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_expensive_warning_code(CODE) CODE
#endif // CGAL_NO_WARNINGS

#if defined(CGAL_NO_WARNINGS) || !defined(CGAL_CHECK_EXACTNESS) || !defined(CGAL_CHECK_EXPENSIVE)
#  define CGAL_expensive_exactness_warning(EX) (static_cast<void>(0))
#  define CGAL_expensive_exactness_warning_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_expensive_exactness_warning_code(CODE)
#else
#  define CGAL_expensive_exactness_warning(EX) \
   ((!CGAL::get_use_assertions() || CGAL::possibly(EX))?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__))
#  define CGAL_expensive_exactness_warning_msg(EX,MSG) \
   ((!CGAL::get_use_assertions() || CGAL::possibly(EX))?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_expensive_exactness_warning_code(CODE) CODE
#endif // CGAL_NO_WARNINGS

// CGAL error
#define CGAL_error_msg(MSG) ::CGAL::assertion_fail( "", __FILE__, __LINE__, MSG )
#define CGAL_error()        ::CGAL::assertion_fail( "", __FILE__, __LINE__ )

} //namespace CGAL

// This comes last as it is dependent on the macros to be defined.
// But the macros need CGAL::possibly().
#include <CGAL/Uncertain.h>

#ifdef CGAL_HEADER_ONLY
#include <CGAL/assertions_impl.h>
#endif // CGAL_HEADER_ONLY

#endif // CGAL_ASSERTIONS_H
