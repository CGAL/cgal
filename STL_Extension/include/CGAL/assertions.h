// Copyright (c) 1999  Utrecht University (The Netherlands),
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
// $URL$
// $Id$
//
//
// Author(s)     : Geert-Jan Giezeman and Sven Schoenherr

#include <CGAL/config.h>

#ifndef CGAL_ASSERTIONS_H
#define CGAL_ASSERTIONS_H

#include <CGAL/Uncertain.h>

CGAL_BEGIN_NAMESPACE

// types
// =====

enum Failure_behaviour { ABORT, EXIT, EXIT_WITH_SUCCESS, CONTINUE,
                         THROW_EXCEPTION };

// function declarations
// =====================
// failure functions
// -----------------
void assertion_fail      ( const char*, const char*, int, const char* = "");
void precondition_fail   ( const char*, const char*, int, const char* = "");
void postcondition_fail  ( const char*, const char*, int, const char* = "");

// warning function
// ----------------
void warning_fail( const char*, const char*, int, const char* = "");


// macro definitions
// =================
// assertions
// ----------

#ifdef NDEBUG
#  define CGAL_NDEBUG
#endif

#ifdef CGAL_NDEBUG
#  define CGAL_NO_ASSERTIONS
#  define CGAL_NO_PRECONDITIONS
#  define CGAL_NO_POSTCONDITIONS
#  define CGAL_NO_WARNINGS
#endif

#if defined(CGAL_NO_ASSERTIONS)
#  define CGAL_assertion(EX) (static_cast<void>(0))
#  define CGAL_assertion_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_assertion_code(CODE)
#else
#  define CGAL_assertion(EX) \
   (CGAL::certainly(EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__))
#  define CGAL_assertion_msg(EX,MSG) \
   (CGAL::certainly(EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_assertion_code(CODE) CODE
#endif // CGAL_NO_ASSERTIONS

#if defined(CGAL_NO_ASSERTIONS) || !defined(CGAL_CHECK_EXACTNESS)
#  define CGAL_exactness_assertion(EX) (static_cast<void>(0))
#  define CGAL_exactness_assertion_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_exactness_assertion_code(CODE)
#else
#  define CGAL_exactness_assertion(EX) \
   (CGAL::certainly(EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__))
#  define CGAL_exactness_assertion_msg(EX,MSG) \
   (CGAL::certainly(EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_exactness_assertion_code(CODE) CODE
#endif // CGAL_NO_ASSERTIONS

#if defined(CGAL_NO_ASSERTIONS) || !defined(CGAL_CHECK_EXPENSIVE)
#  define CGAL_expensive_assertion(EX) (static_cast<void>(0))
#  define CGAL_expensive_assertion_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_expensive_assertion_code(CODE)
#else
#  define CGAL_expensive_assertion(EX) \
   (CGAL::certainly(EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__))
#  define CGAL_expensive_assertion_msg(EX,MSG) \
   (CGAL::certainly(EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_expensive_assertion_code(CODE) CODE
#endif // CGAL_NO_ASSERTIONS

#if defined(CGAL_NO_ASSERTIONS) || !defined(CGAL_CHECK_EXACTNESS) || !defined(CGAL_CHECK_EXPENSIVE)
#  define CGAL_expensive_exactness_assertion(EX) (static_cast<void>(0))
#  define CGAL_expensive_exactness_assertion_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_expensive_exactness_assertion_code(CODE)
#else
#  define CGAL_expensive_exactness_assertion(EX) \
   (CGAL::certainly(EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__))
#  define CGAL_expensive_exactness_assertion_msg(EX,MSG) \
   (CGAL::certainly(EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
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
   (CGAL::certainly(EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_precondition_msg(EX,MSG) \
   (CGAL::certainly(EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_precondition_code(CODE) CODE
#endif // CGAL_NO_PRECONDITIONS

#if defined(CGAL_NO_PRECONDITIONS) || !defined(CGAL_CHECK_EXACTNESS)
#  define CGAL_exactness_precondition(EX) (static_cast<void>(0))
#  define CGAL_exactness_precondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_exactness_precondition_code(CODE)
#else
#  define CGAL_exactness_precondition(EX) \
   (CGAL::certainly(EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_exactness_precondition_msg(EX,MSG) \
   (CGAL::certainly(EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_exactness_precondition_code(CODE) CODE
#endif // CGAL_NO_PRECONDITIONS

#if defined(CGAL_NO_PRECONDITIONS) || !defined(CGAL_CHECK_EXPENSIVE)
#  define CGAL_expensive_precondition(EX) (static_cast<void>(0))
#  define CGAL_expensive_precondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_expensive_precondition_code(CODE)
#else
#  define CGAL_expensive_precondition(EX) \
   (CGAL::certainly(EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_expensive_precondition_msg(EX,MSG) \
   (CGAL::certainly(EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_expensive_precondition_code(CODE) CODE
#endif // CGAL_NO_PRECONDITIONS

#if defined(CGAL_NO_PRECONDITIONS) || !defined(CGAL_CHECK_EXACTNESS) || !defined(CGAL_CHECK_EXPENSIVE) 
#  define CGAL_expensive_exactness_precondition(EX) (static_cast<void>(0))
#  define CGAL_expensive_exactness_precondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_expensive_exactness_precondition_code(CODE)
#else
#  define CGAL_expensive_exactness_precondition(EX) \
   (CGAL::certainly(EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_expensive_exactness_precondition_msg(EX,MSG) \
   (CGAL::certainly(EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, MSG))
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
   (CGAL::certainly(EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_postcondition_msg(EX,MSG) \
   (CGAL::certainly(EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_postcondition_code(CODE) CODE
#endif // CGAL_NO_POSTCONDITIONS

#if defined(CGAL_NO_POSTCONDITIONS) || !defined(CGAL_CHECK_EXACTNESS)
#  define CGAL_exactness_postcondition(EX) (static_cast<void>(0))
#  define CGAL_exactness_postcondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_exactness_postcondition_code(CODE)
#else
#  define CGAL_exactness_postcondition(EX) \
   (CGAL::certainly(EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_exactness_postcondition_msg(EX,MSG) \
   (CGAL::certainly(EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_exactness_postcondition_code(CODE) CODE
#endif // CGAL_NO_POSTCONDITIONS

#if defined(CGAL_NO_POSTCONDITIONS) || !defined(CGAL_CHECK_EXPENSIVE)
#  define CGAL_expensive_postcondition(EX) (static_cast<void>(0))
#  define CGAL_expensive_postcondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_expensive_postcondition_code(CODE)
#else
#  define CGAL_expensive_postcondition(EX) \
   (CGAL::certainly(EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_expensive_postcondition_msg(EX,MSG) \
   (CGAL::certainly(EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_expensive_postcondition_code(CODE) CODE
#endif // CGAL_NO_POSTCONDITIONS

#if defined(CGAL_NO_POSTCONDITIONS) || !defined(CGAL_CHECK_EXACTNESS) || !defined(CGAL_CHECK_EXPENSIVE)
#  define CGAL_expensive_exactness_postcondition(EX) (static_cast<void>(0))
#  define CGAL_expensive_exactness_postcondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_expensive_exactness_postcondition_code(CODE)
#else
#  define CGAL_expensive_exactness_postcondition(EX) \
   (CGAL::certainly(EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_expensive_exactness_postcondition_msg(EX,MSG) \
   (CGAL::certainly(EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, MSG))
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
   (CGAL::certainly(EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__))
#  define CGAL_warning_msg(EX,MSG) \
   (CGAL::certainly(EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_warning_code(CODE) CODE
#endif // CGAL_NO_WARNINGS

#if defined(CGAL_NO_WARNINGS) || !defined(CGAL_CHECK_EXACTNESS)
#  define CGAL_exactness_warning(EX) (static_cast<void>(0))
#  define CGAL_exactness_warning_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_exactness_warning_code(CODE)
#else
#  define CGAL_exactness_warning(EX) \
   (CGAL::certainly(EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__))
#  define CGAL_exactness_warning_msg(EX,MSG) \
   (CGAL::certainly(EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_exactness_warning_code(CODE) CODE
#endif // CGAL_NO_WARNINGS

#if defined(CGAL_NO_WARNINGS) || !defined(CGAL_CHECK_EXPENSIVE)
#  define CGAL_expensive_warning(EX) (static_cast<void>(0))
#  define CGAL_expensive_warning_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_expensive_warning_code(CODE)
#else
#  define CGAL_expensive_warning(EX) \
   (CGAL::certainly(EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__))
#  define CGAL_expensive_warning_msg(EX,MSG) \
   (CGAL::certainly(EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_expensive_warning_code(CODE) CODE
#endif // CGAL_NO_WARNINGS

#if defined(CGAL_NO_WARNINGS) || !defined(CGAL_CHECK_EXACTNESS) || !defined(CGAL_CHECK_EXPENSIVE)
#  define CGAL_expensive_exactness_warning(EX) (static_cast<void>(0))
#  define CGAL_expensive_exactness_warning_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_expensive_exactness_warning_code(CODE)
#else
#  define CGAL_expensive_exactness_warning(EX) \
   (CGAL::certainly(EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__))
#  define CGAL_expensive_exactness_warning_msg(EX,MSG) \
   (CGAL::certainly(EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_expensive_exactness_warning_code(CODE) CODE
#endif // CGAL_NO_WARNINGS

// CGAL error
#define CGAL_error_msg(MSG) ::CGAL::assertion_fail( "", __FILE__, __LINE__, MSG )
#define CGAL_error()        ::CGAL::assertion_fail( "", __FILE__, __LINE__ )

// failure handler declarations
// ==========================
// failure handler
// ---------------
typedef
    void
    (*Failure_function)(
        const char*, const char*, const char*, int, const char*);

Failure_function
set_error_handler( Failure_function handler);

Failure_function
set_warning_handler( Failure_function handler);

// failure behaviour handler
// -------------------------
Failure_behaviour
set_error_behaviour(Failure_behaviour eb);

Failure_behaviour
set_warning_behaviour(Failure_behaviour eb);

CGAL_END_NAMESPACE

#endif // CGAL_ASSERTIONS_H
