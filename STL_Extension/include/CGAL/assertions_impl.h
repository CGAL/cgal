// Copyright (c) 1997
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
// Author(s)     : Geert-Jan Giezeman and Sven Schönherr

#ifdef CGAL_HEADER_ONLY
#define CGAL_INLINE_FUNCTION inline
#else
#define CGAL_INLINE_FUNCTION
#endif

#include <CGAL/config.h>
#include <CGAL/assertions.h>
#include <CGAL/assertions_behaviour.h>
#include <CGAL/exceptions.h>

#include <cstdlib>
#include <iostream>

namespace CGAL {

#ifdef CGAL_HEADER_ONLY

inline Failure_behaviour& get_static_error_behaviour()
{
  static Failure_behaviour _error_behaviour = THROW_EXCEPTION;
  return _error_behaviour;
}
inline Failure_behaviour& get_static_warning_behaviour()
{
  static Failure_behaviour _warning_behaviour = CONTINUE;
  return _warning_behaviour;
}

#else // CGAL_HEADER_ONLY

namespace {

// behaviour variables
// -------------------

Failure_behaviour _error_behaviour   = THROW_EXCEPTION;
Failure_behaviour _warning_behaviour = CONTINUE;

inline Failure_behaviour& get_static_error_behaviour()
{ return _error_behaviour; }
inline Failure_behaviour& get_static_warning_behaviour()
{ return _warning_behaviour; }

} // anonymous namespace

#endif // CGAL_HEADER_ONLY

namespace {

// standard error handlers
// -----------------------
CGAL_INLINE_FUNCTION
void
_standard_error_handler(
        const char* what,
        const char* expr,
        const char* file,
        int         line,
        const char* msg )
{
#if defined(__GNUG__) && !defined(__llvm__)
    // After g++ 3.4, std::terminate defaults to printing to std::cerr itself.
    if (get_static_error_behaviour() == THROW_EXCEPTION)
        return;
#endif
    std::cerr << "CGAL error: " << what << " violation!" << std::endl
         << "Expression : " << expr << std::endl
         << "File       : " << file << std::endl
         << "Line       : " << line << std::endl
         << "Explanation: " << msg << std::endl
         << "Refer to the bug-reporting instructions at https://www.cgal.org/bug_report.html"
         << std::endl;
}


// standard warning handler
// ------------------------
CGAL_INLINE_FUNCTION
void
_standard_warning_handler( const char *,
                          const char* expr,
                          const char* file,
                          int         line,
                          const char* msg )
{
#if defined(__GNUG__) && !defined(__llvm__)
    // After g++ 3.4, std::terminate defaults to printing to std::cerr itself.
    if (get_static_warning_behaviour() == THROW_EXCEPTION)
        return;
#endif
    std::cerr << "CGAL warning: check violation!" << std::endl
         << "Expression : " << expr << std::endl
         << "File       : " << file << std::endl
         << "Line       : " << line << std::endl
         << "Explanation: " << msg << std::endl
         << "Refer to the bug-reporting instructions at https://www.cgal.org/bug_report.html"
         << std::endl;
}

} // anonymous namespace

#ifdef CGAL_HEADER_ONLY

inline Failure_function& get_static_error_handler()
{
  static Failure_function _error_handler = _standard_error_handler;
  return _error_handler;
}
inline Failure_function& get_static_warning_handler()
{
  static Failure_function _warning_handler = _standard_warning_handler;
  return _warning_handler;
}

#else // CGAL_HEADER_ONLY

namespace {

// default handler settings
// ------------------------
Failure_function _error_handler   = _standard_error_handler;
Failure_function _warning_handler = _standard_warning_handler;

inline Failure_function& get_static_error_handler()
{ return _error_handler; }
inline Failure_function& get_static_warning_handler()
{ return _warning_handler; }

} // anonymous namespace

#endif // CGAL_HEADER_ONLY

// failure functions
// -----------------
CGAL_INLINE_FUNCTION
void
assertion_fail( const char* expr,
                const char* file,
                int         line,
                const char* msg)
{
    get_static_error_handler()("assertion", expr, file, line, msg);
    switch (get_static_error_behaviour()) {
    case ABORT:
        std::abort();
    case EXIT:
        std::exit(1);  // EXIT_FAILURE
    case EXIT_WITH_SUCCESS:
        std::exit(0);  // EXIT_SUCCESS
    case CONTINUE: // The CONTINUE case should not be used anymore.
    case THROW_EXCEPTION:
    default:
        throw Assertion_exception("CGAL", expr, file, line, msg);
    }
}

CGAL_INLINE_FUNCTION
void
precondition_fail( const char* expr,
                   const char* file,
                   int         line,
                   const char* msg)
{
    get_static_error_handler()("precondition", expr, file, line, msg);
    switch (get_static_error_behaviour()) {
    case ABORT:
        std::abort();
    case EXIT:
        std::exit(1);  // EXIT_FAILURE
    case EXIT_WITH_SUCCESS:
        std::exit(0);  // EXIT_SUCCESS
    case CONTINUE:
    case THROW_EXCEPTION:
    default:
        throw Precondition_exception("CGAL", expr, file, line, msg);
    }
}

CGAL_INLINE_FUNCTION
void
postcondition_fail(const char* expr,
                   const char* file,
                   int         line,
                   const char* msg)
{
    get_static_error_handler()("postcondition", expr, file, line, msg);
    switch (get_static_error_behaviour()) {
    case ABORT:
        std::abort();
    case EXIT:
        std::exit(1);  // EXIT_FAILURE
    case EXIT_WITH_SUCCESS:
        std::exit(0);  // EXIT_SUCCESS
    case CONTINUE:
    case THROW_EXCEPTION:
    default:
        throw Postcondition_exception("CGAL", expr, file, line, msg);
    }
}


// warning function
// ----------------
CGAL_INLINE_FUNCTION
void
warning_fail( const char* expr,
              const char* file,
              int         line,
              const char* msg)
{
    get_static_warning_handler()("warning", expr, file, line, msg);
    switch (get_static_warning_behaviour()) {
    case ABORT:
        std::abort();
    case EXIT:
        std::exit(1);  // EXIT_FAILURE
    case EXIT_WITH_SUCCESS:
        std::exit(0);  // EXIT_SUCCESS
    case THROW_EXCEPTION:
        throw Warning_exception("CGAL", expr, file, line, msg);
    case CONTINUE:
        ;
    }
}


// error handler set functions
// ---------------------------
CGAL_INLINE_FUNCTION
Failure_function
set_error_handler( Failure_function handler)
{
    Failure_function result = get_static_error_handler();
    get_static_error_handler() = handler;
    return result;
}

CGAL_INLINE_FUNCTION
Failure_function
set_warning_handler( Failure_function handler)
{
    Failure_function result = get_static_warning_handler();
    get_static_warning_handler() = handler;
    return result;
}

CGAL_INLINE_FUNCTION
Failure_behaviour
set_error_behaviour(Failure_behaviour eb)
{
    Failure_behaviour result = get_static_error_behaviour();
    get_static_error_behaviour() = eb;
    return result;
}

CGAL_INLINE_FUNCTION
Failure_behaviour
set_warning_behaviour(Failure_behaviour eb)
{
    Failure_behaviour result = get_static_warning_behaviour();
    get_static_warning_behaviour() = eb;
    return result;
}

} //namespace CGAL
