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
// Author(s)     : Geert-Jan Giezeman, Sven Schoenherr
//
// Generated from script create_assertions.sh


// Note that this header file is intentionnaly not protected with a
// macro (as <cassert>). Calling it a second time with another value
// for NDEBUG for example must make a difference.

#include <CGAL/assertions.h>

// macro definitions
// =================
// assertions
// ----------

#undef CGAL_polygon_assertion
#undef CGAL_polygon_assertion_msg
#undef CGAL_polygon_assertion_code


#ifndef CGAL_POLYGON_ASSERTIONS_H
#define CGAL_POLYGON_ASSERTIONS_H
#if defined(CGAL_POLYGON_NO_ASSERTIONS) || defined(CGAL_NO_ASSERTIONS) \
|| defined(NDEBUG)
namespace CGAL{
inline void set_use_polygon_assertions(bool){}
inline bool get_use_polygon_assertions(){return true;}
}
#else
#include <CGAL/tss.h>
namespace CGAL{
inline bool& get_use_polygon_assertions()
{
  CGAL_STATIC_THREAD_LOCAL_VARIABLE(bool, b, true);
  return b;
}
inline void set_use_polygon_assertions(bool b)
{
  get_use_polygon_assertions() = b;
}
}
#endif
#endif // CGAL_POLYGON_ASSERTIONS_H

#if defined(CGAL_POLYGON_NO_ASSERTIONS) || defined(CGAL_NO_ASSERTIONS) \
  || defined(NDEBUG)
#  define CGAL_polygon_assertion(EX) (static_cast<void>(0))
#  define CGAL_polygon_assertion_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_polygon_assertion_code(CODE)
#else
#  define CGAL_polygon_assertion(EX) \
   ((!get_use_polygon_assertions() || CGAL::possibly(EX))?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__))
#  define CGAL_polygon_assertion_msg(EX,MSG) \
   ((!get_use_polygon_assertions() || CGAL::possibly(EX))?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_polygon_assertion_code(CODE) CODE
#  define CGAL_polygon_assertions 1
#endif // CGAL_POLYGON_NO_ASSERTIONS


#undef CGAL_polygon_exactness_assertion
#undef CGAL_polygon_exactness_assertion_msg
#undef CGAL_polygon_exactness_assertion_code

#if defined(CGAL_POLYGON_NO_ASSERTIONS) || defined(CGAL_NO_ASSERTIONS) \
  || (!defined(CGAL_POLYGON_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || defined(NDEBUG)
#  define CGAL_polygon_exactness_assertion(EX) (static_cast<void>(0))
#  define CGAL_polygon_exactness_assertion_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_polygon_exactness_assertion_code(CODE)
#else
#  define CGAL_polygon_exactness_assertion(EX) \
   ((!get_use_polygon_assertions() || CGAL::possibly(EX))?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__))
#  define CGAL_polygon_exactness_assertion_msg(EX,MSG) \
   ((!get_use_polygon_assertions() || CGAL::possibly(EX))?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_polygon_exactness_assertion_code(CODE) CODE
#  define CGAL_polygon_exactness_assertions 1
#endif // CGAL_POLYGON_NO_ASSERTIONS


#undef CGAL_polygon_expensive_assertion
#undef CGAL_polygon_expensive_assertion_msg
#undef CGAL_polygon_expensive_assertion_code

#if defined(CGAL_POLYGON_NO_ASSERTIONS) \
  || defined(CGAL_NO_ASSERTIONS) \
  || (!defined(CGAL_POLYGON_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_polygon_expensive_assertion(EX) (static_cast<void>(0))
#  define CGAL_polygon_expensive_assertion_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_polygon_expensive_assertion_code(CODE)
#else
#  define CGAL_polygon_expensive_assertion(EX) \
   ((!get_use_polygon_assertions() || CGAL::possibly(EX))?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__))
#  define CGAL_polygon_expensive_assertion_msg(EX,MSG) \
   ((!get_use_polygon_assertions() || CGAL::possibly(EX))?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_polygon_expensive_assertion_code(CODE) CODE
#  define CGAL_polygon_expensive_assertions 1
#endif // CGAL_POLYGON_NO_ASSERTIONS


#undef CGAL_polygon_expensive_exactness_assertion
#undef CGAL_polygon_expensive_exactness_assertion_msg
#undef CGAL_polygon_expensive_exactness_assertion_code

#if defined(CGAL_POLYGON_NO_ASSERTIONS) || defined(CGAL_NO_ASSERTIONS) \
  || (!defined(CGAL_POLYGON_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || (!defined(CGAL_POLYGON_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_polygon_expensive_exactness_assertion(EX) (static_cast<void>(0))
#  define CGAL_polygon_expensive_exactness_assertion_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_polygon_expensive_exactness_assertion_code(CODE)
#else
#  define CGAL_polygon_expensive_exactness_assertion(EX) \
   ((!get_use_polygon_assertions() || CGAL::possibly(EX))?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__))
#  define CGAL_polygon_expensive_exactness_assertion_msg(EX,MSG) \
   ((!get_use_polygon_assertions() || CGAL::possibly(EX))?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_polygon_expensive_exactness_assertion_code(CODE) CODE
#  define CGAL_polygon_expensive_exactness_assertions 1
#endif // CGAL_POLYGON_NO_ASSERTIONS


// preconditions
// -------------

#undef CGAL_polygon_precondition
#undef CGAL_polygon_precondition_msg
#undef CGAL_polygon_precondition_code

#if defined(CGAL_POLYGON_NO_PRECONDITIONS) || defined(CGAL_NO_PRECONDITIONS) \
  || defined(NDEBUG)
#  define CGAL_polygon_precondition(EX) (static_cast<void>(0))
#  define CGAL_polygon_precondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_polygon_precondition_code(CODE)
#else
#  define CGAL_polygon_precondition(EX) \
   ((!get_use_polygon_assertions() || CGAL::possibly(EX))?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_polygon_precondition_msg(EX,MSG) \
   ((!get_use_polygon_assertions() || CGAL::possibly(EX))?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_polygon_precondition_code(CODE) CODE
#  define CGAL_polygon_preconditions 1
#endif // CGAL_POLYGON_NO_PRECONDITIONS


#undef CGAL_polygon_exactness_precondition
#undef CGAL_polygon_exactness_precondition_msg
#undef CGAL_polygon_exactness_precondition_code

#if defined(CGAL_POLYGON_NO_PRECONDITIONS) || defined(CGAL_NO_PRECONDITIONS) \
  || (!defined(CGAL_POLYGON_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || defined(NDEBUG)
#  define CGAL_polygon_exactness_precondition(EX) (static_cast<void>(0))
#  define CGAL_polygon_exactness_precondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_polygon_exactness_precondition_code(CODE)
#else
#  define CGAL_polygon_exactness_precondition(EX) \
   ((!get_use_polygon_assertions() || CGAL::possibly(EX))?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_polygon_exactness_precondition_msg(EX,MSG) \
   ((!get_use_polygon_assertions() || CGAL::possibly(EX))?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_polygon_exactness_precondition_code(CODE) CODE
#  define CGAL_polygon_exactness_preconditions 1
#endif // CGAL_POLYGON_NO_PRECONDITIONS


#undef CGAL_polygon_expensive_precondition
#undef CGAL_polygon_expensive_precondition_msg
#undef CGAL_polygon_expensive_precondition_code

#if defined(CGAL_POLYGON_NO_PRECONDITIONS) || defined(CGAL_NO_PRECONDITIONS) \
  || (!defined(CGAL_POLYGON_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_polygon_expensive_precondition(EX) (static_cast<void>(0))
#  define CGAL_polygon_expensive_precondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_polygon_expensive_precondition_code(CODE)
#else
#  define CGAL_polygon_expensive_precondition(EX) \
   ((!get_use_polygon_assertions() || CGAL::possibly(EX))?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_polygon_expensive_precondition_msg(EX,MSG) \
   ((!get_use_polygon_assertions() || CGAL::possibly(EX))?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_polygon_expensive_precondition_code(CODE) CODE
#  define CGAL_polygon_expensive_preconditions 1
#endif // CGAL_POLYGON_NO_PRECONDITIONS


#undef CGAL_polygon_expensive_exactness_precondition
#undef CGAL_polygon_expensive_exactness_precondition_msg
#undef CGAL_polygon_expensive_exactness_precondition_code

#if defined(CGAL_POLYGON_NO_PRECONDITIONS) || defined(CGAL_NO_PRECONDITIONS) \
  || (!defined(CGAL_POLYGON_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || (!defined(CGAL_POLYGON_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_polygon_expensive_exactness_precondition(EX) (static_cast<void>(0))
#  define CGAL_polygon_expensive_exactness_precondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_polygon_expensive_exactness_precondition_code(CODE)
#else
#  define CGAL_polygon_expensive_exactness_precondition(EX) \
   ((!get_use_polygon_assertions() || CGAL::possibly(EX))?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_polygon_expensive_exactness_precondition_msg(EX,MSG) \
   ((!get_use_polygon_assertions() || CGAL::possibly(EX))?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_polygon_expensive_exactness_precondition_code(CODE) CODE
#  define CGAL_polygon_expensive_exactness_preconditions 1
#endif // CGAL_POLYGON_NO_PRECONDITIONS


// postconditions
// --------------

#undef CGAL_polygon_postcondition
#undef CGAL_polygon_postcondition_msg
#undef CGAL_polygon_postcondition_code

#if defined(CGAL_POLYGON_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
  || defined(NDEBUG)
#  define CGAL_polygon_postcondition(EX) (static_cast<void>(0))
#  define CGAL_polygon_postcondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_polygon_postcondition_code(CODE)
#else
#  define CGAL_polygon_postcondition(EX) \
   ((!get_use_polygon_assertions() || CGAL::possibly(EX))?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_polygon_postcondition_msg(EX,MSG) \
   ((!get_use_polygon_assertions() || CGAL::possibly(EX))?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_polygon_postcondition_code(CODE) CODE
#  define CGAL_polygon_postconditions 1
#endif // CGAL_POLYGON_NO_POSTCONDITIONS


#undef CGAL_polygon_exactness_postcondition
#undef CGAL_polygon_exactness_postcondition_msg
#undef CGAL_polygon_exactness_postcondition_code

#if defined(CGAL_POLYGON_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
  || (!defined(CGAL_POLYGON_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || defined(NDEBUG)
#  define CGAL_polygon_exactness_postcondition(EX) (static_cast<void>(0))
#  define CGAL_polygon_exactness_postcondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_polygon_exactness_postcondition_code(CODE)
#else
#  define CGAL_polygon_exactness_postcondition(EX) \
   ((!get_use_polygon_assertions() || CGAL::possibly(EX))?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_polygon_exactness_postcondition_msg(EX,MSG) \
   ((!get_use_polygon_assertions() || CGAL::possibly(EX))?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_polygon_exactness_postcondition_code(CODE) CODE
#  define CGAL_polygon_exactness_postconditions 1
#endif // CGAL_POLYGON_NO_POSTCONDITIONS


#undef CGAL_polygon_expensive_postcondition
#undef CGAL_polygon_expensive_postcondition_msg
#undef CGAL_polygon_expensive_postcondition_code

#if defined(CGAL_POLYGON_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
  || (!defined(CGAL_POLYGON_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_polygon_expensive_postcondition(EX) (static_cast<void>(0))
#  define CGAL_polygon_expensive_postcondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_polygon_expensive_postcondition_code(CODE)
#else
#  define CGAL_polygon_expensive_postcondition(EX) \
   ((!get_use_polygon_assertions() || CGAL::possibly(EX))?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_polygon_expensive_postcondition_msg(EX,MSG) \
   ((!get_use_polygon_assertions() || CGAL::possibly(EX))?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_polygon_expensive_postcondition_code(CODE) CODE
#  define CGAL_polygon_expensive_postconditions 1
#endif // CGAL_POLYGON_NO_POSTCONDITIONS


#undef CGAL_polygon_expensive_exactness_postcondition
#undef CGAL_polygon_expensive_exactness_postcondition_msg
#undef CGAL_polygon_expensive_exactness_postcondition_code

#if defined(CGAL_POLYGON_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
  || (!defined(CGAL_POLYGON_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || (!defined(CGAL_POLYGON_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_polygon_expensive_exactness_postcondition(EX) (static_cast<void>(0))
#  define CGAL_polygon_expensive_exactness_postcondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_polygon_expensive_exactness_postcondition_code(CODE)
#else
#  define CGAL_polygon_expensive_exactness_postcondition(EX) \
   ((!get_use_polygon_assertions() || CGAL::possibly(EX))?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_polygon_expensive_exactness_postcondition_msg(EX,MSG) \
   ((!get_use_polygon_assertions() || CGAL::possibly(EX))?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_polygon_expensive_exactness_postcondition_code(CODE) CODE
#  define CGAL_polygon_expensive_exactness_postconditions 1
#endif // CGAL_POLYGON_NO_POSTCONDITIONS


// warnings
// --------

#undef CGAL_polygon_warning
#undef CGAL_polygon_warning_msg
#undef CGAL_polygon_warning_code

#if defined(CGAL_POLYGON_NO_WARNINGS) || defined(CGAL_NO_WARNINGS) \
  || defined(NDEBUG)
#  define CGAL_polygon_warning(EX) (static_cast<void>(0))
#  define CGAL_polygon_warning_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_polygon_warning_code(CODE)
#else
#  define CGAL_polygon_warning(EX) \
   ((!get_use_polygon_assertions() || CGAL::possibly(EX))?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__))
#  define CGAL_polygon_warning_msg(EX,MSG) \
   ((!get_use_polygon_assertions() || CGAL::possibly(EX))?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_polygon_warning_code(CODE) CODE
#  define CGAL_polygon_warnings 1
#endif // CGAL_POLYGON_NO_WARNINGS


#undef CGAL_polygon_exactness_warning
#undef CGAL_polygon_exactness_warning_msg
#undef CGAL_polygon_exactness_warning_code

#if defined(CGAL_POLYGON_NO_WARNINGS) || defined(CGAL_NO_WARNINGS) \
  || (!defined(CGAL_POLYGON_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || defined(NDEBUG)
#  define CGAL_polygon_exactness_warning(EX) (static_cast<void>(0))
#  define CGAL_polygon_exactness_warning_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_polygon_exactness_warning_code(CODE)
#else
#  define CGAL_polygon_exactness_warning(EX) \
   ((!get_use_polygon_assertions() || CGAL::possibly(EX))?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__))
#  define CGAL_polygon_exactness_warning_msg(EX,MSG) \
   ((!get_use_polygon_assertions() || CGAL::possibly(EX))?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_polygon_exactness_warning_code(CODE) CODE
#  define CGAL_polygon_exactness_warnings 1
#endif // CGAL_POLYGON_NO_WARNINGS


#undef CGAL_polygon_expensive_warning
#undef CGAL_polygon_expensive_warning_msg
#undef CGAL_polygon_expensive_warning_code

#if defined(CGAL_POLYGON_NO_WARNINGS) || defined(CGAL_NO_WARNINGS) \
  || (!defined(CGAL_POLYGON_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_polygon_expensive_warning(EX) (static_cast<void>(0))
#  define CGAL_polygon_expensive_warning_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_polygon_expensive_warning_code(CODE)
#else
#  define CGAL_polygon_expensive_warning(EX) \
   ((!get_use_polygon_assertions() || CGAL::possibly(EX))?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__))
#  define CGAL_polygon_expensive_warning_msg(EX,MSG) \
   ((!get_use_polygon_assertions() || CGAL::possibly(EX))?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_polygon_expensive_warning_code(CODE) CODE
#  define CGAL_polygon_expensive_warnings 1
#endif // CGAL_POLYGON_NO_WARNINGS


#undef CGAL_polygon_expensive_exactness_warning
#undef CGAL_polygon_expensive_exactness_warning_msg
#undef CGAL_polygon_expensive_exactness_warning_code

#if defined(CGAL_POLYGON_NO_WARNINGS) || defined(CGAL_NO_WARNINGS) \
  || (!defined(CGAL_POLYGON_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || (!defined(CGAL_POLYGON_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_polygon_expensive_exactness_warning(EX) (static_cast<void>(0))
#  define CGAL_polygon_expensive_exactness_warning_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_polygon_expensive_exactness_warning_code(CODE)
#else
#  define CGAL_polygon_expensive_exactness_warning(EX) \
   ((!get_use_polygon_assertions() || CGAL::possibly(EX))?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__))
#  define CGAL_polygon_expensive_exactness_warning_msg(EX,MSG) \
   ((!get_use_polygon_assertions() || CGAL::possibly(EX))?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_polygon_expensive_exactness_warning_code(CODE) CODE
#  define CGAL_polygon_expensive_exactness_warnings 1
#endif // CGAL_POLYGON_NO_WARNINGS
