// Copyright (c) 2005  Tel-Aviv University (Israel).  All rights reserved.
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
// Author(s)     : Geert-Jan Giezeman, Sven Schönherr
//
// Generated from script create_assertions.sh

// macro definitions
// =================
// assertions
// ----------


#if defined(CGAL_MULTISET_NO_ASSERTIONS) || defined(CGAL_NO_ASSERTIONS) \
  || defined(NDEBUG)
#  define CGAL_multiset_assertion(EX) (static_cast<void>(0))
#  define CGAL_multiset_assertion_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_multiset_assertion_code(CODE)
#else
#  define CGAL_multiset_assertion(EX) \
   (CGAL::certainly(EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__))
#  define CGAL_multiset_assertion_msg(EX,MSG) \
   (CGAL::certainly(EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_multiset_assertion_code(CODE) CODE
#  define CGAL_multiset_assertions 1
#endif // CGAL_MULTISET_NO_ASSERTIONS

#if defined(CGAL_MULTISET_NO_ASSERTIONS) || defined(CGAL_NO_ASSERTIONS) \
  || (!defined(CGAL_MULTISET_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || defined(NDEBUG)
#  define CGAL_multiset_exactness_assertion(EX) (static_cast<void>(0))
#  define CGAL_multiset_exactness_assertion_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_multiset_exactness_assertion_code(CODE)
#else
#  define CGAL_multiset_exactness_assertion(EX) \
   (CGAL::certainly(EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__))
#  define CGAL_multiset_exactness_assertion_msg(EX,MSG) \
   (CGAL::certainly(EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_multiset_exactness_assertion_code(CODE) CODE
#  define CGAL_multiset_exactness_assertions 1
#endif // CGAL_MULTISET_NO_ASSERTIONS

#if defined(CGAL_MULTISET_NO_ASSERTIONS) \
  || defined(CGAL_NO_ASSERTIONS) \
  || (!defined(CGAL_MULTISET_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_multiset_expensive_assertion(EX) (static_cast<void>(0))
#  define CGAL_multiset_expensive_assertion_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_multiset_expensive_assertion_code(CODE)
#else
#  define CGAL_multiset_expensive_assertion(EX) \
   (CGAL::certainly(EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__))
#  define CGAL_multiset_expensive_assertion_msg(EX,MSG) \
   (CGAL::certainly(EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_multiset_expensive_assertion_code(CODE) CODE
#  define CGAL_multiset_expensive_assertions 1
#endif // CGAL_MULTISET_NO_ASSERTIONS

#if defined(CGAL_MULTISET_NO_ASSERTIONS) || defined(CGAL_NO_ASSERTIONS) \
  || (!defined(CGAL_MULTISET_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || (!defined(CGAL_MULTISET_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_multiset_expensive_exactness_assertion(EX) (static_cast<void>(0))
#  define CGAL_multiset_expensive_exactness_assertion_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_multiset_expensive_exactness_assertion_code(CODE)
#else
#  define CGAL_multiset_expensive_exactness_assertion(EX) \
   (CGAL::certainly(EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__))
#  define CGAL_multiset_expensive_exactness_assertion_msg(EX,MSG) \
   (CGAL::certainly(EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_multiset_expensive_exactness_assertion_code(CODE) CODE
#  define CGAL_multiset_expensive_exactness_assertions 1
#endif // CGAL_MULTISET_NO_ASSERTIONS


// preconditions
// -------------

#if defined(CGAL_MULTISET_NO_PRECONDITIONS) || defined(CGAL_NO_PRECONDITIONS) \
  || defined(NDEBUG)
#  define CGAL_multiset_precondition(EX) (static_cast<void>(0))
#  define CGAL_multiset_precondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_multiset_precondition_code(CODE)
#else
#  define CGAL_multiset_precondition(EX) \
   (CGAL::certainly(EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_multiset_precondition_msg(EX,MSG) \
   (CGAL::certainly(EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_multiset_precondition_code(CODE) CODE
#  define CGAL_multiset_preconditions 1
#endif // CGAL_MULTISET_NO_PRECONDITIONS

#if defined(CGAL_MULTISET_NO_PRECONDITIONS) || defined(CGAL_NO_PRECONDITIONS) \
  || (!defined(CGAL_MULTISET_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || defined(NDEBUG)
#  define CGAL_multiset_exactness_precondition(EX) (static_cast<void>(0))
#  define CGAL_multiset_exactness_precondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_multiset_exactness_precondition_code(CODE)
#else
#  define CGAL_multiset_exactness_precondition(EX) \
   (CGAL::certainly(EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_multiset_exactness_precondition_msg(EX,MSG) \
   (CGAL::certainly(EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_multiset_exactness_precondition_code(CODE) CODE
#  define CGAL_multiset_exactness_preconditions 1
#endif // CGAL_MULTISET_NO_PRECONDITIONS

#if defined(CGAL_MULTISET_NO_PRECONDITIONS) || defined(CGAL_NO_PRECONDITIONS) \
  || (!defined(CGAL_MULTISET_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_multiset_expensive_precondition(EX) (static_cast<void>(0))
#  define CGAL_multiset_expensive_precondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_multiset_expensive_precondition_code(CODE)
#else
#  define CGAL_multiset_expensive_precondition(EX) \
   (CGAL::certainly(EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_multiset_expensive_precondition_msg(EX,MSG) \
   (CGAL::certainly(EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_multiset_expensive_precondition_code(CODE) CODE
#  define CGAL_multiset_expensive_preconditions 1
#endif // CGAL_MULTISET_NO_PRECONDITIONS

#if defined(CGAL_MULTISET_NO_PRECONDITIONS) || defined(CGAL_NO_PRECONDITIONS) \
  || (!defined(CGAL_MULTISET_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || (!defined(CGAL_MULTISET_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_multiset_expensive_exactness_precondition(EX) (static_cast<void>(0))
#  define CGAL_multiset_expensive_exactness_precondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_multiset_expensive_exactness_precondition_code(CODE)
#else
#  define CGAL_multiset_expensive_exactness_precondition(EX) \
   (CGAL::certainly(EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_multiset_expensive_exactness_precondition_msg(EX,MSG) \
   (CGAL::certainly(EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_multiset_expensive_exactness_precondition_code(CODE) CODE
#  define CGAL_multiset_expensive_exactness_preconditions 1
#endif // CGAL_MULTISET_NO_PRECONDITIONS


// postconditions
// --------------

#if defined(CGAL_MULTISET_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
  || defined(NDEBUG)
#  define CGAL_multiset_postcondition(EX) (static_cast<void>(0))
#  define CGAL_multiset_postcondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_multiset_postcondition_code(CODE)
#else
#  define CGAL_multiset_postcondition(EX) \
   (CGAL::certainly(EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_multiset_postcondition_msg(EX,MSG) \
   (CGAL::certainly(EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_multiset_postcondition_code(CODE) CODE
#  define CGAL_multiset_postconditions 1
#endif // CGAL_MULTISET_NO_POSTCONDITIONS

#if defined(CGAL_MULTISET_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
  || (!defined(CGAL_MULTISET_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || defined(NDEBUG)
#  define CGAL_multiset_exactness_postcondition(EX) (static_cast<void>(0))
#  define CGAL_multiset_exactness_postcondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_multiset_exactness_postcondition_code(CODE)
#else
#  define CGAL_multiset_exactness_postcondition(EX) \
   (CGAL::certainly(EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_multiset_exactness_postcondition_msg(EX,MSG) \
   (CGAL::certainly(EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_multiset_exactness_postcondition_code(CODE) CODE
#  define CGAL_multiset_exactness_postconditions 1
#endif // CGAL_MULTISET_NO_POSTCONDITIONS

#if defined(CGAL_MULTISET_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
  || (!defined(CGAL_MULTISET_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_multiset_expensive_postcondition(EX) (static_cast<void>(0))
#  define CGAL_multiset_expensive_postcondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_multiset_expensive_postcondition_code(CODE)
#else
#  define CGAL_multiset_expensive_postcondition(EX) \
   (CGAL::certainly(EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_multiset_expensive_postcondition_msg(EX,MSG) \
   (CGAL::certainly(EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_multiset_expensive_postcondition_code(CODE) CODE
#  define CGAL_multiset_expensive_postconditions 1
#endif // CGAL_MULTISET_NO_POSTCONDITIONS

#if defined(CGAL_MULTISET_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
  || (!defined(CGAL_MULTISET_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || (!defined(CGAL_MULTISET_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_multiset_expensive_exactness_postcondition(EX) (static_cast<void>(0))
#  define CGAL_multiset_expensive_exactness_postcondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_multiset_expensive_exactness_postcondition_code(CODE)
#else
#  define CGAL_multiset_expensive_exactness_postcondition(EX) \
   (CGAL::certainly(EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_multiset_expensive_exactness_postcondition_msg(EX,MSG) \
   (CGAL::certainly(EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_multiset_expensive_exactness_postcondition_code(CODE) CODE
#  define CGAL_multiset_expensive_exactness_postconditions 1
#endif // CGAL_MULTISET_NO_POSTCONDITIONS


// warnings
// --------

#if defined(CGAL_MULTISET_NO_WARNINGS) || defined(CGAL_NO_WARNINGS) \
  || defined(NDEBUG)
#  define CGAL_multiset_warning(EX) (static_cast<void>(0))
#  define CGAL_multiset_warning_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_multiset_warning_code(CODE)
#else
#  define CGAL_multiset_warning(EX) \
   (CGAL::certainly(EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__))
#  define CGAL_multiset_warning_msg(EX,MSG) \
   (CGAL::certainly(EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_multiset_warning_code(CODE) CODE
#  define CGAL_multiset_warnings 1
#endif // CGAL_MULTISET_NO_WARNINGS

#if defined(CGAL_MULTISET_NO_WARNINGS) || defined(CGAL_NO_WARNINGS) \
  || (!defined(CGAL_MULTISET_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || defined(NDEBUG)
#  define CGAL_multiset_exactness_warning(EX) (static_cast<void>(0))
#  define CGAL_multiset_exactness_warning_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_multiset_exactness_warning_code(CODE)
#else
#  define CGAL_multiset_exactness_warning(EX) \
   (CGAL::certainly(EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__))
#  define CGAL_multiset_exactness_warning_msg(EX,MSG) \
   (CGAL::certainly(EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_multiset_exactness_warning_code(CODE) CODE
#  define CGAL_multiset_exactness_warnings 1
#endif // CGAL_MULTISET_NO_WARNINGS

#if defined(CGAL_MULTISET_NO_WARNINGS) || defined(CGAL_NO_WARNINGS) \
  || (!defined(CGAL_MULTISET_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_multiset_expensive_warning(EX) (static_cast<void>(0))
#  define CGAL_multiset_expensive_warning_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_multiset_expensive_warning_code(CODE)
#else
#  define CGAL_multiset_expensive_warning(EX) \
   (CGAL::certainly(EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__))
#  define CGAL_multiset_expensive_warning_msg(EX,MSG) \
   (CGAL::certainly(EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_multiset_expensive_warning_code(CODE) CODE
#  define CGAL_multiset_expensive_warnings 1
#endif // CGAL_MULTISET_NO_WARNINGS

#if defined(CGAL_MULTISET_NO_WARNINGS) || defined(CGAL_NO_WARNINGS) \
  || (!defined(CGAL_MULTISET_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || (!defined(CGAL_MULTISET_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_multiset_expensive_exactness_warning(EX) (static_cast<void>(0))
#  define CGAL_multiset_expensive_exactness_warning_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_multiset_expensive_exactness_warning_code(CODE)
#else
#  define CGAL_multiset_expensive_exactness_warning(EX) \
   (CGAL::certainly(EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__))
#  define CGAL_multiset_expensive_exactness_warning_msg(EX,MSG) \
   (CGAL::certainly(EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_multiset_expensive_exactness_warning_code(CODE) CODE
#  define CGAL_multiset_expensive_exactness_warnings 1
#endif // CGAL_MULTISET_NO_WARNINGS


