// Copyright (c) 2005  INRIA (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Geert-Jan Giezeman, Sven Schönherr

// Generated from script create_assertions.sh

// macro definitions
// =================
// assertions
// ----------


#if defined(CGAL_PARAMETERIZATION_NO_ASSERTIONS) || defined(CGAL_NO_ASSERTIONS) \
  || defined(NDEBUG)
#  define CGAL_parameterization_assertion(EX) (static_cast<void>(0))
#  define CGAL_parameterization_assertion_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_parameterization_assertion_code(CODE)
#else
#  define CGAL_parameterization_assertion(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_parameterization_assertion_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_parameterization_assertion_code(CODE) CODE
#  define CGAL_parameterization_assertions 1
#endif // CGAL_PARAMETERIZATION_NO_ASSERTIONS

#if defined(CGAL_PARAMETERIZATION_NO_ASSERTIONS) || defined(CGAL_NO_ASSERTIONS) \
  || (!defined(CGAL_PARAMETERIZATION_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || defined(NDEBUG)
#  define CGAL_parameterization_exactness_assertion(EX) (static_cast<void>(0))
#  define CGAL_parameterization_exactness_assertion_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_parameterization_exactness_assertion_code(CODE)
#else
#  define CGAL_parameterization_exactness_assertion(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_parameterization_exactness_assertion_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_parameterization_exactness_assertion_code(CODE) CODE
#  define CGAL_parameterization_exactness_assertions 1
#endif // CGAL_PARAMETERIZATION_NO_ASSERTIONS

#if defined(CGAL_PARAMETERIZATION_NO_ASSERTIONS) \
  || defined(CGAL_NO_ASSERTIONS) \
  || (!defined(CGAL_PARAMETERIZATION_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_parameterization_expensive_assertion(EX) (static_cast<void>(0))
#  define CGAL_parameterization_expensive_assertion_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_parameterization_expensive_assertion_code(CODE)
#else
#  define CGAL_parameterization_expensive_assertion(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_parameterization_expensive_assertion_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_parameterization_expensive_assertion_code(CODE) CODE
#  define CGAL_parameterization_expensive_assertions 1
#endif // CGAL_PARAMETERIZATION_NO_ASSERTIONS

#if defined(CGAL_PARAMETERIZATION_NO_ASSERTIONS) || defined(CGAL_NO_ASSERTIONS) \
  || (!defined(CGAL_PARAMETERIZATION_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || (!defined(CGAL_PARAMETERIZATION_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_parameterization_expensive_exactness_assertion(EX) (static_cast<void>(0))
#  define CGAL_parameterization_expensive_exactness_assertion_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_parameterization_expensive_exactness_assertion_code(CODE)
#else
#  define CGAL_parameterization_expensive_exactness_assertion(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_parameterization_expensive_exactness_assertion_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_parameterization_expensive_exactness_assertion_code(CODE) CODE
#  define CGAL_parameterization_expensive_exactness_assertions 1
#endif // CGAL_PARAMETERIZATION_NO_ASSERTIONS


// preconditions
// -------------

#if defined(CGAL_PARAMETERIZATION_NO_PRECONDITIONS) || defined(CGAL_NO_PRECONDITIONS) \
  || defined(NDEBUG)
#  define CGAL_parameterization_precondition(EX) (static_cast<void>(0))
#  define CGAL_parameterization_precondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_parameterization_precondition_code(CODE)
#else
#  define CGAL_parameterization_precondition(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_parameterization_precondition_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_parameterization_precondition_code(CODE) CODE
#  define CGAL_parameterization_preconditions 1
#endif // CGAL_PARAMETERIZATION_NO_PRECONDITIONS

#if defined(CGAL_PARAMETERIZATION_NO_PRECONDITIONS) || defined(CGAL_NO_PRECONDITIONS) \
  || (!defined(CGAL_PARAMETERIZATION_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || defined(NDEBUG)
#  define CGAL_parameterization_exactness_precondition(EX) (static_cast<void>(0))
#  define CGAL_parameterization_exactness_precondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_parameterization_exactness_precondition_code(CODE)
#else
#  define CGAL_parameterization_exactness_precondition(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_parameterization_exactness_precondition_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_parameterization_exactness_precondition_code(CODE) CODE
#  define CGAL_parameterization_exactness_preconditions 1
#endif // CGAL_PARAMETERIZATION_NO_PRECONDITIONS

#if defined(CGAL_PARAMETERIZATION_NO_PRECONDITIONS) || defined(CGAL_NO_PRECONDITIONS) \
  || (!defined(CGAL_PARAMETERIZATION_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_parameterization_expensive_precondition(EX) (static_cast<void>(0))
#  define CGAL_parameterization_expensive_precondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_parameterization_expensive_precondition_code(CODE)
#else
#  define CGAL_parameterization_expensive_precondition(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_parameterization_expensive_precondition_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_parameterization_expensive_precondition_code(CODE) CODE
#  define CGAL_parameterization_expensive_preconditions 1
#endif // CGAL_PARAMETERIZATION_NO_PRECONDITIONS

#if defined(CGAL_PARAMETERIZATION_NO_PRECONDITIONS) || defined(CGAL_NO_PRECONDITIONS) \
  || (!defined(CGAL_PARAMETERIZATION_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || (!defined(CGAL_PARAMETERIZATION_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_parameterization_expensive_exactness_precondition(EX) (static_cast<void>(0))
#  define CGAL_parameterization_expensive_exactness_precondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_parameterization_expensive_exactness_precondition_code(CODE)
#else
#  define CGAL_parameterization_expensive_exactness_precondition(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_parameterization_expensive_exactness_precondition_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_parameterization_expensive_exactness_precondition_code(CODE) CODE
#  define CGAL_parameterization_expensive_exactness_preconditions 1
#endif // CGAL_PARAMETERIZATION_NO_PRECONDITIONS


// postconditions
// --------------

#if defined(CGAL_PARAMETERIZATION_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
  || defined(NDEBUG)
#  define CGAL_parameterization_postcondition(EX) (static_cast<void>(0))
#  define CGAL_parameterization_postcondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_parameterization_postcondition_code(CODE)
#else
#  define CGAL_parameterization_postcondition(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_parameterization_postcondition_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_parameterization_postcondition_code(CODE) CODE
#  define CGAL_parameterization_postconditions 1
#endif // CGAL_PARAMETERIZATION_NO_POSTCONDITIONS

#if defined(CGAL_PARAMETERIZATION_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
  || (!defined(CGAL_PARAMETERIZATION_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || defined(NDEBUG)
#  define CGAL_parameterization_exactness_postcondition(EX) (static_cast<void>(0))
#  define CGAL_parameterization_exactness_postcondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_parameterization_exactness_postcondition_code(CODE)
#else
#  define CGAL_parameterization_exactness_postcondition(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_parameterization_exactness_postcondition_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_parameterization_exactness_postcondition_code(CODE) CODE
#  define CGAL_parameterization_exactness_postconditions 1
#endif // CGAL_PARAMETERIZATION_NO_POSTCONDITIONS

#if defined(CGAL_PARAMETERIZATION_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
  || (!defined(CGAL_PARAMETERIZATION_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_parameterization_expensive_postcondition(EX) (static_cast<void>(0))
#  define CGAL_parameterization_expensive_postcondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_parameterization_expensive_postcondition_code(CODE)
#else
#  define CGAL_parameterization_expensive_postcondition(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_parameterization_expensive_postcondition_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_parameterization_expensive_postcondition_code(CODE) CODE
#  define CGAL_parameterization_expensive_postconditions 1
#endif // CGAL_PARAMETERIZATION_NO_POSTCONDITIONS

#if defined(CGAL_PARAMETERIZATION_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
  || (!defined(CGAL_PARAMETERIZATION_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || (!defined(CGAL_PARAMETERIZATION_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_parameterization_expensive_exactness_postcondition(EX) (static_cast<void>(0))
#  define CGAL_parameterization_expensive_exactness_postcondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_parameterization_expensive_exactness_postcondition_code(CODE)
#else
#  define CGAL_parameterization_expensive_exactness_postcondition(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_parameterization_expensive_exactness_postcondition_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_parameterization_expensive_exactness_postcondition_code(CODE) CODE
#  define CGAL_parameterization_expensive_exactness_postconditions 1
#endif // CGAL_PARAMETERIZATION_NO_POSTCONDITIONS


// warnings
// --------

#if defined(CGAL_PARAMETERIZATION_NO_WARNINGS) || defined(CGAL_NO_WARNINGS) \
  || defined(NDEBUG)
#  define CGAL_parameterization_warning(EX) (static_cast<void>(0))
#  define CGAL_parameterization_warning_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_parameterization_warning_code(CODE)
#else
#  define CGAL_parameterization_warning(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_parameterization_warning_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_parameterization_warning_code(CODE) CODE
#  define CGAL_parameterization_warnings 1
#endif // CGAL_PARAMETERIZATION_NO_WARNINGS

#if defined(CGAL_PARAMETERIZATION_NO_WARNINGS) || defined(CGAL_NO_WARNINGS) \
  || (!defined(CGAL_PARAMETERIZATION_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || defined(NDEBUG)
#  define CGAL_parameterization_exactness_warning(EX) (static_cast<void>(0))
#  define CGAL_parameterization_exactness_warning_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_parameterization_exactness_warning_code(CODE)
#else
#  define CGAL_parameterization_exactness_warning(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_parameterization_exactness_warning_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_parameterization_exactness_warning_code(CODE) CODE
#  define CGAL_parameterization_exactness_warnings 1
#endif // CGAL_PARAMETERIZATION_NO_WARNINGS

#if defined(CGAL_PARAMETERIZATION_NO_WARNINGS) || defined(CGAL_NO_WARNINGS) \
  || (!defined(CGAL_PARAMETERIZATION_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_parameterization_expensive_warning(EX) (static_cast<void>(0))
#  define CGAL_parameterization_expensive_warning_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_parameterization_expensive_warning_code(CODE)
#else
#  define CGAL_parameterization_expensive_warning(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_parameterization_expensive_warning_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_parameterization_expensive_warning_code(CODE) CODE
#  define CGAL_parameterization_expensive_warnings 1
#endif // CGAL_PARAMETERIZATION_NO_WARNINGS

#if defined(CGAL_PARAMETERIZATION_NO_WARNINGS) || defined(CGAL_NO_WARNINGS) \
  || (!defined(CGAL_PARAMETERIZATION_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || (!defined(CGAL_PARAMETERIZATION_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_parameterization_expensive_exactness_warning(EX) (static_cast<void>(0))
#  define CGAL_parameterization_expensive_exactness_warning_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_parameterization_expensive_exactness_warning_code(CODE)
#else
#  define CGAL_parameterization_expensive_exactness_warning(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_parameterization_expensive_exactness_warning_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_parameterization_expensive_exactness_warning_code(CODE) CODE
#  define CGAL_parameterization_expensive_exactness_warnings 1
#endif // CGAL_PARAMETERIZATION_NO_WARNINGS


