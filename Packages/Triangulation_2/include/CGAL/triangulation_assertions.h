// Copyright (c) 2001  INRIA Sophia-Antipolis (France).
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
// Author(s)     : script by Geert-Jan Giezeman and Sven Schoenherr 



// macro definitions
// =================
// assertions
// ----------


#if defined(CGAL_TRIANGULATION_NO_ASSERTIONS) || defined(CGAL_NO_ASSERTIONS) \
  || defined(NDEBUG)
#  define CGAL_triangulation_assertion(EX) (static_cast<void>(0))
#  define CGAL_triangulation_assertion_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_triangulation_assertion_code(CODE)
#else
#  define CGAL_triangulation_assertion(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_triangulation_assertion_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_triangulation_assertion_code(CODE) CODE
#endif // CGAL_TRIANGULATION_NO_ASSERTIONS

#if defined(CGAL_TRIANGULATION_NO_ASSERTIONS) || defined(CGAL_NO_ASSERTIONS) \
  || (!defined(CGAL_TRIANGULATION_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || defined(NDEBUG)
#  define CGAL_triangulation_exactness_assertion(EX) (static_cast<void>(0))
#  define CGAL_triangulation_exactness_assertion_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_triangulation_exactness_assertion_code(CODE)
#else
#  define CGAL_triangulation_exactness_assertion(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_triangulation_exactness_assertion_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_triangulation_exactness_assertion_code(CODE) CODE
#endif // CGAL_TRIANGULATION_NO_ASSERTIONS

#if defined(CGAL_TRIANGULATION_NO_ASSERTIONS) \
  || defined(CGAL_NO_ASSERTIONS) \
  || (!defined(CGAL_TRIANGULATION_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_triangulation_expensive_assertion(EX) (static_cast<void>(0))
#  define CGAL_triangulation_expensive_assertion_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_triangulation_expensive_assertion_code(CODE)
#else
#  define CGAL_triangulation_expensive_assertion(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_triangulation_expensive_assertion_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_triangulation_expensive_assertion_code(CODE) CODE
#endif // CGAL_TRIANGULATION_NO_ASSERTIONS

#if defined(CGAL_TRIANGULATION_NO_ASSERTIONS) || defined(CGAL_NO_ASSERTIONS) \
  || (!defined(CGAL_TRIANGULATION_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || (!defined(CGAL_TRIANGULATION_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_triangulation_expensive_exactness_assertion(EX) (static_cast<void>(0))
#  define CGAL_triangulation_expensive_exactness_assertion_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_triangulation_expensive_exactness_assertion_code(CODE)
#else
#  define CGAL_triangulation_expensive_exactness_assertion(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_triangulation_expensive_exactness_assertion_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_triangulation_expensive_exactness_assertion_code(CODE) CODE
#endif // CGAL_TRIANGULATION_NO_ASSERTIONS


// preconditions
// -------------

#if defined(CGAL_TRIANGULATION_NO_PRECONDITIONS) || defined(CGAL_NO_PRECONDITIONS) \
  || defined(NDEBUG)
#  define CGAL_triangulation_precondition(EX) (static_cast<void>(0))
#  define CGAL_triangulation_precondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_triangulation_precondition_code(CODE)
#else
#  define CGAL_triangulation_precondition(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_triangulation_precondition_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_triangulation_precondition_code(CODE) CODE
#endif // CGAL_TRIANGULATION_NO_PRECONDITIONS

#if defined(CGAL_TRIANGULATION_NO_PRECONDITIONS) || defined(CGAL_NO_PRECONDITIONS) \
  || (!defined(CGAL_TRIANGULATION_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || defined(NDEBUG)
#  define CGAL_triangulation_exactness_precondition(EX) (static_cast<void>(0))
#  define CGAL_triangulation_exactness_precondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_triangulation_exactness_precondition_code(CODE)
#else
#  define CGAL_triangulation_exactness_precondition(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_triangulation_exactness_precondition_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_triangulation_exactness_precondition_code(CODE) CODE
#endif // CGAL_TRIANGULATION_NO_PRECONDITIONS

#if defined(CGAL_TRIANGULATION_NO_PRECONDITIONS) || defined(CGAL_NO_PRECONDITIONS) \
  || (!defined(CGAL_TRIANGULATION_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_triangulation_expensive_precondition(EX) (static_cast<void>(0))
#  define CGAL_triangulation_expensive_precondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_triangulation_expensive_precondition_code(CODE)
#else
#  define CGAL_triangulation_expensive_precondition(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_triangulation_expensive_precondition_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_triangulation_expensive_precondition_code(CODE) CODE
#endif // CGAL_TRIANGULATION_NO_PRECONDITIONS

#if defined(CGAL_TRIANGULATION_NO_PRECONDITIONS) || defined(CGAL_NO_PRECONDITIONS) \
  || (!defined(CGAL_TRIANGULATION_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || (!defined(CGAL_TRIANGULATION_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_triangulation_expensive_exactness_precondition(EX) (static_cast<void>(0))
#  define CGAL_triangulation_expensive_exactness_precondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_triangulation_expensive_exactness_precondition_code(CODE)
#else
#  define CGAL_triangulation_expensive_exactness_precondition(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_triangulation_expensive_exactness_precondition_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_triangulation_expensive_exactness_precondition_code(CODE) CODE
#endif // CGAL_TRIANGULATION_NO_PRECONDITIONS


// postconditions
// --------------

#if defined(CGAL_TRIANGULATION_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
  || defined(NDEBUG)
#  define CGAL_triangulation_postcondition(EX) (static_cast<void>(0))
#  define CGAL_triangulation_postcondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_triangulation_postcondition_code(CODE)
#else
#  define CGAL_triangulation_postcondition(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_triangulation_postcondition_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_triangulation_postcondition_code(CODE) CODE
#endif // CGAL_TRIANGULATION_NO_POSTCONDITIONS

#if defined(CGAL_TRIANGULATION_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
  || (!defined(CGAL_TRIANGULATION_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || defined(NDEBUG)
#  define CGAL_triangulation_exactness_postcondition(EX) (static_cast<void>(0))
#  define CGAL_triangulation_exactness_postcondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_triangulation_exactness_postcondition_code(CODE)
#else
#  define CGAL_triangulation_exactness_postcondition(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_triangulation_exactness_postcondition_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_triangulation_exactness_postcondition_code(CODE) CODE
#endif // CGAL_TRIANGULATION_NO_POSTCONDITIONS

#if defined(CGAL_TRIANGULATION_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
  || (!defined(CGAL_TRIANGULATION_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_triangulation_expensive_postcondition(EX) (static_cast<void>(0))
#  define CGAL_triangulation_expensive_postcondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_triangulation_expensive_postcondition_code(CODE)
#else
#  define CGAL_triangulation_expensive_postcondition(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_triangulation_expensive_postcondition_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_triangulation_expensive_postcondition_code(CODE) CODE
#endif // CGAL_TRIANGULATION_NO_POSTCONDITIONS

#if defined(CGAL_TRIANGULATION_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
  || (!defined(CGAL_TRIANGULATION_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || (!defined(CGAL_TRIANGULATION_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_triangulation_expensive_exactness_postcondition(EX) (static_cast<void>(0))
#  define CGAL_triangulation_expensive_exactness_postcondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_triangulation_expensive_exactness_postcondition_code(CODE)
#else
#  define CGAL_triangulation_expensive_exactness_postcondition(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_triangulation_expensive_exactness_postcondition_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_triangulation_expensive_exactness_postcondition_code(CODE) CODE
#endif // CGAL_TRIANGULATION_NO_POSTCONDITIONS


// warnings
// --------

#if defined(CGAL_TRIANGULATION_NO_WARNINGS) || defined(CGAL_NO_WARNINGS) \
  || defined(NDEBUG)
#  define CGAL_triangulation_warning(EX) (static_cast<void>(0))
#  define CGAL_triangulation_warning_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_triangulation_warning_code(CODE)
#else
#  define CGAL_triangulation_warning(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_triangulation_warning_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_triangulation_warning_code(CODE) CODE
#endif // CGAL_TRIANGULATION_NO_WARNINGS

#if defined(CGAL_TRIANGULATION_NO_WARNINGS) || defined(CGAL_NO_WARNINGS) \
  || (!defined(CGAL_TRIANGULATION_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || defined(NDEBUG)
#  define CGAL_triangulation_exactness_warning(EX) (static_cast<void>(0))
#  define CGAL_triangulation_exactness_warning_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_triangulation_exactness_warning_code(CODE)
#else
#  define CGAL_triangulation_exactness_warning(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_triangulation_exactness_warning_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_triangulation_exactness_warning_code(CODE) CODE
#endif // CGAL_TRIANGULATION_NO_WARNINGS

#if defined(CGAL_TRIANGULATION_NO_WARNINGS) || defined(CGAL_NO_WARNINGS) \
  || (!defined(CGAL_TRIANGULATION_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_triangulation_expensive_warning(EX) (static_cast<void>(0))
#  define CGAL_triangulation_expensive_warning_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_triangulation_expensive_warning_code(CODE)
#else
#  define CGAL_triangulation_expensive_warning(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_triangulation_expensive_warning_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_triangulation_expensive_warning_code(CODE) CODE
#endif // CGAL_TRIANGULATION_NO_WARNINGS

#if defined(CGAL_TRIANGULATION_NO_WARNINGS) || defined(CGAL_NO_WARNINGS) \
  || (!defined(CGAL_TRIANGULATION_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || (!defined(CGAL_TRIANGULATION_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_triangulation_expensive_exactness_warning(EX) (static_cast<void>(0))
#  define CGAL_triangulation_expensive_exactness_warning_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_triangulation_expensive_exactness_warning_code(CODE)
#else
#  define CGAL_triangulation_expensive_exactness_warning(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_triangulation_expensive_exactness_warning_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_triangulation_expensive_exactness_warning_code(CODE) CODE
#endif // CGAL_TRIANGULATION_NO_WARNINGS


