

// ============================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision:  $
// release_date  : $CGAL_Date:  $
//
// file          : include/CGAL/polygon_assertions.h
// source        : Generated from script create_assertions.sh
// author(s)     : script by Geert-Jan Giezeman and Sven Schönherr 
//
// coordinator   : MPI, Saarbruecken
//
// ============================================================================



// macro definitions
// =================
// assertions
// ----------


#if defined(CGAL_POLYGON_NO_ASSERTIONS) || defined(CGAL_NO_ASSERTIONS) \
  || defined(NDEBUG)
#  define CGAL_polygon_assertion(EX) (static_cast<void>(0))
#  define CGAL_polygon_assertion_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_polygon_assertion_code(CODE)
#else
#  define CGAL_polygon_assertion(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_polygon_assertion_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_polygon_assertion_code(CODE) CODE
#endif // CGAL_POLYGON_NO_ASSERTIONS

#if defined(CGAL_POLYGON_NO_ASSERTIONS) || defined(CGAL_NO_ASSERTIONS) \
  || (!defined(CGAL_POLYGON_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || defined(NDEBUG)
#  define CGAL_polygon_exactness_assertion(EX) (static_cast<void>(0))
#  define CGAL_polygon_exactness_assertion_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_polygon_exactness_assertion_code(CODE)
#else
#  define CGAL_polygon_exactness_assertion(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_polygon_exactness_assertion_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_polygon_exactness_assertion_code(CODE) CODE
#endif // CGAL_POLYGON_NO_ASSERTIONS

#if defined(CGAL_POLYGON_NO_ASSERTIONS) \
  || defined(CGAL_NO_ASSERTIONS) \
  || (!defined(CGAL_POLYGON_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_polygon_expensive_assertion(EX) (static_cast<void>(0))
#  define CGAL_polygon_expensive_assertion_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_polygon_expensive_assertion_code(CODE)
#else
#  define CGAL_polygon_expensive_assertion(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_polygon_expensive_assertion_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_polygon_expensive_assertion_code(CODE) CODE
#endif // CGAL_POLYGON_NO_ASSERTIONS

#if defined(CGAL_POLYGON_NO_ASSERTIONS) || defined(CGAL_NO_ASSERTIONS) \
  || (!defined(CGAL_POLYGON_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || (!defined(CGAL_POLYGON_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_polygon_expensive_exactness_assertion(EX) (static_cast<void>(0))
#  define CGAL_polygon_expensive_exactness_assertion_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_polygon_expensive_exactness_assertion_code(CODE)
#else
#  define CGAL_polygon_expensive_exactness_assertion(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_polygon_expensive_exactness_assertion_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_polygon_expensive_exactness_assertion_code(CODE) CODE
#endif // CGAL_POLYGON_NO_ASSERTIONS


// preconditions
// -------------

#if defined(CGAL_POLYGON_NO_PRECONDITIONS) || defined(CGAL_NO_PRECONDITIONS) \
  || defined(NDEBUG)
#  define CGAL_polygon_precondition(EX) (static_cast<void>(0))
#  define CGAL_polygon_precondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_polygon_precondition_code(CODE)
#else
#  define CGAL_polygon_precondition(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_polygon_precondition_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_polygon_precondition_code(CODE) CODE
#endif // CGAL_POLYGON_NO_PRECONDITIONS

#if defined(CGAL_POLYGON_NO_PRECONDITIONS) || defined(CGAL_NO_PRECONDITIONS) \
  || (!defined(CGAL_POLYGON_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || defined(NDEBUG)
#  define CGAL_polygon_exactness_precondition(EX) (static_cast<void>(0))
#  define CGAL_polygon_exactness_precondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_polygon_exactness_precondition_code(CODE)
#else
#  define CGAL_polygon_exactness_precondition(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_polygon_exactness_precondition_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_polygon_exactness_precondition_code(CODE) CODE
#endif // CGAL_POLYGON_NO_PRECONDITIONS

#if defined(CGAL_POLYGON_NO_PRECONDITIONS) || defined(CGAL_NO_PRECONDITIONS) \
  || (!defined(CGAL_POLYGON_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_polygon_expensive_precondition(EX) (static_cast<void>(0))
#  define CGAL_polygon_expensive_precondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_polygon_expensive_precondition_code(CODE)
#else
#  define CGAL_polygon_expensive_precondition(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_polygon_expensive_precondition_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_polygon_expensive_precondition_code(CODE) CODE
#endif // CGAL_POLYGON_NO_PRECONDITIONS

#if defined(CGAL_POLYGON_NO_PRECONDITIONS) || defined(CGAL_NO_PRECONDITIONS) \
  || (!defined(CGAL_POLYGON_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || (!defined(CGAL_POLYGON_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_polygon_expensive_exactness_precondition(EX) (static_cast<void>(0))
#  define CGAL_polygon_expensive_exactness_precondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_polygon_expensive_exactness_precondition_code(CODE)
#else
#  define CGAL_polygon_expensive_exactness_precondition(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_polygon_expensive_exactness_precondition_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_polygon_expensive_exactness_precondition_code(CODE) CODE
#endif // CGAL_POLYGON_NO_PRECONDITIONS


// postconditions
// --------------

#if defined(CGAL_POLYGON_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
  || defined(NDEBUG)
#  define CGAL_polygon_postcondition(EX) (static_cast<void>(0))
#  define CGAL_polygon_postcondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_polygon_postcondition_code(CODE)
#else
#  define CGAL_polygon_postcondition(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_polygon_postcondition_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_polygon_postcondition_code(CODE) CODE
#endif // CGAL_POLYGON_NO_POSTCONDITIONS

#if defined(CGAL_POLYGON_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
  || (!defined(CGAL_POLYGON_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || defined(NDEBUG)
#  define CGAL_polygon_exactness_postcondition(EX) (static_cast<void>(0))
#  define CGAL_polygon_exactness_postcondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_polygon_exactness_postcondition_code(CODE)
#else
#  define CGAL_polygon_exactness_postcondition(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_polygon_exactness_postcondition_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_polygon_exactness_postcondition_code(CODE) CODE
#endif // CGAL_POLYGON_NO_POSTCONDITIONS

#if defined(CGAL_POLYGON_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
  || (!defined(CGAL_POLYGON_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_polygon_expensive_postcondition(EX) (static_cast<void>(0))
#  define CGAL_polygon_expensive_postcondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_polygon_expensive_postcondition_code(CODE)
#else
#  define CGAL_polygon_expensive_postcondition(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_polygon_expensive_postcondition_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_polygon_expensive_postcondition_code(CODE) CODE
#endif // CGAL_POLYGON_NO_POSTCONDITIONS

#if defined(CGAL_POLYGON_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
  || (!defined(CGAL_POLYGON_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || (!defined(CGAL_POLYGON_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_polygon_expensive_exactness_postcondition(EX) (static_cast<void>(0))
#  define CGAL_polygon_expensive_exactness_postcondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_polygon_expensive_exactness_postcondition_code(CODE)
#else
#  define CGAL_polygon_expensive_exactness_postcondition(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_polygon_expensive_exactness_postcondition_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_polygon_expensive_exactness_postcondition_code(CODE) CODE
#endif // CGAL_POLYGON_NO_POSTCONDITIONS


// warnings
// --------

#if defined(CGAL_POLYGON_NO_WARNINGS) || defined(CGAL_NO_WARNINGS) \
  || defined(NDEBUG)
#  define CGAL_polygon_warning(EX) (static_cast<void>(0))
#  define CGAL_polygon_warning_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_polygon_warning_code(CODE)
#else
#  define CGAL_polygon_warning(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_polygon_warning_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_polygon_warning_code(CODE) CODE
#endif // CGAL_POLYGON_NO_WARNINGS

#if defined(CGAL_POLYGON_NO_WARNINGS) || defined(CGAL_NO_WARNINGS) \
  || (!defined(CGAL_POLYGON_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || defined(NDEBUG)
#  define CGAL_polygon_exactness_warning(EX) (static_cast<void>(0))
#  define CGAL_polygon_exactness_warning_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_polygon_exactness_warning_code(CODE)
#else
#  define CGAL_polygon_exactness_warning(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_polygon_exactness_warning_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_polygon_exactness_warning_code(CODE) CODE
#endif // CGAL_POLYGON_NO_WARNINGS

#if defined(CGAL_POLYGON_NO_WARNINGS) || defined(CGAL_NO_WARNINGS) \
  || (!defined(CGAL_POLYGON_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_polygon_expensive_warning(EX) (static_cast<void>(0))
#  define CGAL_polygon_expensive_warning_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_polygon_expensive_warning_code(CODE)
#else
#  define CGAL_polygon_expensive_warning(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_polygon_expensive_warning_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_polygon_expensive_warning_code(CODE) CODE
#endif // CGAL_POLYGON_NO_WARNINGS

#if defined(CGAL_POLYGON_NO_WARNINGS) || defined(CGAL_NO_WARNINGS) \
  || (!defined(CGAL_POLYGON_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || (!defined(CGAL_POLYGON_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_polygon_expensive_exactness_warning(EX) (static_cast<void>(0))
#  define CGAL_polygon_expensive_exactness_warning_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_polygon_expensive_exactness_warning_code(CODE)
#else
#  define CGAL_polygon_expensive_exactness_warning(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_polygon_expensive_exactness_warning_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_polygon_expensive_exactness_warning_code(CODE) CODE
#endif // CGAL_POLYGON_NO_WARNINGS


