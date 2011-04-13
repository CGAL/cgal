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
// file          : include/CGAL/partition_assertions.h
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


#if defined(CGAL_PARTITION_NO_ASSERTIONS) || defined(CGAL_NO_ASSERTIONS) \
  || defined(NDEBUG)
#  define CGAL_partition_assertion(EX) ((void)0)
#  define CGAL_partition_assertion_msg(EX,MSG) ((void)0)
#  define CGAL_partition_assertion_code(CODE)
#else
#  define CGAL_partition_assertion(EX) \
   ((EX)?((void)0): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_partition_assertion_msg(EX,MSG) \
   ((EX)?((void)0): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_partition_assertion_code(CODE) CODE
#endif // CGAL_PARTITION_NO_ASSERTIONS

#if defined(CGAL_PARTITION_NO_ASSERTIONS) || defined(CGAL_NO_ASSERTIONS) \
  || (!defined(CGAL_PARTITION_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || defined(NDEBUG)
#  define CGAL_partition_exactness_assertion(EX) ((void)0)
#  define CGAL_partition_exactness_assertion_msg(EX,MSG) ((void)0)
#  define CGAL_partition_exactness_assertion_code(CODE)
#else
#  define CGAL_partition_exactness_assertion(EX) \
   ((EX)?((void)0): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_partition_exactness_assertion_msg(EX,MSG) \
   ((EX)?((void)0): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_partition_exactness_assertion_code(CODE) CODE
#endif // CGAL_PARTITION_NO_ASSERTIONS

#if defined(CGAL_PARTITION_NO_ASSERTIONS) \
  || defined(CGAL_NO_ASSERTIONS) \
  || (!defined(CGAL_PARTITION_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_partition_expensive_assertion(EX) ((void)0)
#  define CGAL_partition_expensive_assertion_msg(EX,MSG) ((void)0)
#  define CGAL_partition_expensive_assertion_code(CODE)
#else
#  define CGAL_partition_expensive_assertion(EX) \
   ((EX)?((void)0): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_partition_expensive_assertion_msg(EX,MSG) \
   ((EX)?((void)0): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_partition_expensive_assertion_code(CODE) CODE
#endif // CGAL_PARTITION_NO_ASSERTIONS

#if defined(CGAL_PARTITION_NO_ASSERTIONS) || defined(CGAL_NO_ASSERTIONS) \
  || (!defined(CGAL_PARTITION_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || (!defined(CGAL_PARTITION_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_partition_expensive_exactness_assertion(EX) ((void)0)
#  define CGAL_partition_expensive_exactness_assertion_msg(EX,MSG) ((void)0)
#  define CGAL_partition_expensive_exactness_assertion_code(CODE)
#else
#  define CGAL_partition_expensive_exactness_assertion(EX) \
   ((EX)?((void)0): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_partition_expensive_exactness_assertion_msg(EX,MSG) \
   ((EX)?((void)0): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_partition_expensive_exactness_assertion_code(CODE) CODE
#endif // CGAL_PARTITION_NO_ASSERTIONS


// preconditions
// -------------

#if defined(CGAL_PARTITION_NO_PRECONDITIONS) || defined(CGAL_NO_PRECONDITIONS) \
  || defined(NDEBUG)
#  define CGAL_partition_precondition(EX) ((void)0)
#  define CGAL_partition_precondition_msg(EX,MSG) ((void)0)
#  define CGAL_partition_precondition_code(CODE)
#else
#  define CGAL_partition_precondition(EX) \
   ((EX)?((void)0): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_partition_precondition_msg(EX,MSG) \
   ((EX)?((void)0): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_partition_precondition_code(CODE) CODE
#endif // CGAL_PARTITION_NO_PRECONDITIONS

#if defined(CGAL_PARTITION_NO_PRECONDITIONS) || defined(CGAL_NO_PRECONDITIONS) \
  || (!defined(CGAL_PARTITION_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || defined(NDEBUG)
#  define CGAL_partition_exactness_precondition(EX) ((void)0)
#  define CGAL_partition_exactness_precondition_msg(EX,MSG) ((void)0)
#  define CGAL_partition_exactness_precondition_code(CODE)
#else
#  define CGAL_partition_exactness_precondition(EX) \
   ((EX)?((void)0): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_partition_exactness_precondition_msg(EX,MSG) \
   ((EX)?((void)0): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_partition_exactness_precondition_code(CODE) CODE
#endif // CGAL_PARTITION_NO_PRECONDITIONS

#if defined(CGAL_PARTITION_NO_PRECONDITIONS) || defined(CGAL_NO_PRECONDITIONS) \
  || (!defined(CGAL_PARTITION_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_partition_expensive_precondition(EX) ((void)0)
#  define CGAL_partition_expensive_precondition_msg(EX,MSG) ((void)0)
#  define CGAL_partition_expensive_precondition_code(CODE)
#else
#  define CGAL_partition_expensive_precondition(EX) \
   ((EX)?((void)0): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_partition_expensive_precondition_msg(EX,MSG) \
   ((EX)?((void)0): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_partition_expensive_precondition_code(CODE) CODE
#endif // CGAL_PARTITION_NO_PRECONDITIONS

#if defined(CGAL_PARTITION_NO_PRECONDITIONS) || defined(CGAL_NO_PRECONDITIONS) \
  || (!defined(CGAL_PARTITION_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || (!defined(CGAL_PARTITION_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_partition_expensive_exactness_precondition(EX) ((void)0)
#  define CGAL_partition_expensive_exactness_precondition_msg(EX,MSG) ((void)0)
#  define CGAL_partition_expensive_exactness_precondition_code(CODE)
#else
#  define CGAL_partition_expensive_exactness_precondition(EX) \
   ((EX)?((void)0): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_partition_expensive_exactness_precondition_msg(EX,MSG) \
   ((EX)?((void)0): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_partition_expensive_exactness_precondition_code(CODE) CODE
#endif // CGAL_PARTITION_NO_PRECONDITIONS


// postconditions
// --------------

#if defined(CGAL_PARTITION_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
  || defined(NDEBUG)
#  define CGAL_partition_postcondition(EX) ((void)0)
#  define CGAL_partition_postcondition_msg(EX,MSG) ((void)0)
#  define CGAL_partition_postcondition_code(CODE)
#else
#  define CGAL_partition_postcondition(EX) \
   ((EX)?((void)0): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_partition_postcondition_msg(EX,MSG) \
   ((EX)?((void)0): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_partition_postcondition_code(CODE) CODE
#endif // CGAL_PARTITION_NO_POSTCONDITIONS

#if defined(CGAL_PARTITION_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
  || (!defined(CGAL_PARTITION_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || defined(NDEBUG)
#  define CGAL_partition_exactness_postcondition(EX) ((void)0)
#  define CGAL_partition_exactness_postcondition_msg(EX,MSG) ((void)0)
#  define CGAL_partition_exactness_postcondition_code(CODE)
#else
#  define CGAL_partition_exactness_postcondition(EX) \
   ((EX)?((void)0): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_partition_exactness_postcondition_msg(EX,MSG) \
   ((EX)?((void)0): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_partition_exactness_postcondition_code(CODE) CODE
#endif // CGAL_PARTITION_NO_POSTCONDITIONS

#if defined(CGAL_PARTITION_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
  || (!defined(CGAL_PARTITION_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_partition_expensive_postcondition(EX) ((void)0)
#  define CGAL_partition_expensive_postcondition_msg(EX,MSG) ((void)0)
#  define CGAL_partition_expensive_postcondition_code(CODE)
#else
#  define CGAL_partition_expensive_postcondition(EX) \
   ((EX)?((void)0): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_partition_expensive_postcondition_msg(EX,MSG) \
   ((EX)?((void)0): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_partition_expensive_postcondition_code(CODE) CODE
#endif // CGAL_PARTITION_NO_POSTCONDITIONS

#if defined(CGAL_PARTITION_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
  || (!defined(CGAL_PARTITION_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || (!defined(CGAL_PARTITION_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_partition_expensive_exactness_postcondition(EX) ((void)0)
#  define CGAL_partition_expensive_exactness_postcondition_msg(EX,MSG) ((void)0)
#  define CGAL_partition_expensive_exactness_postcondition_code(CODE)
#else
#  define CGAL_partition_expensive_exactness_postcondition(EX) \
   ((EX)?((void)0): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_partition_expensive_exactness_postcondition_msg(EX,MSG) \
   ((EX)?((void)0): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_partition_expensive_exactness_postcondition_code(CODE) CODE
#endif // CGAL_PARTITION_NO_POSTCONDITIONS


// warnings
// --------

#if defined(CGAL_PARTITION_NO_WARNINGS) || defined(CGAL_NO_WARNINGS) \
  || defined(NDEBUG)
#  define CGAL_partition_warning(EX) ((void)0)
#  define CGAL_partition_warning_msg(EX,MSG) ((void)0)
#  define CGAL_partition_warning_code(CODE)
#else
#  define CGAL_partition_warning(EX) \
   ((EX)?((void)0): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_partition_warning_msg(EX,MSG) \
   ((EX)?((void)0): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_partition_warning_code(CODE) CODE
#endif // CGAL_PARTITION_NO_WARNINGS

#if defined(CGAL_PARTITION_NO_WARNINGS) || defined(CGAL_NO_WARNINGS) \
  || (!defined(CGAL_PARTITION_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || defined(NDEBUG)
#  define CGAL_partition_exactness_warning(EX) ((void)0)
#  define CGAL_partition_exactness_warning_msg(EX,MSG) ((void)0)
#  define CGAL_partition_exactness_warning_code(CODE)
#else
#  define CGAL_partition_exactness_warning(EX) \
   ((EX)?((void)0): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_partition_exactness_warning_msg(EX,MSG) \
   ((EX)?((void)0): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_partition_exactness_warning_code(CODE) CODE
#endif // CGAL_PARTITION_NO_WARNINGS

#if defined(CGAL_PARTITION_NO_WARNINGS) || defined(CGAL_NO_WARNINGS) \
  || (!defined(CGAL_PARTITION_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_partition_expensive_warning(EX) ((void)0)
#  define CGAL_partition_expensive_warning_msg(EX,MSG) ((void)0)
#  define CGAL_partition_expensive_warning_code(CODE)
#else
#  define CGAL_partition_expensive_warning(EX) \
   ((EX)?((void)0): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_partition_expensive_warning_msg(EX,MSG) \
   ((EX)?((void)0): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_partition_expensive_warning_code(CODE) CODE
#endif // CGAL_PARTITION_NO_WARNINGS

#if defined(CGAL_PARTITION_NO_WARNINGS) || defined(CGAL_NO_WARNINGS) \
  || (!defined(CGAL_PARTITION_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || (!defined(CGAL_PARTITION_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_partition_expensive_exactness_warning(EX) ((void)0)
#  define CGAL_partition_expensive_exactness_warning_msg(EX,MSG) ((void)0)
#  define CGAL_partition_expensive_exactness_warning_code(CODE)
#else
#  define CGAL_partition_expensive_exactness_warning(EX) \
   ((EX)?((void)0): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_partition_expensive_exactness_warning_msg(EX,MSG) \
   ((EX)?((void)0): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_partition_expensive_exactness_warning_code(CODE) CODE
#endif // CGAL_PARTITION_NO_WARNINGS


