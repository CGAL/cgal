// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// release       : 
// release_date  : 
//
// file          : include/CGAL/ch_assertions.h
// package       : Convex_hull_2
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : script by Geert-Jan Giezeman and Sven Schönherr 
//
// coordinator   : MPI, Saarbruecken
// ======================================================================




// macro definitions
// =================
// assertions
// ----------


#if defined(CGAL_CH_NO_ASSERTIONS) || defined(CGAL_NO_ASSERTIONS) \
  || defined(NDEBUG)
#  define CGAL_ch_assertion(EX) ((void)0)
#  define CGAL_ch_assertion_msg(EX,MSG) ((void)0)
#  define CGAL_ch_assertion_code(CODE)
#else
#  define CGAL_ch_assertion(EX) \
   ((EX)?((void)0): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_ch_assertion_msg(EX,MSG) \
   ((EX)?((void)0): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_ch_assertion_code(CODE) CODE
#endif // CH_NO_ASSERTIONS

#if defined(CGAL_CH_NO_ASSERTIONS) || defined(CGAL_NO_ASSERTIONS) \
  || (!defined(CGAL_CH_CHECK_EXACTNESS) && !defined(CHECK_EXACTNESS))\
  || defined(NDEBUG)
#  define CGAL_ch_exactness_assertion(EX) ((void)0)
#  define CGAL_ch_exactness_assertion_msg(EX,MSG) ((void)0)
#  define CGAL_ch_exactness_assertion_code(CODE)
#else
#  define CGAL_ch_exactness_assertion(EX) \
   ((EX)?((void)0): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_ch_exactness_assertion_msg(EX,MSG) \
   ((EX)?((void)0): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_ch_exactness_assertion_code(CODE) CODE
#endif // CH_NO_ASSERTIONS

#if defined(CGAL_CH_NO_ASSERTIONS) \
  || defined(CGAL_NO_ASSERTIONS) \
  || (!defined(CGAL_CH_CHECK_EXPENSIVE) && !defined(CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_ch_expensive_assertion(EX) ((void)0)
#  define CGAL_ch_expensive_assertion_msg(EX,MSG) ((void)0)
#  define CGAL_ch_expensive_assertion_code(CODE)
#else
#  define CGAL_ch_expensive_assertion(EX) \
   ((EX)?((void)0): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_ch_expensive_assertion_msg(EX,MSG) \
   ((EX)?((void)0): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_ch_expensive_assertion_code(CODE) CODE
#endif // CH_NO_ASSERTIONS

#if defined(CGAL_CH_NO_ASSERTIONS) || defined(CGAL_NO_ASSERTIONS) \
  || (!defined(CGAL_CH_CHECK_EXACTNESS) && !defined(CHECK_EXACTNESS))\
  || (!defined(CGAL_CH_CHECK_EXPENSIVE) && !defined(CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_ch_expensive_exactness_assertion(EX) ((void)0)
#  define CGAL_ch_expensive_exactness_assertion_msg(EX,MSG) ((void)0)
#  define CGAL_ch_expensive_exactness_assertion_code(CODE)
#else
#  define CGAL_ch_expensive_exactness_assertion(EX) \
   ((EX)?((void)0): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_ch_expensive_exactness_assertion_msg(EX,MSG) \
   ((EX)?((void)0): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_ch_expensive_exactness_assertion_code(CODE) CODE
#endif // CH_NO_ASSERTIONS


// preconditions
// -------------

#if defined(CGAL_CH_NO_PRECONDITIONS) || defined(CGAL_NO_PRECONDITIONS) \
  || defined(NDEBUG)
#  define CGAL_ch_precondition(EX) ((void)0)
#  define CGAL_ch_precondition_msg(EX,MSG) ((void)0)
#  define CGAL_ch_precondition_code(CODE)
#else
#  define CGAL_ch_precondition(EX) \
   ((EX)?((void)0): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_ch_precondition_msg(EX,MSG) \
   ((EX)?((void)0): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_ch_precondition_code(CODE) CODE
#endif // CH_NO_PRECONDITIONS

#if defined(CGAL_CH_NO_PRECONDITIONS) || defined(CGAL_NO_PRECONDITIONS) \
  || (!defined(CGAL_CH_CHECK_EXACTNESS) && !defined(CHECK_EXACTNESS))\
  || defined(NDEBUG)
#  define CGAL_ch_exactness_precondition(EX) ((void)0)
#  define CGAL_ch_exactness_precondition_msg(EX,MSG) ((void)0)
#  define CGAL_ch_exactness_precondition_code(CODE)
#else
#  define CGAL_ch_exactness_precondition(EX) \
   ((EX)?((void)0): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_ch_exactness_precondition_msg(EX,MSG) \
   ((EX)?((void)0): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_ch_exactness_precondition_code(CODE) CODE
#endif // CH_NO_PRECONDITIONS

#if defined(CGAL_CH_NO_PRECONDITIONS) || defined(CGAL_NO_PRECONDITIONS) \
  || (!defined(CGAL_CH_CHECK_EXPENSIVE) && !defined(CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_ch_expensive_precondition(EX) ((void)0)
#  define CGAL_ch_expensive_precondition_msg(EX,MSG) ((void)0)
#  define CGAL_ch_expensive_precondition_code(CODE)
#else
#  define CGAL_ch_expensive_precondition(EX) \
   ((EX)?((void)0): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_ch_expensive_precondition_msg(EX,MSG) \
   ((EX)?((void)0): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_ch_expensive_precondition_code(CODE) CODE
#endif // CH_NO_PRECONDITIONS

#if defined(CGAL_CH_NO_PRECONDITIONS) || defined(CGAL_NO_PRECONDITIONS) \
  || (!defined(CGAL_CH_CHECK_EXACTNESS) && !defined(CHECK_EXACTNESS))\
  || (!defined(CGAL_CH_CHECK_EXPENSIVE) && !defined(CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_ch_expensive_exactness_precondition(EX) ((void)0)
#  define CGAL_ch_expensive_exactness_precondition_msg(EX,MSG) ((void)0)
#  define CGAL_ch_expensive_exactness_precondition_code(CODE)
#else
#  define CGAL_ch_expensive_exactness_precondition(EX) \
   ((EX)?((void)0): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_ch_expensive_exactness_precondition_msg(EX,MSG) \
   ((EX)?((void)0): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_ch_expensive_exactness_precondition_code(CODE) CODE
#endif // CH_NO_PRECONDITIONS


// postconditions
// --------------

#if defined(CGAL_CH_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
  || defined(NDEBUG)
#  define CGAL_ch_postcondition(EX) ((void)0)
#  define CGAL_ch_postcondition_msg(EX,MSG) ((void)0)
#  define CGAL_ch_postcondition_code(CODE)
#else
#  define CGAL_ch_postcondition(EX) \
   ((EX)?((void)0): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_ch_postcondition_msg(EX,MSG) \
   ((EX)?((void)0): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_ch_postcondition_code(CODE) CODE
#endif // CH_NO_POSTCONDITIONS

#if defined(CGAL_CH_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
  || (!defined(CGAL_CH_CHECK_EXACTNESS) && !defined(CHECK_EXACTNESS))\
  || defined(NDEBUG)
#  define CGAL_ch_exactness_postcondition(EX) ((void)0)
#  define CGAL_ch_exactness_postcondition_msg(EX,MSG) ((void)0)
#  define CGAL_ch_exactness_postcondition_code(CODE)
#else
#  define CGAL_ch_exactness_postcondition(EX) \
   ((EX)?((void)0): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_ch_exactness_postcondition_msg(EX,MSG) \
   ((EX)?((void)0): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_ch_exactness_postcondition_code(CODE) CODE
#endif // CH_NO_POSTCONDITIONS

#if defined(CGAL_CH_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
  || (!defined(CGAL_CH_CHECK_EXPENSIVE) && !defined(CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_ch_expensive_postcondition(EX) ((void)0)
#  define CGAL_ch_expensive_postcondition_msg(EX,MSG) ((void)0)
#  define CGAL_ch_expensive_postcondition_code(CODE)
#else
#  define CGAL_ch_expensive_postcondition(EX) \
   ((EX)?((void)0): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_ch_expensive_postcondition_msg(EX,MSG) \
   ((EX)?((void)0): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_ch_expensive_postcondition_code(CODE) CODE
#endif // CH_NO_POSTCONDITIONS

#if defined(CGAL_CH_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
  || (!defined(CGAL_CH_CHECK_EXACTNESS) && !defined(CHECK_EXACTNESS))\
  || (!defined(CGAL_CH_CHECK_EXPENSIVE) && !defined(CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_ch_expensive_exactness_postcondition(EX) ((void)0)
#  define CGAL_ch_expensive_exactness_postcondition_msg(EX,MSG) ((void)0)
#  define CGAL_ch_expensive_exactness_postcondition_code(CODE)
#else
#  define CGAL_ch_expensive_exactness_postcondition(EX) \
   ((EX)?((void)0): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_ch_expensive_exactness_postcondition_msg(EX,MSG) \
   ((EX)?((void)0): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_ch_expensive_exactness_postcondition_code(CODE) CODE
#endif // CH_NO_POSTCONDITIONS


// warnings
// --------

#if defined(CGAL_CH_NO_WARNINGS) || defined(CGAL_NO_WARNINGS) \
  || defined(NDEBUG)
#  define CGAL_ch_warning(EX) ((void)0)
#  define CGAL_ch_warning_msg(EX,MSG) ((void)0)
#  define CGAL_ch_warning_code(CODE)
#else
#  define CGAL_ch_warning(EX) \
   ((EX)?((void)0): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_ch_warning_msg(EX,MSG) \
   ((EX)?((void)0): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_ch_warning_code(CODE) CODE
#endif // CH_NO_WARNINGS

#if defined(CGAL_CH_NO_WARNINGS) || defined(CGAL_NO_WARNINGS) \
  || (!defined(CGAL_CH_CHECK_EXACTNESS) && !defined(CHECK_EXACTNESS))\
  || defined(NDEBUG)
#  define CGAL_ch_exactness_warning(EX) ((void)0)
#  define CGAL_ch_exactness_warning_msg(EX,MSG) ((void)0)
#  define CGAL_ch_exactness_warning_code(CODE)
#else
#  define CGAL_ch_exactness_warning(EX) \
   ((EX)?((void)0): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_ch_exactness_warning_msg(EX,MSG) \
   ((EX)?((void)0): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_ch_exactness_warning_code(CODE) CODE
#endif // CH_NO_WARNINGS

#if defined(CGAL_CH_NO_WARNINGS) || defined(CGAL_NO_WARNINGS) \
  || (!defined(CGAL_CH_CHECK_EXPENSIVE) && !defined(CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_ch_expensive_warning(EX) ((void)0)
#  define CGAL_ch_expensive_warning_msg(EX,MSG) ((void)0)
#  define CGAL_ch_expensive_warning_code(CODE)
#else
#  define CGAL_ch_expensive_warning(EX) \
   ((EX)?((void)0): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_ch_expensive_warning_msg(EX,MSG) \
   ((EX)?((void)0): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_ch_expensive_warning_code(CODE) CODE
#endif // CH_NO_WARNINGS

#if defined(CGAL_CH_NO_WARNINGS) || defined(CGAL_NO_WARNINGS) \
  || (!defined(CGAL_CH_CHECK_EXACTNESS) && !defined(CHECK_EXACTNESS))\
  || (!defined(CGAL_CH_CHECK_EXPENSIVE) && !defined(CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_ch_expensive_exactness_warning(EX) ((void)0)
#  define CGAL_ch_expensive_exactness_warning_msg(EX,MSG) ((void)0)
#  define CGAL_ch_expensive_exactness_warning_code(CODE)
#else
#  define CGAL_ch_expensive_exactness_warning(EX) \
   ((EX)?((void)0): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_ch_expensive_exactness_warning_msg(EX,MSG) \
   ((EX)?((void)0): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_ch_expensive_exactness_warning_code(CODE) CODE
#endif // CH_NO_WARNINGS


