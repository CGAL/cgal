// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : include/CGAL/kernel_assertions.h
// package       : Kernel_basic
// revision      :
// revision_date : $Date$
// author(s)     : script by Geert-Jan Giezeman and Sven Schönherr 
//
// coordinator   : MPI, Saarbruecken
//
// ======================================================================



// macro definitions
// =================
// assertions
// ----------


#if defined(CGAL_KERNEL_NO_ASSERTIONS) || defined(CGAL_NO_ASSERTIONS) \
  || defined(NDEBUG)
#  define CGAL_kernel_assertion(EX) (static_cast<void>(0))
#  define CGAL_kernel_assertion_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_kernel_assertion_code(CODE)
#else
#  define CGAL_kernel_assertion(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_kernel_assertion_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_kernel_assertion_code(CODE) CODE
#endif // CGAL_KERNEL_NO_ASSERTIONS

#if defined(CGAL_KERNEL_NO_ASSERTIONS) || defined(CGAL_NO_ASSERTIONS) \
  || (!defined(CGAL_KERNEL_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || defined(NDEBUG)
#  define CGAL_kernel_exactness_assertion(EX) (static_cast<void>(0))
#  define CGAL_kernel_exactness_assertion_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_kernel_exactness_assertion_code(CODE)
#else
#  define CGAL_kernel_exactness_assertion(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_kernel_exactness_assertion_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_kernel_exactness_assertion_code(CODE) CODE
#endif // CGAL_KERNEL_NO_ASSERTIONS

#if defined(CGAL_KERNEL_NO_ASSERTIONS) \
  || defined(CGAL_NO_ASSERTIONS) \
  || (!defined(CGAL_KERNEL_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_kernel_expensive_assertion(EX) (static_cast<void>(0))
#  define CGAL_kernel_expensive_assertion_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_kernel_expensive_assertion_code(CODE)
#else
#  define CGAL_kernel_expensive_assertion(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_kernel_expensive_assertion_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_kernel_expensive_assertion_code(CODE) CODE
#endif // CGAL_KERNEL_NO_ASSERTIONS

#if defined(CGAL_KERNEL_NO_ASSERTIONS) || defined(CGAL_NO_ASSERTIONS) \
  || (!defined(CGAL_KERNEL_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || (!defined(CGAL_KERNEL_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_kernel_expensive_exactness_assertion(EX) (static_cast<void>(0))
#  define CGAL_kernel_expensive_exactness_assertion_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_kernel_expensive_exactness_assertion_code(CODE)
#else
#  define CGAL_kernel_expensive_exactness_assertion(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_kernel_expensive_exactness_assertion_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_kernel_expensive_exactness_assertion_code(CODE) CODE
#endif // CGAL_KERNEL_NO_ASSERTIONS


// preconditions
// -------------

#if defined(CGAL_KERNEL_NO_PRECONDITIONS) || defined(CGAL_NO_PRECONDITIONS) \
  || defined(NDEBUG)
#  define CGAL_kernel_precondition(EX) (static_cast<void>(0))
#  define CGAL_kernel_precondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_kernel_precondition_code(CODE)
#else
#  define CGAL_kernel_precondition(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_kernel_precondition_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_kernel_precondition_code(CODE) CODE
#endif // CGAL_KERNEL_NO_PRECONDITIONS

#if defined(CGAL_KERNEL_NO_PRECONDITIONS) || defined(CGAL_NO_PRECONDITIONS) \
  || (!defined(CGAL_KERNEL_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || defined(NDEBUG)
#  define CGAL_kernel_exactness_precondition(EX) (static_cast<void>(0))
#  define CGAL_kernel_exactness_precondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_kernel_exactness_precondition_code(CODE)
#else
#  define CGAL_kernel_exactness_precondition(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_kernel_exactness_precondition_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_kernel_exactness_precondition_code(CODE) CODE
#endif // CGAL_KERNEL_NO_PRECONDITIONS

#if defined(CGAL_KERNEL_NO_PRECONDITIONS) || defined(CGAL_NO_PRECONDITIONS) \
  || (!defined(CGAL_KERNEL_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_kernel_expensive_precondition(EX) (static_cast<void>(0))
#  define CGAL_kernel_expensive_precondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_kernel_expensive_precondition_code(CODE)
#else
#  define CGAL_kernel_expensive_precondition(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_kernel_expensive_precondition_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_kernel_expensive_precondition_code(CODE) CODE
#endif // CGAL_KERNEL_NO_PRECONDITIONS

#if defined(CGAL_KERNEL_NO_PRECONDITIONS) || defined(CGAL_NO_PRECONDITIONS) \
  || (!defined(CGAL_KERNEL_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || (!defined(CGAL_KERNEL_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_kernel_expensive_exactness_precondition(EX) (static_cast<void>(0))
#  define CGAL_kernel_expensive_exactness_precondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_kernel_expensive_exactness_precondition_code(CODE)
#else
#  define CGAL_kernel_expensive_exactness_precondition(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_kernel_expensive_exactness_precondition_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_kernel_expensive_exactness_precondition_code(CODE) CODE
#endif // CGAL_KERNEL_NO_PRECONDITIONS


// postconditions
// --------------

#if defined(CGAL_KERNEL_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
  || defined(NDEBUG)
#  define CGAL_kernel_postcondition(EX) (static_cast<void>(0))
#  define CGAL_kernel_postcondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_kernel_postcondition_code(CODE)
#else
#  define CGAL_kernel_postcondition(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_kernel_postcondition_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_kernel_postcondition_code(CODE) CODE
#endif // CGAL_KERNEL_NO_POSTCONDITIONS

#if defined(CGAL_KERNEL_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
  || (!defined(CGAL_KERNEL_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || defined(NDEBUG)
#  define CGAL_kernel_exactness_postcondition(EX) (static_cast<void>(0))
#  define CGAL_kernel_exactness_postcondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_kernel_exactness_postcondition_code(CODE)
#else
#  define CGAL_kernel_exactness_postcondition(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_kernel_exactness_postcondition_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_kernel_exactness_postcondition_code(CODE) CODE
#endif // CGAL_KERNEL_NO_POSTCONDITIONS

#if defined(CGAL_KERNEL_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
  || (!defined(CGAL_KERNEL_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_kernel_expensive_postcondition(EX) (static_cast<void>(0))
#  define CGAL_kernel_expensive_postcondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_kernel_expensive_postcondition_code(CODE)
#else
#  define CGAL_kernel_expensive_postcondition(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_kernel_expensive_postcondition_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_kernel_expensive_postcondition_code(CODE) CODE
#endif // CGAL_KERNEL_NO_POSTCONDITIONS

#if defined(CGAL_KERNEL_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
  || (!defined(CGAL_KERNEL_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || (!defined(CGAL_KERNEL_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_kernel_expensive_exactness_postcondition(EX) (static_cast<void>(0))
#  define CGAL_kernel_expensive_exactness_postcondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_kernel_expensive_exactness_postcondition_code(CODE)
#else
#  define CGAL_kernel_expensive_exactness_postcondition(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_kernel_expensive_exactness_postcondition_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_kernel_expensive_exactness_postcondition_code(CODE) CODE
#endif // CGAL_KERNEL_NO_POSTCONDITIONS


// warnings
// --------

#if defined(CGAL_KERNEL_NO_WARNINGS) || defined(CGAL_NO_WARNINGS) \
  || defined(NDEBUG)
#  define CGAL_kernel_warning(EX) (static_cast<void>(0))
#  define CGAL_kernel_warning_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_kernel_warning_code(CODE)
#else
#  define CGAL_kernel_warning(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_kernel_warning_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_kernel_warning_code(CODE) CODE
#endif // CGAL_KERNEL_NO_WARNINGS

#if defined(CGAL_KERNEL_NO_WARNINGS) || defined(CGAL_NO_WARNINGS) \
  || (!defined(CGAL_KERNEL_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || defined(NDEBUG)
#  define CGAL_kernel_exactness_warning(EX) (static_cast<void>(0))
#  define CGAL_kernel_exactness_warning_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_kernel_exactness_warning_code(CODE)
#else
#  define CGAL_kernel_exactness_warning(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_kernel_exactness_warning_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_kernel_exactness_warning_code(CODE) CODE
#endif // CGAL_KERNEL_NO_WARNINGS

#if defined(CGAL_KERNEL_NO_WARNINGS) || defined(CGAL_NO_WARNINGS) \
  || (!defined(CGAL_KERNEL_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_kernel_expensive_warning(EX) (static_cast<void>(0))
#  define CGAL_kernel_expensive_warning_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_kernel_expensive_warning_code(CODE)
#else
#  define CGAL_kernel_expensive_warning(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_kernel_expensive_warning_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_kernel_expensive_warning_code(CODE) CODE
#endif // CGAL_KERNEL_NO_WARNINGS

#if defined(CGAL_KERNEL_NO_WARNINGS) || defined(CGAL_NO_WARNINGS) \
  || (!defined(CGAL_KERNEL_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || (!defined(CGAL_KERNEL_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_kernel_expensive_exactness_warning(EX) (static_cast<void>(0))
#  define CGAL_kernel_expensive_exactness_warning_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_kernel_expensive_exactness_warning_code(CODE)
#else
#  define CGAL_kernel_expensive_exactness_warning(EX) \
   ((EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_kernel_expensive_exactness_warning_msg(EX,MSG) \
   ((EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_kernel_expensive_exactness_warning_code(CODE) CODE
#endif // CGAL_KERNEL_NO_WARNINGS


