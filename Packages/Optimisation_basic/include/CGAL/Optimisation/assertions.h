// ============================================================================
//
// Copyright (c) 1997-2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-wip $
// release_date  : $CGAL_Date$
//
// file          : include/CGAL/Optimisation/assertions.h
// package       : Optimisation_basic 3.8.3 (20 Mar 2001)
// chapter       : Geometric Optimisation
//
// source        : web/Optimisation_basic.aw
// revision      : 5.9
// revision_date : 2001/03/20 16:54:14
//
// author(s)     : Geert-Jan Giezeman, Sven Schönherr
// maintainer    : Sven Schönherr <sven@inf.ethz.ch>
// coordinator   : ETH Zürich (Bernd Gärtner <gaertner@inf.ethz.ch>)
//
// implementation: assertion macros for optimisation algorithms
// ============================================================================

#ifndef CGAL_OPTIMISATION_ASSERTIONS_H
#define CGAL_OPTIMISATION_ASSERTIONS_H

// macro definitions
// =================

// assertions
// ----------
#if (    defined( CGAL_OPTIMISATION_NO_ASSERTIONS) \
      || defined( CGAL_NO_ASSERTIONS) || defined( NDEBUG))
#  define  CGAL_optimisation_assertion(EX)         ((void)0)
#  define  CGAL_optimisation_assertion_msg(EX,MSG) ((void)0)
#  define  CGAL_optimisation_assertion_code(CODE)
#else
#  define  CGAL_optimisation_assertion(EX) \
     ((EX)?((void)0): ::CGAL::assertion_fail( # EX ,__FILE__,__LINE__,0))
#  define  CGAL_optimisation_assertion_msg(EX,MSG) \
     ((EX)?((void)0): ::CGAL::assertion_fail( # EX ,__FILE__,__LINE__,MSG))
#  define  CGAL_optimisation_assertion_code(CODE) CODE
#endif // optimisation assertions

#if (    ! (    defined( CGAL_OPTIMISATION_CHECK_EXACTNESS) \
             || defined( CGAL_CHECK_EXACTNESS)              ) \
      || defined( CGAL_OPTIMISATION_NO_ASSERTIONS) \
      || defined( CGAL_NO_ASSERTIONS) || defined( NDEBUG))
#  define  CGAL_optimisation_exactness_assertion(EX)         ((void)0)
#  define  CGAL_optimisation_exactness_assertion_msg(EX,MSG) ((void)0)
#  define  CGAL_optimisation_exactness_assertion_code(CODE)
#else
#  define  CGAL_optimisation_exactness_assertion(EX) \
     ((EX)?((void)0): ::CGAL::assertion_fail( # EX ,__FILE__,__LINE__,0))
#  define  CGAL_optimisation_exactness_assertion_msg(EX,MSG) \
     ((EX)?((void)0): ::CGAL::assertion_fail( # EX ,__FILE__,__LINE__,MSG))
#  define  CGAL_optimisation_exactness_assertion_code(CODE) CODE
#endif // optimisation exactness assertions

#if (    ! (    defined( CGAL_OPTIMISATION_CHECK_EXPENSIVE) \
             || defined( CGAL_CHECK_EXPENSIVE)              ) \
      || defined( CGAL_OPTIMISATION_NO_ASSERTIONS) \
      || defined( CGAL_NO_ASSERTIONS) || defined( NDEBUG))
#  define  CGAL_optimisation_expensive_assertion(EX)         ((void)0)
#  define  CGAL_optimisation_expensive_assertion_msg(EX,MSG) ((void)0)
#  define  CGAL_optimisation_expensive_assertion_code(CODE)
#else
#  define  CGAL_optimisation_expensive_assertion(EX) \
     ((EX)?((void)0): ::CGAL::assertion_fail( # EX ,__FILE__,__LINE__,0))
#  define  CGAL_optimisation_expensive_assertion_msg(EX,MSG) \
     ((EX)?((void)0): ::CGAL::assertion_fail( # EX ,__FILE__,__LINE__,MSG))
#  define  CGAL_optimisation_expensive_assertion_code(CODE) CODE
#endif // optimisation expensive assertions

#if (    ! (    defined( CGAL_OPTIMISATION_CHECK_EXACTNESS) \
             || defined( CGAL_OPTIMISATION_CHECK_EXPENSIVE) \
             || defined( CGAL_CHECK_EXACTNESS)              \
             || defined( CGAL_CHECK_EXPENSIVE)              ) \
      || defined( CGAL_OPTIMISATION_NO_ASSERTIONS) \
      || defined( CGAL_NO_ASSERTIONS) || defined( NDEBUG))
#  define  CGAL_optimisation_expensive_exactness_assertion(EX) \
                                                                  ((void)0)
#  define  CGAL_optimisation_expensive_exactness_assertion_msg(EX,MSG) \
                                                                  ((void)0)
#  define  CGAL_optimisation_expensive_exactness_assertion_code(CODE)
#else
#  define  CGAL_optimisation_expensive_exactness_assertion(EX) \
     ((EX)?((void)0): ::CGAL::assertion_fail( # EX ,__FILE__,__LINE__,0))
#  define  CGAL_optimisation_expensive_exactness_assertion_msg(EX,MSG) \
     ((EX)?((void)0): ::CGAL::assertion_fail( # EX ,__FILE__,__LINE__,MSG))
#  define  CGAL_optimisation_expensive_exactness_assertion_code(CODE) CODE
#endif // optimisation expensive exactness assertions



// preconditions
// -------------
#if (    defined( CGAL_OPTIMISATION_NO_PRECONDITIONS) \
      || defined( CGAL_NO_PRECONDITIONS) || defined( NDEBUG))
#  define  CGAL_optimisation_precondition(EX)         ((void)0)
#  define  CGAL_optimisation_precondition_msg(EX,MSG) ((void)0)
#  define  CGAL_optimisation_precondition_code(CODE)
#else
#  define  CGAL_optimisation_precondition(EX) \
     ((EX)?((void)0): ::CGAL::precondition_fail( # EX ,__FILE__,__LINE__,0))
#  define  CGAL_optimisation_precondition_msg(EX,MSG) \
     ((EX)?((void)0): ::CGAL::precondition_fail( # EX ,__FILE__,__LINE__,MSG))
#  define  CGAL_optimisation_precondition_code(CODE) CODE
#endif // optimisation preconditions

#if (    ! (    defined( CGAL_OPTIMISATION_CHECK_EXACTNESS) \
             || defined( CGAL_CHECK_EXACTNESS)              ) \
      || defined( CGAL_OPTIMISATION_NO_PRECONDITIONS) \
      || defined( CGAL_NO_PRECONDITIONS) || defined( NDEBUG))
#  define  CGAL_optimisation_exactness_precondition(EX)         ((void)0)
#  define  CGAL_optimisation_exactness_precondition_msg(EX,MSG) ((void)0)
#  define  CGAL_optimisation_exactness_precondition_code(CODE)
#else
#  define  CGAL_optimisation_exactness_precondition(EX) \
     ((EX)?((void)0): ::CGAL::precondition_fail( # EX ,__FILE__,__LINE__,0))
#  define  CGAL_optimisation_exactness_precondition_msg(EX,MSG) \
     ((EX)?((void)0): ::CGAL::precondition_fail( # EX ,__FILE__,__LINE__,MSG))
#  define  CGAL_optimisation_exactness_precondition_code(CODE) CODE
#endif // optimisation exactness preconditions

#if (    ! (    defined( CGAL_OPTIMISATION_CHECK_EXPENSIVE) \
             || defined( CGAL_CHECK_EXPENSIVE)              ) \
      || defined( CGAL_OPTIMISATION_NO_PRECONDITIONS) \
      || defined( CGAL_NO_PRECONDITIONS) || defined( NDEBUG))
#  define  CGAL_optimisation_expensive_precondition(EX)         ((void)0)
#  define  CGAL_optimisation_expensive_precondition_msg(EX,MSG) ((void)0)
#  define  CGAL_optimisation_expensive_precondition_code(CODE)
#else
#  define  CGAL_optimisation_expensive_precondition(EX) \
     ((EX)?((void)0): ::CGAL::precondition_fail( # EX ,__FILE__,__LINE__,0))
#  define  CGAL_optimisation_expensive_precondition_msg(EX,MSG) \
     ((EX)?((void)0): ::CGAL::precondition_fail( # EX ,__FILE__,__LINE__,MSG))
#  define  CGAL_optimisation_expensive_precondition_code(CODE) CODE
#endif // optimisation expensive preconditions

#if (    ! (    defined( CGAL_OPTIMISATION_CHECK_EXACTNESS) \
             || defined( CGAL_OPTIMISATION_CHECK_EXPENSIVE) \
             || defined( CGAL_CHECK_EXACTNESS)              \
             || defined( CGAL_CHECK_EXPENSIVE)              ) \
      || defined( CGAL_OPTIMISATION_NO_PRECONDITIONS) \
      || defined( CGAL_NO_PRECONDITIONS) || defined( NDEBUG))
#  define  CGAL_optimisation_expensive_exactness_precondition(EX) \
                                                                  ((void)0)
#  define  CGAL_optimisation_expensive_exactness_precondition_msg(EX,MSG) \
                                                                  ((void)0)
#  define  CGAL_optimisation_expensive_exactness_precondition_code(CODE)
#else
#  define  CGAL_optimisation_expensive_exactness_precondition(EX) \
     ((EX)?((void)0): ::CGAL::precondition_fail( # EX ,__FILE__,__LINE__,0))
#  define  CGAL_optimisation_expensive_exactness_precondition_msg(EX,MSG) \
     ((EX)?((void)0): ::CGAL::precondition_fail( # EX ,__FILE__,__LINE__,MSG))
#  define  CGAL_optimisation_expensive_exactness_precondition_code(CODE) CODE
#endif // optimisation expensive exactness preconditions



// postconditions
// --------------
#if (    defined( CGAL_OPTIMISATION_NO_POSTCONDITIONS) \
      || defined( CGAL_NO_POSTCONDITIONS) || defined( NDEBUG))
#  define  CGAL_optimisation_postcondition(EX)         ((void)0)
#  define  CGAL_optimisation_postcondition_msg(EX,MSG) ((void)0)
#  define  CGAL_optimisation_postcondition_code(CODE)
#else
#  define  CGAL_optimisation_postcondition(EX) \
     ((EX)?((void)0): ::CGAL::postcondition_fail( # EX ,__FILE__,__LINE__,0))
#  define  CGAL_optimisation_postcondition_msg(EX,MSG) \
     ((EX)?((void)0): ::CGAL::postcondition_fail( # EX ,__FILE__,__LINE__,MSG))
#  define  CGAL_optimisation_postcondition_code(CODE) CODE
#endif // optimisation postconditions

#if (    ! (    defined( CGAL_OPTIMISATION_CHECK_EXACTNESS) \
             || defined( CGAL_CHECK_EXACTNESS)              ) \
      || defined( CGAL_OPTIMISATION_NO_POSTCONDITIONS) \
      || defined( CGAL_NO_POSTCONDITIONS) || defined( NDEBUG))
#  define  CGAL_optimisation_exactness_postcondition(EX)         ((void)0)
#  define  CGAL_optimisation_exactness_postcondition_msg(EX,MSG) ((void)0)
#  define  CGAL_optimisation_exactness_postcondition_code(CODE)
#else
#  define  CGAL_optimisation_exactness_postcondition(EX) \
     ((EX)?((void)0): ::CGAL::postcondition_fail( # EX ,__FILE__,__LINE__,0))
#  define  CGAL_optimisation_exactness_postcondition_msg(EX,MSG) \
     ((EX)?((void)0): ::CGAL::postcondition_fail( # EX ,__FILE__,__LINE__,MSG))
#  define  CGAL_optimisation_exactness_postcondition_code(CODE) CODE
#endif // optimisation exactness postconditions

#if (    ! (    defined( CGAL_OPTIMISATION_CHECK_EXPENSIVE) \
             || defined( CGAL_CHECK_EXPENSIVE)              ) \
      || defined( CGAL_OPTIMISATION_NO_POSTCONDITIONS) \
      || defined( CGAL_NO_POSTCONDITIONS) || defined( NDEBUG))
#  define  CGAL_optimisation_expensive_postcondition(EX)         ((void)0)
#  define  CGAL_optimisation_expensive_postcondition_msg(EX,MSG) ((void)0)
#  define  CGAL_optimisation_expensive_postcondition_code(CODE)
#else
#  define  CGAL_optimisation_expensive_postcondition(EX) \
     ((EX)?((void)0): ::CGAL::postcondition_fail( # EX ,__FILE__,__LINE__,0))
#  define  CGAL_optimisation_expensive_postcondition_msg(EX,MSG) \
     ((EX)?((void)0): ::CGAL::postcondition_fail( # EX ,__FILE__,__LINE__,MSG))
#  define  CGAL_optimisation_expensive_postcondition_code(CODE) CODE
#endif // optimisation expensive postconditions

#if (    ! (    defined( CGAL_OPTIMISATION_CHECK_EXACTNESS) \
             || defined( CGAL_OPTIMISATION_CHECK_EXPENSIVE) \
             || defined( CGAL_CHECK_EXACTNESS)              \
             || defined( CGAL_CHECK_EXPENSIVE)              ) \
      || defined( CGAL_OPTIMISATION_NO_POSTCONDITIONS) \
      || defined( CGAL_NO_POSTCONDITIONS) || defined( NDEBUG))
#  define  CGAL_optimisation_expensive_exactness_postcondition(EX) \
                                                                  ((void)0)
#  define  CGAL_optimisation_expensive_exactness_postcondition_msg(EX,MSG) \
                                                                  ((void)0)
#  define  CGAL_optimisation_expensive_exactness_postcondition_code(CODE)
#else
#  define  CGAL_optimisation_expensive_exactness_postcondition(EX) \
     ((EX)?((void)0): ::CGAL::postcondition_fail( # EX ,__FILE__,__LINE__,0))
#  define  CGAL_optimisation_expensive_exactness_postcondition_msg(EX,MSG) \
     ((EX)?((void)0): ::CGAL::postcondition_fail( # EX ,__FILE__,__LINE__,MSG))
#  define  CGAL_optimisation_expensive_exactness_postcondition_code(CODE) CODE
#endif // optimisation expensive exactness postconditions



// warnings
// --------
#if (    defined( CGAL_OPTIMISATION_NO_WARNINGS) \
      || defined( CGAL_NO_WARNINGS) || defined( NDEBUG))
#  define  CGAL_optimisation_warning(EX)         ((void)0)
#  define  CGAL_optimisation_warning_msg(EX,MSG) ((void)0)
#  define  CGAL_optimisation_warning_code(CODE)
#else
#  define  CGAL_optimisation_warning(EX) \
     ((EX)?((void)0): ::CGAL::warning_fail( # EX ,__FILE__,__LINE__,0))
#  define  CGAL_optimisation_warning_msg(EX,MSG) \
     ((EX)?((void)0): ::CGAL::warning_fail( # EX ,__FILE__,__LINE__,MSG))
#  define  CGAL_optimisation_warning_code(CODE) CODE
#endif // optimisation warnings

#if (    ! (    defined( CGAL_OPTIMISATION_CHECK_EXACTNESS) \
             || defined( CGAL_CHECK_EXACTNESS)              ) \
      || defined( CGAL_OPTIMISATION_NO_WARNINGS) \
      || defined( CGAL_NO_WARNINGS) || defined( NDEBUG))
#  define  CGAL_optimisation_exactness_warning(EX)         ((void)0)
#  define  CGAL_optimisation_exactness_warning_msg(EX,MSG) ((void)0)
#  define  CGAL_optimisation_exactness_warning_code(CODE)
#else
#  define  CGAL_optimisation_exactness_warning(EX) \
     ((EX)?((void)0): ::CGAL::warning_fail( # EX ,__FILE__,__LINE__,0))
#  define  CGAL_optimisation_exactness_warning_msg(EX,MSG) \
     ((EX)?((void)0): ::CGAL::warning_fail( # EX ,__FILE__,__LINE__,MSG))
#  define  CGAL_optimisation_exactness_warning_code(CODE) CODE
#endif // optimisation exactness warnings

#if (    ! (    defined( CGAL_OPTIMISATION_CHECK_EXPENSIVE) \
             || defined( CGAL_CHECK_EXPENSIVE)              ) \
      || defined( CGAL_OPTIMISATION_NO_WARNINGS) \
      || defined( CGAL_NO_WARNINGS) || defined( NDEBUG))
#  define  CGAL_optimisation_expensive_warning(EX)         ((void)0)
#  define  CGAL_optimisation_expensive_warning_msg(EX,MSG) ((void)0)
#  define  CGAL_optimisation_expensive_warning_code(CODE)
#else
#  define  CGAL_optimisation_expensive_warning(EX) \
     ((EX)?((void)0): ::CGAL::warning_fail( # EX ,__FILE__,__LINE__,0))
#  define  CGAL_optimisation_expensive_warning_msg(EX,MSG) \
     ((EX)?((void)0): ::CGAL::warning_fail( # EX ,__FILE__,__LINE__,MSG))
#  define  CGAL_optimisation_expensive_warning_code(CODE) CODE
#endif // optimisation expensive warnings

#if (    ! (    defined( CGAL_OPTIMISATION_CHECK_EXACTNESS) \
             || defined( CGAL_OPTIMISATION_CHECK_EXPENSIVE) \
             || defined( CGAL_CHECK_EXACTNESS)              \
             || defined( CGAL_CHECK_EXPENSIVE)              ) \
      || defined( CGAL_OPTIMISATION_NO_WARNINGS) \
      || defined( CGAL_NO_WARNINGS) || defined( NDEBUG))
#  define  CGAL_optimisation_expensive_exactness_warning(EX) \
                                                                  ((void)0)
#  define  CGAL_optimisation_expensive_exactness_warning_msg(EX,MSG) \
                                                                  ((void)0)
#  define  CGAL_optimisation_expensive_exactness_warning_code(CODE)
#else
#  define  CGAL_optimisation_expensive_exactness_warning(EX) \
     ((EX)?((void)0): ::CGAL::warning_fail( # EX ,__FILE__,__LINE__,0))
#  define  CGAL_optimisation_expensive_exactness_warning_msg(EX,MSG) \
     ((EX)?((void)0): ::CGAL::warning_fail( # EX ,__FILE__,__LINE__,MSG))
#  define  CGAL_optimisation_expensive_exactness_warning_code(CODE) CODE
#endif // optimisation expensive exactness warnings



#endif // CGAL_OPTIMISATION_ASSERTIONS_H

// ===== EOF ==================================================================
