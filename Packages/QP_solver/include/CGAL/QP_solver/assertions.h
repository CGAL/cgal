// ============================================================================
//
// Copyright (c) 1997-2004 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-I $
// release_date  : $CGAL_Date$
//
// file          : include/CGAL/QP_engine/assertions.h
// package       : $CGAL_Package: QP_engine $
// chapter       : Quadratic Programming Engine
//
// revision      : 3.0alpha
// revision_date : 2004/06
//
// author(s)     : Sven Schönherr <sven@inf.fu-berlin.de>
// coordinator   : ETH Zürich (Bernd Gärtner <gaertner@inf.ethz.ch>)
//
// implementation: assertion macros for the QP engine
// ======================================================================

#ifndef CGAL_QPE_ASSERTIONS_H
#define CGAL_QPE_ASSERTIONS_H

// =================
// macro definitions
// =================

// ----------
// assertions
// ----------
#if (    defined( CGAL_QPE_NO_ASSERTIONS) \
      || defined( CGAL_NO_ASSERTIONS) || defined( NDEBUG))
#  define  CGAL_qpe_assertion(EX)         ((void)0)
#  define  CGAL_qpe_assertion_msg(EX,MSG) ((void)0)
#  define  CGAL_qpe_assertion_code(CODE)
#else
#  define  CGAL_qpe_assertion(EX) \
     ((EX)?((void)0): ::CGAL::assertion_fail( # EX ,__FILE__,__LINE__,0))
#  define  CGAL_qpe_assertion_msg(EX,MSG) \
     ((EX)?((void)0): ::CGAL::assertion_fail( # EX ,__FILE__,__LINE__,MSG))
#  define  CGAL_qpe_assertion_code(CODE) CODE
#endif // qpe assertions

#if (    ! (    defined( CGAL_QPE_CHECK_EXACTNESS) \
             || defined( CGAL_CHECK_EXACTNESS)              ) \
      || defined( CGAL_QPE_NO_ASSERTIONS) \
      || defined( CGAL_NO_ASSERTIONS) || defined( NDEBUG))
#  define  CGAL_qpe_exactness_assertion(EX)         ((void)0)
#  define  CGAL_qpe_exactness_assertion_msg(EX,MSG) ((void)0)
#  define  CGAL_qpe_exactness_assertion_code(CODE)
#else
#  define  CGAL_qpe_exactness_assertion(EX) \
     ((EX)?((void)0): ::CGAL::assertion_fail( # EX ,__FILE__,__LINE__,0))
#  define  CGAL_qpe_exactness_assertion_msg(EX,MSG) \
     ((EX)?((void)0): ::CGAL::assertion_fail( # EX ,__FILE__,__LINE__,MSG))
#  define  CGAL_qpe_exactness_assertion_code(CODE) CODE
#endif // qpe exactness assertions

#if (    ! (    defined( CGAL_QPE_CHECK_EXPENSIVE) \
             || defined( CGAL_CHECK_EXPENSIVE)              ) \
      || defined( CGAL_QPE_NO_ASSERTIONS) \
      || defined( CGAL_NO_ASSERTIONS) || defined( NDEBUG))
#  define  CGAL_qpe_expensive_assertion(EX)         ((void)0)
#  define  CGAL_qpe_expensive_assertion_msg(EX,MSG) ((void)0)
#  define  CGAL_qpe_expensive_assertion_code(CODE)
#else
#  define  CGAL_qpe_expensive_assertion(EX) \
     ((EX)?((void)0): ::CGAL::assertion_fail( # EX ,__FILE__,__LINE__,0))
#  define  CGAL_qpe_expensive_assertion_msg(EX,MSG) \
     ((EX)?((void)0): ::CGAL::assertion_fail( # EX ,__FILE__,__LINE__,MSG))
#  define  CGAL_qpe_expensive_assertion_code(CODE) CODE
#endif // qpe expensive assertions

#if (    ! (    defined( CGAL_QPE_CHECK_EXACTNESS) \
             || defined( CGAL_QPE_CHECK_EXPENSIVE) \
             || defined( CGAL_CHECK_EXACTNESS)              \
             || defined( CGAL_CHECK_EXPENSIVE)              ) \
      || defined( CGAL_QPE_NO_ASSERTIONS) \
      || defined( CGAL_NO_ASSERTIONS) || defined( NDEBUG))
#  define  CGAL_qpe_expensive_exactness_assertion(EX) \
                                                                  ((void)0)
#  define  CGAL_qpe_expensive_exactness_assertion_msg(EX,MSG) \
                                                                  ((void)0)
#  define  CGAL_qpe_expensive_exactness_assertion_code(CODE)
#else
#  define  CGAL_qpe_expensive_exactness_assertion(EX) \
     ((EX)?((void)0): ::CGAL::assertion_fail( # EX ,__FILE__,__LINE__,0))
#  define  CGAL_qpe_expensive_exactness_assertion_msg(EX,MSG) \
     ((EX)?((void)0): ::CGAL::assertion_fail( # EX ,__FILE__,__LINE__,MSG))
#  define  CGAL_qpe_expensive_exactness_assertion_code(CODE) CODE
#endif // qpe expensive exactness assertions


// -------------
// preconditions
// -------------
#if (    defined( CGAL_QPE_NO_PRECONDITIONS) \
      || defined( CGAL_NO_PRECONDITIONS) || defined( NDEBUG))
#  define  CGAL_qpe_precondition(EX)         ((void)0)
#  define  CGAL_qpe_precondition_msg(EX,MSG) ((void)0)
#  define  CGAL_qpe_precondition_code(CODE)
#else
#  define  CGAL_qpe_precondition(EX) \
     ((EX)?((void)0): ::CGAL::precondition_fail( # EX ,__FILE__,__LINE__,0))
#  define  CGAL_qpe_precondition_msg(EX,MSG) \
     ((EX)?((void)0): ::CGAL::precondition_fail( # EX ,__FILE__,__LINE__,MSG))
#  define  CGAL_qpe_precondition_code(CODE) CODE
#endif // qpe preconditions

#if (    ! (    defined( CGAL_QPE_CHECK_EXACTNESS) \
             || defined( CGAL_CHECK_EXACTNESS)              ) \
      || defined( CGAL_QPE_NO_PRECONDITIONS) \
      || defined( CGAL_NO_PRECONDITIONS) || defined( NDEBUG))
#  define  CGAL_qpe_exactness_precondition(EX)         ((void)0)
#  define  CGAL_qpe_exactness_precondition_msg(EX,MSG) ((void)0)
#  define  CGAL_qpe_exactness_precondition_code(CODE)
#else
#  define  CGAL_qpe_exactness_precondition(EX) \
     ((EX)?((void)0): ::CGAL::precondition_fail( # EX ,__FILE__,__LINE__,0))
#  define  CGAL_qpe_exactness_precondition_msg(EX,MSG) \
     ((EX)?((void)0): ::CGAL::precondition_fail( # EX ,__FILE__,__LINE__,MSG))
#  define  CGAL_qpe_exactness_precondition_code(CODE) CODE
#endif // qpe exactness preconditions

#if (    ! (    defined( CGAL_QPE_CHECK_EXPENSIVE) \
             || defined( CGAL_CHECK_EXPENSIVE)              ) \
      || defined( CGAL_QPE_NO_PRECONDITIONS) \
      || defined( CGAL_NO_PRECONDITIONS) || defined( NDEBUG))
#  define  CGAL_qpe_expensive_precondition(EX)         ((void)0)
#  define  CGAL_qpe_expensive_precondition_msg(EX,MSG) ((void)0)
#  define  CGAL_qpe_expensive_precondition_code(CODE)
#else
#  define  CGAL_qpe_expensive_precondition(EX) \
     ((EX)?((void)0): ::CGAL::precondition_fail( # EX ,__FILE__,__LINE__,0))
#  define  CGAL_qpe_expensive_precondition_msg(EX,MSG) \
     ((EX)?((void)0): ::CGAL::precondition_fail( # EX ,__FILE__,__LINE__,MSG))
#  define  CGAL_qpe_expensive_precondition_code(CODE) CODE
#endif // qpe expensive preconditions

#if (    ! (    defined( CGAL_QPE_CHECK_EXACTNESS) \
             || defined( CGAL_QPE_CHECK_EXPENSIVE) \
             || defined( CGAL_CHECK_EXACTNESS)              \
             || defined( CGAL_CHECK_EXPENSIVE)              ) \
      || defined( CGAL_QPE_NO_PRECONDITIONS) \
      || defined( CGAL_NO_PRECONDITIONS) || defined( NDEBUG))
#  define  CGAL_qpe_expensive_exactness_precondition(EX) \
                                                                  ((void)0)
#  define  CGAL_qpe_expensive_exactness_precondition_msg(EX,MSG) \
                                                                  ((void)0)
#  define  CGAL_qpe_expensive_exactness_precondition_code(CODE)
#else
#  define  CGAL_qpe_expensive_exactness_precondition(EX) \
     ((EX)?((void)0): ::CGAL::precondition_fail( # EX ,__FILE__,__LINE__,0))
#  define  CGAL_qpe_expensive_exactness_precondition_msg(EX,MSG) \
     ((EX)?((void)0): ::CGAL::precondition_fail( # EX ,__FILE__,__LINE__,MSG))
#  define  CGAL_qpe_expensive_exactness_precondition_code(CODE) CODE
#endif // qpe expensive exactness preconditions


// --------------
// postconditions
// --------------
#if (    defined( CGAL_QPE_NO_POSTCONDITIONS) \
      || defined( CGAL_NO_POSTCONDITIONS) || defined( NDEBUG))
#  define  CGAL_qpe_postcondition(EX)         ((void)0)
#  define  CGAL_qpe_postcondition_msg(EX,MSG) ((void)0)
#  define  CGAL_qpe_postcondition_code(CODE)
#else
#  define  CGAL_qpe_postcondition(EX) \
     ((EX)?((void)0): ::CGAL::postcondition_fail( # EX ,__FILE__,__LINE__,0))
#  define  CGAL_qpe_postcondition_msg(EX,MSG) \
     ((EX)?((void)0): ::CGAL::postcondition_fail( # EX ,__FILE__,__LINE__,MSG))
#  define  CGAL_qpe_postcondition_code(CODE) CODE
#endif // qpe postconditions

#if (    ! (    defined( CGAL_QPE_CHECK_EXACTNESS) \
             || defined( CGAL_CHECK_EXACTNESS)              ) \
      || defined( CGAL_QPE_NO_POSTCONDITIONS) \
      || defined( CGAL_NO_POSTCONDITIONS) || defined( NDEBUG))
#  define  CGAL_qpe_exactness_postcondition(EX)         ((void)0)
#  define  CGAL_qpe_exactness_postcondition_msg(EX,MSG) ((void)0)
#  define  CGAL_qpe_exactness_postcondition_code(CODE)
#else
#  define  CGAL_qpe_exactness_postcondition(EX) \
     ((EX)?((void)0): ::CGAL::postcondition_fail( # EX ,__FILE__,__LINE__,0))
#  define  CGAL_qpe_exactness_postcondition_msg(EX,MSG) \
     ((EX)?((void)0): ::CGAL::postcondition_fail( # EX ,__FILE__,__LINE__,MSG))
#  define  CGAL_qpe_exactness_postcondition_code(CODE) CODE
#endif // qpe exactness postconditions

#if (    ! (    defined( CGAL_QPE_CHECK_EXPENSIVE) \
             || defined( CGAL_CHECK_EXPENSIVE)              ) \
      || defined( CGAL_QPE_NO_POSTCONDITIONS) \
      || defined( CGAL_NO_POSTCONDITIONS) || defined( NDEBUG))
#  define  CGAL_qpe_expensive_postcondition(EX)         ((void)0)
#  define  CGAL_qpe_expensive_postcondition_msg(EX,MSG) ((void)0)
#  define  CGAL_qpe_expensive_postcondition_code(CODE)
#else
#  define  CGAL_qpe_expensive_postcondition(EX) \
     ((EX)?((void)0): ::CGAL::postcondition_fail( # EX ,__FILE__,__LINE__,0))
#  define  CGAL_qpe_expensive_postcondition_msg(EX,MSG) \
     ((EX)?((void)0): ::CGAL::postcondition_fail( # EX ,__FILE__,__LINE__,MSG))
#  define  CGAL_qpe_expensive_postcondition_code(CODE) CODE
#endif // qpe expensive postconditions

#if (    ! (    defined( CGAL_QPE_CHECK_EXACTNESS) \
             || defined( CGAL_QPE_CHECK_EXPENSIVE) \
             || defined( CGAL_CHECK_EXACTNESS)              \
             || defined( CGAL_CHECK_EXPENSIVE)              ) \
      || defined( CGAL_QPE_NO_POSTCONDITIONS) \
      || defined( CGAL_NO_POSTCONDITIONS) || defined( NDEBUG))
#  define  CGAL_qpe_expensive_exactness_postcondition(EX) \
                                                                  ((void)0)
#  define  CGAL_qpe_expensive_exactness_postcondition_msg(EX,MSG) \
                                                                  ((void)0)
#  define  CGAL_qpe_expensive_exactness_postcondition_code(CODE)
#else
#  define  CGAL_qpe_expensive_exactness_postcondition(EX) \
     ((EX)?((void)0): ::CGAL::postcondition_fail( # EX ,__FILE__,__LINE__,0))
#  define  CGAL_qpe_expensive_exactness_postcondition_msg(EX,MSG) \
     ((EX)?((void)0): ::CGAL::postcondition_fail( # EX ,__FILE__,__LINE__,MSG))
#  define  CGAL_qpe_expensive_exactness_postcondition_code(CODE) CODE
#endif // qpe expensive exactness postconditions


// --------
// warnings
// --------
#if (    defined( CGAL_QPE_NO_WARNINGS) \
      || defined( CGAL_NO_WARNINGS) || defined( NDEBUG))
#  define  CGAL_qpe_warning(EX)         ((void)0)
#  define  CGAL_qpe_warning_msg(EX,MSG) ((void)0)
#  define  CGAL_qpe_warning_code(CODE)
#else
#  define  CGAL_qpe_warning(EX) \
     ((EX)?((void)0): ::CGAL::warning_fail( # EX ,__FILE__,__LINE__,0))
#  define  CGAL_qpe_warning_msg(EX,MSG) \
     ((EX)?((void)0): ::CGAL::warning_fail( # EX ,__FILE__,__LINE__,MSG))
#  define  CGAL_qpe_warning_code(CODE) CODE
#endif // qpe warnings

#if (    ! (    defined( CGAL_QPE_CHECK_EXACTNESS) \
             || defined( CGAL_CHECK_EXACTNESS)              ) \
      || defined( CGAL_QPE_NO_WARNINGS) \
      || defined( CGAL_NO_WARNINGS) || defined( NDEBUG))
#  define  CGAL_qpe_exactness_warning(EX)         ((void)0)
#  define  CGAL_qpe_exactness_warning_msg(EX,MSG) ((void)0)
#  define  CGAL_qpe_exactness_warning_code(CODE)
#else
#  define  CGAL_qpe_exactness_warning(EX) \
     ((EX)?((void)0): ::CGAL::warning_fail( # EX ,__FILE__,__LINE__,0))
#  define  CGAL_qpe_exactness_warning_msg(EX,MSG) \
     ((EX)?((void)0): ::CGAL::warning_fail( # EX ,__FILE__,__LINE__,MSG))
#  define  CGAL_qpe_exactness_warning_code(CODE) CODE
#endif // qpe exactness warnings

#if (    ! (    defined( CGAL_QPE_CHECK_EXPENSIVE) \
             || defined( CGAL_CHECK_EXPENSIVE)              ) \
      || defined( CGAL_QPE_NO_WARNINGS) \
      || defined( CGAL_NO_WARNINGS) || defined( NDEBUG))
#  define  CGAL_qpe_expensive_warning(EX)         ((void)0)
#  define  CGAL_qpe_expensive_warning_msg(EX,MSG) ((void)0)
#  define  CGAL_qpe_expensive_warning_code(CODE)
#else
#  define  CGAL_qpe_expensive_warning(EX) \
     ((EX)?((void)0): ::CGAL::warning_fail( # EX ,__FILE__,__LINE__,0))
#  define  CGAL_qpe_expensive_warning_msg(EX,MSG) \
     ((EX)?((void)0): ::CGAL::warning_fail( # EX ,__FILE__,__LINE__,MSG))
#  define  CGAL_qpe_expensive_warning_code(CODE) CODE
#endif // qpe expensive warnings

#if (    ! (    defined( CGAL_QPE_CHECK_EXACTNESS) \
             || defined( CGAL_QPE_CHECK_EXPENSIVE) \
             || defined( CGAL_CHECK_EXACTNESS)              \
             || defined( CGAL_CHECK_EXPENSIVE)              ) \
      || defined( CGAL_QPE_NO_WARNINGS) \
      || defined( CGAL_NO_WARNINGS) || defined( NDEBUG))
#  define  CGAL_qpe_expensive_exactness_warning(EX) \
                                                                  ((void)0)
#  define  CGAL_qpe_expensive_exactness_warning_msg(EX,MSG) \
                                                                  ((void)0)
#  define  CGAL_qpe_expensive_exactness_warning_code(CODE)
#else
#  define  CGAL_qpe_expensive_exactness_warning(EX) \
     ((EX)?((void)0): ::CGAL::warning_fail( # EX ,__FILE__,__LINE__,0))
#  define  CGAL_qpe_expensive_exactness_warning_msg(EX,MSG) \
     ((EX)?((void)0): ::CGAL::warning_fail( # EX ,__FILE__,__LINE__,MSG))
#  define  CGAL_qpe_expensive_exactness_warning_code(CODE) CODE
#endif // qpe expensive exactness warnings


#endif // CGAL_QPE_ASSERTIONS_H

// ===== EOF ==================================================================
