// Copyright (c) 1997-2007  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
//
// Author(s)     : Sven Schoenherr
//                 Bernd Gaertner <gaertner@inf.ethz.ch>
//                 Franz Wessendorp 
//                 Kaspar Fischer

#ifndef CGAL_QP_ASSERTIONS_H
#define CGAL_QP_ASSERTIONS_H

#include <CGAL/license/QP_solver.h>


// =================
// macro definitions
// =================

// ----------
// assertions
// ----------
#if (    defined( CGAL_QP_NO_ASSERTIONS) \
      || defined( CGAL_NO_ASSERTIONS) || defined( NDEBUG))
#  define  CGAL_qpe_assertion(EX)         ((void)0)
#  define  CGAL_qpe_assertion_msg(EX,MSG) ((void)0)
#  define  CGAL_qpe_assertion_code(CODE)
#else
#  define  CGAL_qpe_assertion(EX) \
     (CGAL::possibly(EX)?((void)0): ::CGAL::assertion_fail( # EX ,__FILE__,__LINE__))
#  define  CGAL_qpe_assertion_msg(EX,MSG) \
     (CGAL::possibly(EX)?((void)0): ::CGAL::assertion_fail( # EX ,__FILE__,__LINE__,MSG))
#  define  CGAL_qpe_assertion_code(CODE) CODE
#endif // qpe assertions

#if (    ! (    defined( CGAL_QP_CHECK_EXACTNESS) \
             || defined( CGAL_CHECK_EXACTNESS)              ) \
      || defined( CGAL_QP_NO_ASSERTIONS) \
      || defined( CGAL_NO_ASSERTIONS) || defined( NDEBUG))
#  define  CGAL_qpe_exactness_assertion(EX)         ((void)0)
#  define  CGAL_qpe_exactness_assertion_msg(EX,MSG) ((void)0)
#  define  CGAL_qpe_exactness_assertion_code(CODE)
#else
#  define  CGAL_qpe_exactness_assertion(EX) \
     (CGAL::possibly(EX)?((void)0): ::CGAL::assertion_fail( # EX ,__FILE__,__LINE__))
#  define  CGAL_qpe_exactness_assertion_msg(EX,MSG) \
     (CGAL::possibly(EX)?((void)0): ::CGAL::assertion_fail( # EX ,__FILE__,__LINE__,MSG))
#  define  CGAL_qpe_exactness_assertion_code(CODE) CODE
#endif // qpe exactness assertions

#if (    ! (    defined( CGAL_QP_CHECK_EXPENSIVE) \
             || defined( CGAL_CHECK_EXPENSIVE)              ) \
      || defined( CGAL_QP_NO_ASSERTIONS) \
      || defined( CGAL_NO_ASSERTIONS) || defined( NDEBUG))
#  define  CGAL_qpe_expensive_assertion(EX)         ((void)0)
#  define  CGAL_qpe_expensive_assertion_msg(EX,MSG) ((void)0)
#  define  CGAL_qpe_expensive_assertion_code(CODE)
#else
#  define  CGAL_qpe_expensive_assertion(EX) \
     (CGAL::possibly(EX)?((void)0): ::CGAL::assertion_fail( # EX ,__FILE__,__LINE__))
#  define  CGAL_qpe_expensive_assertion_msg(EX,MSG) \
     (CGAL::possibly(EX)?((void)0): ::CGAL::assertion_fail( # EX ,__FILE__,__LINE__,MSG))
#  define  CGAL_qpe_expensive_assertion_code(CODE) CODE
#endif // qpe expensive assertions

#if (    ! (    defined( CGAL_QP_CHECK_EXACTNESS) \
             || defined( CGAL_QP_CHECK_EXPENSIVE) \
             || defined( CGAL_CHECK_EXACTNESS)              \
             || defined( CGAL_CHECK_EXPENSIVE)              ) \
      || defined( CGAL_QP_NO_ASSERTIONS) \
      || defined( CGAL_NO_ASSERTIONS) || defined( NDEBUG))
#  define  CGAL_qpe_expensive_exactness_assertion(EX) \
                                                                  ((void)0)
#  define  CGAL_qpe_expensive_exactness_assertion_msg(EX,MSG) \
                                                                  ((void)0)
#  define  CGAL_qpe_expensive_exactness_assertion_code(CODE)
#else
#  define  CGAL_qpe_expensive_exactness_assertion(EX) \
     (CGAL::possibly(EX)?((void)0): ::CGAL::assertion_fail( # EX ,__FILE__,__LINE__))
#  define  CGAL_qpe_expensive_exactness_assertion_msg(EX,MSG) \
     (CGAL::possibly(EX)?((void)0): ::CGAL::assertion_fail( # EX ,__FILE__,__LINE__,MSG))
#  define  CGAL_qpe_expensive_exactness_assertion_code(CODE) CODE
#endif // qpe expensive exactness assertions


// -------------
// preconditions
// -------------
#if (    defined( CGAL_QP_NO_PRECONDITIONS) \
      || defined( CGAL_NO_PRECONDITIONS) || defined( NDEBUG))
#  define  CGAL_qpe_precondition(EX)         ((void)0)
#  define  CGAL_qpe_precondition_msg(EX,MSG) ((void)0)
#  define  CGAL_qpe_precondition_code(CODE)
#else
#  define  CGAL_qpe_precondition(EX) \
     (CGAL::possibly(EX)?((void)0): ::CGAL::precondition_fail( # EX ,__FILE__,__LINE__))
#  define  CGAL_qpe_precondition_msg(EX,MSG) \
     (CGAL::possibly(EX)?((void)0): ::CGAL::precondition_fail( # EX ,__FILE__,__LINE__,MSG))
#  define  CGAL_qpe_precondition_code(CODE) CODE
#endif // qpe preconditions

#if (    ! (    defined( CGAL_QP_CHECK_EXACTNESS) \
             || defined( CGAL_CHECK_EXACTNESS)              ) \
      || defined( CGAL_QP_NO_PRECONDITIONS) \
      || defined( CGAL_NO_PRECONDITIONS) || defined( NDEBUG))
#  define  CGAL_qpe_exactness_precondition(EX)         ((void)0)
#  define  CGAL_qpe_exactness_precondition_msg(EX,MSG) ((void)0)
#  define  CGAL_qpe_exactness_precondition_code(CODE)
#else
#  define  CGAL_qpe_exactness_precondition(EX) \
     (CGAL::possibly(EX)?((void)0): ::CGAL::precondition_fail( # EX ,__FILE__,__LINE__))
#  define  CGAL_qpe_exactness_precondition_msg(EX,MSG) \
     (CGAL::possibly(EX)?((void)0): ::CGAL::precondition_fail( # EX ,__FILE__,__LINE__,MSG))
#  define  CGAL_qpe_exactness_precondition_code(CODE) CODE
#endif // qpe exactness preconditions

#if (    ! (    defined( CGAL_QP_CHECK_EXPENSIVE) \
             || defined( CGAL_CHECK_EXPENSIVE)              ) \
      || defined( CGAL_QP_NO_PRECONDITIONS) \
      || defined( CGAL_NO_PRECONDITIONS) || defined( NDEBUG))
#  define  CGAL_qpe_expensive_precondition(EX)         ((void)0)
#  define  CGAL_qpe_expensive_precondition_msg(EX,MSG) ((void)0)
#  define  CGAL_qpe_expensive_precondition_code(CODE)
#else
#  define  CGAL_qpe_expensive_precondition(EX) \
     (CGAL::possibly(EX)?((void)0): ::CGAL::precondition_fail( # EX ,__FILE__,__LINE__))
#  define  CGAL_qpe_expensive_precondition_msg(EX,MSG) \
     (CGAL::possibly(EX)?((void)0): ::CGAL::precondition_fail( # EX ,__FILE__,__LINE__,MSG))
#  define  CGAL_qpe_expensive_precondition_code(CODE) CODE
#endif // qpe expensive preconditions

#if (    ! (    defined( CGAL_QP_CHECK_EXACTNESS) \
             || defined( CGAL_QP_CHECK_EXPENSIVE) \
             || defined( CGAL_CHECK_EXACTNESS)              \
             || defined( CGAL_CHECK_EXPENSIVE)              ) \
      || defined( CGAL_QP_NO_PRECONDITIONS) \
      || defined( CGAL_NO_PRECONDITIONS) || defined( NDEBUG))
#  define  CGAL_qpe_expensive_exactness_precondition(EX) \
                                                                  ((void)0)
#  define  CGAL_qpe_expensive_exactness_precondition_msg(EX,MSG) \
                                                                  ((void)0)
#  define  CGAL_qpe_expensive_exactness_precondition_code(CODE)
#else
#  define  CGAL_qpe_expensive_exactness_precondition(EX) \
     (CGAL::possibly(EX)?((void)0): ::CGAL::precondition_fail( # EX ,__FILE__,__LINE__))
#  define  CGAL_qpe_expensive_exactness_precondition_msg(EX,MSG) \
     (CGAL::possibly(EX)?((void)0): ::CGAL::precondition_fail( # EX ,__FILE__,__LINE__,MSG))
#  define  CGAL_qpe_expensive_exactness_precondition_code(CODE) CODE
#endif // qpe expensive exactness preconditions


// --------------
// postconditions
// --------------
#if (    defined( CGAL_QP_NO_POSTCONDITIONS) \
      || defined( CGAL_NO_POSTCONDITIONS) || defined( NDEBUG))
#  define  CGAL_qpe_postcondition(EX)         ((void)0)
#  define  CGAL_qpe_postcondition_msg(EX,MSG) ((void)0)
#  define  CGAL_qpe_postcondition_code(CODE)
#else
#  define  CGAL_qpe_postcondition(EX) \
     (CGAL::possibly(EX)?((void)0): ::CGAL::postcondition_fail( # EX ,__FILE__,__LINE__))
#  define  CGAL_qpe_postcondition_msg(EX,MSG) \
     (CGAL::possibly(EX)?((void)0): ::CGAL::postcondition_fail( # EX ,__FILE__,__LINE__,MSG))
#  define  CGAL_qpe_postcondition_code(CODE) CODE
#endif // qpe postconditions

#if (    ! (    defined( CGAL_QP_CHECK_EXACTNESS) \
             || defined( CGAL_CHECK_EXACTNESS)              ) \
      || defined( CGAL_QP_NO_POSTCONDITIONS) \
      || defined( CGAL_NO_POSTCONDITIONS) || defined( NDEBUG))
#  define  CGAL_qpe_exactness_postcondition(EX)         ((void)0)
#  define  CGAL_qpe_exactness_postcondition_msg(EX,MSG) ((void)0)
#  define  CGAL_qpe_exactness_postcondition_code(CODE)
#else
#  define  CGAL_qpe_exactness_postcondition(EX) \
     (CGAL::possibly(EX)?((void)0): ::CGAL::postcondition_fail( # EX ,__FILE__,__LINE__))
#  define  CGAL_qpe_exactness_postcondition_msg(EX,MSG) \
     (CGAL::possibly(EX)?((void)0): ::CGAL::postcondition_fail( # EX ,__FILE__,__LINE__,MSG))
#  define  CGAL_qpe_exactness_postcondition_code(CODE) CODE
#endif // qpe exactness postconditions

#if (    ! (    defined( CGAL_QP_CHECK_EXPENSIVE) \
             || defined( CGAL_CHECK_EXPENSIVE)              ) \
      || defined( CGAL_QP_NO_POSTCONDITIONS) \
      || defined( CGAL_NO_POSTCONDITIONS) || defined( NDEBUG))
#  define  CGAL_qpe_expensive_postcondition(EX)         ((void)0)
#  define  CGAL_qpe_expensive_postcondition_msg(EX,MSG) ((void)0)
#  define  CGAL_qpe_expensive_postcondition_code(CODE)
#else
#  define  CGAL_qpe_expensive_postcondition(EX) \
     (CGAL::possibly(EX)?((void)0): ::CGAL::postcondition_fail( # EX ,__FILE__,__LINE__))
#  define  CGAL_qpe_expensive_postcondition_msg(EX,MSG) \
     (CGAL::possibly(EX)?((void)0): ::CGAL::postcondition_fail( # EX ,__FILE__,__LINE__,MSG))
#  define  CGAL_qpe_expensive_postcondition_code(CODE) CODE
#endif // qpe expensive postconditions

#if (    ! (    defined( CGAL_QP_CHECK_EXACTNESS) \
             || defined( CGAL_QP_CHECK_EXPENSIVE) \
             || defined( CGAL_CHECK_EXACTNESS)              \
             || defined( CGAL_CHECK_EXPENSIVE)              ) \
      || defined( CGAL_QP_NO_POSTCONDITIONS) \
      || defined( CGAL_NO_POSTCONDITIONS) || defined( NDEBUG))
#  define  CGAL_qpe_expensive_exactness_postcondition(EX) \
                                                                  ((void)0)
#  define  CGAL_qpe_expensive_exactness_postcondition_msg(EX,MSG) \
                                                                  ((void)0)
#  define  CGAL_qpe_expensive_exactness_postcondition_code(CODE)
#else
#  define  CGAL_qpe_expensive_exactness_postcondition(EX) \
     (CGAL::possibly(EX)?((void)0): ::CGAL::postcondition_fail( # EX ,__FILE__,__LINE__))
#  define  CGAL_qpe_expensive_exactness_postcondition_msg(EX,MSG) \
     (CGAL::possibly(EX)?((void)0): ::CGAL::postcondition_fail( # EX ,__FILE__,__LINE__,MSG))
#  define  CGAL_qpe_expensive_exactness_postcondition_code(CODE) CODE
#endif // qpe expensive exactness postconditions


// --------
// warnings
// --------
#if (    defined( CGAL_QP_NO_WARNINGS) \
      || defined( CGAL_NO_WARNINGS) || defined( NDEBUG))
#  define  CGAL_qpe_warning(EX)         ((void)0)
#  define  CGAL_qpe_warning_msg(EX,MSG) ((void)0)
#  define  CGAL_qpe_warning_code(CODE)
#else
#  define  CGAL_qpe_warning(EX) \
     (CGAL::possibly(EX)?((void)0): ::CGAL::warning_fail( # EX ,__FILE__,__LINE__))
#  define  CGAL_qpe_warning_msg(EX,MSG) \
     (CGAL::possibly(EX)?((void)0): ::CGAL::warning_fail( # EX ,__FILE__,__LINE__,MSG))
#  define  CGAL_qpe_warning_code(CODE) CODE
#endif // qpe warnings

#if (    ! (    defined( CGAL_QP_CHECK_EXACTNESS) \
             || defined( CGAL_CHECK_EXACTNESS)              ) \
      || defined( CGAL_QP_NO_WARNINGS) \
      || defined( CGAL_NO_WARNINGS) || defined( NDEBUG))
#  define  CGAL_qpe_exactness_warning(EX)         ((void)0)
#  define  CGAL_qpe_exactness_warning_msg(EX,MSG) ((void)0)
#  define  CGAL_qpe_exactness_warning_code(CODE)
#else
#  define  CGAL_qpe_exactness_warning(EX) \
     (CGAL::possibly(EX)?((void)0): ::CGAL::warning_fail( # EX ,__FILE__,__LINE__))
#  define  CGAL_qpe_exactness_warning_msg(EX,MSG) \
     (CGAL::possibly(EX)?((void)0): ::CGAL::warning_fail( # EX ,__FILE__,__LINE__,MSG))
#  define  CGAL_qpe_exactness_warning_code(CODE) CODE
#endif // qpe exactness warnings

#if (    ! (    defined( CGAL_QP_CHECK_EXPENSIVE) \
             || defined( CGAL_CHECK_EXPENSIVE)              ) \
      || defined( CGAL_QP_NO_WARNINGS) \
      || defined( CGAL_NO_WARNINGS) || defined( NDEBUG))
#  define  CGAL_qpe_expensive_warning(EX)         ((void)0)
#  define  CGAL_qpe_expensive_warning_msg(EX,MSG) ((void)0)
#  define  CGAL_qpe_expensive_warning_code(CODE)
#else
#  define  CGAL_qpe_expensive_warning(EX) \
     (CGAL::possibly(EX)?((void)0): ::CGAL::warning_fail( # EX ,__FILE__,__LINE__))
#  define  CGAL_qpe_expensive_warning_msg(EX,MSG) \
     (CGAL::possibly(EX)?((void)0): ::CGAL::warning_fail( # EX ,__FILE__,__LINE__,MSG))
#  define  CGAL_qpe_expensive_warning_code(CODE) CODE
#endif // qpe expensive warnings

#if (    ! (    defined( CGAL_QP_CHECK_EXACTNESS) \
             || defined( CGAL_QP_CHECK_EXPENSIVE) \
             || defined( CGAL_CHECK_EXACTNESS)              \
             || defined( CGAL_CHECK_EXPENSIVE)              ) \
      || defined( CGAL_QP_NO_WARNINGS) \
      || defined( CGAL_NO_WARNINGS) || defined( NDEBUG))
#  define  CGAL_qpe_expensive_exactness_warning(EX) \
                                                                  ((void)0)
#  define  CGAL_qpe_expensive_exactness_warning_msg(EX,MSG) \
                                                                  ((void)0)
#  define  CGAL_qpe_expensive_exactness_warning_code(CODE)
#else
#  define  CGAL_qpe_expensive_exactness_warning(EX) \
     (CGAL::possibly(EX)?((void)0): ::CGAL::warning_fail( # EX ,__FILE__,__LINE__))
#  define  CGAL_qpe_expensive_exactness_warning_msg(EX,MSG) \
     (CGAL::possibly(EX)?((void)0): ::CGAL::warning_fail( # EX ,__FILE__,__LINE__,MSG))
#  define  CGAL_qpe_expensive_exactness_warning_code(CODE) CODE
#endif // qpe expensive exactness warnings


#endif // CGAL_QP_ASSERTIONS_H

// ===== EOF ==================================================================
