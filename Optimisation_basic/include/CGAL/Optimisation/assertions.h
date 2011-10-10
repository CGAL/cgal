// Copyright (c) 1997-2001  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
// Author(s)     : Geert-Jan Giezeman, Sven Schoenherr <sven@inf.ethz.ch>

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
#  undef   CGAL_OPTIMISATION_ASSERTION_TAG
#else
#  define  CGAL_optimisation_assertion(EX) \
     (CGAL::possibly(EX)?((void)0): ::CGAL::assertion_fail( # EX ,__FILE__,__LINE__))
#  define  CGAL_optimisation_assertion_msg(EX,MSG) \
     (CGAL::possibly(EX)?((void)0): ::CGAL::assertion_fail( # EX ,__FILE__,__LINE__,MSG))
#  define  CGAL_optimisation_assertion_code(CODE) CODE
#  define  CGAL_OPTIMISATION_ASSERTION_TAG 1
#endif // optimisation assertions

#if (    ! (    defined( CGAL_OPTIMISATION_CHECK_EXACTNESS) \
             || defined( CGAL_CHECK_EXACTNESS)              ) \
      || defined( CGAL_OPTIMISATION_NO_ASSERTIONS) \
      || defined( CGAL_NO_ASSERTIONS) || defined( NDEBUG))
#  define  CGAL_optimisation_exactness_assertion(EX)         ((void)0)
#  define  CGAL_optimisation_exactness_assertion_msg(EX,MSG) ((void)0)
#  define  CGAL_optimisation_exactness_assertion_code(CODE)
#  undef   CGAL_OPTIMISATION_EXACTNESS_ASSERTION_TAG
#else
#  define  CGAL_optimisation_exactness_assertion(EX) \
     (CGAL::possibly(EX)?((void)0): ::CGAL::assertion_fail( # EX ,__FILE__,__LINE__))
#  define  CGAL_optimisation_exactness_assertion_msg(EX,MSG) \
     (CGAL::possibly(EX)?((void)0): ::CGAL::assertion_fail( # EX ,__FILE__,__LINE__,MSG))
#  define  CGAL_optimisation_exactness_assertion_code(CODE) CODE
#  define  CGAL_OPTIMISATION_EXACTNESS_ASSERTION_TAG 1
#endif // optimisation exactness assertions

#if (    ! (    defined( CGAL_OPTIMISATION_CHECK_EXPENSIVE) \
             || defined( CGAL_CHECK_EXPENSIVE)              ) \
      || defined( CGAL_OPTIMISATION_NO_ASSERTIONS) \
      || defined( CGAL_NO_ASSERTIONS) || defined( NDEBUG))
#  define  CGAL_optimisation_expensive_assertion(EX)         ((void)0)
#  define  CGAL_optimisation_expensive_assertion_msg(EX,MSG) ((void)0)
#  define  CGAL_optimisation_expensive_assertion_code(CODE)
#  undef   CGAL_OPTIMISATION_EXPENSIVE_ASSERTION_TAG
#else
#  define  CGAL_optimisation_expensive_assertion(EX) \
     (CGAL::possibly(EX)?((void)0): ::CGAL::assertion_fail( # EX ,__FILE__,__LINE__))
#  define  CGAL_optimisation_expensive_assertion_msg(EX,MSG) \
     (CGAL::possibly(EX)?((void)0): ::CGAL::assertion_fail( # EX ,__FILE__,__LINE__,MSG))
#  define  CGAL_optimisation_expensive_assertion_code(CODE) CODE
#  define  CGAL_OPTIMISATION_EXPENSIVE_ASSERTION_TAG 1
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
#  undef   CGAL_OPTIMISATION_EXPENSIVE_EXACTNESS_ASSERTION_TAG
#else
#  define  CGAL_optimisation_expensive_exactness_assertion(EX) \
     (CGAL::possibly(EX)?((void)0): ::CGAL::assertion_fail( # EX ,__FILE__,__LINE__))
#  define  CGAL_optimisation_expensive_exactness_assertion_msg(EX,MSG) \
     (CGAL::possibly(EX)?((void)0): ::CGAL::assertion_fail( # EX ,__FILE__,__LINE__,MSG))
#  define  CGAL_optimisation_expensive_exactness_assertion_code(CODE) CODE
#  define  CGAL_OPTIMISATION_EXPENSIVE_EXACTNESS_ASSERTION_TAG 1
#endif // optimisation expensive exactness assertions



// preconditions
// -------------
#if (    defined( CGAL_OPTIMISATION_NO_PRECONDITIONS) \
      || defined( CGAL_NO_PRECONDITIONS) || defined( NDEBUG))
#  define  CGAL_optimisation_precondition(EX)         ((void)0)
#  define  CGAL_optimisation_precondition_msg(EX,MSG) ((void)0)
#  define  CGAL_optimisation_precondition_code(CODE)
#  undef   CGAL_OPTIMISATION_PRECONDITION_TAG
#else
#  define  CGAL_optimisation_precondition(EX) \
     (CGAL::possibly(EX)?((void)0): ::CGAL::precondition_fail( # EX ,__FILE__,__LINE__))
#  define  CGAL_optimisation_precondition_msg(EX,MSG) \
     (CGAL::possibly(EX)?((void)0): ::CGAL::precondition_fail( # EX ,__FILE__,__LINE__,MSG))
#  define  CGAL_optimisation_precondition_code(CODE) CODE
#  define  CGAL_OPTIMISATION_PRECONDITION_TAG 1
#endif // optimisation preconditions

#if (    ! (    defined( CGAL_OPTIMISATION_CHECK_EXACTNESS) \
             || defined( CGAL_CHECK_EXACTNESS)              ) \
      || defined( CGAL_OPTIMISATION_NO_PRECONDITIONS) \
      || defined( CGAL_NO_PRECONDITIONS) || defined( NDEBUG))
#  define  CGAL_optimisation_exactness_precondition(EX)         ((void)0)
#  define  CGAL_optimisation_exactness_precondition_msg(EX,MSG) ((void)0)
#  define  CGAL_optimisation_exactness_precondition_code(CODE)
#  undef   CGAL_OPTIMISATION_EXACTNESS_PRECONDITION_TAG
#else
#  define  CGAL_optimisation_exactness_precondition(EX) \
     (CGAL::possibly(EX)?((void)0): ::CGAL::precondition_fail( # EX ,__FILE__,__LINE__))
#  define  CGAL_optimisation_exactness_precondition_msg(EX,MSG) \
     (CGAL::possibly(EX)?((void)0): ::CGAL::precondition_fail( # EX ,__FILE__,__LINE__,MSG))
#  define  CGAL_optimisation_exactness_precondition_code(CODE) CODE
#  define  CGAL_OPTIMISATION_EXACTNESS_PRECONDITION_TAG 1
#endif // optimisation exactness preconditions

#if (    ! (    defined( CGAL_OPTIMISATION_CHECK_EXPENSIVE) \
             || defined( CGAL_CHECK_EXPENSIVE)              ) \
      || defined( CGAL_OPTIMISATION_NO_PRECONDITIONS) \
      || defined( CGAL_NO_PRECONDITIONS) || defined( NDEBUG))
#  define  CGAL_optimisation_expensive_precondition(EX)         ((void)0)
#  define  CGAL_optimisation_expensive_precondition_msg(EX,MSG) ((void)0)
#  define  CGAL_optimisation_expensive_precondition_code(CODE)
#  undef   CGAL_OPTIMISATION_EXPENSIVE_PRECONDITION_TAG
#else
#  define  CGAL_optimisation_expensive_precondition(EX) \
     (CGAL::possibly(EX)?((void)0): ::CGAL::precondition_fail( # EX ,__FILE__,__LINE__))
#  define  CGAL_optimisation_expensive_precondition_msg(EX,MSG) \
     (CGAL::possibly(EX)?((void)0): ::CGAL::precondition_fail( # EX ,__FILE__,__LINE__,MSG))
#  define  CGAL_optimisation_expensive_precondition_code(CODE) CODE
#  define  CGAL_OPTIMISATION_EXPENSIVE_PRECONDITION_TAG 1
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
#  undef   CGAL_OPTIMISATION_EXPENSIVE_EXACTNESS_PRECONDITION_TAG
#else
#  define  CGAL_optimisation_expensive_exactness_precondition(EX) \
     (CGAL::possibly(EX)?((void)0): ::CGAL::precondition_fail( # EX ,__FILE__,__LINE__))
#  define  CGAL_optimisation_expensive_exactness_precondition_msg(EX,MSG) \
     (CGAL::possibly(EX)?((void)0): ::CGAL::precondition_fail( # EX ,__FILE__,__LINE__,MSG))
#  define  CGAL_optimisation_expensive_exactness_precondition_code(CODE) CODE
#  define  CGAL_OPTIMISATION_EXPENSIVE_EXACTNESS_PRECONDITION_TAG 1
#endif // optimisation expensive exactness preconditions



// postconditions
// --------------
#if (    defined( CGAL_OPTIMISATION_NO_POSTCONDITIONS) \
      || defined( CGAL_NO_POSTCONDITIONS) || defined( NDEBUG))
#  define  CGAL_optimisation_postcondition(EX)         ((void)0)
#  define  CGAL_optimisation_postcondition_msg(EX,MSG) ((void)0)
#  define  CGAL_optimisation_postcondition_code(CODE)
#  undef   CGAL_OPTIMISATION_POSTCONDITION_TAG
#else
#  define  CGAL_optimisation_postcondition(EX) \
     (CGAL::possibly(EX)?((void)0): ::CGAL::postcondition_fail( # EX ,__FILE__,__LINE__))
#  define  CGAL_optimisation_postcondition_msg(EX,MSG) \
     (CGAL::possibly(EX)?((void)0): ::CGAL::postcondition_fail( # EX ,__FILE__,__LINE__,MSG))
#  define  CGAL_optimisation_postcondition_code(CODE) CODE
#  define  CGAL_OPTIMISATION_POSTCONDITION_TAG 1
#endif // optimisation postconditions

#if (    ! (    defined( CGAL_OPTIMISATION_CHECK_EXACTNESS) \
             || defined( CGAL_CHECK_EXACTNESS)              ) \
      || defined( CGAL_OPTIMISATION_NO_POSTCONDITIONS) \
      || defined( CGAL_NO_POSTCONDITIONS) || defined( NDEBUG))
#  define  CGAL_optimisation_exactness_postcondition(EX)         ((void)0)
#  define  CGAL_optimisation_exactness_postcondition_msg(EX,MSG) ((void)0)
#  define  CGAL_optimisation_exactness_postcondition_code(CODE)
#  undef   CGAL_OPTIMISATION_EXACTNESS_POSTCONDITION_TAG
#else
#  define  CGAL_optimisation_exactness_postcondition(EX) \
     (CGAL::possibly(EX)?((void)0): ::CGAL::postcondition_fail( # EX ,__FILE__,__LINE__))
#  define  CGAL_optimisation_exactness_postcondition_msg(EX,MSG) \
     (CGAL::possibly(EX)?((void)0): ::CGAL::postcondition_fail( # EX ,__FILE__,__LINE__,MSG))
#  define  CGAL_optimisation_exactness_postcondition_code(CODE) CODE
#  define  CGAL_OPTIMISATION_EXACTNESS_POSTCONDITION_TAG 1
#endif // optimisation exactness postconditions

#if (    ! (    defined( CGAL_OPTIMISATION_CHECK_EXPENSIVE) \
             || defined( CGAL_CHECK_EXPENSIVE)              ) \
      || defined( CGAL_OPTIMISATION_NO_POSTCONDITIONS) \
      || defined( CGAL_NO_POSTCONDITIONS) || defined( NDEBUG))
#  define  CGAL_optimisation_expensive_postcondition(EX)         ((void)0)
#  define  CGAL_optimisation_expensive_postcondition_msg(EX,MSG) ((void)0)
#  define  CGAL_optimisation_expensive_postcondition_code(CODE)
#  undef   CGAL_OPTIMISATION_EXPENSIVE_POSTCONDITION_TAG
#else
#  define  CGAL_optimisation_expensive_postcondition(EX) \
     (CGAL::possibly(EX)?((void)0): ::CGAL::postcondition_fail( # EX ,__FILE__,__LINE__))
#  define  CGAL_optimisation_expensive_postcondition_msg(EX,MSG) \
     (CGAL::possibly(EX)?((void)0): ::CGAL::postcondition_fail( # EX ,__FILE__,__LINE__,MSG))
#  define  CGAL_optimisation_expensive_postcondition_code(CODE) CODE
#  define  CGAL_OPTIMISATION_EXPENSIVE_POSTCONDITION_TAG 1
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
#  undef   CGAL_OPTIMISATION_EXPENSIVE_EXACTNESS_POSTCONDITION_TAG
#else
#  define  CGAL_optimisation_expensive_exactness_postcondition(EX) \
     (CGAL::possibly(EX)?((void)0): ::CGAL::postcondition_fail( # EX ,__FILE__,__LINE__))
#  define  CGAL_optimisation_expensive_exactness_postcondition_msg(EX,MSG) \
     (CGAL::possibly(EX)?((void)0): ::CGAL::postcondition_fail( # EX ,__FILE__,__LINE__,MSG))
#  define  CGAL_optimisation_expensive_exactness_postcondition_code(CODE) CODE
#  define  CGAL_OPTIMISATION_EXPENSIVE_EXACTNESS_POSTCONDITION_TAG 1
#endif // optimisation expensive exactness postconditions



// warnings
// --------
#if (    defined( CGAL_OPTIMISATION_NO_WARNINGS) \
      || defined( CGAL_NO_WARNINGS) || defined( NDEBUG))
#  define  CGAL_optimisation_warning(EX)         ((void)0)
#  define  CGAL_optimisation_warning_msg(EX,MSG) ((void)0)
#  define  CGAL_optimisation_warning_code(CODE)
#  undef   CGAL_OPTIMISATION_WARNING_TAG
#else
#  define  CGAL_optimisation_warning(EX) \
     (CGAL::possibly(EX)?((void)0): ::CGAL::warning_fail( # EX ,__FILE__,__LINE__))
#  define  CGAL_optimisation_warning_msg(EX,MSG) \
     (CGAL::possibly(EX)?((void)0): ::CGAL::warning_fail( # EX ,__FILE__,__LINE__,MSG))
#  define  CGAL_optimisation_warning_code(CODE) CODE
#  define  CGAL_OPTIMISATION_WARNING_TAG 1
#endif // optimisation warnings

#if (    ! (    defined( CGAL_OPTIMISATION_CHECK_EXACTNESS) \
             || defined( CGAL_CHECK_EXACTNESS)              ) \
      || defined( CGAL_OPTIMISATION_NO_WARNINGS) \
      || defined( CGAL_NO_WARNINGS) || defined( NDEBUG))
#  define  CGAL_optimisation_exactness_warning(EX)         ((void)0)
#  define  CGAL_optimisation_exactness_warning_msg(EX,MSG) ((void)0)
#  define  CGAL_optimisation_exactness_warning_code(CODE)
#  undef   CGAL_OPTIMISATION_EXACTNESS_WARNING_TAG
#else
#  define  CGAL_optimisation_exactness_warning(EX) \
     (CGAL::possibly(EX)?((void)0): ::CGAL::warning_fail( # EX ,__FILE__,__LINE__))
#  define  CGAL_optimisation_exactness_warning_msg(EX,MSG) \
     (CGAL::possibly(EX)?((void)0): ::CGAL::warning_fail( # EX ,__FILE__,__LINE__,MSG))
#  define  CGAL_optimisation_exactness_warning_code(CODE) CODE
#  define  CGAL_OPTIMISATION_EXACTNESS_WARNING_TAG 1
#endif // optimisation exactness warnings

#if (    ! (    defined( CGAL_OPTIMISATION_CHECK_EXPENSIVE) \
             || defined( CGAL_CHECK_EXPENSIVE)              ) \
      || defined( CGAL_OPTIMISATION_NO_WARNINGS) \
      || defined( CGAL_NO_WARNINGS) || defined( NDEBUG))
#  define  CGAL_optimisation_expensive_warning(EX)         ((void)0)
#  define  CGAL_optimisation_expensive_warning_msg(EX,MSG) ((void)0)
#  define  CGAL_optimisation_expensive_warning_code(CODE)
#  undef   CGAL_OPTIMISATION_EXPENSIVE_WARNING_TAG
#else
#  define  CGAL_optimisation_expensive_warning(EX) \
     (CGAL::possibly(EX)?((void)0): ::CGAL::warning_fail( # EX ,__FILE__,__LINE__))
#  define  CGAL_optimisation_expensive_warning_msg(EX,MSG) \
     (CGAL::possibly(EX)?((void)0): ::CGAL::warning_fail( # EX ,__FILE__,__LINE__,MSG))
#  define  CGAL_optimisation_expensive_warning_code(CODE) CODE
#  define  CGAL_OPTIMISATION_EXPENSIVE_WARNING_TAG 1
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
#  undef   CGAL_OPTIMISATION_EXPENSIVE_EXACTNESS_WARNING_TAG
#else
#  define  CGAL_optimisation_expensive_exactness_warning(EX) \
     (CGAL::possibly(EX)?((void)0): ::CGAL::warning_fail( # EX ,__FILE__,__LINE__))
#  define  CGAL_optimisation_expensive_exactness_warning_msg(EX,MSG) \
     (CGAL::possibly(EX)?((void)0): ::CGAL::warning_fail( # EX ,__FILE__,__LINE__,MSG))
#  define  CGAL_optimisation_expensive_exactness_warning_code(CODE) CODE
#  define  CGAL_OPTIMISATION_EXPENSIVE_EXACTNESS_WARNING_TAG 1
#endif // optimisation expensive exactness warnings



#endif // CGAL_OPTIMISATION_ASSERTIONS_H

// ===== EOF ==================================================================
