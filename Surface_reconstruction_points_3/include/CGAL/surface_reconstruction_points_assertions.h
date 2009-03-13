// Copyright (c) 2007  INRIA (France).
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
// $URL$
// $Id$
// 
//
// Author(s)     : Geert-Jan Giezeman, Sven Schönherr, Laurent Saboret
//
// Generated from script create_assertions.sh

/// @file surface_reconstruction_points_assertions.h
/// Define checking macros for the Surface_reconstruction_points_3 package

#include <CGAL/assertions.h>

#include <stdio.h>


// macro definitions
// =================

// assertions
// ----------

#if defined(CGAL_SURFACE_RECONSTRUCTION_POINTS_NO_ASSERTIONS) || defined(CGAL_NO_ASSERTIONS) \
  || defined(CGAL_NDEBUG)
#  define CGAL_surface_reconstruction_points_assertion(EX) (static_cast<void>(0))
#  define CGAL_surface_reconstruction_points_assertion_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_surface_reconstruction_points_assertion_code(CODE)
#else
#  define CGAL_surface_reconstruction_points_assertion(EX) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__))
#  define CGAL_surface_reconstruction_points_assertion_msg(EX,MSG) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_surface_reconstruction_points_assertion_code(CODE) CODE
#  define CGAL_surface_reconstruction_points_points_assertions 1
#endif // CGAL_SURFACE_RECONSTRUCTION_POINTS_NO_ASSERTIONS

#if defined(CGAL_SURFACE_RECONSTRUCTION_POINTS_NO_ASSERTIONS) || defined(CGAL_NO_ASSERTIONS) \
  || (!defined(CGAL_SURFACE_RECONSTRUCTION_POINTS_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || defined(CGAL_NDEBUG)
#  define CGAL_surface_reconstruction_points_exactness_assertion(EX) (static_cast<void>(0))
#  define CGAL_surface_reconstruction_points_exactness_assertion_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_surface_reconstruction_points_exactness_assertion_code(CODE)
#else
#  define CGAL_surface_reconstruction_points_exactness_assertion(EX) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__))
#  define CGAL_surface_reconstruction_points_exactness_assertion_msg(EX,MSG) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_surface_reconstruction_points_exactness_assertion_code(CODE) CODE
#  define CGAL_surface_reconstruction_points_exactness_assertions 1
#endif // CGAL_SURFACE_RECONSTRUCTION_POINTS_NO_ASSERTIONS

#if defined(CGAL_SURFACE_RECONSTRUCTION_POINTS_NO_ASSERTIONS) \
  || defined(CGAL_NO_ASSERTIONS) \
  || (!defined(CGAL_SURFACE_RECONSTRUCTION_POINTS_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(CGAL_NDEBUG)
#  define CGAL_surface_reconstruction_points_expensive_assertion(EX) (static_cast<void>(0))
#  define CGAL_surface_reconstruction_points_expensive_assertion_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_surface_reconstruction_points_expensive_assertion_code(CODE)
#else
#  define CGAL_surface_reconstruction_points_expensive_assertion(EX) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__))
#  define CGAL_surface_reconstruction_points_expensive_assertion_msg(EX,MSG) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_surface_reconstruction_points_expensive_assertion_code(CODE) CODE
#  define CGAL_surface_reconstruction_points_expensive_assertions 1
#endif // CGAL_SURFACE_RECONSTRUCTION_POINTS_NO_ASSERTIONS

#if defined(CGAL_SURFACE_RECONSTRUCTION_POINTS_NO_ASSERTIONS) || defined(CGAL_NO_ASSERTIONS) \
  || (!defined(CGAL_SURFACE_RECONSTRUCTION_POINTS_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || (!defined(CGAL_SURFACE_RECONSTRUCTION_POINTS_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(CGAL_NDEBUG)
#  define CGAL_surface_reconstruction_points_expensive_exactness_assertion(EX) (static_cast<void>(0))
#  define CGAL_surface_reconstruction_points_expensive_exactness_assertion_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_surface_reconstruction_points_expensive_exactness_assertion_code(CODE)
#else
#  define CGAL_surface_reconstruction_points_expensive_exactness_assertion(EX) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__))
#  define CGAL_surface_reconstruction_points_expensive_exactness_assertion_msg(EX,MSG) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_surface_reconstruction_points_expensive_exactness_assertion_code(CODE) CODE
#  define CGAL_surface_reconstruction_points_expensive_exactness_assertions 1
#endif // CGAL_SURFACE_RECONSTRUCTION_POINTS_NO_ASSERTIONS


// preconditions
// -------------

#if defined(CGAL_SURFACE_RECONSTRUCTION_POINTS_NO_PRECONDITIONS) || defined(CGAL_NO_PRECONDITIONS) \
  || defined(CGAL_NDEBUG)
#  define CGAL_surface_reconstruction_points_precondition(EX) (static_cast<void>(0))
#  define CGAL_surface_reconstruction_points_precondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_surface_reconstruction_points_precondition_code(CODE)
#else
#  define CGAL_surface_reconstruction_points_precondition(EX) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_surface_reconstruction_points_precondition_msg(EX,MSG) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_surface_reconstruction_points_precondition_code(CODE) CODE
#  define CGAL_surface_reconstruction_points_preconditions 1
#endif // CGAL_SURFACE_RECONSTRUCTION_POINTS_NO_PRECONDITIONS

#if defined(CGAL_SURFACE_RECONSTRUCTION_POINTS_NO_PRECONDITIONS) || defined(CGAL_NO_PRECONDITIONS) \
  || (!defined(CGAL_SURFACE_RECONSTRUCTION_POINTS_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || defined(CGAL_NDEBUG)
#  define CGAL_surface_reconstruction_points_exactness_precondition(EX) (static_cast<void>(0))
#  define CGAL_surface_reconstruction_points_exactness_precondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_surface_reconstruction_points_exactness_precondition_code(CODE)
#else
#  define CGAL_surface_reconstruction_points_exactness_precondition(EX) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_surface_reconstruction_points_exactness_precondition_msg(EX,MSG) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_surface_reconstruction_points_exactness_precondition_code(CODE) CODE
#  define CGAL_surface_reconstruction_points_exactness_preconditions 1
#endif // CGAL_SURFACE_RECONSTRUCTION_POINTS_NO_PRECONDITIONS

#if defined(CGAL_SURFACE_RECONSTRUCTION_POINTS_NO_PRECONDITIONS) || defined(CGAL_NO_PRECONDITIONS) \
  || (!defined(CGAL_SURFACE_RECONSTRUCTION_POINTS_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(CGAL_NDEBUG)
#  define CGAL_surface_reconstruction_points_expensive_precondition(EX) (static_cast<void>(0))
#  define CGAL_surface_reconstruction_points_expensive_precondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_surface_reconstruction_points_expensive_precondition_code(CODE)
#else
#  define CGAL_surface_reconstruction_points_expensive_precondition(EX) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_surface_reconstruction_points_expensive_precondition_msg(EX,MSG) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_surface_reconstruction_points_expensive_precondition_code(CODE) CODE
#  define CGAL_surface_reconstruction_points_expensive_preconditions 1
#endif // CGAL_SURFACE_RECONSTRUCTION_POINTS_NO_PRECONDITIONS

#if defined(CGAL_SURFACE_RECONSTRUCTION_POINTS_NO_PRECONDITIONS) || defined(CGAL_NO_PRECONDITIONS) \
  || (!defined(CGAL_SURFACE_RECONSTRUCTION_POINTS_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || (!defined(CGAL_SURFACE_RECONSTRUCTION_POINTS_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(CGAL_NDEBUG)
#  define CGAL_surface_reconstruction_points_expensive_exactness_precondition(EX) (static_cast<void>(0))
#  define CGAL_surface_reconstruction_points_expensive_exactness_precondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_surface_reconstruction_points_expensive_exactness_precondition_code(CODE)
#else
#  define CGAL_surface_reconstruction_points_expensive_exactness_precondition(EX) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_surface_reconstruction_points_expensive_exactness_precondition_msg(EX,MSG) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_surface_reconstruction_points_expensive_exactness_precondition_code(CODE) CODE
#  define CGAL_surface_reconstruction_points_expensive_exactness_preconditions 1
#endif // CGAL_SURFACE_RECONSTRUCTION_POINTS_NO_PRECONDITIONS


// postconditions
// --------------

#if defined(CGAL_SURFACE_RECONSTRUCTION_POINTS_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
  || defined(CGAL_NDEBUG)
#  define CGAL_surface_reconstruction_points_postcondition(EX) (static_cast<void>(0))
#  define CGAL_surface_reconstruction_points_postcondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_surface_reconstruction_points_postcondition_code(CODE)
#else
#  define CGAL_surface_reconstruction_points_postcondition(EX) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_surface_reconstruction_points_postcondition_msg(EX,MSG) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_surface_reconstruction_points_postcondition_code(CODE) CODE
#  define CGAL_surface_reconstruction_points_postconditions 1
#endif // CGAL_SURFACE_RECONSTRUCTION_POINTS_NO_POSTCONDITIONS

#if defined(CGAL_SURFACE_RECONSTRUCTION_POINTS_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
  || (!defined(CGAL_SURFACE_RECONSTRUCTION_POINTS_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || defined(CGAL_NDEBUG)
#  define CGAL_surface_reconstruction_points_exactness_postcondition(EX) (static_cast<void>(0))
#  define CGAL_surface_reconstruction_points_exactness_postcondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_surface_reconstruction_points_exactness_postcondition_code(CODE)
#else
#  define CGAL_surface_reconstruction_points_exactness_postcondition(EX) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_surface_reconstruction_points_exactness_postcondition_msg(EX,MSG) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_surface_reconstruction_points_exactness_postcondition_code(CODE) CODE
#  define CGAL_surface_reconstruction_points_exactness_postconditions 1
#endif // CGAL_SURFACE_RECONSTRUCTION_POINTS_NO_POSTCONDITIONS

#if defined(CGAL_SURFACE_RECONSTRUCTION_POINTS_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
  || (!defined(CGAL_SURFACE_RECONSTRUCTION_POINTS_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(CGAL_NDEBUG)
#  define CGAL_surface_reconstruction_points_expensive_postcondition(EX) (static_cast<void>(0))
#  define CGAL_surface_reconstruction_points_expensive_postcondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_surface_reconstruction_points_expensive_postcondition_code(CODE)
#else
#  define CGAL_surface_reconstruction_points_expensive_postcondition(EX) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_surface_reconstruction_points_expensive_postcondition_msg(EX,MSG) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_surface_reconstruction_points_expensive_postcondition_code(CODE) CODE
#  define CGAL_surface_reconstruction_points_expensive_postconditions 1
#endif // CGAL_SURFACE_RECONSTRUCTION_POINTS_NO_POSTCONDITIONS

#if defined(CGAL_SURFACE_RECONSTRUCTION_POINTS_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
  || (!defined(CGAL_SURFACE_RECONSTRUCTION_POINTS_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || (!defined(CGAL_SURFACE_RECONSTRUCTION_POINTS_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(CGAL_NDEBUG)
#  define CGAL_surface_reconstruction_points_expensive_exactness_postcondition(EX) (static_cast<void>(0))
#  define CGAL_surface_reconstruction_points_expensive_exactness_postcondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_surface_reconstruction_points_expensive_exactness_postcondition_code(CODE)
#else
#  define CGAL_surface_reconstruction_points_expensive_exactness_postcondition(EX) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_surface_reconstruction_points_expensive_exactness_postcondition_msg(EX,MSG) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_surface_reconstruction_points_expensive_exactness_postcondition_code(CODE) CODE
#  define CGAL_surface_reconstruction_points_expensive_exactness_postconditions 1
#endif // CGAL_SURFACE_RECONSTRUCTION_POINTS_NO_POSTCONDITIONS


// warnings
// --------

#if defined(CGAL_SURFACE_RECONSTRUCTION_POINTS_NO_WARNINGS) || defined(CGAL_NO_WARNINGS) \
  || defined(CGAL_NDEBUG)
#  define CGAL_surface_reconstruction_points_warning(EX) (static_cast<void>(0))
#  define CGAL_surface_reconstruction_points_warning_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_surface_reconstruction_points_warning_code(CODE)
#else
#  define CGAL_surface_reconstruction_points_warning(EX) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__))
#  define CGAL_surface_reconstruction_points_warning_msg(EX,MSG) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_surface_reconstruction_points_warning_code(CODE) CODE
#  define CGAL_surface_reconstruction_points_warnings 1
#endif // CGAL_SURFACE_RECONSTRUCTION_POINTS_NO_WARNINGS

#if defined(CGAL_SURFACE_RECONSTRUCTION_POINTS_NO_WARNINGS) || defined(CGAL_NO_WARNINGS) \
  || (!defined(CGAL_SURFACE_RECONSTRUCTION_POINTS_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || defined(CGAL_NDEBUG)
#  define CGAL_surface_reconstruction_points_exactness_warning(EX) (static_cast<void>(0))
#  define CGAL_surface_reconstruction_points_exactness_warning_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_surface_reconstruction_points_exactness_warning_code(CODE)
#else
#  define CGAL_surface_reconstruction_points_exactness_warning(EX) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__))
#  define CGAL_surface_reconstruction_points_exactness_warning_msg(EX,MSG) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_surface_reconstruction_points_exactness_warning_code(CODE) CODE
#  define CGAL_surface_reconstruction_points_exactness_warnings 1
#endif // CGAL_SURFACE_RECONSTRUCTION_POINTS_NO_WARNINGS

#if defined(CGAL_SURFACE_RECONSTRUCTION_POINTS_NO_WARNINGS) || defined(CGAL_NO_WARNINGS) \
  || (!defined(CGAL_SURFACE_RECONSTRUCTION_POINTS_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(CGAL_NDEBUG)
#  define CGAL_surface_reconstruction_points_expensive_warning(EX) (static_cast<void>(0))
#  define CGAL_surface_reconstruction_points_expensive_warning_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_surface_reconstruction_points_expensive_warning_code(CODE)
#else
#  define CGAL_surface_reconstruction_points_expensive_warning(EX) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__))
#  define CGAL_surface_reconstruction_points_expensive_warning_msg(EX,MSG) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_surface_reconstruction_points_expensive_warning_code(CODE) CODE
#  define CGAL_surface_reconstruction_points_expensive_warnings 1
#endif // CGAL_SURFACE_RECONSTRUCTION_POINTS_NO_WARNINGS

#if defined(CGAL_SURFACE_RECONSTRUCTION_POINTS_NO_WARNINGS) || defined(CGAL_NO_WARNINGS) \
  || (!defined(CGAL_SURFACE_RECONSTRUCTION_POINTS_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || (!defined(CGAL_SURFACE_RECONSTRUCTION_POINTS_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(CGAL_NDEBUG)
#  define CGAL_surface_reconstruction_points_expensive_exactness_warning(EX) (static_cast<void>(0))
#  define CGAL_surface_reconstruction_points_expensive_exactness_warning_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_surface_reconstruction_points_expensive_exactness_warning_code(CODE)
#else
#  define CGAL_surface_reconstruction_points_expensive_exactness_warning(EX) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__))
#  define CGAL_surface_reconstruction_points_expensive_exactness_warning_msg(EX,MSG) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_surface_reconstruction_points_expensive_exactness_warning_code(CODE) CODE
#  define CGAL_surface_reconstruction_points_expensive_exactness_warnings 1
#endif // CGAL_SURFACE_RECONSTRUCTION_POINTS_NO_WARNINGS


// Traces
// ------

#ifdef DEBUG_TRACE
  #define CGAL_TRACE  printf
#else
  #define CGAL_TRACE  if (false) printf
#endif

#ifdef DEBUG_TRACE
  #define CGAL_TRACE_STREAM  std::cerr
#else
  #define CGAL_TRACE_STREAM  if (false) std::cerr
#endif


