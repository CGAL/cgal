// Copyright (c) 1997  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Geert-Jan Giezeman, Sven Schoenherr, Gabriele Neyer
//
// Generated from script create_assertions.sh


// Note that this header file is intentionnaly not protected with a
// macro (as <cassert>). Calling it a second time with another value
// for NDEBUG for example must make a difference.

#include <CGAL/assertions.h>

// macro definitions
// =================
// assertions
// ----------

#undef CGAL_Tree_assertion
#undef CGAL_Tree_assertion_msg
#undef CGAL_Tree_assertion_code

#if defined(CGAL_TREE_NO_ASSERTIONS) || defined(CGAL_NO_ASSERTIONS) \
  || defined(NDEBUG)
#  define CGAL_Tree_assertion(EX) (static_cast<void>(0))

#include <CGAL/license/SearchStructures.h>

#  define CGAL_Tree_assertion_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_Tree_assertion_code(CODE)
#else
#  define CGAL_Tree_assertion(EX) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__))
#  define CGAL_Tree_assertion_msg(EX,MSG) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_Tree_assertion_code(CODE) CODE
#  define CGAL_Tree_assertions 1
#endif // CGAL_TREE_NO_ASSERTIONS


#undef CGAL_Tree_exactness_assertion
#undef CGAL_Tree_exactness_assertion_msg
#undef CGAL_Tree_exactness_assertion_code

#if defined(CGAL_TREE_NO_ASSERTIONS) || defined(CGAL_NO_ASSERTIONS) \
  || (!defined(CGAL_TREE_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || defined(NDEBUG)
#  define CGAL_Tree_exactness_assertion(EX) (static_cast<void>(0))
#  define CGAL_Tree_exactness_assertion_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_Tree_exactness_assertion_code(CODE)
#else
#  define CGAL_Tree_exactness_assertion(EX) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__))
#  define CGAL_Tree_exactness_assertion_msg(EX,MSG) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_Tree_exactness_assertion_code(CODE) CODE
#  define CGAL_Tree_exactness_assertions 1
#endif // CGAL_TREE_NO_ASSERTIONS


#undef CGAL_Tree_expensive_assertion
#undef CGAL_Tree_expensive_assertion_msg
#undef CGAL_Tree_expensive_assertion_code

#if defined(CGAL_TREE_NO_ASSERTIONS) \
  || defined(CGAL_NO_ASSERTIONS) \
  || (!defined(CGAL_TREE_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_Tree_expensive_assertion(EX) (static_cast<void>(0))
#  define CGAL_Tree_expensive_assertion_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_Tree_expensive_assertion_code(CODE)
#else
#  define CGAL_Tree_expensive_assertion(EX) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__))
#  define CGAL_Tree_expensive_assertion_msg(EX,MSG) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_Tree_expensive_assertion_code(CODE) CODE
#  define CGAL_Tree_expensive_assertions 1
#endif // CGAL_TREE_NO_ASSERTIONS


#undef CGAL_Tree_expensive_exactness_assertion
#undef CGAL_Tree_expensive_exactness_assertion_msg
#undef CGAL_Tree_expensive_exactness_assertion_code

#if defined(CGAL_TREE_NO_ASSERTIONS) || defined(CGAL_NO_ASSERTIONS) \
  || (!defined(CGAL_TREE_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || (!defined(CGAL_TREE_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_Tree_expensive_exactness_assertion(EX) (static_cast<void>(0))
#  define CGAL_Tree_expensive_exactness_assertion_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_Tree_expensive_exactness_assertion_code(CODE)
#else
#  define CGAL_Tree_expensive_exactness_assertion(EX) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__))
#  define CGAL_Tree_expensive_exactness_assertion_msg(EX,MSG) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_Tree_expensive_exactness_assertion_code(CODE) CODE
#  define CGAL_Tree_expensive_exactness_assertions 1
#endif // CGAL_TREE_NO_ASSERTIONS


// preconditions
// -------------

#undef CGAL_Tree_precondition
#undef CGAL_Tree_precondition_msg
#undef CGAL_Tree_precondition_code

#if defined(CGAL_TREE_NO_PRECONDITIONS) || defined(CGAL_NO_PRECONDITIONS) \
  || defined(NDEBUG)
#  define CGAL_Tree_precondition(EX) (static_cast<void>(0))
#  define CGAL_Tree_precondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_Tree_precondition_code(CODE)
#else
#  define CGAL_Tree_precondition(EX) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_Tree_precondition_msg(EX,MSG) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_Tree_precondition_code(CODE) CODE
#  define CGAL_Tree_preconditions 1
#endif // CGAL_TREE_NO_PRECONDITIONS


#undef CGAL_Tree_exactness_precondition
#undef CGAL_Tree_exactness_precondition_msg
#undef CGAL_Tree_exactness_precondition_code

#if defined(CGAL_TREE_NO_PRECONDITIONS) || defined(CGAL_NO_PRECONDITIONS) \
  || (!defined(CGAL_TREE_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || defined(NDEBUG)
#  define CGAL_Tree_exactness_precondition(EX) (static_cast<void>(0))
#  define CGAL_Tree_exactness_precondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_Tree_exactness_precondition_code(CODE)
#else
#  define CGAL_Tree_exactness_precondition(EX) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_Tree_exactness_precondition_msg(EX,MSG) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_Tree_exactness_precondition_code(CODE) CODE
#  define CGAL_Tree_exactness_preconditions 1
#endif // CGAL_TREE_NO_PRECONDITIONS


#undef CGAL_Tree_expensive_precondition
#undef CGAL_Tree_expensive_precondition_msg
#undef CGAL_Tree_expensive_precondition_code

#if defined(CGAL_TREE_NO_PRECONDITIONS) || defined(CGAL_NO_PRECONDITIONS) \
  || (!defined(CGAL_TREE_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_Tree_expensive_precondition(EX) (static_cast<void>(0))
#  define CGAL_Tree_expensive_precondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_Tree_expensive_precondition_code(CODE)
#else
#  define CGAL_Tree_expensive_precondition(EX) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_Tree_expensive_precondition_msg(EX,MSG) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_Tree_expensive_precondition_code(CODE) CODE
#  define CGAL_Tree_expensive_preconditions 1
#endif // CGAL_TREE_NO_PRECONDITIONS


#undef CGAL_Tree_expensive_exactness_precondition
#undef CGAL_Tree_expensive_exactness_precondition_msg
#undef CGAL_Tree_expensive_exactness_precondition_code

#if defined(CGAL_TREE_NO_PRECONDITIONS) || defined(CGAL_NO_PRECONDITIONS) \
  || (!defined(CGAL_TREE_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || (!defined(CGAL_TREE_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_Tree_expensive_exactness_precondition(EX) (static_cast<void>(0))
#  define CGAL_Tree_expensive_exactness_precondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_Tree_expensive_exactness_precondition_code(CODE)
#else
#  define CGAL_Tree_expensive_exactness_precondition(EX) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_Tree_expensive_exactness_precondition_msg(EX,MSG) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_Tree_expensive_exactness_precondition_code(CODE) CODE
#  define CGAL_Tree_expensive_exactness_preconditions 1
#endif // CGAL_TREE_NO_PRECONDITIONS


// postconditions
// --------------

#undef CGAL_Tree_postcondition
#undef CGAL_Tree_postcondition_msg
#undef CGAL_Tree_postcondition_code

#if defined(CGAL_TREE_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
  || defined(NDEBUG)
#  define CGAL_Tree_postcondition(EX) (static_cast<void>(0))
#  define CGAL_Tree_postcondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_Tree_postcondition_code(CODE)
#else
#  define CGAL_Tree_postcondition(EX) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_Tree_postcondition_msg(EX,MSG) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_Tree_postcondition_code(CODE) CODE
#  define CGAL_Tree_postconditions 1
#endif // CGAL_TREE_NO_POSTCONDITIONS


#undef CGAL_Tree_exactness_postcondition
#undef CGAL_Tree_exactness_postcondition_msg
#undef CGAL_Tree_exactness_postcondition_code

#if defined(CGAL_TREE_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
  || (!defined(CGAL_TREE_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || defined(NDEBUG)
#  define CGAL_Tree_exactness_postcondition(EX) (static_cast<void>(0))
#  define CGAL_Tree_exactness_postcondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_Tree_exactness_postcondition_code(CODE)
#else
#  define CGAL_Tree_exactness_postcondition(EX) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_Tree_exactness_postcondition_msg(EX,MSG) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_Tree_exactness_postcondition_code(CODE) CODE
#  define CGAL_Tree_exactness_postconditions 1
#endif // CGAL_TREE_NO_POSTCONDITIONS


#undef CGAL_Tree_expensive_postcondition
#undef CGAL_Tree_expensive_postcondition_msg
#undef CGAL_Tree_expensive_postcondition_code

#if defined(CGAL_TREE_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
  || (!defined(CGAL_TREE_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_Tree_expensive_postcondition(EX) (static_cast<void>(0))
#  define CGAL_Tree_expensive_postcondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_Tree_expensive_postcondition_code(CODE)
#else
#  define CGAL_Tree_expensive_postcondition(EX) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_Tree_expensive_postcondition_msg(EX,MSG) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_Tree_expensive_postcondition_code(CODE) CODE
#  define CGAL_Tree_expensive_postconditions 1
#endif // CGAL_TREE_NO_POSTCONDITIONS


#undef CGAL_Tree_expensive_exactness_postcondition
#undef CGAL_Tree_expensive_exactness_postcondition_msg
#undef CGAL_Tree_expensive_exactness_postcondition_code

#if defined(CGAL_TREE_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
  || (!defined(CGAL_TREE_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || (!defined(CGAL_TREE_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_Tree_expensive_exactness_postcondition(EX) (static_cast<void>(0))
#  define CGAL_Tree_expensive_exactness_postcondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_Tree_expensive_exactness_postcondition_code(CODE)
#else
#  define CGAL_Tree_expensive_exactness_postcondition(EX) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_Tree_expensive_exactness_postcondition_msg(EX,MSG) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_Tree_expensive_exactness_postcondition_code(CODE) CODE
#  define CGAL_Tree_expensive_exactness_postconditions 1
#endif // CGAL_TREE_NO_POSTCONDITIONS


// warnings
// --------

#undef CGAL_Tree_warning
#undef CGAL_Tree_warning_msg
#undef CGAL_Tree_warning_code

#if defined(CGAL_TREE_NO_WARNINGS) || defined(CGAL_NO_WARNINGS) \
  || defined(NDEBUG)
#  define CGAL_Tree_warning(EX) (static_cast<void>(0))
#  define CGAL_Tree_warning_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_Tree_warning_code(CODE)
#else
#  define CGAL_Tree_warning(EX) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__))
#  define CGAL_Tree_warning_msg(EX,MSG) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_Tree_warning_code(CODE) CODE
#  define CGAL_Tree_warnings 1
#endif // CGAL_TREE_NO_WARNINGS


#undef CGAL_Tree_exactness_warning
#undef CGAL_Tree_exactness_warning_msg
#undef CGAL_Tree_exactness_warning_code

#if defined(CGAL_TREE_NO_WARNINGS) || defined(CGAL_NO_WARNINGS) \
  || (!defined(CGAL_TREE_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || defined(NDEBUG)
#  define CGAL_Tree_exactness_warning(EX) (static_cast<void>(0))
#  define CGAL_Tree_exactness_warning_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_Tree_exactness_warning_code(CODE)
#else
#  define CGAL_Tree_exactness_warning(EX) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__))
#  define CGAL_Tree_exactness_warning_msg(EX,MSG) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_Tree_exactness_warning_code(CODE) CODE
#  define CGAL_Tree_exactness_warnings 1
#endif // CGAL_TREE_NO_WARNINGS


#undef CGAL_Tree_expensive_warning
#undef CGAL_Tree_expensive_warning_msg
#undef CGAL_Tree_expensive_warning_code

#if defined(CGAL_TREE_NO_WARNINGS) || defined(CGAL_NO_WARNINGS) \
  || (!defined(CGAL_TREE_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_Tree_expensive_warning(EX) (static_cast<void>(0))
#  define CGAL_Tree_expensive_warning_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_Tree_expensive_warning_code(CODE)
#else
#  define CGAL_Tree_expensive_warning(EX) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__))
#  define CGAL_Tree_expensive_warning_msg(EX,MSG) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_Tree_expensive_warning_code(CODE) CODE
#  define CGAL_Tree_expensive_warnings 1
#endif // CGAL_TREE_NO_WARNINGS


#undef CGAL_Tree_expensive_exactness_warning
#undef CGAL_Tree_expensive_exactness_warning_msg
#undef CGAL_Tree_expensive_exactness_warning_code

#if defined(CGAL_TREE_NO_WARNINGS) || defined(CGAL_NO_WARNINGS) \
  || (!defined(CGAL_TREE_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || (!defined(CGAL_TREE_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_Tree_expensive_exactness_warning(EX) (static_cast<void>(0))
#  define CGAL_Tree_expensive_exactness_warning_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_Tree_expensive_exactness_warning_code(CODE)
#else
#  define CGAL_Tree_expensive_exactness_warning(EX) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__))
#  define CGAL_Tree_expensive_exactness_warning_msg(EX,MSG) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_Tree_expensive_exactness_warning_code(CODE) CODE
#  define CGAL_Tree_expensive_exactness_warnings 1
#endif // CGAL_TREE_NO_WARNINGS
