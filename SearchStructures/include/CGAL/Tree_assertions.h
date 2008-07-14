// Copyright (c) 1997  ETH Zurich (Switzerland).
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
// Author(s)     : Gabriele Neyer


#ifndef CGAL_TREE_ASSERTIONS_H
#define CGAL_TREE_ASSERTIONS_H

#  include <CGAL/assertions.h>


// macro definitions
// =================
// assertions
// ----------


#if defined(CGAL_TREE_NO_ASSERTIONS) || defined(CGAL_NO_ASSERTIONS) \
  || defined(NDEBUG)
#  define CGAL_Tree_assertion(EX) ((void)0)
#  define CGAL_Tree_assertion_msg(EX,MSG) ((void)0)
#  define CGAL_Tree_assertion_code(CODE)
#else
#  define CGAL_Tree_assertion(EX) \
   (CGAL::possibly(EX)?((void)0):assertion_fail( # EX , __FILE__, __LINE__))
#  define CGAL_Tree_assertion_msg(EX,MSG) \
   (CGAL::possibly(EX)?((void)0):assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_Tree_assertion_code(CODE) CODE
#endif // CGAL_TREE_NO_ASSERTIONS

#if defined(CGAL_TREE_NO_ASSERTIONS) || defined(CGAL_NO_ASSERTIONS) \
  || (!defined(CGAL_TREE_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || defined(NDEBUG)
#  define CGAL_Tree_exactness_assertion(EX) ((void)0)
#  define CGAL_Tree_exactness_assertion_msg(EX,MSG) ((void)0)
#  define CGAL_Tree_exactness_assertion_code(CODE)
#else
#  define CGAL_Tree_exactness_assertion(EX) \
   (CGAL::possibly(EX)?((void)0):assertion_fail( # EX , __FILE__, __LINE__))
#  define CGAL_Tree_exactness_assertion_msg(EX,MSG) \
   (CGAL::possibly(EX)?((void)0):assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_Tree_exactness_assertion_code(CODE) CODE
#endif // CGAL_TREE_NO_ASSERTIONS

#if defined(CGAL_TREE_NO_ASSERTIONS) \
  || defined(CGAL_NO_ASSERTIONS) \
  || (!defined(CGAL_TREE_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_Tree_expensive_assertion(EX) ((void)0)
#  define CGAL_Tree_expensive_assertion_msg(EX,MSG) ((void)0)
#  define CGAL_Tree_expensive_assertion_code(CODE)
#else
#  define CGAL_Tree_expensive_assertion(EX) \
   (CGAL::possibly(EX)?((void)0):assertion_fail( # EX , __FILE__, __LINE__))
#  define CGAL_Tree_expensive_assertion_msg(EX,MSG) \
   (CGAL::possibly(EX)?((void)0):assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_Tree_expensive_assertion_code(CODE) CODE
#endif // CGAL_TREE_NO_ASSERTIONS

#if defined(CGAL_TREE_NO_ASSERTIONS) || defined(CGAL_NO_ASSERTIONS) \
  || (!defined(CGAL_TREE_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || (!defined(CGAL_TREE_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_Tree_expensive_exactness_assertion(EX) ((void)0)
#  define CGAL_Tree_expensive_exactness_assertion_msg(EX,MSG) ((void)0)
#  define CGAL_Tree_expensive_exactness_assertion_code(CODE)
#else
#  define CGAL_Tree_expensive_exactness_assertion(EX) \
   (CGAL::possibly(EX)?((void)0):assertion_fail( # EX , __FILE__, __LINE__))
#  define CGAL_Tree_expensive_exactness_assertion_msg(EX,MSG) \
   (CGAL::possibly(EX)?((void)0):assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_Tree_expensive_exactness_assertion_code(CODE) CODE
#endif // CGAL_TREE_NO_ASSERTIONS


// preconditions
// -------------

#if defined(CGAL_TREE_NO_PRECONDITIONS) || defined(CGAL_NO_PRECONDITIONS) \
  || defined(NDEBUG)
#  define CGAL_Tree_precondition(EX) ((void)0)
#  define CGAL_Tree_precondition_msg(EX,MSG) ((void)0)
#  define CGAL_Tree_precondition_code(CODE)
#else
#  define CGAL_Tree_precondition(EX) \
   (CGAL::possibly(EX)?((void)0):precondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_Tree_precondition_msg(EX,MSG) \
   (CGAL::possibly(EX)?((void)0):precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_Tree_precondition_code(CODE) CODE
#endif // CGAL_TREE_NO_PRECONDITIONS

#if defined(CGAL_TREE_NO_PRECONDITIONS) || defined(CGAL_NO_PRECONDITIONS) \
  || (!defined(CGAL_TREE_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || defined(NDEBUG)
#  define CGAL_Tree_exactness_precondition(EX) ((void)0)
#  define CGAL_Tree_exactness_precondition_msg(EX,MSG) ((void)0)
#  define CGAL_Tree_exactness_precondition_code(CODE)
#else
#  define CGAL_Tree_exactness_precondition(EX) \
   (CGAL::possibly(EX)?((void)0):precondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_Tree_exactness_precondition_msg(EX,MSG) \
   (CGAL::possibly(EX)?((void)0):precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_Tree_exactness_precondition_code(CODE) CODE
#endif // CGAL_TREE_NO_PRECONDITIONS

#if defined(CGAL_TREE_NO_PRECONDITIONS) || defined(CGAL_NO_PRECONDITIONS) \
  || (!defined(CGAL_TREE_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_Tree_expensive_precondition(EX) ((void)0)
#  define CGAL_Tree_expensive_precondition_msg(EX,MSG) ((void)0)
#  define CGAL_Tree_expensive_precondition_code(CODE)
#else
#  define CGAL_Tree_expensive_precondition(EX) \
   (CGAL::possibly(EX)?((void)0):precondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_Tree_expensive_precondition_msg(EX,MSG) \
   (CGAL::possibly(EX)?((void)0):precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_Tree_expensive_precondition_code(CODE) CODE
#endif // CGAL_TREE_NO_PRECONDITIONS

#if defined(CGAL_TREE_NO_PRECONDITIONS) || defined(CGAL_NO_PRECONDITIONS) \
  || (!defined(CGAL_TREE_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || (!defined(CGAL_TREE_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_Tree_expensive_exactness_precondition(EX) ((void)0)
#  define CGAL_Tree_expensive_exactness_precondition_msg(EX,MSG) ((void)0)
#  define CGAL_Tree_expensive_exactness_precondition_code(CODE)
#else
#  define CGAL_Tree_expensive_exactness_precondition(EX) \
   (CGAL::possibly(EX)?((void)0):precondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_Tree_expensive_exactness_precondition_msg(EX,MSG) \
   (CGAL::possibly(EX)?((void)0):precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_Tree_expensive_exactness_precondition_code(CODE) CODE
#endif // CGAL_TREE_NO_PRECONDITIONS


// postconditions
// --------------

#if defined(CGAL_TREE_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
  || defined(NDEBUG)
#  define CGAL_Tree_postcondition(EX) ((void)0)
#  define CGAL_Tree_postcondition_msg(EX,MSG) ((void)0)
#  define CGAL_Tree_postcondition_code(CODE)
#else
#  define CGAL_Tree_postcondition(EX) \
   (CGAL::possibly(EX)?((void)0):postcondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_Tree_postcondition_msg(EX,MSG) \
   (CGAL::possibly(EX)?((void)0):postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_Tree_postcondition_code(CODE) CODE
#endif // CGAL_TREE_NO_POSTCONDITIONS

#if defined(CGAL_TREE_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
  || (!defined(CGAL_TREE_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || defined(NDEBUG)
#  define CGAL_Tree_exactness_postcondition(EX) ((void)0)
#  define CGAL_Tree_exactness_postcondition_msg(EX,MSG) ((void)0)
#  define CGAL_Tree_exactness_postcondition_code(CODE)
#else
#  define CGAL_Tree_exactness_postcondition(EX) \
   (CGAL::possibly(EX)?((void)0):postcondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_Tree_exactness_postcondition_msg(EX,MSG) \
   (CGAL::possibly(EX)?((void)0):postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_Tree_exactness_postcondition_code(CODE) CODE
#endif // CGAL_TREE_NO_POSTCONDITIONS

#if defined(CGAL_TREE_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
  || (!defined(CGAL_TREE_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_Tree_expensive_postcondition(EX) ((void)0)
#  define CGAL_Tree_expensive_postcondition_msg(EX,MSG) ((void)0)
#  define CGAL_Tree_expensive_postcondition_code(CODE)
#else
#  define CGAL_Tree_expensive_postcondition(EX) \
   (CGAL::possibly(EX)?((void)0):postcondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_Tree_expensive_postcondition_msg(EX,MSG) \
   (CGAL::possibly(EX)?((void)0):postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_Tree_expensive_postcondition_code(CODE) CODE
#endif // CGAL_TREE_NO_POSTCONDITIONS

#if defined(CGAL_TREE_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
  || (!defined(CGAL_TREE_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || (!defined(CGAL_TREE_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_Tree_expensive_exactness_postcondition(EX) ((void)0)
#  define CGAL_Tree_expensive_exactness_postcondition_msg(EX,MSG) ((void)0)
#  define CGAL_Tree_expensive_exactness_postcondition_code(CODE)
#else
#  define CGAL_Tree_expensive_exactness_postcondition(EX) \
   (CGAL::possibly(EX)?((void)0):postcondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_Tree_expensive_exactness_postcondition_msg(EX,MSG) \
   (CGAL::possibly(EX)?((void)0):postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_Tree_expensive_exactness_postcondition_code(CODE) CODE
#endif // CGAL_TREE_NO_POSTCONDITIONS


// warnings
// --------

#if defined(CGAL_TREE_NO_WARNINGS) || defined(CGAL_NO_WARNINGS) \
  || defined(NDEBUG)
#  define CGAL_Tree_warning(EX) ((void)0)
#  define CGAL_Tree_warning_msg(EX,MSG) ((void)0)
#  define CGAL_Tree_warning_code(CODE)
#else
#  define CGAL_Tree_warning(EX) \
   (CGAL::possibly(EX)?((void)0):warning_fail( # EX , __FILE__, __LINE__))
#  define CGAL_Tree_warning_msg(EX,MSG) \
   (CGAL::possibly(EX)?((void)0):warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_Tree_warning_code(CODE) CODE
#endif // CGAL_TREE_NO_WARNINGS

#if defined(CGAL_TREE_NO_WARNINGS) || defined(CGAL_NO_WARNINGS) \
  || (!defined(CGAL_TREE_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || defined(NDEBUG)
#  define CGAL_Tree_exactness_warning(EX) ((void)0)
#  define CGAL_Tree_exactness_warning_msg(EX,MSG) ((void)0)
#  define CGAL_Tree_exactness_warning_code(CODE)
#else
#  define CGAL_Tree_exactness_warning(EX) \
   (CGAL::possibly(EX)?((void)0):warning_fail( # EX , __FILE__, __LINE__))
#  define CGAL_Tree_exactness_warning_msg(EX,MSG) \
   (CGAL::possibly(EX)?((void)0):warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_Tree_exactness_warning_code(CODE) CODE
#endif // CGAL_TREE_NO_WARNINGS

#if defined(CGAL_TREE_NO_WARNINGS) || defined(CGAL_NO_WARNINGS) \
  || (!defined(CGAL_TREE_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_Tree_expensive_warning(EX) ((void)0)
#  define CGAL_Tree_expensive_warning_msg(EX,MSG) ((void)0)
#  define CGAL_Tree_expensive_warning_code(CODE)
#else
#  define CGAL_Tree_expensive_warning(EX) \
   (CGAL::possibly(EX)?((void)0):warning_fail( # EX , __FILE__, __LINE__))
#  define CGAL_Tree_expensive_warning_msg(EX,MSG) \
   (CGAL::possibly(EX)?((void)0):warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_Tree_expensive_warning_code(CODE) CODE
#endif // CGAL_TREE_NO_WARNINGS

#if defined(CGAL_TREE_NO_WARNINGS) || defined(CGAL_NO_WARNINGS) \
  || (!defined(CGAL_TREE_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || (!defined(CGAL_TREE_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_Tree_expensive_exactness_warning(EX) ((void)0)
#  define CGAL_Tree_expensive_exactness_warning_msg(EX,MSG) ((void)0)
#  define CGAL_Tree_expensive_exactness_warning_code(CODE)
#else
#  define CGAL_Tree_expensive_exactness_warning(EX) \
   (CGAL::possibly(EX)?((void)0):warning_fail( # EX , __FILE__, __LINE__))
#  define CGAL_Tree_expensive_exactness_warning_msg(EX,MSG) \
   (CGAL::possibly(EX)?((void)0):warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_Tree_expensive_exactness_warning_code(CODE) CODE
#endif // CGAL_TREE_NO_WARNINGS


#endif
