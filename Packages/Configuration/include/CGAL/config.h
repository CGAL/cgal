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
// release       : $CGAL_Revision: CGAL-0.9-I-05 $
// release_date  : $CGAL_Date: 1997/12/17 $
//
// file          : include/CGAL/config.h
// source        :
// revision      : 1.11
// revision_date : 30 Mar 1998
// author(s)     : Wieger Wesselink <wieger@cs.ruu.nl>
//                 Michael Hoffmann
//
// coordinator   : Utrecht University
//
// ============================================================================

#ifndef CGAL_CONFIG_H
#define CGAL_CONFIG_H

#define CGAL_VERSION 1.1-I-01
#define CGAL_VERSION_NR 1001100001

//----------------------------------------------------------------------//
//             no namespaces for MIPS 7.2.1
//----------------------------------------------------------------------//

#if defined(__sgi) && !defined(__GNUC__) && defined(_COMPILER_VERSION)
#if (_COMPILER_VERSION >= 721) && defined(_NAMESPACES)
#define __STL_NO_NAMESPACES
#include <stl_config.h>
#undef __STL_USE_NAMESPACES
#endif
#endif

//----------------------------------------------------------------------//
//             include platform specific workaround flags (CGAL_CFG_...)
//----------------------------------------------------------------------//

#include <CGAL/compiler_config.h>

//----------------------------------------------------------------------//
//             do some post processing for the flags
//----------------------------------------------------------------------//

#ifdef CGAL_CFG_NO_TYPENAME
#  define typename
#endif

#ifdef CGAL_CFG_NO_EXPLICIT
#define explicit
#endif

#ifdef CGAL_CFG_NO_NAMESPACE
#  define CGAL_USING_NAMESPACE_STD
#  define CGAL_STD
#else
#  define CGAL_USING_NAMESPACE_STD using namespace std;
#  define CGAL_STD std
#  ifndef CGAL_USE_NAMESPACE
#    define CGAL_USE_NAMESPACE 1
#  endif
#endif

#if CGAL_USE_NAMESPACE
#  define CGAL_NAMESPACE_BEGIN namespace CGAL {
#  define CGAL_NAMESPACE_END }
#else
#  define CGAL_NAMESPACE_BEGIN
#  define CGAL_NAMESPACE_END
#endif

#ifdef CGAL_CFG_NO_MUTABLE
#  define mutable
#endif

// unset the flag CGAL_CFG_NO_SCOPE_MEMBER_FUNCTION_PARAMETERS if it is
// just the missing typename keyword that gives problems
#ifdef CGAL_CFG_NO_SCOPE_MEMBER_FUNCTION_PARAMETERS
#  ifndef CGAL_CFG_NO_SCOPE_MEMBER_FUNCTION_PARAMETERS_NO_TYPENAME
#    undef CGAL_CFG_NO_SCOPE_MEMBER_FUNCTION_PARAMETERS
#  endif
#endif

// unset the flag CGAL_CFG_RETURN_TYPE_BUG_1 if it is
// just the missing typename keyword that gives problems
#ifdef CGAL_CFG_RETURN_TYPE_BUG_1
#  ifndef CGAL_CFG_RETURN_TYPE_BUG_1_NO_TYPENAME
#   undef CGAL_CFG_RETURN_TYPE_BUG_1
#  endif
#endif

// unset the flag CGAL_CFG_INCOMPLETE_TYPE_BUG_4 if it is
// just the missing typename keyword that gives problems
#ifdef CGAL_CFG_INCOMPLETE_TYPE_BUG_4
#  ifndef CGAL_CFG_INCOMPLETE_TYPE_BUG_4_NO_TYPENAME
#    undef CGAL_CFG_INCOMPLETE_TYPE_BUG_4
#  endif
#endif

#ifdef CGAL_CFG_NO_TEMPLATE_FRIEND_DISTINCTION
#  define CGAL_NULL_TMPL_ARGS
#else
#  define CGAL_NULL_TMPL_ARGS <>
#endif

#ifdef CGAL_CFG_NO_EXPLICIT_CLASS_TEMPLATE_SPECIALISATION
#  define CGAL_TEMPLATE_NULL
#else
#  define CGAL_TEMPLATE_NULL template <>
#endif

//----------------------------------------------------------------------//
//             include separate workaround files
//----------------------------------------------------------------------//

#include <CGAL/workaround_return_type.h>
#include <CGAL/workaround_casts.h>
//#include <CGAL/workaround_stl.h>

//----------------------------------------------------------------------//
//             definition of type bool
//----------------------------------------------------------------------//

// if there is no built-in bool then we borrow the definition from STL
#ifdef CGAL_CFG_NO_BUILTIN_BOOL
#  include <pair.h>
#endif

//----------------------------------------------------------------------//
//             select old or new style headers
//----------------------------------------------------------------------//


#ifndef CGAL_USE_NEWSTYLE_HEADERS
#  ifndef CGAL_CFG_NO_STANDARD_HEADERS
#    ifndef CGAL_NO_NEWSTYLE_HEADERS
#      define CGAL_USE_NEWSTYLE_HEADERS
#    endif // ! CGAL_NO_NEWSTYLE_HEADERS
#  endif // ! CGAL_CFG_NO_STANDARD_HEADERS
#endif // ! CGAL_USE_NEWSTYLE_HEADERS

#endif // CGAL_CONFIG_H

