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
// release       : $CGAL_Revision: CGAL-2.0-I-1 $
// release_date  : $CGAL_Date: 1999/02/09 $
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
// ======================================================================

#ifndef CGAL_CONFIG_H
#define CGAL_CONFIG_H

#define CGAL_VERSION 1.1-I-01
#define CGAL_VERSION_NR 1001100001

#define CGAL_CFG_NO_ADVANCED_KERNEL 1

//----------------------------------------------------------------------//
//             STLport fix for MSVC
//----------------------------------------------------------------------//


#ifdef _MSC_VER
#   define CGAL_SCOPE
#   define CGAL_LIMITED_ITERATOR_TRAITS_SUPPORT 1
#   include <stl_config.h>
#   include <stl_iterator_base.h>
#else  // not _MSC_VER
#   define CGAL_SCOPE CGAL::
#   define CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(a)
#endif // _MSC_VER


//----------------------------------------------------------------------//
//             include platform specific workaround flags (CGAL_CFG_...)
//----------------------------------------------------------------------//

#include <CGAL/compiler_config.h>

//----------------------------------------------------------------------//
//             do some post processing for the flags
//----------------------------------------------------------------------//


#ifdef CGAL_CFG_TYPENAME_BUG
#   define CGAL_TYPENAME_MSVC_NULL
#else
#   define CGAL_TYPENAME_MSVC_NULL typename
#endif


#ifdef CGAL_CFG_NO_NAMESPACE
#  define CGAL_USING_NAMESPACE_STD
#  define CGAL_STD
#  define CGAL std
#else
#  define CGAL_USING_NAMESPACE_STD using namespace std;
#  define CGAL_STD std
#  ifndef CGAL_USE_NAMESPACE
#    define CGAL_USE_NAMESPACE 1
#  endif
#endif

#if CGAL_USE_NAMESPACE
#  define CGAL_BEGIN_NAMESPACE namespace CGAL {
#  define CGAL_END_NAMESPACE }
#else
#  define CGAL_BEGIN_NAMESPACE
#  define CGAL_END_NAMESPACE
#endif

#ifdef CGAL_CFG_NO_MUTABLE
#  define mutable
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


#ifdef CGAL_CFG_NO_STDC_NAMESPACE
#define CGAL_CLIB_STD
#else
#define CGAL_CLIB_STD std
#endif

//----------------------------------------------------------------------//
//             include separate workaround files
//----------------------------------------------------------------------//

#ifdef _MSC_VER
#  include <CGAL/MSVC_standard_header_fixes.h>
#endif
#if defined(__BORLANDC__) && __BORLANDC__ > 0x520
#include <CGAL/Borland_fixes.h>
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

