// Copyright (c) 1997-2004  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
// Author(s)     : Wieger Wesselink 
//                 Michael Hoffmann <hoffmann@inf.ethz.ch>
//                 Sylvain Pion

#ifndef CGAL_CONFIG_H
#define CGAL_CONFIG_H


#ifdef CGAL_INCLUDE_WINDOWS_DOT_H
// Mimic users including this file which defines min max macros
// and other names leading to name clashes
#include <windows.h>
#endif

// The following header file defines among other things  BOOST_PREVENT_MACRO_SUBSTITUTION 
#include <boost/config.hpp>

#include <CGAL/version.h>

//----------------------------------------------------------------------//
//             include platform specific workaround flags (CGAL_CFG_...)
//----------------------------------------------------------------------//

#include <CGAL/compiler_config.h>

//----------------------------------------------------------------------//
//        auto-link the CGAL library on platforms that support it
//----------------------------------------------------------------------//
#include <CGAL/auto_link/CGAL.h>

//----------------------------------------------------------------------//
//             do some post processing for the flags
//----------------------------------------------------------------------//

#ifdef CGAL_CFG_NO_STL
#  error "This compiler does not have a working STL"
#endif

// This macro computes the version number from an x.y.z release number.
// It only works for public releases.
#define CGAL_VERSION_NUMBER(x,y,z) (1000001 + 10000*x + 100*y + 10*z) * 1000

#define CGAL_BEGIN_NAMESPACE namespace CGAL {
#define CGAL_END_NAMESPACE }

#ifndef CGAL_CFG_NO_CPP0X_LONG_LONG
#  define CGAL_USE_LONG_LONG
#endif


#ifndef CGAL_CFG_TYPENAME_BEFORE_DEFAULT_ARGUMENT_BUG
#  define CGAL_TYPENAME_DEFAULT_ARG typename
#else
#  define CGAL_TYPENAME_DEFAULT_ARG
#endif


#ifdef CGAL_CFG_NO_CPP0X_DELETED_AND_DEFAULT_FUNCTIONS
#  define CGAL_DELETED
#  define CGAL_EXPLICITLY_DEFAULTED
#else
#  define CGAL_DELETED = delete
#  define CGAL_EXPLICITLY_DEFAULTED = default
#endif


// Big endian or little endian machine.
// ====================================

#if defined (__GLIBC__)
#  include <endian.h>
#  if (__BYTE_ORDER == __LITTLE_ENDIAN)
#    define CGAL_LITTLE_ENDIAN
#  elif (__BYTE_ORDER == __BIG_ENDIAN)
#    define CGAL_BIG_ENDIAN
#  else
#    error Unknown endianness
#  endif
#elif defined(__sparc) || defined(__sparc__) \
   || defined(_POWER) || defined(__powerpc__) \
   || defined(__ppc__) || defined(__hppa) \
   || defined(_MIPSEB) || defined(_POWER) \
   || defined(__s390__)
#  define CGAL_BIG_ENDIAN
#elif defined(__i386__) || defined(__alpha__) \
   || defined(__x86_64) || defined(__x86_64__) \
   || defined(__ia64) || defined(__ia64__) \
   || defined(_M_IX86) || defined(_M_IA64) \
   || defined(_M_ALPHA) || defined(_WIN64)
#  define CGAL_LITTLE_ENDIAN
#else
#  error Unknown endianness
#endif


// Symbolic constants to tailor inlining. Inlining Policy.
// =======================================================
#ifndef CGAL_MEDIUM_INLINE
#  define CGAL_MEDIUM_INLINE inline
#endif

#ifndef CGAL_LARGE_INLINE
#  define CGAL_LARGE_INLINE
#endif

#ifndef CGAL_HUGE_INLINE
#  define CGAL_HUGE_INLINE
#endif


//----------------------------------------------------------------------//
// SunPRO specific.
//----------------------------------------------------------------------//
#ifdef __SUNPRO_CC
#  include <iterator>
#  ifdef _RWSTD_NO_CLASS_PARTIAL_SPEC
#    warning "CGAL does not support SunPRO with the old Rogue Wave STL: use STLPort."
#  endif

// Sun CC has an issue with templates that means overloading
// Qualified_result_of does not work so well.
#  define CGAL_CFG_DONT_OVERLOAD_TOO_MUCH 1

#endif

#ifdef __SUNPRO_CC
// SunPRO 5.9 emits warnings "The variable tag has not yet been assigned a value"
// even for empty "tag" variables.  No way to write a config/testfile for this.
#  define CGAL_SUNPRO_INITIALIZE(C) C
#else
#  define CGAL_SUNPRO_INITIALIZE(C)
#endif


//-------------------------------------------------------------------//
// When the global min and max are no longer defined (as macros) 
// because of NOMINMAX flag definition, we define our own global 
// min/max functions to make the Microsoft headers compile. (afxtempl.h)
// Users that does not want the global min/max 
// should define CGAL_NOMINMAX
//-------------------------------------------------------------------//
#include <algorithm>
#if defined NOMINMAX && !defined CGAL_NOMINMAX
using std::min;
using std::max;
#endif

//-------------------------------------------------------------------//
// Is CORE usable ?
#ifdef CGAL_USE_GMP
#  define CGAL_USE_CORE CGAL_USE_GMP
#endif


//-------------------------------------------------------------------//
// Is Geomview usable ?
#if !defined(_MSC_VER) && !defined(__MINGW32__)
#  define CGAL_USE_GEOMVIEW
#endif


//-------------------------------------------------------------------//
// Compilers provide different macros to access the current function name
#ifdef _MSC_VER
#  define CGAL_PRETTY_FUNCTION __FUNCSIG__
#elif defined __GNUG__
#  define CGAL_PRETTY_FUNCTION __PRETTY_FUNCTION__
#else 
#  define CGAL_PRETTY_FUNCTION __func__
// with sunpro, this requires -features=extensions
#endif


// Macro to trigger deprecation warnings
#ifdef CGAL_NO_DEPRECATION_WARNINGS
#  define CGAL_DEPRECATED
#elif defined (__GNUC__) && (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 1))
#  define CGAL_DEPRECATED __attribute__((__deprecated__))
#elif defined (_MSC_VER) && (_MSC_VER > 1300)
#  define CGAL_DEPRECATED __declspec(deprecated)
#else
#  define CGAL_DEPRECATED
#endif


// Macro to specify a noreturn attribute.
#ifdef __GNUG__
#  define CGAL_NORETURN  __attribute__ ((__noreturn__))
#else
#  define CGAL_NORETURN
#endif


// If CGAL_HAS_THREADS is not defined, then CGAL code assumes
// it can do any thread-unsafe things (like using static variables).
#if !defined CGAL_HAS_THREADS && !defined CGAL_HAS_NO_THREADS
  #if defined BOOST_HAS_THREADS || defined _OPENMP
  #  define CGAL_HAS_THREADS
  #endif
#endif


CGAL_BEGIN_NAMESPACE

// Typedef for the type of NULL.
typedef const void * Nullptr_t;   // Anticipate C++0x's std::nullptr_t

CGAL_END_NAMESPACE

#endif // CGAL_CONFIG_H
