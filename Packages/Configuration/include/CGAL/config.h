// Copyright (c) 1997-2003  Utrecht University (The Netherlands),
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Wieger Wesselink 
//                 Michael Hoffmann <hoffmann@inf.ethz.ch>

#ifndef CGAL_CONFIG_H
#define CGAL_CONFIG_H

#define CGAL_VERSION 2.4-I-65
#define CGAL_VERSION_NR 1002004065


//----------------------------------------------------------------------//
//             STLport fix for MSVC
//----------------------------------------------------------------------//


#if defined( _MSC_VER) && (_MSC_VER <=1300)
#   if ! defined(__INTEL_COMPILER)
#     define CGAL_LIMITED_ITERATOR_TRAITS_SUPPORT 1
#   endif
#else
#   define CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(a)
#endif // _MSC_VER


//----------------------------------------------------------------------//
//             include platform specific workaround flags (CGAL_CFG_...)
//----------------------------------------------------------------------//

#include <CGAL/compiler_config.h>

//----------------------------------------------------------------------//
//             do some post processing for the flags
//----------------------------------------------------------------------//


// Used to depend on config macros.
#define CGAL_TYPENAME_MSVC_NULL typename
#define CGAL_TEMPLATE_NULL      template <>

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

#ifdef CGAL_CFG_VC7_PRIVATE_TYPE_BUG
#  define CGAL_VC7_BUG_PROTECTED protected:
#else
#  define CGAL_VC7_BUG_PROTECTED
#endif

#ifdef CGAL_CFG_MATCHING_BUG_2
#   define CGAL_MSVC_DUMMY_ARGUMENT , int dummy=1
#else
#   define CGAL_MSVC_DUMMY_ARGUMENT
#endif

#ifdef CGAL_CFG_NO_TEMPLATE_FRIEND_DISTINCTION
#  define CGAL_NULL_TMPL_ARGS
#else
#  define CGAL_NULL_TMPL_ARGS <>
#endif

#ifdef CGAL_CFG_NO_STDC_NAMESPACE
#  define CGAL_CLIB_STD
#else
#  define CGAL_CLIB_STD std
#endif

#ifndef CGAL_CFG_NO_LONG_LONG
#  define CGAL_USE_LONG_LONG
#endif

// FIXME: what is the problem with Sun 5.5? MATCHING_BUG_4 is not
// triggered, but there are runtime errors, e.g., with Distance_3,
// that do not appear when using the wrapper...
#if defined(CGAL_CFG_MATCHING_BUG_4) || \
  (defined(__sun) && defined(__SUNPRO_CC))
#  define CGAL_WRAP(K) Matching_bug_wrapper<K>
#  include <CGAL/Matching_bug_wrapper.h>
#else
#  define CGAL_WRAP(K) K
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
#if defined(__sun) && defined(__SUNPRO_CC)
#include <CGAL/Sun_fixes.h>
#endif

//--------------------------------------------------------------------//
// This addresses a bug in VC++ 7.0 that (re)defines min(a, b)
// and max(a, b) in windows.h and windef.h 
//-------------------------------------------------------------------//

#ifdef _MSC_VER
#  define NOMINMAX 1
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



//-----------------------------------------------------------------------//
// the MSVC 6.0 and 7.0 compilers cannot deal with function overloading
// very well, so we have to use specific templating here with the CGAL
// Polyhedron_3 type in its two different forms (one that is swallowed by
// MSVC6 and the other by MSVC 7.0). 
//----------------------------------------------------------------------//

#if defined(_MSC_VER) && ! defined(__INTEL_COMPILER) && (_MSC_VER < 1310)
#  define CGAL_CFG_FUNCTION_OVERLOAD_BUG
#endif

#endif // CGAL_CONFIG_H
