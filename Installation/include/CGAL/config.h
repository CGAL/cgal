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

#ifndef CGAL_CFG_NO_LONG_LONG
#  define CGAL_USE_LONG_LONG
#endif

#ifdef CGAL_CFG_NO_LONG_DOUBLE_IO
#include <iostream>
namespace std {
  template < typename _CharT, typename _Traits >
  inline basic_ostream<_CharT, _Traits> &
  operator<<(basic_ostream<_CharT, _Traits> & os, const long double &ld)
  {
      return os << (double) ld;
  }
}
#endif


#ifndef CGAL_CFG_TYPENAME_BEFORE_DEFAULT_ARGUMENT_BUG
#  define CGAL_TYPENAME_DEFAULT_ARG typename
#else
#  define CGAL_TYPENAME_DEFAULT_ARG
#endif


#ifndef CGAL_CFG_SUNPRO_RWSTD
#  define CGAL_reverse_iterator(T) std::reverse_iterator< T >
#else
#  define CGAL_reverse_iterator(T) std::reverse_iterator< T , \
                                   typename T::iterator_category , \
                                   typename T::value_type , \
                                   typename T::reference , \
                                   typename T::pointer , \
                                   typename T::difference_type >
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


#ifndef CGAL_USE_LEDA
#  define CGAL_USE_CGAL_WINDOW
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
//             include separate workaround files
//----------------------------------------------------------------------//

#if defined(__BORLANDC__) && __BORLANDC__ > 0x520
#  include <CGAL/Borland_fixes.h>
#elif defined(__SUNPRO_CC)
#  include <CGAL/Sun_fixes.h>
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
#if !defined(__BORLANDC__) && !defined(_MSC_VER) && !defined(__MINGW32__)
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


// SunPRO's STL is missing the std::vector constructor from an iterator range.
#ifdef CGAL_CFG_MISSING_TEMPLATE_VECTOR_CONSTRUCTORS_BUG
#include <vector>
#include <list>
#include <iterator>
#include <algorithm>
namespace CGAL { namespace CGALi {
template < typename Iterator >
inline
std::vector<typename std::iterator_traits<Iterator>::value_type>
make_vector(Iterator begin, Iterator end)
{
  std::vector<typename std::iterator_traits<Iterator>::value_type> v;
  std::copy(begin, end, std::back_inserter(v));
  return v;
}
template < typename Iterator >
inline
std::list<typename std::iterator_traits<Iterator>::value_type>
make_list(Iterator begin, Iterator end)
{
  std::list<typename std::iterator_traits<Iterator>::value_type> v;
  std::copy(begin, end, std::back_inserter(v));
  return v;
}
} }
#  define CGAL_make_vector(begin, end) (CGAL::CGALi::make_vector(begin, end))
#  define CGAL_make_list(begin, end) (CGAL::CGALi::make_list(begin, end))
#else
#  define CGAL_make_vector(begin, end) (begin, end)
#  define CGAL_make_list(begin, end) (begin, end)
#endif

#if (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 1))
#define CGAL_DEPRECATED  __attribute__((deprecated))
#elif _MSC_VER > 1300
#define CGAL_DEPRECATED __declspec(deprecated)
#else
#define CGAL_DEPRECATED
#endif



#endif // CGAL_CONFIG_H
