// Copyright (c) 1997-2013
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
// Author(s)     : Wieger Wesselink 
//                 Michael Hoffmann <hoffmann@inf.ethz.ch>
//                 Sylvain Pion
//                 Laurent Rineau

#ifndef CGAL_CONFIG_H
#define CGAL_CONFIG_H

#ifdef CGAL_HEADER_ONLY
#  define CGAL_NO_AUTOLINK 1
#endif

// Workaround for a bug in Boost, that checks WIN64 instead of _WIN64
//   https://svn.boost.org/trac/boost/ticket/5519
#if defined(_WIN64) && ! defined(WIN64)
#  define WIN64
#endif

#ifdef CGAL_INCLUDE_WINDOWS_DOT_H
// Mimic users including this file which defines min max macros
// and other names leading to name clashes
#include <windows.h>
#endif

#if defined(CGAL_TEST_SUITE) && defined(NDEBUG)
#  error The test-suite needs no NDEBUG defined
#endif // CGAL_TEST_SUITE and NDEBUG

// See [[Small features/Visual_Leak_Detector]] in CGAL developers wiki
// See also: http://vld.codeplex.com/
#if defined(CGAL_ENABLE_VLD)
#  include <vld.h>
#endif // CGAL_ENABLE_VLD

// Workaround to the following bug:
// https://bugreports.qt-project.org/browse/QTBUG-22829
#ifdef Q_MOC_RUN
// When Qt moc runs on CGAL files, do not process
// <boost/type_traits/detail/has_binary_operator.hpp>
#  define BOOST_TT_HAS_OPERATOR_HPP_INCLUDED
#  define BOOST_TT_HAS_BIT_AND_HPP_INCLUDED
#  define BOOST_TT_HAS_BIT_AND_ASSIGN_HPP_INCLUDED
#  define BOOST_TT_HAS_BIT_OR_HPP_INCLUDED
#  define BOOST_TT_HAS_BIT_OR_ASSIGN_HPP_INCLUDED
#  define BOOST_TT_HAS_BIT_XOR_HPP_INCLUDED
#  define BOOST_TT_HAS_BIT_XOR_ASSIGN_HPP_INCLUDED
#  define BOOST_TT_HAS_DIVIDES_HPP_INCLUDED
#  define BOOST_TT_HAS_DIVIDES_ASSIGN_HPP_INCLUDED
#  define BOOST_TT_HAS_EQUAL_TO_HPP_INCLUDED
#  define BOOST_TT_HAS_GREATER_HPP_INCLUDED
#  define BOOST_TT_HAS_GREATER_EQUAL_HPP_INCLUDED
#  define BOOST_TT_HAS_LEFT_SHIFT_HPP_INCLUDED
#  define BOOST_TT_HAS_LEFT_SHIFT_ASSIGN_HPP_INCLUDED
#  define BOOST_TT_HAS_LESS_HPP_INCLUDED
#  define BOOST_TT_HAS_LESS_EQUAL_HPP_INCLUDED
#  define BOOST_TT_HAS_LOGICAL_AND_HPP_INCLUDED
#  define BOOST_TT_HAS_LOGICAL_OR_HPP_INCLUDED
#  define BOOST_TT_HAS_MINUS_HPP_INCLUDED
#  define BOOST_TT_HAS_MINUS_ASSIGN_HPP_INCLUDED
#  define BOOST_TT_HAS_MODULUS_HPP_INCLUDED
#  define BOOST_TT_HAS_MODULUS_ASSIGN_HPP_INCLUDED
#  define BOOST_TT_HAS_MULTIPLIES_HPP_INCLUDED
#  define BOOST_TT_HAS_MULTIPLIES_ASSIGN_HPP_INCLUDED
#  define BOOST_TT_HAS_NOT_EQUAL_TO_HPP_INCLUDED
#  define BOOST_TT_HAS_PLUS_HPP_INCLUDED
#  define BOOST_TT_HAS_PLUS_ASSIGN_HPP_INCLUDED
#  define BOOST_TT_HAS_RIGHT_SHIFT_HPP_INCLUDED
#  define BOOST_TT_HAS_RIGHT_SHIFT_ASSIGN_HPP_INCLUDED
// do not include <boost/random.hpp> either
// it includes <boost/type_traits/has_binary_operator.hpp>
#  define BOOST_RANDOM_HPP
// <boost/type_traits/detail/has_prefix_operator.hpp> fails as well
#  define BOOST_TT_HAS_COMPLEMENT_HPP_INCLUDED
#  define BOOST_TT_HAS_DEREFERENCE_HPP_INCLUDED
#  define BOOST_TT_HAS_LOGICAL_NOT_HPP_INCLUDED
#  define BOOST_TT_HAS_NEGATE_HPP_INCLUDED
#  define BOOST_TT_HAS_PRE_DECREMENT_HPP_INCLUDED
#  define BOOST_TT_HAS_PRE_INCREMENT_HPP_INCLUDED
#  define BOOST_TT_HAS_UNARY_MINUS_HPP_INCLUDED
#  define BOOST_TT_HAS_UNARY_PLUS_HPP_INCLUDED
// <boost/type_traits/detail/has_postfix_operator.hpp> fails as well
#  define BOOST_TT_HAS_POST_DECREMENT_HPP_INCLUDED
#  define BOOST_TT_HAS_POST_INCREMENT_HPP_INCLUDED
#endif

// The following header file defines among other things  BOOST_PREVENT_MACRO_SUBSTITUTION 
#include <boost/config.hpp>
#include <boost/version.hpp>

// bug-fix for g++-5.x and Boost.Config<1.57
//    https://svn.boost.org/trac/boost/ticket/10500
#if BOOST_VERSION < 105700 && BOOST_GCC < 60000 && \
  ! defined(__GXX_EXPERIMENTAL_CXX0X__) && defined(BOOST_HAS_VARIADIC_TMPL)
#  undef BOOST_HAS_VARIADIC_TMPL
#  define BOOST_NO_CXX11_VARIADIC_TEMPLATES
#endif

// workaround for the bug https://svn.boost.org/trac10/ticket/12534
// That bug was introduced in Boost 1.62 and fixed in 1.63.
#if BOOST_VERSION >= 106200 && BOOSTS_VERSION < 106300
#  include <boost/container/flat_map.hpp>
#endif

// Hack: Boost<1.55 does not detect correctly the C++11 features of ICC.
// We declare by hand two features that we need (variadic templates and
// rvalue references).
#if defined(__INTEL_COMPILER) && defined(__GXX_EXPERIMENTAL_CXX0X__)
#  undef BOOST_NO_VARIADIC_TEMPLATES
#  undef BOOST_NO_CXX11_VARIADIC_TEMPLATES
#  undef BOOST_NO_RVALUE_REFERENCES
#  undef BOOST_NO_CXX11_RVALUE_REFERENCES
#endif

#include <CGAL/version.h>

//----------------------------------------------------------------------//
//  platform specific workaround flags (CGAL_CFG_...)
//----------------------------------------------------------------------//

#include <CGAL/compiler_config.h>

//----------------------------------------------------------------------//
//  Support for DLL on Windows (CGAL_EXPORT macro)
//----------------------------------------------------------------------//

#include <CGAL/export/CGAL.h>

//----------------------------------------------------------------------//
//  Use an implementation of fabs with sse2 on Windows
//----------------------------------------------------------------------//

#if (_M_IX86_FP >= 2) || defined(_M_X64)
#define CGAL_USE_SSE2_FABS
#endif

//----------------------------------------------------------------------//
//  Detect features at compile-time. Some macros have only been
//  introduced as of Boost 1.40. In that case, we simply say that the
//  feature is not available, even if that is wrong.
//  ----------------------------------------------------------------------//

#if defined(BOOST_NO_CXX11_RANGE_BASED_FOR) || BOOST_VERSION < 105000
#define CGAL_CFG_NO_CPP0X_RANGE_BASED_FOR 1
#endif
#if defined(BOOST_NO_0X_HDR_ARRAY) || \
    defined(BOOST_NO_CXX11_HDR_ARRAY) || BOOST_VERSION < 104000
#define CGAL_CFG_NO_CPP0X_ARRAY 1
#endif
#if defined(BOOST_NO_0X_HDR_UNORDERED_SET) || \
    defined(BOOST_NO_0X_HDR_UNORDERED_MAP) || \
    defined(BOOST_NO_CXX11_HDR_UNORDERED_SET) || \
    defined(BOOST_NO_CXX11_HDR_UNORDERED_MAP)
#define CGAL_CFG_NO_CPP0X_UNORDERED 1
#endif
#if defined(BOOST_NO_DECLTYPE) || \
    defined(BOOST_NO_CXX11_DECLTYPE) || (BOOST_VERSION < 103600)
#define CGAL_CFG_NO_CPP0X_DECLTYPE 1
#endif
#if defined(BOOST_NO_DELETED_FUNCTIONS) || \
    defined(BOOST_NO_DEFAULTED_FUNCTIONS) || \
    defined(BOOST_NO_CXX11_DELETED_FUNCTIONS) || \
    defined(BOOST_NO_CXX11_DEFAULTED_FUNCTIONS) || (BOOST_VERSION < 103600)
#define CGAL_CFG_NO_CPP0X_DELETED_AND_DEFAULT_FUNCTIONS 1
#endif
#if defined(BOOST_NO_FUNCTION_TEMPLATE_DEFAULT_ARGS) || \
    defined(BOOST_NO_CXX11_FUNCTION_TEMPLATE_DEFAULT_ARGS) || \
    (BOOST_VERSION < 104100)
#define CGAL_CFG_NO_CPP0X_DEFAULT_TEMPLATE_ARGUMENTS_FOR_FUNCTION_TEMPLATES 1
#endif
#if defined(BOOST_NO_INITIALIZER_LISTS) || \
    defined(BOOST_NO_CXX11_HDR_INITIALIZER_LIST) || (BOOST_VERSION < 103900)
#define CGAL_CFG_NO_CPP0X_INITIALIZER_LISTS 1
#endif
#if defined(BOOST_MSVC)
#define CGAL_CFG_NO_CPP0X_ISFINITE 1 // used in <CGAL/CORE/Filter.h>
#endif
#if defined(BOOST_NO_LONG_LONG) || (BOOST_VERSION < 103600)
#define CGAL_CFG_NO_CPP0X_LONG_LONG 1
#endif
#if defined(BOOST_NO_LAMBDAS) || \
    defined(BOOST_NO_CXX11_LAMBDAS) || BOOST_VERSION < 104000
#define CGAL_CFG_NO_CPP0X_LAMBDAS 1
#endif
#if defined(BOOST_NO_RVALUE_REFERENCES) || \
    defined(BOOST_NO_CXX11_RVALUE_REFERENCES) || (BOOST_VERSION < 103600)
#define CGAL_CFG_NO_CPP0X_RVALUE_REFERENCE 1
#endif
#if defined(BOOST_NO_STATIC_ASSERT) || \
    defined(BOOST_NO_CXX11_STATIC_ASSERT) || (BOOST_VERSION < 103600)
#define CGAL_CFG_NO_CPP0X_STATIC_ASSERT 1
#endif
#if defined(BOOST_NO_0X_HDR_TUPLE) || \
    defined(BOOST_NO_CXX11_HDR_TUPLE) || (BOOST_VERSION < 104000)
#define CGAL_CFG_NO_CPP0X_TUPLE 1
#endif
#if defined(BOOST_NO_VARIADIC_TEMPLATES) || \
    defined(BOOST_NO_CXX11_VARIADIC_TEMPLATES) || (BOOST_VERSION < 103600)
#define CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES 1
#endif
// never use TR1
#define CGAL_CFG_NO_TR1_ARRAY 1
// never use TR1
#define CGAL_CFG_NO_TR1_TUPLE 1
#if !defined(__GNUC__) || defined(__INTEL_COMPILER)
#define CGAL_CFG_NO_STATEMENT_EXPRESSIONS 1
#endif
#if defined(BOOST_NO_CXX11_UNIFIED_INITIALIZATION_SYNTAX) || (BOOST_VERSION < 105100) || _MSC_VER==1800
#define CGAL_CFG_NO_CPP0X_UNIFIED_INITIALIZATION_SYNTAX
#endif
#if __cplusplus < 201103L && !(_MSC_VER >= 1600)
#define CGAL_CFG_NO_CPP0X_COPY_N 1
#define CGAL_CFG_NO_CPP0X_NEXT_PREV 1
#endif
#if defined(BOOST_NO_EXPLICIT_CONVERSION_OPERATIONS) \
    || defined(BOOST_NO_EXPLICIT_CONVERSION_OPERATORS) \
    || defined(BOOST_NO_CXX11_EXPLICIT_CONVERSION_OPERATORS) \
    || (BOOST_VERSION < 103600)
#define CGAL_CFG_NO_CPP0X_EXPLICIT_CONVERSION_OPERATORS 1
#endif

// Some random list to let us write C++11 without thinking about
// each feature we are using.
#if __cplusplus >= 201103L && \
    !defined CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES && \
    !defined CGAL_CFG_NO_CPP0X_RVALUE_REFERENCE && \
    !defined CGAL_CFG_NO_CPP0X_EXPLICIT_CONVERSION_OPERATORS && \
    !defined CGAL_CFG_NO_CPP0X_TUPLE && \
    !defined CGAL_CFG_NO_CPP0X_UNIFIED_INITIALIZATION_SYNTAX && \
    !defined CGAL_CFG_NO_CPP0X_STATIC_ASSERT && \
    !defined CGAL_CFG_NO_CPP0X_DECLTYPE && \
    !defined CGAL_CFG_NO_CPP0X_DELETED_AND_DEFAULT_FUNCTIONS && \
    !defined CGAL_CFG_NO_CPP0X_DEFAULT_TEMPLATE_ARGUMENTS_FOR_FUNCTION_TEMPLATES
#define CGAL_CXX11
#endif

#if defined(BOOST_NO_CXX11_HDR_FUNCTIONAL) || BOOST_VERSION < 105000
#define CGAL_CFG_NO_STD_HASH 1
#endif

//----------------------------------------------------------------------//
//  auto-link the CGAL library on platforms that support it
//----------------------------------------------------------------------//

#include <CGAL/auto_link/CGAL.h>

//----------------------------------------------------------------------//
//  do some post processing for the flags
//----------------------------------------------------------------------//

#ifdef CGAL_CFG_NO_STL
#  error "This compiler does not have a working STL"
#endif

// This macro computes the version number from an x.y.z release number.
// It only works for public releases.
#define CGAL_VERSION_NUMBER(x,y,z) (1000001 + 10000*x + 100*y + 10*z) * 1000

#ifndef CGAL_NO_DEPRECATED_CODE
#define CGAL_BEGIN_NAMESPACE  namespace CGAL { 
#define CGAL_END_NAMESPACE }
#endif


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

#if (BOOST_VERSION >= 105500)
#  include <boost/predef.h>
#  if BOOST_ENDIAN_BIG_BYTE
#    define CGAL_BIG_ENDIAN
#  elif BOOST_ENDIAN_LITTLE_BYTE
#    define CGAL_LITTLE_ENDIAN
#  else
#    error Unknown endianness
#  endif
#elif defined (__GLIBC__)
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
#    error "CGAL does not support SunPRO with the old Rogue Wave STL: use STLPort."
#  endif
#endif

#ifdef __SUNPRO_CC
// SunPRO 5.9 emits warnings "The variable tag has not yet been assigned a value"
// even for empty "tag" variables.  No way to write a config/testfile for this.
#  define CGAL_SUNPRO_INITIALIZE(C) C
#else
#  define CGAL_SUNPRO_INITIALIZE(C)
#endif

//----------------------------------------------------------------------//
// MacOSX specific.
//----------------------------------------------------------------------//

#ifdef __APPLE__
#  if defined(__GNUG__) && (__GNUG__ == 4) && (__GNUC_MINOR__ == 0) \
   && defined(__OPTIMIZE__) && !defined(CGAL_NO_WARNING_FOR_MACOSX_GCC_4_0_BUG)
#    warning "Your configuration may exhibit run-time errors in CGAL code"
#    warning "This appears with g++ 4.0 on MacOSX when optimizing"
#    warning "You can disable this warning using -DCGAL_NO_WARNING_FOR_MACOSX_GCC_4_0_BUG"
#    warning "For more information, see http://www.cgal.org/FAQ.html#mac_optimization_bug"
#  endif
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

// Macro to detect GCC versions.
// It evaluates to 0 if the compiler is not GCC. Be careful that the Intel
// compilers on Linux, and the LLVM/clang compiler both define GCC version
// macros.
#define CGAL_GCC_VERSION (__GNUC__ * 10000       \
                          + __GNUC_MINOR__ * 100 \
                          + __GNUC_PATCHLEVEL__)

// Macros to detect features of clang. We define them for the other
// compilers.
// See http://clang.llvm.org/docs/LanguageExtensions.html
// See also http://en.cppreference.com/w/cpp/experimental/feature_test
#ifndef __has_feature
  #define __has_feature(x) 0  // Compatibility with non-clang compilers.
#endif
#ifndef __has_include
  #define __has_include(x) 0  // Compatibility with non-clang compilers.
#endif
#ifndef __has_extension
  #define __has_extension __has_feature // Compatibility with pre-3.0 compilers.
#endif
#ifndef __has_builtin
  #define __has_builtin(x) 0  // Compatibility with non-clang compilers.
#endif
#ifndef __has_attribute
  #define __has_attribute(x) 0  // Compatibility with non-clang compilers.
#endif
#ifndef __has_cpp_attribute
  #define __has_cpp_attribute(x) 0  // Compatibility with non-supporting compilers.
#endif
#ifndef __has_warning
  #define __has_warning(x) 0  // Compatibility with non-clang compilers.
#endif

// Macro to specify a 'unused' attribute.
#if defined(__GNUG__) || __has_attribute(__unused__)
#  define CGAL_UNUSED  __attribute__ ((__unused__))
#else
#  define CGAL_UNUSED
#endif

// Macro to trigger deprecation warnings
#ifdef CGAL_NO_DEPRECATION_WARNINGS
#  define CGAL_DEPRECATED
#  define CGAL_DEPRECATED_UNUSED CGAL_UNUSED
#elif defined(__GNUC__) || __has_attribute(__deprecated__)
#  define CGAL_DEPRECATED __attribute__((__deprecated__))
#if __has_attribute(__unused__)
#  define CGAL_DEPRECATED_UNUSED __attribute__((__deprecated__, __unused__))
#else
#  define CGAL_DEPRECATED_UNUSED __attribute__((__deprecated__))
#endif
#elif defined (_MSC_VER) && (_MSC_VER > 1300)
#  define CGAL_DEPRECATED __declspec(deprecated)
#  define CGAL_DEPRECATED_UNUSED __declspec(deprecated)
#else
#  define CGAL_DEPRECATED
#  define CGAL_DEPRECATED_UNUSED
#endif

// Macro to trigger deprecation warnings with a custom message
#ifdef CGAL_NO_DEPRECATION_WARNINGS
#  define CGAL_DEPRECATED_MSG(msg)
#elif defined(__GNUC__) || __has_attribute(__deprecated__)
#  if BOOST_GCC >= 40500 || __has_attribute(__deprecated__)
#    define CGAL_DEPRECATED_MSG(msg) __attribute__ ((deprecated(msg)))
#  else
#    define CGAL_DEPRECATED_MSG(msg) __attribute__ ((deprecated))
#  endif
#elif defined (_MSC_VER) && (_MSC_VER > 1300)
#  define CGAL_DEPRECATED_MSG(msg) __declspec(deprecated(msg))
#elif defined(__clang__)
#  define CGAL_DEPRECATED_MSG(msg) __attribute__ ((deprecated(msg)))
#else
#  define CGAL_DEPRECATED_MSG(msg)
#endif

// Macro to specify a 'noreturn' attribute.
#if defined(__GNUG__) || __has_attribute(__noreturn__)
#  define CGAL_NORETURN  __attribute__ ((__noreturn__))
#else
#  define CGAL_NORETURN
#endif

// Macro CGAL_ASSUME
// Call a builtin of the compiler to pass a hint to the compiler
#if __has_builtin(__builtin_unreachable) || (CGAL_GCC_VERSION >= 40500 && !__STRICT_ANSI__)
// From g++ 4.5, there exists a __builtin_unreachable()
// Also in LLVM/clang
#  define CGAL_ASSUME(EX) if(!(EX)) { __builtin_unreachable(); }
#elif defined(_MSC_VER)
// MSVC has __assume
#  define CGAL_ASSUME(EX) __assume(EX)
#endif
// If CGAL_ASSUME is not defined, then CGAL_assume and CGAL_assume_code are
// defined differently, in <CGAL/assertions.h>

// If CGAL_HAS_THREADS is not defined, then CGAL code assumes
// it can do any thread-unsafe things (like using static variables).
#if !defined CGAL_HAS_THREADS && !defined CGAL_HAS_NO_THREADS
#  if defined BOOST_HAS_THREADS || defined _OPENMP
#    define CGAL_HAS_THREADS
#  endif
#endif

#if __has_feature(cxx_thread_local) || \
    ( (__GNUC__ * 100 + __GNUC_MINOR__) >= 408 && __cplusplus >= 201103L ) || \
    ( _MSC_VER >= 1900 )
// see also Installation/config/support/CGAL_test_cpp_version.cpp
#define CGAL_CAN_USE_CXX11_THREAD_LOCAL
#endif

#if ( BOOST_VERSION >= 105000 && ! defined(BOOST_NO_CXX11_HDR_MUTEX) ) || \
    (__has_include(<mutex>) && __cplusplus >= 201103L ) | \
    ( (__GNUC__ * 100 + __GNUC_MINOR__) >= 408 && __cplusplus >= 201103L ) || \
    ( _MSC_VER >= 1700 )
#define CGAL_CAN_USE_CXX11_MUTEX
#endif

#if ( BOOST_VERSION >= 105600 && ! defined(BOOST_NO_CXX11_HDR_ATOMIC) ) || \
    (__has_include(<atomic>) && __cplusplus >= 201103L ) || \
    ( (__GNUC__ * 100 + __GNUC_MINOR__) >= 408 && __cplusplus >= 201103L ) || \
    ( _MSC_VER >= 1700 )
#define CGAL_CAN_USE_CXX11_ATOMIC
#endif

// Support for LEDA with threads
//   Not that, if CGAL_HAS_THREADS is defined, and you want to use LEDA,
//   you must link with a version of LEDA libraries that support threads.
#if defined(CGAL_HAS_THREADS) && CGAL_USE_LEDA
#  define LEDA_MULTI_THREAD 1
#endif
// Support for LEDA_numbers on Windows
#define LEDA_NUMBERS_DLL 1

// Helper macros to disable macros
#if defined(__clang__) || (BOOST_GCC >= 40600)
#  define CGAL_PRAGMA_DIAG_PUSH _Pragma("GCC diagnostic push")
#  define CGAL_PRAGMA_DIAG_POP  _Pragma("GCC diagnostic pop")
#else
#  define CGAL_PRAGMA_DIAG_PUSH
#  define CGAL_PRAGMA_DIAG_POP
#endif

namespace CGAL {

// Typedef for the type of NULL.
typedef const void * Nullptr_t;   // Anticipate C++0x's std::nullptr_t

} //namespace CGAL

//Support for c++11 noexcept
#if BOOST_VERSION > 104600 && !defined(BOOST_NO_CXX11_NOEXCEPT) && !defined(BOOST_NO_NOEXCEPT)
#define CGAL_NOEXCEPT(x) noexcept(x)
#else
#define CGAL_NOEXCEPT(x)
#endif

// The fallthrough attribute
// See for clang:
//   http://clang.llvm.org/docs/AttributeReference.html#statement-attributes
// See for gcc:
//   https://gcc.gnu.org/onlinedocs/gcc/Warning-Options.html
#if __has_cpp_attribute(fallthrough)
#  define CGAL_FALLTHROUGH [[fallthrough]]
#elif __has_cpp_attribute(gnu::fallthrough)
#  define CGAL_FALLTHROUGH [[gnu::fallthrough]]
#elif __has_cpp_attribute(clang::fallthrough)
#  define CGAL_FALLTHROUGH [[clang::fallthrough]]
#elif __has_attribute(fallthrough) && ! __clang__
#  define CGAL_FALLTHROUGH __attribute__ ((fallthrough))
#else
#  define CGAL_FALLTHROUGH while(false){}
#endif

// https://svn.boost.org/trac/boost/ticket/2839
#if defined(BOOST_MSVC) && BOOST_VERSION < 105600
#define CGAL_CFG_BOOST_VARIANT_SWAP_BUG 1
#endif

#ifndef CGAL_NO_ASSERTIONS
#  define CGAL_NO_ASSERTIONS_BOOL false
#else
#  define CGAL_NO_ASSERTIONS_BOOL true
#endif

#if defined( __INTEL_COMPILER)
#define CGAL_ADDITIONAL_VARIANT_FOR_ICL ,int
#else
#define CGAL_ADDITIONAL_VARIANT_FOR_ICL
#endif

#if !defined CGAL_EIGEN3_ENABLED && \
    !defined CGAL_EIGEN3_DISABLED && \
    __has_include(<Eigen/Jacobi>)
#  define CGAL_EIGEN3_ENABLED 1
#endif

#define CGAL_STRINGIZE_HELPER(x) #x
#define CGAL_STRINGIZE(x) CGAL_STRINGIZE_HELPER(x)

/// Macro `CGAL_WARNING`.
/// Must be used with `#pragma`, this way:
///
///     #pragma CGAL_WARNING(This line should trigger a warning)
///
/// @{
#ifdef BOOST_MSVC
#  define CGAL_WARNING(desc) message(__FILE__ "(" CGAL_STRINGIZE(__LINE__) ") : warning: " desc)
#else // not BOOST_MSVC
#  define CGAL_WARNING(desc) message( "warning: " desc)
#endif // not BOOST_MSVC
/// @}

/// Macro `CGAL_pragma_warning`.
#ifdef BOOST_MSVC
#  define CGAL_pragma_warning(desc) __pragma(CGAL_WARNING(desc))
#else // not BOOST_MSVC
#  define CGAL_pragma_warning(desc) _Pragma(CGAL_STRINGIZE(CGAL_WARNING(desc)))
#endif // not BOOST_MSVC
/// @}
#include <CGAL/license/lgpl.h>


#endif // CGAL_CONFIG_H
