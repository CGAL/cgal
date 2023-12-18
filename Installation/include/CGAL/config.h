// Copyright (c) 1997-2013
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
//
// Author(s)     : Wieger Wesselink
//                 Michael Hoffmann <hoffmann@inf.ethz.ch>
//                 Sylvain Pion
//                 Laurent Rineau

#ifndef CGAL_CONFIG_H
#define CGAL_CONFIG_H

// CGAL is header-only by default since CGAL-5.0.
#if !defined(CGAL_HEADER_ONLY) && ! CGAL_NOT_HEADER_ONLY
#  define CGAL_HEADER_ONLY 1
#endif

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
// See also: https://kinddragon.github.io/vld/
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
//work around for moc bug : https://bugreports.qt.io/browse/QTBUG-80990
#if defined(CGAL_LINKED_WITH_TBB)
#undef CGAL_LINKED_WITH_TBB
#endif
#endif

// The following header file defines among other things  BOOST_PREVENT_MACRO_SUBSTITUTION
#include <boost/config.hpp>
#include <boost/version.hpp>

#include <CGAL/version.h>
#include <CGAL/version_checker.h>

//----------------------------------------------------------------------//
//  platform specific workaround flags (CGAL_CFG_...)
//----------------------------------------------------------------------//

#if CGAL_HEADER_ONLY
#  include <CGAL/Installation/internal/enable_third_party_libraries.h>
#else
#  include <CGAL/compiler_config.h>
#endif

#if BOOST_MSVC && CGAL_TEST_SUITE
#  include <CGAL/Testsuite/vc_debug_hook.h>
#endif

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

// Same for C++17
#if !(__cplusplus >= 201703L || _MSVC_LANG >= 201703L)
#error "CGAL requires C++ 17"
#endif
// Same for C++20
#if __cplusplus >= 202002L || _MSVC_LANG >= 202002L
#  define CGAL_CXX20 1
#endif


//----------------------------------------------------------------------//
//  As std::unary_function and std::binary_function are deprecated
//  we use internally equivalent class templates from
// <CGAL/functional.h>.
//----------------------------------------------------------------------//

#include <CGAL/functional.h>

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

#ifndef CGAL_CFG_TYPENAME_BEFORE_DEFAULT_ARGUMENT_BUG
#  define CGAL_TYPENAME_DEFAULT_ARG typename
#else
#  define CGAL_TYPENAME_DEFAULT_ARG
#endif

// Big endian or little endian machine.
// ====================================

#include <boost/predef.h>
#if BOOST_ENDIAN_BIG_BYTE
#  define CGAL_BIG_ENDIAN
#elif BOOST_ENDIAN_LITTLE_BYTE
#  define CGAL_LITTLE_ENDIAN
#endif

#if ! defined(CGAL_LITTLE_ENDIAN) && ! defined(CGAL_BIG_ENDIAN)
#  ifdef CGAL_DEFAULT_IS_LITTLE_ENDIAN
#    if CGAL_DEFAULT_IS_LITTLE_ENDIAN
#      define CGAL_LITTLE_ENDIAN
#    else
#      define CGAL_BIG_ENDIAN
#    endif
#  else
#    error Unknown endianness: Define CGAL_DEFAULT_IS_LITTLE_ENDIAN to 1 for little endian and to 0 for big endian.
#  endif
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
#    warning "For more information, see https://www.cgal.org/FAQ.html#mac_optimization_bug"
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
// See https://clang.llvm.org/docs/LanguageExtensions.html
// See also https://en.cppreference.com/w/cpp/experimental/feature_test
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
#if __has_cpp_attribute(maybe_unused)
#  define CGAL_UNUSED [[maybe_unused]]
#elif defined(__GNUG__) || __has_attribute(__unused__) // [[maybe_unused]] is C++17
#  define CGAL_UNUSED __attribute__ ((__unused__))
#else
#  define CGAL_UNUSED
#endif

// Macro to trigger deprecation warnings
#ifdef CGAL_NO_DEPRECATION_WARNINGS
#  define CGAL_DEPRECATED
#  define CGAL_DEPRECATED_MSG(msg)
#  define CGAL_DEPRECATED_UNUSED CGAL_UNUSED
#else
#  define CGAL_DEPRECATED [[deprecated]]
#  define CGAL_DEPRECATED_MSG(msg) [[deprecated(msg)]]
#  define CGAL_DEPRECATED_UNUSED [[deprecated]] CGAL_UNUSED
#endif

// Macro to specify a 'noreturn' attribute.
// (This macro existed in CGAL before we switched to C++11. Let's keep
// the macro defined for backward-compatibility. That cannot harm.)
#define CGAL_NORETURN  [[noreturn]]

// Macro to specify [[no_unique_address]] if supported
#if _MSC_VER >= 1929 && _MSVC_LANG >= 202002L
// see https://devblogs.microsoft.com/cppblog/msvc-cpp20-and-the-std-cpp20-switch/#c20-no_unique_address
#  define CGAL_NO_UNIQUE_ADDRESS [[msvc::no_unique_address]]
#elif __has_cpp_attribute(no_unique_address)
#  define CGAL_NO_UNIQUE_ADDRESS [[no_unique_address]]
#else
#  define CGAL_NO_UNIQUE_ADDRESS
#endif

// Macro CGAL_ASSUME and CGAL_UNREACHABLE
// Call a builtin of the compiler to pass a hint to the compiler
#if __has_builtin(__builtin_unreachable) || (CGAL_GCC_VERSION > 0 && !__STRICT_ANSI__)
// From g++ 4.5, there exists a __builtin_unreachable()
// Also in LLVM/clang
#  define CGAL_ASSUME(EX) if(!(EX)) { __builtin_unreachable(); }
#  define CGAL_UNREACHABLE() __builtin_unreachable()
#elif defined(_MSC_VER)
// MSVC has __assume
#  define CGAL_ASSUME(EX) __assume(EX)
#  define CGAL_UNREACHABLE() __assume(0)
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

#ifndef CGAL_HAS_THREADS
  namespace CGAL { inline bool is_currently_single_threaded(){ return true; } }
#elif __has_include(<sys/single_threaded.h>)
#  include <sys/single_threaded.h>
  namespace CGAL { inline bool is_currently_single_threaded(){ return ::__libc_single_threaded; } }
#else
  /* This is the conservative default */
  namespace CGAL { inline bool is_currently_single_threaded(){ return false; } }
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
#if defined(__clang__) || defined(BOOST_GCC)
#  define CGAL_PRAGMA_DIAG_PUSH _Pragma("GCC diagnostic push")
#  define CGAL_PRAGMA_DIAG_POP  _Pragma("GCC diagnostic pop")
#else
#  define CGAL_PRAGMA_DIAG_PUSH
#  define CGAL_PRAGMA_DIAG_POP
#endif

//
// Compatibility with CGAL-4.14.
#ifndef CGAL_NO_DEPRECATED_CODE
//
// That is temporary, and will be replaced by a namespace alias, as
// soon as we can remove cpp11::result_of, and <CGAL/atomic.h> and
// <CGAL/thread.h>.
//
#  include <iterator>
#  include <array>
#  include <utility>
#  include <type_traits>
#  include <unordered_set>
#  include <unordered_map>
#  include <functional>
#  include <thread>
#  include <chrono>
#  include <atomic>
//
namespace CGAL {
//
  namespace cpp11 {
    using std::next;
    using std::prev;
    using std::copy_n;
    using std::array;
    using std::function;
    using std::tuple;
    using std::make_tuple;
    using std::tie;
    using std::get;
    using std::tuple_size;
    using std::tuple_element;
    using std::is_enum;
    using std::unordered_set;
    using std::unordered_map;
    using std::atomic;
    using std::memory_order_relaxed;
    using std::memory_order_consume;
    using std::memory_order_acquire;
    using std::memory_order_release;
    using std::memory_order_acq_rel;
    using std::memory_order_seq_cst;
    using std::atomic_thread_fence;
    using std::thread;

  }
//
  namespace cpp0x = cpp11;
  using cpp11::array;
  using cpp11::copy_n;
} // end of the temporary compatibility with CGAL-4.14
#endif // CGAL_NO_DEPRECATED_CODE
namespace CGAL {

// Typedef for the type of nullptr.
typedef const void * Nullptr_t;   // Anticipate C++0x's std::nullptr_t
namespace cpp11{
#if CGAL_CXX20 || __cpp_lib_is_invocable>=201703L
    template<typename Signature> class result_of;
    template<typename F, typename... Args>
    class result_of<F(Args...)> : public std::invoke_result<F, Args...> { };
#else
    using std::result_of;
#endif
}//namespace cpp11
} //namespace CGAL

// The fallthrough attribute
// See for clang:
//   https://clang.llvm.org/docs/AttributeReference.html#statement-attributes
// See for gcc:
//   https://gcc.gnu.org/onlinedocs/gcc/Warning-Options.html
#if __cplusplus > 201402L && __has_cpp_attribute(fallthrough)
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
///     #pragma CGAL_WARNING("This line should trigger a warning")
///
/// @{
#ifdef BOOST_MSVC
#  define CGAL_WARNING(desc) message(__FILE__ "(" CGAL_STRINGIZE(__LINE__) ") : warning: " desc)
#else // not BOOST_MSVC
#  define CGAL_WARNING(desc) message( "warning: " desc)
#endif // not BOOST_MSVC
/// @}

/// Macro `CGAL_pragma_warning`.
/// @{
#ifdef BOOST_MSVC
#  define CGAL_pragma_warning(desc) __pragma(CGAL_WARNING(desc))
#else // not BOOST_MSVC
#  define CGAL_pragma_warning(desc) _Pragma(CGAL_STRINGIZE(CGAL_WARNING(desc)))
#endif // not BOOST_MSVC
/// @}
#include <CGAL/license/lgpl.h>

//----------------------------------------------------------------------//
//  Function to define data directory
//----------------------------------------------------------------------//
#include <cstdlib>
#include <string>
#include <fstream>
#include <iostream>

namespace CGAL {

// Returns filename prefixed by the directory of CGAL containing data.
// This directory is either defined in the environment variable CGAL_DATA_DIR,
// otherwise it is taken from the constant CGAL_DATA_DIR (defined in CMake),
// otherwise it is empty (and thus returns filename unmodified).
inline std::string data_file_path(const std::string& filename)
{
  const char* cgal_dir=nullptr;

#ifdef _MSC_VER
  char* cgal_dir_windows=nullptr;
  _dupenv_s( &cgal_dir_windows, nullptr, "CGAL_DATA_DIR");
  if (cgal_dir_windows!=nullptr)
  { cgal_dir=cgal_dir_windows; }
#else
  cgal_dir=getenv("CGAL_DATA_DIR");
#endif

#ifdef CGAL_DATA_DIR
 if (cgal_dir==nullptr)
 { cgal_dir=CGAL_DATA_DIR; }
#endif

 std::string cgal_dir_string;
 if (cgal_dir!=nullptr)
 { cgal_dir_string=std::string(cgal_dir); }

 std::string res=cgal_dir_string;
 if (!res.empty() && res.back()!='/')
 { res+=std::string("/"); }
 res+=filename;

 // Test if the file exists, write a warning otherwise
 std::ifstream f(res);
 if (!f)
 {
   std::cerr<<"[WARNING] file "<<res<<" does not exist or cannot be read "
            <<"(CGAL_DATA_DIR='"<<cgal_dir_string<<"')."<<std::endl;
 }

#ifdef _MSC_VER
 if (cgal_dir_windows!=nullptr)
 { free(cgal_dir_windows); }
#endif

 return res;
}

} // end namespace CGAL


#if BOOST_VERSION < 107900

// Workaround for an accidental enable if of Eigen::Matrix in the
// boost::multiprecision::cpp_int constructor for some versions of
// boost

namespace Eigen{
  template <class A, int B, int C, int D, int E, int F>
  class Matrix;
  template <class A, int B, class C>
  class Ref;

  template <class A, class B, int C>
  class Product;

  template<typename BinaryOp, typename Lhs, typename Rhs>  class CwiseBinaryOp;

}

namespace boost {
    namespace multiprecision {
        namespace detail {
            template <typename T>
            struct is_byte_container;


            template <class A, int B, int C, int D, int E, int F>
            struct is_byte_container< Eigen::Matrix<A, B, C, D, E, F>>
            {
                static const bool value = false;
            };

            template <class A, int B, class C>
            struct is_byte_container< Eigen::Ref<A, B, C>>
            {
                static const bool value = false;
            };

            template <class A, class B, int C>
            struct is_byte_container< Eigen::Product<A, B, C>>
            {
                static const bool value = false;
            };

            template <class A, class B, class C>
            struct is_byte_container< Eigen::CwiseBinaryOp<A, B, C>>
            {
                static const bool value = false;
            };

        }
    }
}

#endif // BOOST_VERSION < 107900

#endif // CGAL_CONFIG_H
