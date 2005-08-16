//  (C) Copyright John Maddock 2001 - 2003. 
//  (C) Copyright Darin Adler 2001 - 2002. 
//  (C) Copyright Peter Dimov 2001. 
//  (C) Copyright Aleksey Gurtovoy 2002. 
//  (C) Copyright David Abrahams 2002 - 2003. 
//  (C) Copyright Beman Dawes 2002 - 2003. 
//  Use, modification and distribution are subject to the 
//  Boost Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org for most recent version.

//  Microsoft Visual C++ compiler setup:

#define BOOST_MSVC _MSC_VER

// turn off the warnings before we #include anything
#pragma warning( disable : 4503 ) // warning: decorated name length exceeded

#if _MSC_VER < 1300  // 1200 == VC++ 6.0, 1201 == EVC4.2
#pragma warning( disable : 4786 ) // ident trunc to '255' chars in debug info
#  define BOOST_NO_DEPENDENT_TYPES_IN_TEMPLATE_VALUE_PARAMETERS
#  define BOOST_NO_VOID_RETURNS
#  define BOOST_NO_EXCEPTION_STD_NAMESPACE
   // disable min/max macro defines on vc6:
   //
#endif

#if (_MSC_VER <= 1300)  // 1300 == VC++ 7.0

#if !defined(_MSC_EXTENSIONS) && !defined(BOOST_NO_DEPENDENT_TYPES_IN_TEMPLATE_VALUE_PARAMETERS)      // VC7 bug with /Za
#  define BOOST_NO_DEPENDENT_TYPES_IN_TEMPLATE_VALUE_PARAMETERS
#endif

#  define BOOST_NO_EXPLICIT_FUNCTION_TEMPLATE_ARGUMENTS
#  define BOOST_NO_INCLASS_MEMBER_INITIALIZATION
#  define BOOST_NO_PRIVATE_IN_AGGREGATE
#  define BOOST_NO_ARGUMENT_DEPENDENT_LOOKUP
#  define BOOST_NO_INTEGRAL_INT64_T
#  define BOOST_NO_DEDUCED_TYPENAME
#  define BOOST_NO_USING_DECLARATION_OVERLOADS_FROM_TYPENAME_BASE

//    VC++ 6/7 has member templates but they have numerous problems including
//    cases of silent failure, so for safety we define:
#  define BOOST_NO_MEMBER_TEMPLATES
//    For VC++ experts wishing to attempt workarounds, we define:
#  define BOOST_MSVC6_MEMBER_TEMPLATES

#  define BOOST_NO_MEMBER_TEMPLATE_FRIENDS
#  define BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION
#  define BOOST_NO_CV_VOID_SPECIALIZATIONS
#  define BOOST_NO_FUNCTION_TEMPLATE_ORDERING
#  define BOOST_NO_USING_TEMPLATE
#  define BOOST_NO_SWPRINTF
#  define BOOST_NO_TEMPLATE_TEMPLATES
#  define BOOST_NO_SFINAE
#  define BOOST_NO_POINTER_TO_MEMBER_TEMPLATE_PARAMETERS
#  define BOOST_NO_IS_ABSTRACT
#  if (_MSC_VER > 1200)
#     define BOOST_NO_MEMBER_FUNCTION_SPECIALIZATIONS
#  endif

#endif

#if _MSC_VER < 1310 // 1310 == VC++ 7.1
#  define BOOST_NO_SWPRINTF
#endif

#if _MSC_VER <= 1400  // 1400 == VC++ 8.0
#  define BOOST_NO_MEMBER_TEMPLATE_FRIENDS
#endif

#ifndef _NATIVE_WCHAR_T_DEFINED
#  define BOOST_NO_INTRINSIC_WCHAR_T
#endif

#ifdef _WIN32_WCE
#  define BOOST_NO_THREADEX
#  define BOOST_NO_GETSYSTEMTIMEASFILETIME
#endif

//   
// check for exception handling support:   
#ifndef _CPPUNWIND   
#  define BOOST_NO_EXCEPTIONS   
#endif 

//
// __int64 support:
//
#if (_MSC_VER >= 1200)
#   define BOOST_HAS_MS_INT64
#endif
#if (_MSC_VER >= 1310) && defined(_MSC_EXTENSIONS)
#   define BOOST_HAS_LONG_LONG
#endif
//
// disable Win32 API's if compiler extentions are
// turned off:
//
#ifndef _MSC_EXTENSIONS
#  define BOOST_DISABLE_WIN32
#endif

//
// all versions support __declspec:
//
#define BOOST_HAS_DECLSPEC
//
// prefix and suffix headers:
//
#ifndef BOOST_ABI_PREFIX
#  define BOOST_ABI_PREFIX "boost/config/abi/msvc_prefix.hpp"
#endif
#ifndef BOOST_ABI_SUFFIX
#  define BOOST_ABI_SUFFIX "boost/config/abi/msvc_suffix.hpp"
#endif

# if _MSC_VER == 1200
#   define BOOST_COMPILER_VERSION 6.0
# elif _MSC_VER == 1300
#   define BOOST_COMPILER_VERSION 7.0
# elif _MSC_VER == 1310
#   define BOOST_COMPILER_VERSION 7.1
# elif _MSC_VER == 1400
#   define BOOST_COMPILER_VERSION 8.0
# else
#   define BOOST_COMPILER_VERSION _MSC_VER
# endif

#define BOOST_COMPILER "Microsoft Visual C++ version " BOOST_STRINGIZE(BOOST_COMPILER_VERSION)

//
// versions check:
// we don't support Visual C++ prior to version 6:
#if _MSC_VER < 1200
#error "Compiler not supported or configured - please reconfigure"
#endif
//
// last known and checked version is 1310:
#if (_MSC_VER > 1400)
#  if defined(BOOST_ASSERT_CONFIG)
#     error "Unknown compiler version - please run the configure tests and report the results"
#  else
#     pragma message("Unknown compiler version - please run the configure tests and report the results")
#  endif
#endif
