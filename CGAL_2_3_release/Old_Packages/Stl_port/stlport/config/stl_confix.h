/*
 * Copyright (c) 1999 
 * Boris Fomitchev
 *
 * This material is provided "as is", with absolutely no warranty expressed
 * or implied. Any use is at your own risk.
 *
 * Permission to use or copy this software for any purpose is hereby granted 
 * without fee, provided the above notices are retained on all copies.
 * Permission to modify the code and to distribute modified code is granted,
 * provided the above notices are retained, and a notice that the code was
 * modified is included with the above copyright notice.
 *
 */

/*
 * STLport configuration file
 * It is internal STLport header - DO NOT include it directly
 * Purpose of this file : to define STLport settings that depend on
 * compiler flags or can be otherwise missed by "configure"
 *
 */

#ifndef __STL_CONFIX_H
# define __STL_CONFIX_H

// some fixes to configuration, either "configure"d or
// hardcoded. This also includes modifications
// of STLport switches depending on compiler flags 

// Microsoft Visual C++ 4.0, 4.1, 4.2, 5.0
// distinguish real MSC
#  if defined(_MSC_VER) && !defined(__MWERKS__) && !defined (__BORLANDC__)
#   define __STL_MSVC _MSC_VER
#  endif


// sunpro
# if defined (__SUNPRO_CC)
#  define __STL_NATIVE_INCLUDE_PATH ../CC 

// #  ifndef __STL_HAS_NO_NEW_C_HEADERS
// #   define __STL_HAS_NO_NEW_C_HEADERS
// #  endif

// bit-field problem ?
// #  undef __STL_NEED_MUTABLE
// #  define __STL_NEED_MUTABLE
#  if (__SUNPRO_CC >=0x500)
// #   define __STL_NATIVE_HEADER(header) <../CC/##header##>
#   define __STL_NATIVE_C_HEADER(header) <../CC/##header##.SUNWCCh>
#   define __STL_NATIVE_CPP_C_HEADER(header) <../CC/##header##.SUNWCCh>
#   if !defined (__STLPORT_NEW_IOSTREAMS)
#    undef  __STL_USE_OWN_NAMESPACE
#    define __STL_USE_OWN_NAMESPACE 1
#   endif
#  else
#   define wint_t __wint_t 
// famous CC 4.2 bug
#   define __STL_INLINE_STRING_LITERAL_BUG 1
// /usr/include
#   define __STL_NATIVE_C_INCLUDE_PATH ../include
#  endif
# endif

// msvc

# if defined(__STL_MSVC)

// fix possible "configure" bugs for VC++ and define things
// not settable with "configure"

#  define __STL_NATIVE_HEADER(x) <../include/##x>
#  define __STL_NATIVE_C_HEADER(x) <../include/##x>
#  define __STL_NATIVE_CPP_C_HEADER(x) <../include/##x>

#  define __STL_EXPORT_DECLSPEC __declspec(dllexport)
#  define __STL_IMPORT_DECLSPEC __declspec(dllimport)
#  define __STL_IMPORT_TEMPLATE_KEYWORD  extern
#  if (__STL_MSVC > 0x420)
#   define __STL_CLASS_EXPORT_DECLSPEC __declspec(dllexport)
#   define __STL_CLASS_IMPORT_DECLSPEC __declspec(dllimport)
#  endif

#  define __STL_NO_TYPENAME_IN_TEMPLATE_HEADER
// fails to properly resolve call to sin() from within sin()
#  define __STL_SAME_FUNCTION_NAME_RESOLUTION_BUG

#if __STL_MSVC > 1000

#pragma warning ( disable : 4355)	// ignore warning "this used in base member initializer list.
// #pragma warning ( disable : 4251 )	// ignore template classes being exported in .dll's
#pragma warning ( disable : 4284 )	// ignore arrow operator warning
#pragma warning ( disable : 4231 )      // extern before spesialization

#endif

# if (__STL_MSVC < 1100)

#   undef  __STL_DEFAULT_TYPE_PARAM
#   undef  __STL_HAS_NO_NAMESPACES
#   undef  __STL_NO_AT_MEMBER_FUNCTION
#   undef  __STL_NO_MEMBER_TEMPLATES
#   undef  __STL_NO_MEMBER_TEMPLATE_CLASSES
#   define  __STL_HAS_NO_NAMESPACES 1
#   define  __STL_NO_AT_MEMBER_FUNCTION 1
#   define  __STL_NO_MEMBER_TEMPLATES
#   define  __STL_NO_MEMBER_TEMPLATE_CLASSES

# elif (__STL_MSVC > 1100)
    typedef char __stl_char;
#   define __STL_DEFAULTCHAR __stl_char
# endif /* (__STL_MSVC < 1100 ) */

# define __STL_MSVC50_COMPATIBILITY   1

# endif /* __STL_MSVC */

// list of compilers to use SGI string
# if (defined (__sgi) || defined (__GNUC__)) 
# if !defined (__STL_USE_SGI_STRING)
#  define __STL_USE_SGI_STRING         1
# endif
#  undef  __STL_USE_OWN_NAMESPACE 
# endif

# if defined(__hpux) && !defined(__GNUC__)
#  define __STL_USING_BASE_MEMBER 
// HP aCC with +noeh
#  ifdef __HPACC_NOEH
#   define __STL_HAS_NO_EXCEPTIONS 1
#  endif
# endif

# if defined (__GNUC__)

// strict ANSI prohibits "long long" ( gcc)
#  if defined ( __STRICT_ANSI__ )
#    undef __STL_LONG_LONG
#  endif

#   ifndef __EXCEPTIONS
#     undef  __STL_HAS_NO_EXCEPTIONS
#     define __STL_HAS_NO_EXCEPTIONS  1
#   endif

#if (__GNUC_MINOR__ < 8)
# define __STL_NATIVE_INCLUDE_PATH ../g++-include
#else
# define __STL_NATIVE_INCLUDE_PATH ../g++
#endif

# define __STL_NATIVE_CPP_C_INCLUDE_PATH __STL_NATIVE_INCLUDE_PATH
# define __STL_NATIVE_C_INCLUDE_PATH ../include

#  ifndef __HONOR_STD
#   define __STL_VENDOR_GLOBAL_EXCEPT_STD 1
#  endif

# endif /* __GNUC */

# if defined(__MINGW32__)
#   define __STL_NO_DRAND48
# endif

#if defined (__APOGEE__)
# undef  __STL_LINK_TIME_INSTANTIATION
# define __STL_LINK_TIME_INSTANTIATION 1
# define __STL_TYPENAME_ON_RETURN_TYPE typename
# ifdef __STDLIB
#  undef __STL_HAS_NO_NEW_IOSTREAMS
#  undef __STL_NO_NEW_C_HEADERS
#  undef __STL_NO_NEW_NEW_HEADER
#  undef __STL_NO_BAD_ALLOC
#  undef __STL_LONG_LONG
# else
#  undef  __STL_NO_EXCEPTION_SPEC
#  define __STL_NO_EXCEPTION_SPEC 1
# endif
#endif

// Modena C++ library 

#if defined (__MWERKS__) || defined (__KCC)
# include <mcompile.h>
# define __STL_USE_MSIPL 1
# if defined (__KCC) || (defined(__MSL_CPP__) && \
       ( (__MSL_CPP__ >= 0x5000 && defined( _MSL_NO_MESSAGE_FACET )) || \
	 (__MSL_CPP__ < 0x5000 && defined( MSIPL_NL_TYPES )))  \
	 )
#  define __STL_NO_NATIVE_MESSAGE_FACET 1
# endif
#endif

// common switches for EDG front-end
# if defined (__EDG_SWITCHES)
#   if !(defined(_TYPENAME) || defined (_TYPENAME_IS_KEYWORD))
#     undef  __STL_NEED_TYPENAME
#     define __STL_NEED_TYPENAME 1
#   endif
#   if !defined(_WCHAR_T_IS_KEYWORD)
#     undef __STL_NO_WCHAR_T 
#     define __STL_NO_WCHAR_T 1
#   endif
#   ifndef _PARTIAL_SPECIALIZATION_OF_CLASS_TEMPLATES
#     undef __STL_NO_CLASS_PARTIAL_SPECIALIZATION
#     define __STL_NO_CLASS_PARTIAL_SPECIALIZATION 1
#   endif
#   ifndef _MEMBER_TEMPLATES
#     undef __STL_NO_MEMBER_TEMPLATES
#     define __STL_NO_MEMBER_TEMPLATES 1
#     undef __STL_NO_MEMBER_TEMPLATE_CLASSES
#     define __STL_NO_MEMBER_TEMPLATE_CLASSES 1
#   endif
#   if !defined(_MEMBER_TEMPLATE_KEYWORD)
#     undef  __STL_NO_MEMBER_TEMPLATE_KEYWORD
#     define __STL_NO_MEMBER_TEMPLATE_KEYWORD 1
#   endif
#   if !defined (__EXCEPTIONS) && ! defined (_EXCEPTIONS)
#     undef  __STL_HAS_NO_EXCEPTIONS
#     define __STL_HAS_NO_EXCEPTIONS
#   endif
#   undef __EDG_SWITCHES
# endif /* EDG */
#endif
