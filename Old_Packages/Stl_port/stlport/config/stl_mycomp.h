/*
 * Copyright (c) 1997
 * Moscow Center for SPARC Technology
 *
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
 * Purpose of this file :
 *
 * A list of COMPILER-SPECIFIC portion of STLport settings.
 * This file is provided to help in manulal configuration
 * of STLport. This file is being included by stlcomp.h 
 * when STLport is unable to identify your compiler.
 * Please remove the error diagnostic below before adjusting 
 * macros.
 * 
 */
# ifndef __STLPORT_MYCOMP_H
#  define  __STLPORT_MYCOMP_H

// Error diagnostic 

# error "Your compiler version is not recognized by STLport. Please edit <config/stl_mycomp.h>

//==========================================================

// the values choosen here as defaults try to give
// maximum functionality on the most conservative settings

// Mostly correct guess, change it for Alpha (and other environments
// that has 64-bit "long")
// #  define __STL_UINT32_T unsigned long

// Disables wchar_t functinality
// #  define __STL_NO_WCHAR_T  1

// Define if wchar_t is not a unique type, and is actually a typedef to unsigned short. 
// #  define __STL_WCHAR_T_IS_USHORT 1

// Uncomment if long long is available
// #  define __STL_LONG_LONG 1

// Uncomment if long double is not available
// #  define __STL_NO_LONG_DOUBLE 1

// Uncomment this if your compiler does not support "typename" keyword
// #  define __STL_NEED_TYPENAME 1

// Uncomment this if your compiler does not support "mutable" keyword
// #  define __STL_NEED_MUTABLE 1

// Uncomment this if your compiler does not support "explicit" keyword
// #  define __STL_NEED_EXPLICIT 1

// Uncomment if new-style-casts like const_cast<> are not available
// #  define __STL_NO_NEW_STYLE_CASTS 1

// Uncomment this if your compiler does not have "bool" type
// #  define  __STL_NO_BOOL 1

// Uncomment this if your compiler does not have "bool" type, but has "bool" keyword reserved
// #  define  __STL_DONT_USE_BOOL_TYPEDEF 1

// Uncomment this if your compiler does not have "bool" type, but defines "bool" in <yvals.h>
// #  define  __STL_YVALS_H 1

// Uncomment this if your compiler has limited or no default template arguments for classes
// #  define __STL_LIMITED_DEFAULT_TEMPLATES 1

// Uncomment this if your compiler support only complete (not dependent on other parameters)
// types as default parameters for class templates
// #  define __STL_DEFAULT_TYPE_PARAM 1

// Uncomment this if your compiler has problem with not-type
// default template parameters
// #  define __STL_NO_DEFAULT_NON_TYPE_PARAM 1

// Define if compiler has
// trouble with functions getting non-type-parameterized classes as parameters
// #  define __STL_NON_TYPE_TMPL_PARAM_BUG 1

// Uncomment this if your compiler lacks static data members.
// Uncomment next line if your compiler supports __attribute__((weak))
// #  define __STL_NO_STATIC_TEMPLATE_DATA 1
// #  define __STL_WEAK_ATTRIBUTE 1

// Uncomment this if your compiler does not support namespaces 
// #  define __STL_HAS_NO_NAMESPACES 1

// Uncomment if "using" keyword does not work with template types 
// # define __STL_BROKEN_USING_DIRECTIVE 1

// Uncomment this if your compiler does not support exceptions
// #  define __STL_HAS_NO_EXCEPTIONS 1

// Uncomment this if your compiler does not support exception specifications
// #  define __STL_NO_EXCEPTION_SPEC

// Define this if your compiler requires return statement after throw()
// # define __STL_THROW_RETURN_BUG 1

// Header <new> that comes with the compiler 
// does not define bad_alloc exception
// #  define __STL_NO_BAD_ALLOC  1

// Uncomment if member template methods are not available
// #  define __STL_NO_MEMBER_TEMPLATES   1

// Uncomment if member template classes are not available
// #  define __STL_NO_MEMBER_TEMPLATE_CLASSES   1

// Uncomment if no "template" keyword should be used with member template classes
// #  define __STL_NO_MEMBER_TEMPLATE_KEYWORD   1

// Uncomment if friend member templates are not available
// #  define __STL_NO_FRIEND_TEMPLATES   1

// Compiler does not accept friend declaration qualified with namespace name.
// #  define __STL_NO_QUALIFIED_FRIENDS 1

// Uncomment if partial specialization is not available
// #  define __STL_NO_CLASS_PARTIAL_SPECIALIZATION 1

// Define if class being partially specialized require full name (template parameters)
// of itself for method declarations
// #  define __STL_PARTIAL_SPEC_NEEDS_TEMPLATE_ARGS

// partial specialization has bugs that prevent you from
// using new-style reverse_iterator
// #  define __STL_PARTIAL_SPECIALIZATION_BUG

// Compiler has problems specializing members of partially 
// specialized class
// #  define __STL_MEMBER_SPECIALIZATION_BUG

// Uncomment if partial order of template functions is not available
// #  define __STL_NO_FUNC_PARTIAL_ORDERING 1

// Uncomment if specialization of methods is not allowed
// #  define __STL_NO_METHOD_SPECIALIZATION  1

// Uncomment if full  specialization does not use partial spec. syntax : template <> struct ....
// #  define __STL_NO_PARTIAL_SPECIALIZATION_SYNTAX  1

// Uncomment if compiler does not support explicit template arguments for functions
// # define __STL_NO_EXPLICIT_FUNCTION_TMPL_ARGS

// Uncomment if "__typetraits is being instaniated automatically by the compiler (SGI only ?)
// #  define __STL_AUTOMATIC_TYPE_TRAITS 1

// Uncomment this if your compiler can't inline while(), for()
// #  define __STL_LOOP_INLINE_PROBLEMS 1

// Define if the compiler fails to match a template function argument of base
// #  define __STL_BASE_MATCH_BUG          1

// Define if the compiler fails to match a template function argument of base
// (non-template)
//#  define  __STL_NONTEMPL_BASE_MATCH_BUG 1

// Define if the compiler rejects outline method definition 
// explicitly taking nested types/typedefs
// #  define __STL_NESTED_TYPE_PARAM_BUG   1

// Compiler requires typename keyword on outline method definition 
// explicitly taking nested types/typedefs
// #define  __STL_TYPENAME_ON_RETURN_TYPE

// Define if the baseclass typedefs not visible from outside
// #  define __STL_BASE_TYPEDEF_OUTSIDE_BUG 1

// if your compiler have serious problems with typedefs, try this one
// #  define __STL_BASE_TYPEDEF_BUG          1

// Uncomment if getting errors compiling mem_fun* adaptors
// #  define __STL_MEMBER_POINTER_PARAM_BUG 1

// #  define __STL_UNINITIALIZABLE_PRIVATE  1

// Defined if the compiler
// has trouble instantiating static array members with dimension defined as enum
// # define __STL_STATIC_ARRAY_BUG

// * __STL_STATIC_CONST_INIT_BUG: defined if the compiler can't handle a
//   constant-initializer in the declaration of a static const data member
//   of integer type.  (See section 9.4.2, paragraph 4, of the C++ standard.)
// # define __STL_STATIC_CONST_INIT_BUG

// Define if default constructor for builtin integer type fails to initialize it to 0
// #  define __STL_DEFAULT_CONSTRUCTOR_BUG    1

// Defined if constructor
// required to explicitly call member's default constructors for const objects
// #  define __STL_CONST_CONSTRUCTOR_BUG    1

// Defined if the compiler has trouble calling POD-types constructors/destructors
// #  define __STL_TRIVIAL_CONSTRUCTOR_BUG    1
// #  define __STL_TRIVIAL_DESTRUCTOR_BUG    1

// Define if having problems specializing maps/sets with
// key type being const 
// #  define __STL_MULTI_CONST_TEMPLATE_ARG_BUG

// Uncomment this to disable -> operators on all iterators
// #  define   __SGI_STL_NO_ARROW_OPERATOR 1

// Uncomment this to disble at() member functions for containers
// #  define   __STL_NO_AT_MEMBER_FUNCTION 1

// Uncomment if native new-style iostreams are not available
// #define    __STL_HAS_NO_NEW_IOSTREAMS	1

// Define this if compiler lacks <exception> header
// #  define __STL_NO_EXCEPTION_HEADER 1

// Uncomment this if your C library has lrand48() function
// #  define __STL_RAND48 1

// Uncomment if native new-style C library headers lile <cstddef>, etc are not available.
// #   define __STL_HAS_NO_NEW_C_HEADERS 1

// uncomment if new-style headers <new> is available
// #  define __STL_HAS_NEW_NEW_HEADER 1

// uncomment this if <iostream> and other STD headers put their stuff in ::namespace,
// not std::
// #  define __STL_VENDOR_GLOBAL_STD

// uncomment this if <cstdio> and the like put stuff in ::namespace,
// not std::
// #  define __STL_VENDOR_GLOBAL_CSTD

// Edit relative path below (or put full path) to get native 
// compiler headers included. Default is "../include".
// C headers may reside in different directory, so separate macro is provided.
// Hint : never install STLport in the directory that ends with "include"
// # define __STL_NATIVE_INCLUDE_PATH ../include
// # define __STL_NATIVE_C_INCLUDE_PATH ../include
// # define __STL_NATIVE_CPP_C_INCLUDE_PATH ../include

// This macro constructs header path from directory and name.
// You may change it if your compiler does not understand "/". 
// #  define __STL_MAKE_HEADER(path, header) <path/header>

// This macro constructs native include header path from include path and name.
// You may have do define it if experirncing problems with preprocessor
// # define __STL_NATIVE_HEADER(header) __STL_MAKE_HEADER(__STL_NATIVE_INCLUDE_PATH,header)

// Same for C headers
// __STL_NATIVE_C_HEADER(header)

//==========================================================
# endif
