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

#ifndef __SGI_STL_STLCONF_H
# define __SGI_STL_STLCONF_H

# undef __AUTO_CONFIGURED

//==========================================================
// Getting proper values of autoconf flags
// if you ran 'configure', __AUTO_CONFIGURED is set to 1 and
// specific compiler features will be used.
// Otherwise, the <stlcomp.h> header will be included for per-version
// features recognition.
//==========================================================
# if defined (__AUTO_CONFIGURED)
// auto-configured section

# undef __STL_NO_EXCEPTIONS
# undef __STL_NO_NAMESPACES
# undef __STL_NO_RELOPS_NAMESPACE
# undef __STL_NO_NEW_NEW_HEADER 

# undef __STL_NO_NEW_IOSTREAMS

// select threads strategy
# undef _PTHREADS
# undef _NOTHREADS

// select SGI-style alloc instead of allocator<T>
# undef __STL_USE_SGI_ALLOCATORS

// select allocation method you like
# undef __STL_USE_MALLOC
# undef __STL_USE_NEWALLOC

// this one is not mandatory, just enabled
# undef __STL_USE_DEFALLOC

// define __STL_USE_ABBREVS if your linker has trouble with long 
// external symbols
# undef __STL_USE_ABBREVS


// unsigned 32-bit integer type
#  define __STL_UINT32_T unsigned
#  undef __STL_NO_BOOL
#  undef __STL_DONT_USE_BOOL_TYPEDEF
#  undef __STL_YVALS_H
#  undef __STL_LIMITED_DEFAULT_TEMPLATES
#  undef __STL_DEFAULT_TYPE_PARAM
#  undef __STL_NO_STATIC_TEMPLATE_DATA
#  undef __STL_RAND48
#  undef __STL_LOOP_INLINE_PROBLEMS

#  undef __STL_HAS_NO_NAMESPACES

#  undef __STL_NEED_TYPENAME
#  undef __STL_NEED_EXPLICIT
#  undef __STL_HAS_NO_EXCEPTIONS
#  undef __STL_NO_EXCEPTION_SPEC
#  undef __STL_WEAK_ATTRIBUTE
#  undef __STL_BASE_MATCH_BUG
#  undef __STL_NONTEMPL_BASE_MATCH_BUG
#  undef __STL_NESTED_TYPE_PARAM_BUG
#  undef __SGI_STL_NO_ARROW_OPERATOR
#  undef __STL_UNINITIALIZABLE_PRIVATE
#  undef __STL_BASE_TYPEDEF_BUG
#  undef __STL_BASE_TYPEDEF_OUTSIDE_BUG
#  undef __STL_CONST_CONSTRUCTOR_BUG

#  undef __STL_NO_NEW_STYLE_CASTS
#  undef __STL_NO_WCHAR_T
#  undef __STL_WCHAR_T_IS_USHORT
#  undef __STL_LONG_LONG
#  undef __STL_NO_LONG_DOUBLE
#  undef __STL_NEED_MUTABLE
#  undef __STL_NO_PARTIAL_SPECIALIZATION_SYNTAX
#  undef __STL_NO_BAD_ALLOC
#  undef __STL_DEBUG_ALLOC
#  undef __STL_NO_MEMBER_TEMPLATES
#  undef __STL_NO_MEMBER_TEMPLATE_CLASSES
#  undef __STL_NO_MEMBER_TEMPLATE_KEYWORD
#  undef __STL_NO_FRIEND_TEMPLATES
#  undef __STL_NO_QUALIFIED_FRIENDS
#  undef __STL_NO_CLASS_PARTIAL_SPECIALIZATION
#  undef __STL_NO_FUNCTION_TMPL_PARTIAL_ORDER
#  undef __STL_AUTOMATIC_TYPE_TRAITS
#  undef __STL_MEMBER_POINTER_PARAM_BUG
#  undef __STL_NON_TYPE_TMPL_PARAM_BUG
#  undef __STL_NO_DEFAULT_NON_TYPE_PARAM
#  undef __STL_NO_METHOD_SPECIALIZATION
#  undef __STL_STATIC_ARRAY_BUG
#  undef __STL_STATIC_CONST_INIT_BUG
#  undef __STL_TRIVIAL_CONSTRUCTOR_BUG
#  undef __STL_TRIVIAL_DESTRUCTOR_BUG
#  undef __STL_BROKEN_USING_DIRECTIVE
#  undef __STL_NO_EXPLICIT_FUNCTION_TMPL_ARGS
#  undef __STL_NO_EXCEPTION_HEADER
#  undef __STL_DEFAULT_CONSTRUCTOR_BUG

#  undef __STL_HAS_NO_NEW_IOSTREAMS
#  undef __STL_HAS_NO_NEW_C_HEADERS 
#  undef __STL_STATIC_CONST_INIT_BUG
// new ones
#  undef __STL_THROW_RETURN_BUG
// unimp
#  undef __STL_LINK_TIME_INSTANTIATION
#  undef __STL_PARTIAL_SPEC_NEEDS_TEMPLATE_ARGS
// unimp
#  undef __STL_NO_TEMPLATE_CONVERSIONS
# endif /* AUTO_CONFIGURED */

//==========================================================


#endif /* __STLCONF_H */
