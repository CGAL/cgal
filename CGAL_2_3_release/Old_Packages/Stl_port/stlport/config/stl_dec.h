
# define __STL_NATIVE_HEADER(header) <../cxx/##header>
# define __STL_NATIVE_C_HEADER(x) <../include/##x>
# define __STL_NATIVE_CPP_C_HEADER(header) <../cxx/##header>


#if (__DECCXX_VER < 60000000)

// automatic template instantiation does not
// work with namespaces ;(
# define __STL_HAS_NO_NAMESPACES 1

# define __STL_NO_WCHAR_T  1
# define __STL_NEED_EXPLICIT  1

# define __STL_NO_BOOL  1
# define __STL_NEED_TYPENAME 1
# define __STL_NO_NEW_STYLE_CASTS 1
# define __STL_NEED_MUTABLE 1
# define __STL_NO_BAD_ALLOC 1

# define __STL_NO_NEW_NEW_HEADER 1 
# define __STL_NO_PARTIAL_SPECIALIZATION_SYNTAX 1

# define __STL_NO_MEMBER_TEMPLATES 1
# define __STL_NO_MEMBER_TEMPLATE_CLASSES 1
# define __STL_NO_MEMBER_TEMPLATE_KEYWORD 1
# define __STL_NO_FRIEND_TEMPLATES 1
# define __STL_NO_QUALIFIED_FRIENDS 1
# define __STL_NO_CLASS_PARTIAL_SPECIALIZATION 1
# define __STL_NO_FUNCTION_TMPL_PARTIAL_ORDER 1
# define __STL_NON_TYPE_TMPL_PARAM_BUG 1
# define __STL_BROKEN_USING_DIRECTIVE 1
# define __STL_NO_EXCEPTION_HEADER 1
# define __STL_DEFAULT_CONSTRUCTOR_BUG 1

#endif


#ifdef __NO_USE_STD_IOSTREAM
#  define __STL_HAS_NO_NEW_IOSTREAMS 1
# else
// default is to use new iostreams, anyway
# ifndef __USE_STD_IOSTREAM
#  define __USE_STD_IOSTREAM
# endif
#endif

// # if !defined (__STLPORT_NEW_IOSTREAMS) && ! defined (__STL_DONT_REDEFINE_STD) \
//  && ! defined (__STL_REDEFINE_STD)
// # undef __PRAGMA_ENVIRONMENT
//   #  define __STL_DONT_REDEFINE_STD
// # endif

//# ifndef __STD_STRICT_ANSI_ERRORS
//# endif

#ifndef __EXCEPTIONS
# define __STL_HAS_NO_EXCEPTIONS 1
#endif

# ifdef __IMPLICIT_INCLUDE_ENABLED

#ifndef __STLPORT_IOSTREAMS
// implicit include introduces conflicts
// between stlport and native lib.
# undef __IMPLICIT_INCLUDE_ENABLED
#endif

// but, works with ours ;).
#  define __STL_LINK_TIME_INSTANTIATION 1

# endif

# if defined (__IMPLICIT_USING_STD) && !defined (__NO_USE_STD_IOSTREAM)
// we should ban that !
#  error "STLport won't work with new iostreams and std:: being implicitly included. Please use -std strict_ansi[_errors] or specify __NO_USE_STD_IOSTREAM"
# endif

# if !(defined (__STD_STRICT_ANSI) || defined (__STD_STRICT_ANSI_ERRORS))
// we want to enforce it
#  define __STL_LONG_LONG 1
# endif

// unsigned 32-bit integer type
#  define __STL_UINT32_T unsigned int
#  define __STL_RAND48 1

# define __STL_TYPENAME_ON_RETURN_TYPE typename

#  define __STL_NO_EXPLICIT_FUNCTION_TMPL_ARGS 1

# if (__DECCXX_VER <= 60200000)
#  define __STL_HAS_NO_NEW_C_HEADERS 1 
# endif



