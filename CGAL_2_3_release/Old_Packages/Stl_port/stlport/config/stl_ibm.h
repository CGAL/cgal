// STLport configuration file
// It is internal STLport header - DO NOT include it directly

// AIX xlC 3.1 , 3.0.1 ==0x301
// Visual Age C++ 3.x
// OS-390 C++
// fbp : should be version-specific!

#  define __STL_NO_BOOL 1
#  define __STL_DONT_USE_BOOL_TYPEDEF 1
#  define __STL_LIMITED_DEFAULT_TEMPLATES 1
#  define __STL_HAS_NO_NAMESPACES 1

#  define __STL_NEED_TYPENAME 1
#  define __STL_NEED_EXPLICIT 1

// #  define __STL_HAS_NO_EXCEPTIONS 1
// #  define __STL_NO_EXCEPTION_SPEC 1

#  define __SGI_STL_NO_ARROW_OPERATOR 1
#  define __STL_NO_PROXY_ARROW_OPERATOR 1
#  define __STL_NO_NEW_STYLE_CASTS 1
#  define __STL_NO_WCHAR_T 1
#  define __STL_NEED_MUTABLE 1
#  define __STL_NO_PARTIAL_SPECIALIZATION_SYNTAX 1
#  define __STL_NO_BAD_ALLOC 1

#  define __STL_NO_MEMBER_TEMPLATES 1
#  define __STL_NO_MEMBER_TEMPLATE_CLASSES 1
#  define __STL_NO_MEMBER_TEMPLATE_KEYWORD 1
#  define __STL_NO_FRIEND_TEMPLATES 1
#  define __STL_NO_QUALIFIED_FRIENDS 1
#  define __STL_NO_CLASS_PARTIAL_SPECIALIZATION 1
#  define __STL_NO_FUNCTION_TMPL_PARTIAL_ORDER 1

#  define __STL_NO_DEFAULT_NON_TYPE_PARAM 1
#  define __STL_NO_METHOD_SPECIALIZATION 1
#  define __STL_NO_EXPLICIT_FUNCTION_TMPL_ARGS 1
#  define __STL_NO_EXCEPTION_HEADER 1

#  define __STL_HAS_NO_NEW_IOSTREAMS 1
#  define __STL_HAS_NO_NEW_C_HEADERS 1
#  define __STL_NO_NEW_NEW_HEADER 1

#  define __STL_STATIC_CONST_INIT_BUG 1
// #  define __STL_THROW_RETURN_BUG 1

// I believe this is true for all IBM compilers
// __TEMPINC__ is set when /Ft+ option is used
#  ifdef __TEMPINC__
#    define __STL_LINK_TIME_INSTANTIATION 1
#  endif

#  define __STL_NO_TEMPLATE_CONVERSIONS 1
#  define __STL_UNINITIALIZABLE_PRIVATE 1
#  define __STL_BASE_TYPEDEF_OUTSIDE_BUG 1
#  define __STL_STATIC_ARRAY_BUG 1

// AIX xlC, Visual Age 3.0 for OS/2 and MS 
#  define __STL_TRIVIAL_DESTRUCTOR_BUG  1


#  if ( defined (__MULTI__) && defined (__WINDOWS__))
#   define  __STL_WIN32THREADS 1          // Only Visual Age 3.5 for Windows
#  endif

#  define __STL_NON_TYPE_TMPL_PARAM_BUG 1
#  define __STL_NONTEMPL_BASE_MATCH_BUG 1

#if __IBMCPP__ < 350
#  define __STL_LONG_DOUBLE 1
#endif

#if __IBMCPP__ <= 350
#  define __STL_NEED_UNREACHABLE_RETURN 1
#endif

#if __IBMCPP__ >= 350
#  define __STL_LONG_LONG 1
#endif


