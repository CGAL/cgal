// STLport configuration file
// It is internal STLport header - DO NOT include it directly


#  define __STL_NO_FUNCTION_TMPL_PARTIAL_ORDER 1
#  define __STL_NO_CLASS_PARTIAL_SPECIALIZATION 1
#  define __STL_NO_MEMBER_TEMPLATE_KEYWORD 1
#  define __STL_NO_MEMBER_TEMPLATES 1
#  define __STL_NO_FRIEND_TEMPLATES 1
#  define __STL_NO_MEMBER_TEMPLATE_CLASSES 1


#  define __STL_LIMITED_DEFAULT_TEMPLATES 1
#  define __STL_HAS_NO_NAMESPACES 1
#  define __STL_NEED_TYPENAME 1
#  define __STL_NO_WCHAR_T 1
#  define __STL_NO_PARTIAL_SPECIALIZATION_SYNTAX 1
#  define __STL_NO_EXPLICIT_FUNCTION_TMPL_ARGS 1

#  define __STL_STATIC_CONST_INIT_BUG 1
#  define __STL_THROW_RETURN_BUG 1
#  define __STL_NO_TEMPLATE_CONVERSIONS 1

#  define __STL_BASE_TYPEDEF_OUTSIDE_BUG 1

#  define __STL_HAS_NO_NEW_IOSTREAMS 1
#  define __STL_HAS_NO_NEW_C_HEADERS 1
#  define __STL_NO_NEW_NEW_HEADER 1

#  define __STL_NO_DEFAULT_NON_TYPE_PARAM 1
#  define __STL_NON_TYPE_TMPL_PARAM_BUG 1
#  define __STL_NONTEMPL_BASE_MATCH_BUG
#  define __STL_NO_EXCEPTION_HEADER 1
#  define __STL_NO_BAD_ALLOC 1
#  define __SGI_STL_NO_ARROW_OPERATOR 1
#  define __STL_NESTED_TYPE_PARAM_BUG 1

#  if (__WATCOM_CPLUSPLUS__ < 1100 )
#   define __STL_NO_BOOL 1
#   define __STL_NEED_EXPLICIT 1
#   define __STL_NEED_MUTABLE 1
#   define __STL_NO_NEW_STYLE_CASTS 1
#  endif

// Get rid of Watcom's min and max macros 
#undef min 
#undef max

// for switches (-xs,  -xss,  -xst)
//
#if !(defined (__SW_XS) || defined (__SW_XSS) || defined(__SW_XST))
#    define __STL_HAS_NO_EXCEPTIONS 1
# endif
