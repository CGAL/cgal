// STLport configuration file
// It is internal STLport header - DO NOT include it directly

// AS/400 C++ config

#  define __STL_NO_BOOL
#  define __STL_LIMITED_DEFAULT_TEMPLATES

#  define __STL_HAS_NO_NAMESPACES
#  define __STL_NEED_TYPENAME
#  define __STL_NEED_EXPLICIT
#  define __STL_HAS_NO_EXCEPTIONS
#  define __STL_NO_EXCEPTION_SPEC
#  define __SGI_STL_NO_ARROW_OPERATOR
#  define __STL_NO_NEW_STYLE_CASTS

#  define __STL_NEED_MUTABLE
#  define __STL_NO_PARTIAL_SPECIALIZATION_SYNTAX
#  define __STL_NO_BAD_ALLOC
#  define __STL_NO_MEMBER_TEMPLATES
#  define __STL_NO_MEMBER_TEMPLATE_CLASSES
#  define __STL_NO_MEMBER_TEMPLATE_KEYWORD
#  define __STL_NO_FRIEND_TEMPLATES
#  define __STL_NO_QUALIFIED_FRIENDS
#  define __STL_NO_CLASS_PARTIAL_SPECIALIZATION
#  define __STL_NO_FUNCTION_TMPL_PARTIAL_ORDER

#  define __STL_NO_METHOD_SPECIALIZATION
#  define __STL_NO_EXPLICIT_FUNCTION_TMPL_ARGS
#  define __STL_NO_EXCEPTION_HEADER

#  define __STL_HAS_NO_NEW_IOSTREAMS
#  define __STL_HAS_NO_NEW_C_HEADERS 

#  define __STL_STATIC_CONST_INIT_BUG
#  define __STL_THROW_RETURN_BUG
#  define __STL_LINK_TIME_INSTANTIATION
#  define __STL_NO_TEMPLATE_CONVERSIONS

#  define __STL_UNINITIALIZABLE_PRIVATE 1
#  define __STL_STATIC_ARRAY_BUG 1
#  define __STL_NON_TYPE_TMPL_PARAM_BUG 1
#  define __STL_TRIVIAL_DESTRUCTOR_BUG  1

#  if defined(_LONG_LONG)
#    define __STL_LONG_LONG 1
#  endif
// #  define __STL_LONG_DOUBLE 1
#  if defined(_PTHREADS)
#    define _MULTI_THREADED
#  endif
// fbp : to fix __partition() problem
# define __STL_NONTEMPL_BASE_MATCH_BUG 1
