// STLport configuration file
// It is internal STLport header - DO NOT include it directly

// HP compilers


# if __cplusplus >= 199707L
// it is aCC
// use the most conservative configuration as the base
///#  define __STL_DEFAULT_TEMPLATE_PARAM 1
#  define __STL_RAND48 1


#  define __STL_PARTIAL_SPECIALIZATION_BUG 1

#  define __STL_NO_QUALIFIED_FRIENDS       1

#  define __STL_LONG_LONG 1
// #  define __STL_LONG_DOUBLE 1

#  define __STL_CLASS_PARTIAL_SPECIALIZATION 1

// aCC bug ? need explicit args on constructors of partial specialized
// classes
#  define __STL_PARTIAL_SPEC_NEEDS_TEMPLATE_ARGS 1

// maybe present in later versions ?
#  define __STL_NO_MEMBER_TEMPLATE_CLASSES 1
#  define __STL_NO_MEMBER_TEMPLATE_KEYWORD 1

#  define __STL_FUNC_PARTIAL_ORDERING 1
#  define __STL_METHOD_SPECIALIZATION 1
// ?? fbp: really needed ?
#  define __STL_STATIC_ARRAY_BUG 1

// <exception> and stuff is in global namespace
# define __STL_VENDOR_GLOBAL_STD

# define __STL_TYPENAME_ON_RETURN_TYPE typename

#  else
// it is HP's old cfront-based compiler.

#  define __STL_NO_BOOL 1
#  define __STL_DONT_USE_BOOL_TYPEDEF 1
#  define __STL_NO_SIGNED_BUILTINS

#  define __STL_LIMITED_DEFAULT_TEMPLATES 1
#  define __STL_DEFAULT_TYPE_PARAM 1
#  define __STL_NO_STATIC_TEMPLATE_DATA 1

#  define __STL_HAS_NO_NAMESPACES 1

#  define __STL_NEED_TYPENAME 1
#  define __STL_NEED_EXPLICIT 1
#  define __STL_HAS_NO_EXCEPTIONS 1
#  define __STL_NO_EXCEPTION_SPEC 1

#  define __SGI_STL_NO_ARROW_OPERATOR 1
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

#  define __STL_STATIC_CONST_INIT_BUG 1
#  define __STL_THROW_RETURN_BUG 1
#  define __STL_LINK_TIME_INSTANTIATION 1
#  define __STL_NO_TEMPLATE_CONVERSIONS 1

#  endif /* cfront */
