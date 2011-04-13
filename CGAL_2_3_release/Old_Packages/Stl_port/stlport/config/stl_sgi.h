// STLport configuration file
// It is internal STLport header - DO NOT include it directly

#  include <standards.h>

#  define __EDG_SWITCHES

// any version ???
#  define __STL_AUTOMATIC_TYPE_TRAITS 1

#  define __STL_USE_SGI_STRING         1

#  define __STL_HAS_NO_NEW_C_HEADERS 1

// #  if !(_COMPILER_VERSION >= 700)
#   define __STL_NO_NEW_NEW_HEADER 1
// #  endif

#   if !defined(_BOOL)
#     define __STL_NO_BOOL
#   endif
#   if defined(_MIPS_SIM) && _MIPS_SIM == _ABIO32
#     define __STL_STATIC_CONST_INIT_BUG
#   endif

#   if defined(_MIPS_SIM) && _MIPS_SIM != _MIPS_SIM_ABIO32
#     define __STL_DEFAULT_CONSTRUCTOR_BUG
#   endif
#   if !((_COMPILER_VERSION >= 730) && defined(_MIPS_SIM) && _MIPS_SIM != _ABIO32)
#     define __STL_NO_MEMBER_TEMPLATE_KEYWORD
#   endif
#   if !defined(_STANDARD_C_PLUS_PLUS)
#     define __STL_NO_EXPLICIT_FUNCTION_TMPL_ARGS
#   endif
#   if !((_COMPILER_VERSION >= 721) && defined(_NAMESPACES))
#     define __STL_HAS_NO_NAMESPACES
#   endif 
#   if (_COMPILER_VERSION < 721)
#     define __STL_NO_EXCEPTION_HEADER
#   endif
#   if _COMPILER_VERSION < 730 || !defined(_STANDARD_C_PLUS_PLUS) || !defined(_NAMESPACES)
#     define __STL_NO_BAD_ALLOC
#   endif
#   if defined(_LONGLONG) && defined(_SGIAPI) && _SGIAPI
#     define __STL_LONG_LONG
#   endif
#   if !(_COMPILER_VERSION >= 730 && defined(_STANDARD_C_PLUS_PLUS))
#     define __STL_HAS_NO_NEW_IOSTREAMS
#   endif
#   if !(_COMPILER_VERSION >= 730 && defined(_STANDARD_C_PLUS_PLUS))
#     define __STL_NO_AT_MEMBER_FUNCTION
#   endif
// #   if !(_COMPILER_VERSION >= 730 && defined(_STANDARD_C_PLUS_PLUS))
#   if !(_COMPILER_VERSION >= 721 && defined(_STANDARD_C_PLUS_PLUS))
#     define __STL_NO_TEMPLATE_CONVERSIONS
#   endif
#   if !((_COMPILER_VERSION >= 730) && defined(_MIPS_SIM) && _MIPS_SIM != _ABIO32)
#     define __STL_NO_FUNCTION_TMPL_PARTIAL_ORDER
#   endif
