// STLport configuration file
// It is internal STLport header - DO NOT include it directly

// COMO 4.238 with MSVC

# define __STL_UINT32_T unsigned int

// this one is true only with MS
# if defined (_MSC_VER)
#  define __STL_WCHAR_T_IS_USHORT 1
#  if _MSC_VER <= 1200
#   define __STL_VENDOR_GLOBAL_CSTD
#  endif
#  if _MSC_VER < 1100
#   define __STL_NO_BAD_ALLOC 1
#   define __STL_NO_EXCEPTION_HEADER 1
#   define __STL_HAS_NO_NEW_C_HEADERS 1 
#   define __STL_NO_NEW_NEW_HEADER 1
#   define __STL_NO_NEW_IOSTREAMS 1
#  endif
# endif

// # define __STL_NO_FUNCTION_TMPL_PARTIAL_ORDER 1
// # define __STL_NON_TYPE_TMPL_PARAM_BUG 1
// # define __STL_NO_EXPLICIT_FUNCTION_TMPL_ARGS 1


# define __EDG_SWITCHES


