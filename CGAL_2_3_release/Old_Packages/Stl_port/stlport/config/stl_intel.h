// STLport configuration file
// It is internal STLport header - DO NOT include it directly

// two levels of macros do not work good with icl.
#   define __STL_NATIVE_HEADER(header)    <../include/##header> 
#   define __STL_NATIVE_C_HEADER(header)    <../include/##header> 
#   define __STL_NATIVE_CPP_C_HEADER(header)    <../include/##header> 

// Intel compiler 4.0
#if (__ICL >= 400)

// #   define __STL_NO_EXPLICIT_FUNCTION_TMPL_ARGS 1
// #   define __STL_NO_MEMBER_TEMPLATE_KEYWORD 1

#   define __STL_LONG_LONG 1
// #   define __STL_LONG_DOUBLE            1
#   define __STL_WCHAR_T_IS_USHORT      1
#   define __STL_DEFAULT_CONSTRUCTOR_BUG 1

#   pragma warning(disable:4786)

// <stdio> and the like still put stuff in ::namespace
// up to version 6 
#    define __STL_VENDOR_GLOBAL_CSTD 1

#   if ( _MSC_VER<=1010 )
#    define __STL_DONT_USE_BOOL_TYPEDEF 1
#    define __STL_NO_BAD_ALLOC
#    define __STL_HAS_NO_NEW_IOSTREAMS 1
#    define __STL_NO_NEW_NEW_HEADER 1
#  else
// #   define __STL_YVALS_H 1
#   endif /* 1010 */

#   if (_MSC_VER < 1100)  // MSVC 5.0
#    define __STL_NO_BOOL
#    define __STL_NEED_EXPLICIT      1
// up to 4.2, library is in global namespace
#    define __STL_VENDOR_GLOBAL_STD
#   endif /* 1100 */

#   if _MSC_VER < 1200 /* VC++ 6.0 */
#     define __STL_NON_TYPE_TMPL_PARAM_BUG 1 
#   endif

// these switches depend on compiler flags
#   ifndef _CPPUNWIND
#     define __STL_HAS_NO_EXCEPTIONS
#   endif

#   if defined ( _MT ) && !defined (_NOTHREADS) && !defined (_REENTRANT)
#     define _REENTRANT 1
#   endif

#else /* __ICL >=400 */
// this should work for older versions
# include <config/stl_msvc.h>
#endif /* __ICL >=400 */
