// STLport configuration file
// It is internal STLport header - DO NOT include it directly
// Microsoft Visual C++ 4.0, 4.1, 4.2, 5.0

// Common features for VC++ 4.0 and higher
#   define __STL_NO_FUNCTION_TMPL_PARTIAL_ORDER 1

// CGAL_HACK
// #   define __STL_NO_CLASS_PARTIAL_SPECIALIZATION 1
#define CGAL_LIMITED_ITERATOR_TRAITS_SUPPORT 1

#   define __STL_NO_EXPLICIT_FUNCTION_TMPL_ARGS 1
#   define __STL_NO_MEMBER_TEMPLATE_KEYWORD 1
#   define __STL_NO_MEMBER_TEMPLATE_CLASSES 1
#   define __STL_NO_FRIEND_TEMPLATES

// #   define __STL_LONG_DOUBLE            1
#   define __STL_WCHAR_T_IS_USHORT      1
// up to VC6 ?
#   define __STL_NO_QUALIFIED_FRIENDS    1
#   define __STL_DEFAULT_CONSTRUCTOR_BUG 1
#   define __STL_STATIC_CONST_INIT_BUG   1

// <NBulteau@jouve.fr> : suppressed "truncated debug info" warning
#   pragma warning(disable:4786)

// <stdio> and the like still put stuff in ::namespace
// up to version 6 
#  define __STL_VENDOR_GLOBAL_CSTD

#    define __STL_DONT_USE_BOOL_TYPEDEF 1

#   if ( _MSC_VER<=1010 )
// "bool" is reserved in MSVC 4.1 while <yvals.h> absent, so :
// #    define __STL_USE_ABBREVS           1
#    define __STL_NO_BAD_ALLOC
#    define __STL_HAS_NO_NEW_C_HEADERS 1
#    define __STL_NO_NEW_NEW_HEADER 1
#    define __STL_HAS_NO_NEW_IOSTREAMS 1
#   else
// VC++ 4.2 and higher
#    define __STL_YVALS_H 1
#   endif /* 1010 */

#   if (_MSC_VER >= 1100)  // MSVC 5.0
// these work, as long they are inline
#    define __STL_INLINE_MEMBER_TEMPLATES 1
#   else
#    define __STL_NO_BOOL            1
#    define __STL_NEED_TYPENAME      1
#    define __STL_NEED_EXPLICIT      1
#    define __STL_NEED_MUTABLE       1
#    define __STL_NO_PARTIAL_SPECIALIZATION_SYNTAX
#    define __STL_LIMITED_DEFAULT_TEMPLATES 1
#    define __STL_NO_MEMBER_TEMPLATES 1

// up to 4.2, library is in global namespace
// #    define __STL_NESTED_TYPE_PARAM_BUG
#    define __STL_VENDOR_GLOBAL_STD
#    define __STL_NONTEMPL_BASE_MATCH_BUG 1
#    define __STL_BROKEN_USING_DIRECTIVE  1
#    define __SGI_STL_NO_ARROW_OPERATOR
#    define __STL_NO_EXCEPTION_SPEC 1
#   endif /* 1100 */

#   if _MSC_VER < 1200 /* VC++ 5.0  (not 6.0, apparently DVP)*/
#    define __STL_NON_TYPE_TMPL_PARAM_BUG 1 
#    define __STL_THROW_RETURN_BUG 1
#   endif

// these switches depend on compiler flags
#   ifndef _CPPUNWIND
#     define __STL_HAS_NO_EXCEPTIONS 1
#   endif

#   if defined ( _MT ) && !defined (_NOTHREADS) && !defined (_REENTRANT)
#     define _REENTRANT 1
#   endif

// If we are under Windows CE, include appropriate config

# ifdef UNDER_CE
#   include <config/stl_wince.h>
# endif


