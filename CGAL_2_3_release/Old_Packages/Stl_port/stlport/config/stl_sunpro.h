// STLport configuration file
// It is internal STLport header - DO NOT include it directly

#  define __STL_LONG_LONG    1
// on solaris 2.x only ?
// #  define __STL_LONG_DOUBLE  1
#  define __STL_RAND48 1
#  define __STL_LINK_TIME_INSTANTIATION 1
#  define __STL_USING_BASE_MEMBER


# if (__SUNPRO_CC < 0x600)
#  define __STL_NO_MEMBER_TEMPLATES 1
#  define __STL_NO_MEMBER_TEMPLATE_CLASSES 1
#  define __STL_NO_MEMBER_TEMPLATE_KEYWORD 1
#  define __STL_NO_FRIEND_TEMPLATES 1
#  define __STL_NO_QUALIFIED_FRIENDS 1
#  define __STL_NO_CLASS_PARTIAL_SPECIALIZATION 1
#  define __STL_NO_FUNCTION_TMPL_PARTIAL_ORDER 1
#  define __STL_NON_TYPE_TMPL_PARAM_BUG 1
#  define __STL_STATIC_ARRAY_BUG 1
// #  define __STL_HAS_NO_NEW_C_HEADERS 1 
// no partial , just for explicit one
#  define __STL_PARTIAL_SPEC_NEEDS_TEMPLATE_ARGS
#  define __STL_STATIC_CONST_INIT_BUG 1
#  define __STL_NO_EXPLICIT_FUNCTION_TMPL_ARGS
# endif

// They are present actually, but have problems.
#  if ( __SUNPRO_CC <= 0x500 )
#   define __STL_HAS_NO_NEW_C_HEADERS 1 
#  endif

#  if ( __SUNPRO_CC < 0x500 )
// 4.2 does not like it
#  undef __STL_PARTIAL_SPEC_NEEDS_TEMPLATE_ARGS
#  define __STL_NONTEMPL_BASE_MATCH_BUG 1
#  define __STL_LIMITED_DEFAULT_TEMPLATES 1
#  define __STL_HAS_NO_NEW_IOSTREAMS 1
#  define __STL_NO_NEW_NEW_HEADER 1
#  define __STL_NO_BOOL 1
#  define __STL_HAS_NO_NAMESPACES 1
#  define __STL_NEED_TYPENAME 1
#  define __STL_NEED_EXPLICIT 1
#  define __STL_NEED_MUTABLE  1


#  define __STL_NO_EXCEPTION_HEADER 1
#  define __STL_NO_BAD_ALLOC 1
#  define __STL_UNINITIALIZABLE_PRIVATE 1
#  define __STL_NO_BAD_ALLOC 1
#  define __SGI_STL_NO_ARROW_OPERATOR 1
#  define __STL_DEFAULT_CONSTRUCTOR_BUG 1
#  define __STL_GLOBAL_NESTED_RETURN_TYPE_PARAM_BUG 1
#  define __STL_NO_EXPLICIT_FUNCTION_TMPL_ARGS


#   if ( __SUNPRO_CC < 0x420 )
#    define __STL_NO_PARTIAL_SPECIALIZATION_SYNTAX 1
#    define __STL_NO_NEW_STYLE_CASTS 1
#    define __STL_NO_METHOD_SPECIALIZATION 1

#    if ( __SUNPRO_CC > 0x401 )
#     if (__SUNPRO_CC==0x410)
#      define __STL_BASE_TYPEDEF_OUTSIDE_BUG  1
#     endif
#    else
   // SUNPro C++ 4.0.1
#     define __STL_BASE_MATCH_BUG          1
#     define __STL_BASE_TYPEDEF_BUG        1
#      if ( __SUNPRO_CC < 0x401 )
        __GIVE_UP_WITH_STL(SUNPRO_401)
#      endif
#    endif /* 4.0.1 */
#   endif /* 4.2 */

#  endif /* <  5.0 */

