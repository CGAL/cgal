// STLport configuration file
// It is internal STLport header - DO NOT include it directly

// g++ 2.7.x and above 

#   define __STL_LONG_LONG    1
// #   define __STL_LONG_DOUBLE  1

#   define __STL_HAS_NO_NEW_IOSTREAMS     1
#   define __STL_VENDOR_GLOBAL_CSTD       1
#   define __STL_NO_NATIVE_MBSTATE_T      1
#   define __STL_NO_NATIVE_WIDE_FUNCTIONS 1

// gcc fails to initialize builtin types in expr. like this : new(p) char(); 

# define __STL_DEFAULT_CONSTRUCTOR_BUG 1
// this should always work

#   if (__GNUC_MINOR__ < 90) /* egcs 1.1 */
#     define __STL_NO_TEMPLATE_CONVERSIONS
// #     define __STL_NO_MEMBER_TEMPLATES 1
#     define __STL_NO_MEMBER_TEMPLATE_CLASSES 1
#     define __STL_NO_FUNCTION_TMPL_PARTIAL_ORDER 1
#     define __STL_NO_FRIEND_TEMPLATES 1
#     define __STL_HAS_NO_NAMESPACES 1
//  DJGPP doesn't seem to implement it in 2.8.x
#    ifdef DJGPP
#     define  __STL_NO_STATIC_TEMPLATE_DATA 1
#    endif
#   endif

#  if __GNUC__ <= 2 && __GNUC_MINOR__ <= 7 && ! defined (__CYGWIN32__)

// Will it work with 2.6 ? I doubt it.
#   if ( __GNUC_MINOR__ < 6 )
    __GIVE_UP_WITH_STL(GCC_272);
#   endif

# define  __STL_LIMITED_DEFAULT_TEMPLATES 1
# define  __STL_DEFAULT_TYPE_PARAM 1

# define  __STL_NO_BAD_ALLOC
# define  __SGI_STL_NO_ARROW_OPERATOR 1
# ifndef __STL_NO_STATIC_TEMPLATE_DATA
#  define  __STL_NO_STATIC_TEMPLATE_DATA
# endif
# define  __STL_NO_MEMBER_TEMPLATES 1
# define  __STL_NO_CLASS_PARTIAL_SPECIALIZATION 1
# define  __STL_NO_METHOD_SPECIALIZATION 1

#  if !defined (__CYGWIN32__) 
#   define __STL_NESTED_TYPE_PARAM_BUG   1
#   define __STL_BASE_MATCH_BUG       1
//  unused operators are required (forward)
#   define  __STL_CONST_CONSTRUCTOR_BUG 
#   define __STL_NO_DEFAULT_NON_TYPE_PARAM
#  endif

#   define __STL_NO_PARTIAL_SPECIALIZATION_SYNTAX 1
#   define __STL_NO_EXPLICIT_FUNCTION_TMPL_ARGS 1
#   define __STL_NO_EXCEPTION_HEADER 1

#  else /* ! <= 2.7.* */

#  endif /* ! <= 2.7.* */


// static template data members workaround strategy for gcc tries
// to use weak symbols.
// if you don't want to use that, #define __STL_WEAK_ATTRIBUTE=0 ( you'll
// have to put "#define __PUT_STATIC_DATA_MEMBERS_HERE" line in one of your
// compilation unit ( or CFLAGS for it ) _before_ including any STL header ).
#   if defined (__STL_NO_STATIC_TEMPLATE_DATA) && ! defined (__STL_WEAK_ATTRIBUTE )
// systems using GNU ld or format that supports weak symbols
// may use "weak" attribute
// Linux & Solaris ( x86 & SPARC ) are being auto-recognized here
#    if defined(__STL_GNU_LD) || defined(__ELF__) || \
     (( defined (__SVR4) || defined ( __svr4__ )) && \
      ( defined (sun) || defined ( __sun__ )))
#     define __STL_WEAK_ATTRIBUTE 1
#    endif
#   endif /* __STL_WEAK_ATTRIBUTE */


