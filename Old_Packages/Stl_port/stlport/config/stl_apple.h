// STLport configuration file
// It is internal STLport header - DO NOT include it directly

// common configuration settings for Apple MPW MrCpp / SCpp
#ifdef qMacApp
# ifndef __CONDITIONALMACROS__ // skip including ConditionalMacros_AC.h if ConditionalMacros.h is already included
# include "CoreSwitches_AC.h"
# include "ConditionalMacros_AC.h"
# include "Types_AC.h"
# define __STL_FILE__ _FILE_AC
# define __STL_DEBUG_MESSAGE      
# define __stl_debug_message ProgramBreak_AC
# include <ConditionalMacros.h>
# endif
# include <Types.h>
#else
# include <ConditionalMacros.h>
# include <Types.h>
#endif

#define __STL_UINT32_T UInt32
typedef int wint_t;

#ifndef TYPE_BOOL
# error <ConditionalMacros.h> must be included. (TYPE_BOOL)
#endif
#if !TYPE_BOOL
# define __STL_NO_BOOL
# define __STL_DONT_USE_BOOL_TYPEDEF
#endif

#ifndef TYPE_LONGLONG
# error <ConditionalMacros.h> must be included. (TYPE_LONGLONG)
#endif
#if TYPE_LONGLONG
# define __STL_LONG_LONG
#endif

#if !__option(exceptions)
# define __STL_HAS_NO_EXCEPTIONS
#endif

#define __STL_DEBUG_MESSAGE_POST DebugStr("\pSTL diagnosis issued. See 'stderr' for detail.");
#define __STL_ASSERT_MSG_TRAILER " "
#ifndef __STL_NATIVE_INCLUDE_PATH
# define __STL_NATIVE_INCLUDE_PATH ::CIncludes // expects the alias to {CIncludes} under the same folder as {STL}
#endif
# if !defined(__STL_MAKE_HEADER)
#  define __STL_MAKE_HEADER(path, header) <path:header> // Mac uses ":" for directory delimiter
# endif


// Apple MPW SCpp
#if defined (__SC__)
# define __STL_NO_EXCEPTION_SPEC					// known limitation
# define __STL_NO_MUTABLE							// known limitation
# define __STL_HAS_NO_NAMESPACES					// known limitation
# define __STL_NO_TYPENAME							// known limitation
# define __STL_NO_EXPLICIT							// known limitation

# define __STL_NO_BAD_ALLOC							// known limitation
# define __STL_HAS_NO_NEW_C_HEADERS					// known limitation
# define __STL_NO_NEW_NEW_HEADER					// known limitation
# define __STL_HAS_NO_NEW_IOSTREAMS					// known limitation

# define __STL_LIMITED_DEFAULT_TEMPLATES			// known limitation
# define __STL_NO_MEMBER_TEMPLATES					// known limitation
# define __STL_NO_MEMBER_TEMPLATE_CLASSES			// known limitation
# define __STL_NO_FRIEND_TEMPLATES					// known limitation
# define __STL_NO_PARTIAL_SPECIALIZATION_SYNTAX		// known limitation
# define __STL_NO_CLASS_PARTIAL_SPECIALIZATION		// known limitation
# define __STL_NO_FUNCTION_TMPL_PARTIAL_ORDER		// known limitation

# define __STL_NATIVE_HEADER(header)	<::CIncludes:##header##>	// since SCpp has problem expanding symbols recursively
# define __STL_NATIVE_C_HEADER(header)	<::CIncludes:##header##>
# define __STL_NATIVE_CPP_C_HEADER(header)	<::CIncludes:##header##>

# define __STL_GLOBAL_NESTED_RETURN_TYPE_PARAM_BUG
# define __STL_DEFAULT_PARAM_CONSTRUCTOR_BUG
# define __STL_BOGUS_TEMPLATE_TYPE_MATCHING_BUG
# define __STL_MPW_EXTRA_CONST const  

# define __STL_NON_TYPE_TMPL_PARAM_BUG
# define __STL_THROW_RETURN_BUG
# define __SGI_STL_NO_ARROW_OPERATOR
# define __STL_NO_PROXY_ARROW_OPERATOR
# define __STL_NO_EXPLICIT_FUNCTION_TMPL_ARGS

# define __STL_USE_ABBREVS							// known limitation
// Note : this line is dangerous, since it remanes the public class.
// Only use it if experience major difficulties with symbol length 
// during compilation.
// # define allocator _Al
#endif // defined (__SC__)



// Apple MPW MrCpp 4.1.0
#if defined (__MRC__)
# define __STL_NO_TYPENAME							// known limitation
# define __STL_BROKEN_USING_DIRECTIVE				// known limitation

# define __STL_NO_BAD_ALLOC							// known limitation
# define __STL_HAS_NO_NEW_C_HEADERS					// known limitation
# define __STL_NO_NEW_NEW_HEADER					// known limitation
# define __STL_HAS_NO_NEW_IOSTREAMS					// known limitation

# define __STL_LIMITED_DEFAULT_TEMPLATES			// known limitation
# define __STL_NO_MEMBER_TEMPLATES					// known limitation
# define __STL_NO_MEMBER_TEMPLATE_CLASSES			// known limitation
# define __STL_NO_FRIEND_TEMPLATES					// known limitation
# define __STL_NO_PARTIAL_SPECIALIZATION_SYNTAX		// known limitation
# define __STL_NO_CLASS_PARTIAL_SPECIALIZATION		// known limitation
# define __STL_NO_FUNCTION_TMPL_PARTIAL_ORDER		// known limitation

# define __STL_HAS_NO_NAMESPACES					//*TY 08/01/1999 - still active bug under MrCpp 4.1.0a8
# define __STL_NO_EXPLICIT 							//*TY 08/01/1999 - still active bug under MrCpp 4.1.0a8
# define __STL_NON_TYPE_TMPL_PARAM_BUG				//*TY 08/01/1999 - still active bug under MrCpp 4.1.0a8
# define __STL_THROW_RETURN_BUG						//*TY 08/01/1999 - still active bug under MrCpp 4.1.0a8
# define __SGI_STL_NO_ARROW_OPERATOR 				//*TY 08/01/1999 - still active bug under MrCpp 4.1.0a8
# define __STL_NO_PROXY_ARROW_OPERATOR				//*TY 08/01/1999 - still active bug under MrCpp 4.1.0a8
# define __STL_NO_EXPLICIT_FUNCTION_TMPL_ARGS		//*TY 08/01/1999 - still active bug under MrCpp 4.1.0a8
#endif // defined (__MRC__)



// temporal. workaround for stl_config.h problem
#ifdef __STL_NO_MUTABLE
# define __STL_NEED_MUTABLE		//*TY 08/01/1999 - should correct stl_config.h, stl_iterator.h, stldebug.c
#endif
#ifdef __STL_NO_TYPENAME
# define __STL_NEED_TYPENAME	//*TY 08/01/1999 - should correct stl_config.h
#endif
#ifdef __STL_NO_EXPLICIT
# define __STL_NEED_EXPLICIT 	//*TY 08/01/1999 - should correct stl_config.h
#endif




