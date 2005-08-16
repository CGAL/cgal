/*=============================================================================
    Boost.Wave: A Standard compliant C++ preprocessor library

    Global application configuration
    
    http://www.boost.org/

    Copyright (c) 2001-2005 Hartmut Kaiser. Distributed under the Boost
    Software License, Version 1.0. (See accompanying file
    LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
=============================================================================*/

#if !defined(WAVE_CONFIG_HPP_F143F90A_A63F_4B27_AC41_9CA4F14F538D_INCLUDED)
#define WAVE_CONFIG_HPP_F143F90A_A63F_4B27_AC41_9CA4F14F538D_INCLUDED

#include <boost/config.hpp>
#include <boost/detail/workaround.hpp>
#include <boost/version.hpp>
#include <boost/spirit/version.hpp>
#include <boost/wave/wave_version.hpp>

///////////////////////////////////////////////////////////////////////////////
//  Define the maximal include nesting depth allowed. If this value isn't 
//  defined it defaults to 1024
//
//  To define a new initial include nesting depth uncomment the following and 
//  supply a new integer value.
//
#if !defined(BOOST_WAVE_MAX_INCLUDE_LEVEL_DEPTH)
#define BOOST_WAVE_MAX_INCLUDE_LEVEL_DEPTH 1024
#endif

///////////////////////////////////////////////////////////////////////////////
//  Decide, whether to support variadics and placemarkers
//
//  To implement support variadics and placemarkers uncomment the following
//
#if !defined(BOOST_WAVE_SUPPORT_VARIADICS_PLACEMARKERS)
#define BOOST_WAVE_SUPPORT_VARIADICS_PLACEMARKERS 1
#endif

///////////////////////////////////////////////////////////////////////////////
//  Decide, whether to implement a #warning directive as an extension to the 
//  C++ Standard (same as #error, but emits a warning, not an error)
//
//  To implement #warning directives, uncomment the following
//
#if !defined(BOOST_WAVE_SUPPORT_WARNING_DIRECTIVE)
#define BOOST_WAVE_SUPPORT_WARNING_DIRECTIVE 1
#endif

///////////////////////////////////////////////////////////////////////////////
//  Decide, whether to implement #pragma once 
//
//  To implement #pragma once, uncomment the following
//
#if !defined(BOOST_WAVE_SUPPORT_PRAGMA_ONCE)
#define BOOST_WAVE_SUPPORT_PRAGMA_ONCE 1
#endif

///////////////////////////////////////////////////////////////////////////////
//  Decide, whether to implement #include_next
//  Please note, that this is an extension to the C++ Standard.
//
//  To implement #include_next, uncomment the following
//
#if !defined(BOOST_WAVE_SUPPORT_INCLUDE_NEXT)
#define BOOST_WAVE_SUPPORT_INCLUDE_NEXT 1
#endif

///////////////////////////////////////////////////////////////////////////////
//  Undefine the following, to enable some MS specific language extensions:
//  __int8, __int16, __int32, __int64, __based, __declspec, __cdecl, 
//  __fastcall, __stdcall, __try, __except, __finally, __leave, __inline,
//  __asm, #region, #endregion
//
//  Note: By default this is enabled for Windows based systems, otherwise it's 
//        disabled.
#if !defined(BOOST_WAVE_SUPPORT_MS_EXTENSIONS)
#if defined(BOOST_WINDOWS)
#define BOOST_WAVE_SUPPORT_MS_EXTENSIONS 1
#else
#define BOOST_WAVE_SUPPORT_MS_EXTENSIONS 0
#endif
#endif

///////////////////////////////////////////////////////////////////////////////
//  Allow the message body of the #error and #warning directives to be 
//  preprocessed before the diagnostic is issued.
//
//  Uncommenting the following will preprocess the message bodies of #error and
//  #warning messages before the error (warning) is issued
//
#if !defined(BOOST_WAVE_PREPROCESS_ERROR_MESSAGE_BODY)
#define BOOST_WAVE_PREPROCESS_ERROR_MESSAGE_BODY 1
#endif

///////////////////////////////////////////////////////////////////////////////
//  Allow the #pragma directives to be returned to the caller (optionally after 
//  preprocessing the body) 
//
//  Uncommenting the following will skip #pragma directives, so that the caller
//  will not see them.
//
#if !defined(BOOST_WAVE_EMIT_PRAGMA_DIRECTIVES)
#define BOOST_WAVE_EMIT_PRAGMA_DIRECTIVES 1
#endif

///////////////////////////////////////////////////////////////////////////////
//  Allow the body of a #pragma directive to be preprocessed before the 
//  directive is returned to the caller.
//
//  Uncommenting the following will preprocess the bodies of #pragma directives
//
#if !defined(BOOST_WAVE_PREPROCESS_PRAGMA_BODY)
#define BOOST_WAVE_PREPROCESS_PRAGMA_BODY 1
#endif

///////////////////////////////////////////////////////////////////////////////
//  Allow to define macros with the command line syntax (-DMACRO(x)=definition)
//
//  Uncommenting the following will enable the possibility to define macros
//  based on the command line syntax
//
#if !defined(BOOST_WAVE_ENABLE_COMMANDLINE_MACROS)
#define BOOST_WAVE_ENABLE_COMMANDLINE_MACROS 1
#endif

///////////////////////////////////////////////////////////////////////////////
//  Define the string type to be used to store the token values and the file 
//  names inside a file_position template class
//
#if !defined(BOOST_WAVE_STRINGTYPE)

#if BOOST_WORKAROUND(BOOST_MSVC, <= 1300) || \
    BOOST_WORKAROUND(__MWERKS__, < 0x3200)
// VC7 isn't able to compile the flex_string class, fallback to std::string 
// CW upto 8.3 chokes as well *sigh*
#define BOOST_WAVE_STRINGTYPE std::string

#else
// use the following, if you have a fast std::allocator<char>
#define BOOST_WAVE_STRINGTYPE boost::wave::util::flex_string< \
        char, std::char_traits<char>, std::allocator<char>, \
        boost::wave::util::CowString</*char, */\
            boost::wave::util::AllocatorStringStorage<char> \
        > \
    > \
    /**/
    
/* #define BOOST_WAVE_STRINGTYPE boost::wave::util::flex_string< \
        char, std::char_traits<char>, boost::fast_pool_allocator<char>, \
        boost::wave::util::CowString<char, \
            boost::wave::util::AllocatorStringStorage<char, \
              boost::fast_pool_allocator<char> \
            > \
        > \
    > \
*/    /**/
    
//  This include is needed for the flex_string class used in the 
//  BOOST_WAVE_STRINGTYPE above.
#include <boost/wave/util/flex_string.hpp>

//  This include is needed for the boost::fast_allocator class used in the 
//  BOOST_WAVE_STRINGTYPE above.
//  Configure Boost.Pool thread support (for now: no thread support at all)
//#define BOOST_NO_MT
//#include <boost/pool/pool_alloc.hpp>

// Use the following, if you want to incorporate Maxim Yegorushkin's
// const_string library (http://sourceforge.net/projects/conststring/), which
// may be even faster, than using the flex_string class from above
//#define BOOST_WAVE_STRINGTYPE boost::const_string<char>
//
//#include <boost/const_string/const_string.hpp>
//#include <boost/const_string/io.hpp>
//#include <boost/const_string/concatenation.hpp>

#endif // BOOST_WORKAROUND(_MSC_VER, <= 1300)
#endif

///////////////////////////////////////////////////////////////////////////////
//  Uncomment the following, if you need debug output, the 
//  BOOST_SPIRIT_DEBUG_FLAGS_CPP constants below help to fine control the 
//  amount of the generated debug output.
//#define BOOST_SPIRIT_DEBUG

///////////////////////////////////////////////////////////////////////////////
//  Debug flags for the Wave library, possible flags spcified below.
//
//  Note: These flags take effect only if the BOOST_SPIRIT_DEBUG constant
//        above is defined as well.
#define BOOST_SPIRIT_DEBUG_FLAGS_CPP_GRAMMAR            0x0001
#define BOOST_SPIRIT_DEBUG_FLAGS_TIME_CONVERSION        0x0002
#define BOOST_SPIRIT_DEBUG_FLAGS_CPP_EXPR_GRAMMAR       0x0004
#define BOOST_SPIRIT_DEBUG_FLAGS_INTLIT_GRAMMAR         0x0008
#define BOOST_SPIRIT_DEBUG_FLAGS_CHLIT_GRAMMAR          0x0010
#define BOOST_SPIRIT_DEBUG_FLAGS_DEFINED_GRAMMAR        0x0020
#define BOOST_SPIRIT_DEBUG_FLAGS_PREDEF_MACROS_GRAMMAR  0x0040

#if !defined(BOOST_SPIRIT_DEBUG_FLAGS_CPP)
#define BOOST_SPIRIT_DEBUG_FLAGS_CPP    0    // default is no debugging
#endif 

///////////////////////////////////////////////////////////////////////////////
//
//  For all recognized preprocessor statements the output parse trees 
//  formatted as xml are printed. The formatted parse trees are streamed to the 
//  std::ostream defined by the WAVE_DUMP_PARSE_TREE_OUT constant.
//
//  Uncomment the following, if you want to see these parse trees. 
//
#if !defined(BOOST_WAVE_DUMP_PARSE_TREE)
#define BOOST_WAVE_DUMP_PARSE_TREE 0
#endif
#if BOOST_WAVE_DUMP_PARSE_TREE != 0 && !defined(BOOST_WAVE_DUMP_PARSE_TREE_OUT)
#define BOOST_WAVE_DUMP_PARSE_TREE_OUT std::cerr
#endif

///////////////////////////////////////////////////////////////////////////////
//
//  For all #if and #elif directives the preprocessed expressions are printed.
//  These expressions are streamed to the std::ostream defined by the 
//  BOOST_WAVE_DUMP_CONDITIONAL_EXPRESSIONS_OUT constant.
//
//  Uncomment the following, if you want to see the preprocessed expressions
//
#if !defined(BOOST_WAVE_DUMP_CONDITIONAL_EXPRESSIONS)
#define BOOST_WAVE_DUMP_CONDITIONAL_EXPRESSIONS 0
#endif
#if BOOST_WAVE_DUMP_CONDITIONAL_EXPRESSIONS != 0 && \
   !defined(BOOST_WAVE_DUMP_CONDITIONAL_EXPRESSIONS_OUT)
#define BOOST_WAVE_DUMP_CONDITIONAL_EXPRESSIONS_OUT std::cerr
#endif

///////////////////////////////////////////////////////////////////////////////
//  Decide, whether to use the separate compilation model for the instantiation 
//  of the C++ lexer objects.
//
//  If this is defined, you should explicitly instantiate the C++ lexer
//  template with the correct parameters in a separate compilation unit of
//  your program (see the file instantiate_re2c_lexer.cpp). 
//
//  To use the lexer inclusion model, uncomment the following 
//
#if !defined(BOOST_WAVE_SEPARATE_LEXER_INSTANTIATION)
#define BOOST_WAVE_SEPARATE_LEXER_INSTANTIATION 1
#endif

///////////////////////////////////////////////////////////////////////////////
//  Decide, whether to use the separate compilation model for the instantiation 
//  of the grammar objects.
//
//  If this is defined, you should explicitly instantiate the grammar
//  templates with the correct parameters in a separate compilation unit of
//  your program (see the files instantiate_cpp_grammar.cpp et.al.). 
//
//  To use the grammar inclusion model, uncomment the following 
//
#if !defined(BOOST_WAVE_SEPARATE_GRAMMAR_INSTANTIATION)
#define BOOST_WAVE_SEPARATE_GRAMMAR_INSTANTIATION 1
#endif

///////////////////////////////////////////////////////////////////////////////
//  configure Boost.Pool thread support (for now: no thread support at all)
#if !defined(BOOST_NO_MT)
#define BOOST_NO_MT
#endif // !defined(BOOST_NO_MT)

//#if !defined(BOOST_DISABLE_THREADS)
//#define BOOST_DISABLE_THREADS
//#endif // !defined(BOOST_DISABLE_THREADS)

///////////////////////////////////////////////////////////////////////////////
//  Wave needs at least 4 parameters for phoenix actors
#define PHOENIX_LIMIT 6

///////////////////////////////////////////////////////////////////////////////
//  Auto library naming
#if BOOST_VERSION >= 103100   
// auto link features work beginning from Boost V1.31.0
#if !defined(BOOST_WAVE_SOURCE) && !defined(BOOST_ALL_NO_LIB) && \
    !defined(BOOST_WAVE_NO_LIB)

#define BOOST_LIB_NAME boost_wave

// tell the auto-link code to select a dll when required:
#if defined(BOOST_ALL_DYN_LINK) || defined(BOOST_WAVE_DYN_LINK)
#define BOOST_DYN_LINK
#endif

#include <boost/config/auto_link.hpp>

#endif  // auto-linking disabled
#endif  // BOOST_VERSION

#endif // !defined(WAVE_CONFIG_HPP_F143F90A_A63F_4B27_AC41_9CA4F14F538D_INCLUDED)
