/*
 *
 * Copyright (c) 1998-2002
 * Dr John Maddock
 *
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 *
 */

 /*
  *   LOCATION:    see http://www.boost.org for most recent version.
  *   FILE         config.hpp
  *   VERSION      see <boost/version.hpp>
  *   DESCRIPTION: regex extended config setup.
  */

#ifndef BOOST_REGEX_CONFIG_HPP
#define BOOST_REGEX_CONFIG_HPP
/*
 Borland C++ Fix/error check
 this has to go *before* we include any std lib headers:
*/
#if defined(__BORLANDC__)
#  include <boost/regex/config/borland.hpp>
#endif

/*****************************************************************************
 *
 *  Include all the headers we need here:
 *
 ****************************************************************************/

#ifdef __cplusplus

#  ifndef BOOST_REGEX_USER_CONFIG
#     define BOOST_REGEX_USER_CONFIG <boost/regex/user.hpp>
#  endif

#  include BOOST_REGEX_USER_CONFIG

#  include <cstdlib>
#  include <cstddef>
#  include <cstdio>
#  include <clocale>
#  include <cassert>
#  include <string>
#  include <stdexcept>
#  include <iterator>
#  include <iosfwd>
#  include <vector>
#  include <boost/config.hpp>
#  include <boost/cstdint.hpp>
#  include <boost/regex/config/allocator.hpp>
#  include <boost/regex/config/cstring.hpp>
#  include <boost/throw_exception.hpp>
#  include <boost/scoped_ptr.hpp>
#  ifndef BOOST_NO_STD_LOCALE
#     include <locale>
#  endif
#else
   /*
   * C build,
   * don't include <boost/config.hpp> because that may
   * do C++ specific things in future...
   */
#  include <stdlib.h>
#  include <stddef.h>
#  ifdef _MSC_VER
#     define BOOST_MSVC _MSC_VER
#  endif
#endif

/*****************************************************************************
 *
 *  Boilerplate regex config options:
 *
 ****************************************************************************/

/* Obsolete macro, use BOOST_VERSION instead: */
#define BOOST_RE_VERSION 320

/* fix: */
#if defined(_UNICODE) && !defined(UNICODE)
#define UNICODE
#endif

/*
* If there isn't good enough wide character support then there will
* be no wide character regular expressions:
*/
#if (defined(BOOST_NO_CWCHAR) || defined(BOOST_NO_CWCTYPE) || defined(BOOST_NO_STD_WSTRING))
#  if !defined(BOOST_NO_WREGEX)
#     define BOOST_NO_WREGEX
#  endif
#else
#  if defined(__sgi) && (defined(__SGI_STL_PORT) || defined(_STLPORT_VERSION))
      /* STLPort on IRIX is misconfigured: <cwctype> does not compile
      * as a temporary fix include <wctype.h> instead and prevent inclusion
      * of STLPort version of <cwctype> */
#     include <wctype.h>
#     define __STLPORT_CWCTYPE
#     define _STLP_CWCTYPE
#  endif

#ifdef __cplusplus
#  include <boost/regex/config/cwchar.hpp>
#endif

#endif

/*
* If Win32 support has been disabled for boost in general, then
* it is for regex in particular:
*/
#if defined(BOOST_DISABLE_WIN32) && !defined(BOOST_REGEX_NO_W32)
#  define BOOST_REGEX_NO_W32
#endif

/* some versions of gcc can't merge template instances: */
#if defined(__CYGWIN__)
#  define BOOST_REGEX_NO_TEMPLATE_SWITCH_MERGE
#endif

/* fix problems with bool as a macro,
* this probably doesn't affect any current compilers: */
#if defined(bool) || defined(true) || defined(false)
#  define BOOST_REGEX_NO_BOOL
#endif

/* We don't make our templates external if the compiler
 can't handle it: */
#if (defined(BOOST_NO_MEMBER_FUNCTION_SPECIALIZATIONS) || defined(__HP_aCC) || defined(__MWERKS__) || defined(__COMO__) || defined(BOOST_INTEL))\
   && !defined(BOOST_MSVC) && !defined(__BORLANDC__)
#  define BOOST_REGEX_NO_EXTERNAL_TEMPLATES
#endif

/* disable our own file-iterators and mapfiles if we can't
 support them: */
#if !defined(BOOST_HAS_DIRENT_H) && !(defined(_WIN32) && !defined(BOOST_REGEX_NO_W32))
#  define BOOST_REGEX_NO_FILEITER
#endif

#ifdef __cplusplus
#ifndef MB_CUR_MAX
// yuk!
// better make a conservative guess!
#define MB_CUR_MAX 10
#endif

namespace boost{ namespace re_detail{
#ifdef BOOST_NO_STD_DISTANCE
template <class T>
std::ptrdiff_t distance(const T& x, const T& y)
{ return y - x; }
#else
using std::distance;
#endif
}}


#ifdef BOOST_REGEX_NO_BOOL
#  define BOOST_REGEX_MAKE_BOOL(x) static_cast<bool>((x) ? true : false)
#else
#  ifdef BOOST_MSVC
      // warning suppression with VC6:
#     pragma warning(disable: 4800)
#  endif
#  define BOOST_REGEX_MAKE_BOOL(x) static_cast<bool>(x)
#endif
#endif /* __cplusplus */

/* backwards compatibitity: */
#if defined(BOOST_RE_NO_LIB)
#  define BOOST_REGEX_NO_LIB
#endif

#if defined(__GNUC__) && (defined(_WIN32) || defined(__CYGWIN__))
// gcc on win32 has problems merging switch statements in templates:
#  define BOOST_REGEX_NO_TEMPLATE_SWITCH_MERGE
// gcc on win32 has problems if you include <windows.h>
// (sporadically generates bad code).
#  define BOOST_REGEX_USE_C_LOCALE
#  define BOOST_REGEX_NO_W32
#endif
#if defined(__COMO__) && !defined(BOOST_REGEX_NO_W32) && !defined(_MSC_EXTENSIONS)
#  define BOOST_REGEX_NO_W32
#endif

/*****************************************************************************
 *
 *  Wide character workarounds:
 *
 ****************************************************************************/

#ifdef __cplusplus
#if defined(BOOST_MSVC) && (BOOST_MSVC >= 1300) && !defined(BOOST_REGEX_V3) && !(defined(__SGI_STL_PORT) || defined(_STLPORT_VERSION))
#  define BOOST_REGEX_HAS_SHORT_WCHAR_T
namespace boost{ typedef __wchar_t regex_wchar_type; }
#else
namespace boost{ typedef wchar_t regex_wchar_type; }
#endif
#endif


/*****************************************************************************
 *
 *  Set up dll import/export options:
 *
 ****************************************************************************/

#if defined(BOOST_HAS_DECLSPEC) && (defined(BOOST_REGEX_DYN_LINK) || defined(BOOST_ALL_DYN_LINK)) && !defined(BOOST_REGEX_STATIC_LINK)
#  if defined(BOOST_REGEX_SOURCE)
#     define BOOST_REGEX_DECL __declspec(dllexport)
#     define BOOST_REGEX_BUILD_DLL
#  else
#     define BOOST_REGEX_DECL __declspec(dllimport)
#  endif
#endif

#ifndef BOOST_REGEX_DECL
#  define BOOST_REGEX_DECL
#endif

#if !defined(BOOST_REGEX_NO_LIB) && !defined(BOOST_REGEX_SOURCE) && !defined(BOOST_ALL_NO_LIB) && defined(__cplusplus)
#  define BOOST_LIB_NAME boost_regex
#  if defined(BOOST_REGEX_DYN_LINK) || defined(BOOST_ALL_DYN_LINK)
#     define BOOST_DYN_LINK
#  endif
#ifdef BOOST_REGEX_DIAG
#  define BOOST_LIB_DIAGNOSTIC
#endif
#  include <boost/config/auto_link.hpp>
#endif

/*****************************************************************************
 *
 *  Set up function call type:
 *
 ****************************************************************************/

#if defined(BOOST_MSVC) && (BOOST_MSVC >= 1200) && defined(_MSC_EXTENSIONS)
#if defined(_DEBUG) || defined(__MSVC_RUNTIME_CHECKS)
#  define BOOST_REGEX_CALL __cdecl
#else
#  define BOOST_REGEX_CALL __fastcall
#endif
#  define BOOST_REGEX_CCALL __cdecl
#endif

#if defined(__BORLANDC__) && !defined(BOOST_DISABLE_WIN32)
#  define BOOST_REGEX_CALL __fastcall
#  define BOOST_REGEX_CCALL __stdcall
#endif

#ifndef BOOST_REGEX_CALL
#  define BOOST_REGEX_CALL
#endif
#ifndef BOOST_REGEX_CCALL
#define BOOST_REGEX_CCALL
#endif

/*****************************************************************************
 *
 *  Set up localisation model:
 *
 ****************************************************************************/

/* backwards compatibility: */
#ifdef BOOST_RE_LOCALE_C
#  define BOOST_REGEX_USE_C_LOCALE
#endif

#ifdef BOOST_RE_LOCALE_CPP
#  define BOOST_REGEX_USE_CPP_LOCALE
#endif

/* Win32 defaults to native Win32 locale: */
#if defined(_WIN32) && !defined(BOOST_REGEX_USE_WIN32_LOCALE) && !defined(BOOST_REGEX_USE_C_LOCALE) && !defined(BOOST_REGEX_USE_CPP_LOCALE) && !defined(BOOST_REGEX_NO_W32)
#  define BOOST_REGEX_USE_WIN32_LOCALE
#endif
/* otherwise use C locale: */
#if !defined(BOOST_REGEX_USE_WIN32_LOCALE) && !defined(BOOST_REGEX_USE_C_LOCALE) && !defined(BOOST_REGEX_USE_CPP_LOCALE)
#  define BOOST_REGEX_USE_C_LOCALE
#endif

#if defined(_WIN32) && !defined(BOOST_REGEX_NO_W32)
#  include <windows.h>
#endif

#ifdef MAXPATH
#  define BOOST_REGEX_MAX_PATH MAXPATH
#elif defined(MAX_PATH)
#  define BOOST_REGEX_MAX_PATH MAX_PATH
#elif defined(FILENAME_MAX)
#  define BOOST_REGEX_MAX_PATH FILENAME_MAX
#else
#  define BOOST_REGEX_MAX_PATH 200
#endif

#ifndef BOOST_REGEX_MAX_STATE_COUNT
#  define BOOST_REGEX_MAX_STATE_COUNT 100000000
#endif


/*****************************************************************************
 *
 *  Error Handling for exception free compilers:
 *
 ****************************************************************************/

#ifdef BOOST_NO_EXCEPTIONS
/*
*  If there are no exceptions then we must report critical-errors
*  the only way we know how; by terminating.
*/
#  define BOOST_REGEX_NOEH_ASSERT(x)\
if(0 == (x))\
{\
   std::string s("Error: critical regex++ failure in: ");\
   s.append(#x);\
   std::runtime_error e(s);\
   boost::throw_exception(e);\
}
#else
/*
*  With exceptions then error handling is taken care of and
*  there is no need for these checks:
*/
#  define BOOST_REGEX_NOEH_ASSERT(x)
#endif

/*****************************************************************************
 *
 *  Debugging / tracing support:
 *
 ****************************************************************************/

#if defined(BOOST_REGEX_DEBUG) && defined(__cplusplus)

#  include <iostream>
using std::cout;
using std::cin;
using std::cerr;
using std::endl;
using std::hex;
using std::dec;

#  ifndef jm_assert
#     define jm_assert(x) assert(x)
#  endif
#  ifndef jm_trace
#     define jm_trace(x) cerr << x << endl;
#  endif
#  ifndef jm_instrument
#     define jm_instrument jm_trace(__FILE__<<"#"<<__LINE__)
#  endif

namespace boost{
   namespace re_detail{
class debug_guard
{
public:
   char g1[32];
   const char* pc;
   char* pnc;
   const char* file;
   int line;
   char g2[32];
   debug_guard(const char* f, int l, const char* p1 = 0, char* p2 = 0);
   ~debug_guard();
};

#  define BOOST_RE_GUARD_STACK boost::re_detail::debug_guard sg(__FILE__, __LINE__);
#  define BOOST_RE_GUARD_GLOBAL(x) const char g1##x[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, }; char g2##x[32]; boost::debug_guard g3##x(__FILE__, __LINE__, g1##x, g2##x);

   } // namespace re_detail
} // namespace boost

#else

#  define jm_assert(x)
#  define jm_trace(x)
#  define BOOST_RE_GUARD_STACK
#  define BOOST_RE_GUARD_GLOBAL(x)
#  ifndef jm_instrument
#     define jm_instrument
#  endif
#endif

/*****************************************************************************
 *
 *  Stack protection under MS Windows:
 *
 ****************************************************************************/

#if !defined(BOOST_REGEX_NO_W32) && !defined(BOOST_REGEX_V3)
#  if(defined(_WIN32) || defined(_WIN64) || defined(_WINCE)) \
        && !defined(__GNUC__) \
        && !(defined(__BORLANDC__) && (__BORLANDC__ >= 0x600)) \
        && !(defined(__MWERKS__) && (__MWERKS__ <= 0x3003))
#     define BOOST_REGEX_HAS_MS_STACK_GUARD
#  endif
#elif defined(BOOST_REGEX_HAS_MS_STACK_GUARD)
#  undef BOOST_REGEX_HAS_MS_STACK_GUARD
#endif

#if defined(__cplusplus) && defined(BOOST_REGEX_HAS_MS_STACK_GUARD)

namespace boost{
namespace re_detail{

BOOST_REGEX_DECL void BOOST_REGEX_CALL reset_stack_guard_page();

}
}

#endif


/*****************************************************************************
 *
 *  Error handling:
 *
 ****************************************************************************/

#if defined(__cplusplus)

namespace boost{
namespace re_detail{

BOOST_REGEX_DECL void BOOST_REGEX_CALL raise_regex_exception(const std::string& s);

template <class traits>
void raise_error(const traits& t, unsigned code)
{
   (void)t;  // warning suppression
   raise_regex_exception(t.error_string(code));
}

}
}

#endif

/*****************************************************************************
 *
 *  Algorithm selection and configuration:
 *
 ****************************************************************************/

#if !defined(BOOST_REGEX_RECURSIVE) && !defined(BOOST_REGEX_NON_RECURSIVE)
#  if defined(BOOST_REGEX_HAS_MS_STACK_GUARD) && !defined(_STLP_DEBUG) && !defined(__STL_DEBUG)
#     define BOOST_REGEX_RECURSIVE
#  else
#     define BOOST_REGEX_NON_RECURSIVE
#  endif
#endif

#ifdef BOOST_REGEX_NON_RECURSIVE
#  ifdef BOOST_REGEX_RECURSIVE
#     error "Can't set both BOOST_REGEX_RECURSIVE and BOOST_REGEX_NON_RECURSIVE"
#  endif
#  ifndef BOOST_REGEX_BLOCKSIZE
#     define BOOST_REGEX_BLOCKSIZE 4096
#  endif
#  if BOOST_REGEX_BLOCKSIZE < 512
#     error "BOOST_REGEX_BLOCKSIZE must be at least 512"
#  endif
#  ifndef BOOST_REGEX_MAX_BLOCKS
#     define BOOST_REGEX_MAX_BLOCKS 1024
#  endif
#  ifdef BOOST_REGEX_HAS_MS_STACK_GUARD
#     undef BOOST_REGEX_HAS_MS_STACK_GUARD
#  endif
#  ifndef BOOST_REGEX_MAX_CACHE_BLOCKS
#     define BOOST_REGEX_MAX_CACHE_BLOCKS 16
#  endif
#endif


/*****************************************************************************
 *
 *  Fix broken compilers that wrongly #define some symbols:
 *
 ****************************************************************************/

#ifdef __cplusplus

// the following may be defined as macros; this is
// incompatable with std::something syntax, we have
// no choice but to undef them?

#ifdef sprintf
#undef sprintf
#endif
#ifdef swprintf
#undef swprintf
#endif
#endif

/*****************************************************************************
 *
 *  Fix broken broken namespace support:
 *
 ****************************************************************************/

#if defined(BOOST_NO_STDC_NAMESPACE) && defined(__cplusplus)

namespace std{
   using ::ptrdiff_t;
   using ::size_t;
   using ::sprintf;
   using ::abs;
   using ::setlocale;
#  ifndef BOOST_NO_WREGEX
#     ifndef BOOST_NO_SWPRINTF
   using ::swprintf;
#     endif
   using ::wcstombs;
   using ::mbstowcs;
#     if !defined(BOOST_NO_STD_LOCALE) && !defined (__STL_NO_NATIVE_MBSTATE_T) && !defined(_STLP_NO_NATIVE_MBSTATE_T)
   using ::mbstate_t;
#     endif
#  endif /* BOOST_NO_WREGEX */
   using ::fseek;
   using ::fread;
   using ::ftell;
   using ::fopen;
   using ::fclose;
   using ::FILE;
#ifdef BOOST_NO_EXCEPTIONS
   using ::fprintf;
   using ::abort;
#endif
}

#endif

/*****************************************************************************
 *
 *  helper functions pointer_construct/pointer_destroy:
 *
 ****************************************************************************/

#ifdef __cplusplus
namespace boost{ namespace re_detail{

#ifdef BOOST_MSVC
#pragma warning (push)
#pragma warning (disable : 4100)
#endif

template <class T>
inline void pointer_destroy(T* p)
{ p->~T(); (void)p; }

#ifdef BOOST_MSVC
#pragma warning (pop)
#endif

template <class T>
inline void pointer_construct(T* p, const T& t)
{ new (p) T(t); }

}} // namespaces
#endif

/*****************************************************************************
 *
 *  helper memory allocation functions:
 *
 ****************************************************************************/

#if defined(__cplusplus) && defined(BOOST_REGEX_NON_RECURSIVE)
namespace boost{ namespace re_detail{

BOOST_REGEX_DECL void* BOOST_REGEX_CALL get_mem_block();
BOOST_REGEX_DECL void BOOST_REGEX_CALL put_mem_block(void*);

}} // namespaces
#endif

/*****************************************************************************
 *
 *  Diagnostics:
 *
 ****************************************************************************/

#ifdef BOOST_REGEX_CONFIG_INFO
BOOST_REGEX_DECL void BOOST_REGEX_CALL print_regex_library_info();
#endif

#if defined(BOOST_REGEX_DIAG)
#  pragma message ("BOOST_REGEX_DECL set as: " BOOST_STRINGIZE(BOOST_REGEX_DECL))
#  pragma message ("BOOST_REGEX_CALL set as: " BOOST_STRINGIZE(BOOST_REGEX_CALL))
#  pragma message ("BOOST_REGEX_CCALL set as: " BOOST_STRINGIZE(BOOST_REGEX_CCALL))
#ifdef BOOST_REGEX_USE_C_LOCALE
#  pragma message ("Using C locale in regex traits class")
#elif BOOST_REGEX_USE_CPP_LOCALE
#  pragma message ("Using C++ locale in regex traits class")
#else
#  pragma message ("Using Win32 locale in regex traits class")
#endif
#ifdef BOOST_REGEX_DYN_LINK
#  pragma message ("Dynamic linking enabled")
#endif
#ifdef BOOST_REGEX_NO_LIB
#  pragma message ("Auto-linking disabled")
#endif
#ifdef BOOST_REGEX_NO_EXTERNAL_TEMPLATES
#  pragma message ("Extern templates disabled")
#endif
#ifdef BOOST_REGEX_V3
#  pragma message ("Using Version 3 regex code")
#else
#  pragma message ("Using Version 4 regex code")
#endif

#endif

#endif




