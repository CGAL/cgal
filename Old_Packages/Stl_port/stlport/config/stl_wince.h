/*
 * File to have Windows CE Toolkit for VC++ 5.0 working with STLport
 * 09 - 03 - 1999
 * Origin : Giuseppe Govi - g.govi@iol.it
 */

#ifndef __STL_WINCE_H
#define __STL_WINCE_H

// this flag is being used by STLport
#   define __STL_WINCE

// tell other parts no iostreams are desired
#   define __STL_NO_IOSTREAMS 1

#     undef __STL_HAS_NO_EXCEPTIONS
#     define __STL_HAS_NO_EXCEPTIONS
#     undef __STL_NO_EXCEPTION_HEADER
#     define __STL_NO_EXCEPTION_HEADER

# ifdef __STL_MSVC
#     pragma warning (disable: 4786)
# endif

#ifdef __STL_WINCE_USE_OUTPUTDEBUGSTRING
#define __STL_WINCE_TRACE(msg)   OutputDebugString(msg)
#else
#define __STL_WINCE_TRACE(msg)   MessageBox(NULL,(msg),NULL,MB_OK)
#endif

#ifndef __THROW_BAD_ALLOC
#define __THROW_BAD_ALLOC __STL_WINCE_TRACE(L"out of memory"); ExitThread(1)
#endif

#ifndef _SIZE_T_DEFINED
typedef unsigned int size_t;
#define _SIZE_T_DEFINED
#endif

#ifndef __PLACEMENT_NEW_INLINE
inline void *__cdecl operator new(size_t, void *_P) { return (_P); }
#define __PLACEMENT_NEW_INLINE
#endif

#ifndef _WCHAR_T_DEFINED
typedef unsigned short wchar_t;
#define _WCHAR_T_DEFINED
#endif

//ptrdiff_t is not defined in Windows CE SDK
#ifndef _PTRDIFF_T_DEFINED
typedef int ptrdiff_t;
#define _PTRDIFF_T_DEFINED
#endif

#endif /* __STL_WCE_H */


