/*
 *
 * Copyright (c) 1994
 * Hewlett-Packard Company
 *
 * Copyright (c) 1996,1997
 * Silicon Graphics Computer Systems, Inc.
 *
 * Copyright (c) 1997
 * Moscow Center for SPARC Technology
 *
 * Copyright (c) 1999 
 * Boris Fomitchev
 *
 * This material is provided "as is", with absolutely no warranty expressed
 * or implied. Any use is at your own risk.
 *
 * Permission to use or copy this software for any purpose is hereby granted 
 * without fee, provided the above notices are retained on all copies.
 * Permission to modify the code and to distribute modified code is granted,
 * provided the above notices are retained, and a notice that the code was
 * modified is included with the above copyright notice.
 *
 */

/* NOTE: This is an internal header file, included by other STL headers.
 *   You should not attempt to use it directly.
 */

#ifndef __SGI_STL_INTERNAL_CONSTRUCT_H
#define __SGI_STL_INTERNAL_CONSTRUCT_H

# if defined (__STL_DEBUG_UNINITIALIZED) && ! defined (__STLPORT_CSTRING)
# include <cstring>
# endif

# ifndef __STLPORT_NEW
#  include <new>
# endif

# ifndef __TYPE_TRAITS_H
#  include <type_traits.h>
# endif

#ifndef __SGI_STL_INTERNAL_ITERATOR_BASE_H
# include <stl_iterator_base.h>
#endif

__STL_BEGIN_NAMESPACE

# ifdef __STL_TRIVIAL_DESTRUCTOR_BUG
template <class _Tp>
inline void __destroy_aux(_Tp* __pointer, __false_type) { __pointer->~_Tp(); }
template <class _Tp>
inline void __destroy_aux(_Tp* __pointer, __true_type) {}
# endif

template <class _Tp>
inline void destroy(_Tp* __pointer) {
# if _MSC_VER >= 1010
  __pointer;
# endif	// _MSC_VER >= 1000
# ifdef __STL_TRIVIAL_DESTRUCTOR_BUG
  typedef typename __type_traits<_Tp>::has_trivial_destructor
          _Trivial_destructor;
  __destroy_aux(__pointer, _Trivial_destructor());
# else
#  if ( defined (__BORLANDC__) && ( __BORLANDC__ < 0x500 ) )
    __pointer->_Tp::~_Tp();
#  else
    __pointer->~_Tp();
#  endif
# endif
# ifdef __STL_DEBUG_UNINITIALIZED
	memset((char*)__pointer, (char)__STL_SHRED_BYTE, sizeof(_Tp));
# endif
}

# if defined (new)
#   define __STL_NEW_REDEFINE new
#   undef new
# endif 

template <class _T1, class _T2>
inline void construct(_T1* __p, const _T2& __value) {
# ifdef __STL_DEBUG_UNINITIALIZED
	memset((char*)__p, (char)__STL_SHRED_BYTE, sizeof(_T1));
# endif
    __STL_PLACEMENT_NEW (__p) _T1(__value);
}

template <class _T1>
inline void construct(_T1* __p) {
# ifdef __STL_DEBUG_UNINITIALIZED
	memset((char*)__p, (char)__STL_SHRED_BYTE, sizeof(_T1));
# endif
  __STL_PLACEMENT_NEW (__p) _T1();
}

# if defined(__STL_NEW_REDEFINE)
# ifdef DEBUG_NEW
#  define new DEBUG_NEW
# endif
#  undef __STL_NEW_REDEFINE
# endif 

template <class _ForwardIterator>
__STL_INLINE_LOOP void
__destroy_aux(_ForwardIterator __first, _ForwardIterator __last, __false_type)
{
  for ( ; __first != __last; ++__first)
    destroy(&*__first);
}

template <class _ForwardIterator> 
inline void __destroy_aux(_ForwardIterator, _ForwardIterator, __true_type) {}

template <class _ForwardIterator, class _Tp>
inline void 
__destroy(_ForwardIterator __first, _ForwardIterator __last, _Tp*)
{
  typedef typename __type_traits<_Tp>::has_trivial_destructor
          _Trivial_destructor;
  __destroy_aux(__first, __last, _Trivial_destructor());
}

template <class _ForwardIterator>
inline void destroy(_ForwardIterator __first, _ForwardIterator __last) {
  __destroy(__first, __last, __VALUE_TYPE(__first));
}

inline void destroy(char*, char*) {}
# ifdef __STL_HAS_WCHAR_T // dwa 8/15/97
inline void destroy(wchar_t*, wchar_t*) {}
inline void destroy(const wchar_t*, const wchar_t*) {}
# endif

__STL_END_NAMESPACE

#endif /* __SGI_STL_INTERNAL_CONSTRUCT_H */

// Local Variables:
// mode:C++
// End:
