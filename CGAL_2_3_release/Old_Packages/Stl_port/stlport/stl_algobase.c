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
#ifndef __SGI_STL_ALGOBASE_C
#define __SGI_STL_ALGOBASE_C

__STL_BEGIN_NAMESPACE

template <class _InputIter1, class _InputIter2>
bool lexicographical_compare(_InputIter1 __first1, _InputIter1 __last1,
                             _InputIter2 __first2, _InputIter2 __last2) {
  for ( ; __first1 != __last1 && __first2 != __last2
        ; ++__first1, ++__first2) {
    if (*__first1 < *__first2)
      return true;
    if (*__first2 < *__first1)
      return false;
  }
  return __first1 == __last1 && __first2 != __last2;
}

template <class _InputIter1, class _InputIter2, class _Compare>
bool lexicographical_compare(_InputIter1 __first1, _InputIter1 __last1,
                             _InputIter2 __first2, _InputIter2 __last2,
                             _Compare __comp) {
  for ( ; __first1 != __last1 && __first2 != __last2
        ; ++__first1, ++__first2) {
    if (__comp(*__first1, *__first2))
      return true;
    if (__comp(*__first2, *__first1))
      return false;
  }
  return __first1 == __last1 && __first2 != __last2;
}

template <class _InputIter1, class _InputIter2>
int __lexicographical_compare_3way(_InputIter1 __first1, _InputIter1 __last1,
                                   _InputIter2 __first2, _InputIter2 __last2)
{
  while (__first1 != __last1 && __first2 != __last2) {
    if (*__first1 < *__first2)
      return -1;
    if (*__first2 < *__first1)
      return 1;
    ++__first1;
    ++__first2;
  }
  if (__first2 == __last2) {
    return !(__first1 == __last1);
  }
  else {
    return -1;
  }
}


template <class _InputIter1, class _InputIter2>
int lexicographical_compare_3way(_InputIter1 __first1, _InputIter1 __last1,
                                 _InputIter2 __first2, _InputIter2 __last2)
{
  return __lexicographical_compare_3way(__first1, __last1, __first2, __last2);
}

__STL_END_NAMESPACE

#endif /* __SGI_STL_ALGOBASE_C */

// Local Variables:
// mode:C++
// End:
