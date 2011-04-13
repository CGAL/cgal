/*
 *
 * Copyright (c) 1994
 * Hewlett-Packard Company
 *
 * Copyright (c) 1996,1997
 * Silicon Graphics Computer Systems, Inc.
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


#ifndef __SGI_STL_INTERNAL_NUMERIC_H
#define __SGI_STL_INTERNAL_NUMERIC_H

#ifndef __STL_CONFIG_H
#include <stl_config.h>
#endif

#if !defined (__STLPORT_DEBUG_H) && (defined  (__STL_DEBUG) || defined (__STL_ASSERTIONS))
# include <stldebug.h>
#endif

__STL_BEGIN_NAMESPACE

template <class _InputIterator, class _Tp>
__STL_INLINE_LOOP
_Tp accumulate(_InputIterator __first, _InputIterator __last, _Tp __init)
{
  __stl_debug_check(__check_range(__first, __last));
  for ( ; __first != __last; ++__first)
    __init = __init + *__first;
  return __init;
}

template <class _InputIterator, class _Tp, class _BinaryOperation>
__STL_INLINE_LOOP
_Tp accumulate(_InputIterator __first, _InputIterator __last, _Tp __init,
               _BinaryOperation __binary_op)
{
  __stl_debug_check(__check_range(__first, __last));
  for ( ; __first != __last; ++__first)
    __init = __binary_op(__init, *__first);
  return __init;
}

template <class _InputIterator1, class _InputIterator2, class _Tp>
__STL_INLINE_LOOP
_Tp inner_product(_InputIterator1 __first1, _InputIterator1 __last1,
                  _InputIterator2 __first2, _Tp __init)
{
  __stl_debug_check(__check_range(__first1, __last1));
  for ( ; __first1 != __last1; ++__first1, ++__first2)
    __init = __init + (*__first1 * *__first2);
  return __init;
}

template <class _InputIterator1, class _InputIterator2, class _Tp,
          class _BinaryOperation1, class _BinaryOperation2>
__STL_INLINE_LOOP
_Tp inner_product(_InputIterator1 __first1, _InputIterator1 __last1,
                  _InputIterator2 __first2, _Tp __init, 
                  _BinaryOperation1 __binary_op1,
                  _BinaryOperation2 __binary_op2)
{
  __stl_debug_check(__check_range(__first1, __last1));
  for ( ; __first1 != __last1; ++__first1, ++__first2)
    __init = __binary_op1(__init, __binary_op2(*__first1, *__first2));
  return __init;
}

template <class _InputIterator, class _OutputIterator>
_OutputIterator 
partial_sum(_InputIterator __first, _InputIterator __last,
            _OutputIterator __result);


template <class _InputIterator, class _OutputIterator, class _BinaryOperation>
_OutputIterator 
partial_sum(_InputIterator __first, _InputIterator __last,
            _OutputIterator __result, _BinaryOperation __binary_op);

template <class _InputIterator, class _OutputIterator>
_OutputIterator
adjacent_difference(_InputIterator __first,
                    _InputIterator __last, _OutputIterator __result);

template <class _InputIterator, class _OutputIterator, class _BinaryOperation>
_OutputIterator 
adjacent_difference(_InputIterator __first, _InputIterator __last,
                    _OutputIterator __result, _BinaryOperation __binary_op);


// Returns __x ** __n, where __n >= 0.  _Note that "multiplication"
// is required to be associative, but not necessarily commutative.
// Alias for the internal name __power.  Note that power is an extension,
// not part of the C++ standard.

template <class _Tp, class _Integer, class _MonoidOperation>
_Tp power(_Tp __x, _Integer __n, _MonoidOperation __opr);

template <class _Tp, class _Integer>
_Tp power(_Tp __x, _Integer __n);

// iota is not part of the C++ standard.  It is an extension.

template <class _ForwardIterator, class _Tp>
__STL_INLINE_LOOP
void 
iota(_ForwardIterator __first, _ForwardIterator __last, _Tp __value)
{
  while (__first != __last)
    *__first++ = __value++;
}

__STL_END_NAMESPACE

# if !defined (__STL_LINK_TIME_INSTANTIATION)
#  include <stl_numeric.c>
# endif

#endif /* __SGI_STL_INTERNAL_NUMERIC_H */

// Local Variables:
// mode:C++
// End:
