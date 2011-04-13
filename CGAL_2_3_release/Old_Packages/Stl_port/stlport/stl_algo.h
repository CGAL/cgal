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

#ifndef __SGI_STL_INTERNAL_ALGO_H
#define __SGI_STL_INTERNAL_ALGO_H


#if defined(__sgi) && !defined(__GNUC__) && (_MIPS_SIM != _MIPS_SIM_ABI32)
#pragma set woff 1209
#endif

# ifndef __SGI_STL_INTERNAL_ALGOBASE_H
#  include <stl_algobase.h>
# endif

# ifndef __SGI_STL_INTERNAL_TEMPBUF_H
#  include <stl_tempbuf.h>
# endif

# ifndef __SGI_STL_INTERNAL_HEAP_H
#  include <stl_heap.h>
# endif

# ifndef __SGI_STL_INTERNAL_ITERATOR_H
#  include <stl_iterator.h>
# endif

__STL_BEGIN_NAMESPACE

// for_each.  Apply a function to every element of a range.
template <class _InputIter, class _Function>
_Function for_each(_InputIter __first, _InputIter __last, _Function __f);

// find and find_if.

template <class _InputIter, class _Tp>
inline _InputIter find(_InputIter __first, _InputIter __last,
                       const _Tp& __val,
                       input_iterator_tag)
{
  __stl_debug_check(__check_range(__first, __last));
  while (__first != __last && !(*__first == __val))
    ++__first;
  return __first;
}

template <class _InputIter, class _Predicate>
inline _InputIter find_if(_InputIter __first, __STL_MPW_EXTRA_CONST _InputIter __last,
                          _Predicate __pred,
                          input_iterator_tag)
{
  __stl_debug_check(__check_range(__first, __last));
  while (__first != __last && !__pred(*__first))
    ++__first;
  return __first;
}

#ifdef __STL_CLASS_PARTIAL_SPECIALIZATION

template <class _RandomAccessIter, class _Tp>
_RandomAccessIter find(_RandomAccessIter __first, _RandomAccessIter __last,
                       const _Tp& __val,
                       random_access_iterator_tag);

template <class _RandomAccessIter, class _Predicate>
_RandomAccessIter find_if(_RandomAccessIter __first, _RandomAccessIter __last,
                          _Predicate __pred,
                          random_access_iterator_tag);

#endif /* __STL_CLASS_PARTIAL_SPECIALIZATION */

// fbp : without partial spec, we essentially have one variant for those.
template <class _InputIter, class _Tp>
inline _InputIter find(_InputIter __first, _InputIter __last,
                       const _Tp& __val)
{
# if defined (__STL_CLASS_PARTIAL_SPECIALIZATION)
  return find(__first, __last, __val, __ITERATOR_CATEGORY(__first));
# else
  return find(__first, __last, __val, input_iterator_tag());
# endif /* __STL_CLASS_PARTIAL_SPECIALIZATION */
}

template <class _InputIter, class _Predicate>
inline _InputIter find_if(_InputIter __first, _InputIter __last,
                          _Predicate __pred) {
# if defined (__STL_CLASS_PARTIAL_SPECIALIZATION)
  return find_if(__first, __last, __pred, __ITERATOR_CATEGORY(__first));
# else
  return find_if(__first, __last, __pred, input_iterator_tag());
# endif /* __STL_CLASS_PARTIAL_SPECIALIZATION */
}

// adjacent_find.

template <class _ForwardIter>
_ForwardIter adjacent_find(_ForwardIter __first, _ForwardIter __last);
template <class _ForwardIter, class _BinaryPredicate>
_ForwardIter adjacent_find(_ForwardIter __first, _ForwardIter __last,
                           _BinaryPredicate __binary_pred);

// count and count_if.  There are two version of each, one whose return type
// type is void and one (present only if we have partial specialization)
// whose return type is iterator_traits<_InputIter>::difference_type.  The
// C++ standard only has the latter version, but the former, which was present
// in the HP STL, is retained for backward compatibility.

template <class _InputIter, class _Tp, class _Size>
void count(_InputIter __first, _InputIter __last, const _Tp& __value,
           _Size& __n);
template <class _InputIter, class _Predicate, class _Size>
void count_if(_InputIter __first, _InputIter __last, _Predicate __pred,
              _Size& __n);


template <class _InputIter, class _Tp>
__STL_DIFFERENCE_TYPE(_InputIter)
count(_InputIter __first, _InputIter __last, const _Tp& __value);

template <class _InputIter, class _Predicate>
__STL_DIFFERENCE_TYPE(_InputIter)
count_if(_InputIter __first, _InputIter __last, _Predicate __pred);

// search.

template <class _ForwardIter1, class _ForwardIter2>
_ForwardIter1 search(_ForwardIter1 __first1, _ForwardIter1 __last1,
                     _ForwardIter2 __first2, _ForwardIter2 __last2);
template <class _ForwardIter1, class _ForwardIter2, class _BinaryPred>
_ForwardIter1 search(_ForwardIter1 __first1, _ForwardIter1 __last1,
                     _ForwardIter2 __first2, _ForwardIter2 __last2,
                     _BinaryPred  __predicate);


// search_n.  Search for __count consecutive copies of __val.

template <class _ForwardIter, class _Integer, class _Tp>
_ForwardIter search_n(_ForwardIter __first, _ForwardIter __last,
                      _Integer __count, const _Tp& __val);

template <class _ForwardIter, class _Integer, class _Tp, class _BinaryPred>
_ForwardIter search_n(_ForwardIter __first, _ForwardIter __last,
                      _Integer __count, const _Tp& __val,
                      _BinaryPred __binary_pred);



// swap_ranges

template <class _ForwardIter1, class _ForwardIter2>
_ForwardIter2 swap_ranges(_ForwardIter1 __first1, _ForwardIter1 __last1,
                          _ForwardIter2 __first2);


// transform
template <class _InputIter, class _OutputIter, class _UnaryOperation>
_OutputIter transform(_InputIter __first, _InputIter __last,
                      _OutputIter __result, _UnaryOperation __opr);

template <class _InputIter1, class _InputIter2, class _OutputIter,
          class _BinaryOperation>
_OutputIter transform(_InputIter1 __first1, _InputIter1 __last1,
                      _InputIter2 __first2, _OutputIter __result,
                      _BinaryOperation __binary_op);

// replace, replace_if, replace_copy, replace_copy_if

template <class _ForwardIter, class _Tp>
void replace(_ForwardIter __first, _ForwardIter __last,
             const _Tp& __old_value, const _Tp& __new_value);
template <class _ForwardIter, class _Predicate, class _Tp>
void replace_if(_ForwardIter __first, _ForwardIter __last,
                _Predicate __pred, const _Tp& __new_value);

template <class _InputIter, class _OutputIter, class _Tp>
_OutputIter replace_copy(_InputIter __first, _InputIter __last,
                         _OutputIter __result,
                         const _Tp& __old_value, const _Tp& __new_value);
template <class Iterator, class _OutputIter, class _Predicate, class _Tp>
_OutputIter replace_copy_if(Iterator __first, Iterator __last,
                            _OutputIter __result,
                            _Predicate __pred, const _Tp& __new_value);


// generate and generate_n

template <class _ForwardIter, class _Generator>
void generate(_ForwardIter __first, _ForwardIter __last, _Generator __gen);

template <class _OutputIter, class _Size, class _Generator>
_OutputIter generate_n(_OutputIter __first, _Size __n, _Generator __gen);

// remove, remove_if, remove_copy, remove_copy_if

template <class _InputIter, class _OutputIter, class _Tp>
_OutputIter remove_copy(_InputIter __first, _InputIter __last,
                        _OutputIter __result, const _Tp& __value);

template <class _InputIter, class _OutputIter, class _Predicate>
_OutputIter remove_copy_if(_InputIter __first, _InputIter __last,
                           _OutputIter __result, _Predicate __pred);

template <class _ForwardIter, class _Tp>
_ForwardIter remove(_ForwardIter __first, _ForwardIter __last,
		       const _Tp& __value);

template <class _ForwardIter, class _Predicate>
_ForwardIter remove_if(_ForwardIter __first, _ForwardIter __last,
			  _Predicate __pred);

// unique and unique_copy

template <class _InputIter, class _OutputIter, class _Tp>
_OutputIter __unique_copy(_InputIter __first, _InputIter __last,
                          _OutputIter __result, _Tp*);


template <class _InputIter, class _OutputIter>
inline _OutputIter __unique_copy(_InputIter __first, _InputIter __last,
                                 _OutputIter __result, 
                                 output_iterator_tag) {
  return __unique_copy(__first, __last, __result, __VALUE_TYPE(__first));
}

template <class _InputIter, class _ForwardIter>
_ForwardIter __unique_copy(_InputIter __first, _InputIter __last,
                           _ForwardIter __result, forward_iterator_tag);

# if defined (__STL_NONTEMPL_BASE_MATCH_BUG)
template <class _InputIter, class _ForwardIter>
inline 
_ForwardIter __unique_copy(_InputIter __first, _InputIter __last,
                           _ForwardIter __result, bidirectional_iterator_tag) {
  return __unique_copy(__first, __last, __result, forward_iterator_tag());
}
template <class _InputIter, class _ForwardIter>
inline
_ForwardIter __unique_copy(_InputIter __first, _InputIter __last,
                           _ForwardIter __result, random_access_iterator_tag) {
  return __unique_copy(__first, __last, __result, forward_iterator_tag());
}
# endif /* __STL_NONTEMPL_BASE_MATCH_BUG */

template <class _InputIter, class _OutputIter>
inline _OutputIter unique_copy(_InputIter __first, _InputIter __last,
                               _OutputIter __result) {
  __stl_debug_check(__check_range(__first, __last));
  if (__first == __last) return __result;
  return __unique_copy(__first, __last, __result,
                       __ITERATOR_CATEGORY(__result));
}

template <class InputIterator, class OutputIterator, class BinaryPredicate,
					    class _Tp>
OutputIterator __unique_copy(InputIterator first, InputIterator last,
                             OutputIterator result,
                             BinaryPredicate binary_pred, _Tp*);


template <class _InputIter, class _OutputIter, class _BinaryPredicate>
inline _OutputIter __unique_copy(_InputIter __first, _InputIter __last,
                                 _OutputIter __result,
                                 _BinaryPredicate __binary_pred,
                                 output_iterator_tag) {
  return __unique_copy(__first, __last, __result, __binary_pred,
                       __VALUE_TYPE(__first));
}

template <class _InputIter, class _ForwardIter, class _BinaryPredicate>
_ForwardIter __unique_copy(_InputIter __first, _InputIter __last,
                           _ForwardIter __result, 
                           _BinaryPredicate __binary_pred,
                           forward_iterator_tag);

# if defined (__STL_NONTEMPL_BASE_MATCH_BUG)
template <class InputIterator, class BidirectionalIterator,
				class BinaryPredicate>
inline BidirectionalIterator __unique_copy(InputIterator first, 
					   InputIterator last,
			            	   BidirectionalIterator result, 
					   BinaryPredicate binary_pred,
				    	   bidirectional_iterator_tag) {
  return __unique_copy(first, last, result, binary_pred,
		       forward_iterator_tag());
}

template <class InputIterator, class RandomAccessIterator,
					     class BinaryPredicate>
inline RandomAccessIterator __unique_copy(InputIterator first, 
					  InputIterator last,
			           	  RandomAccessIterator result, 
					  BinaryPredicate binary_pred,
				   	  random_access_iterator_tag) {
  return __unique_copy(first, last, result, binary_pred, 
		       forward_iterator_tag());
}
# endif /* __STL_NONTEMPL_BASE_MATCH_BUG */

template <class _InputIter, class _OutputIter, class _BinaryPredicate>
inline _OutputIter unique_copy(_InputIter __first, _InputIter __last,
                               _OutputIter __result,
                               _BinaryPredicate __binary_pred) {
  __stl_debug_check(__check_range(__first, __last));
  if (__first == __last) return __result;
  return __unique_copy(__first, __last, __result, __binary_pred,
                       __ITERATOR_CATEGORY(__result));
}

template <class _ForwardIter>
inline _ForwardIter unique(_ForwardIter __first, _ForwardIter __last) {
  __first = adjacent_find(__first, __last);
  return unique_copy(__first, __last, __first);
}

template <class _ForwardIter, class _BinaryPredicate>
inline _ForwardIter unique(_ForwardIter __first, _ForwardIter __last,
                    _BinaryPredicate __binary_pred) {
  __first = adjacent_find(__first, __last, __binary_pred);
  return unique_copy(__first, __last, __first, __binary_pred);
}

// reverse and reverse_copy, and their auxiliary functions

template <class _BidirectionalIter>
__STL_INLINE_LOOP
void __reverse(_BidirectionalIter __first, _BidirectionalIter __last, 
               bidirectional_iterator_tag) {
  while (true)
    if (__first == __last || __first == --__last)
      return;
    else
      iter_swap(__first++, __last);
}


template <class _RandomAccessIter>
__STL_INLINE_LOOP
void __reverse(_RandomAccessIter __first, _RandomAccessIter __last,
               random_access_iterator_tag) {
  for (; __first < __last; ++__first) iter_swap(__first, --__last);
}

template <class _BidirectionalIter>
inline void reverse(_BidirectionalIter __first, _BidirectionalIter __last) {
  __stl_debug_check(__check_range(__first, __last));
  __reverse(__first, __last, __ITERATOR_CATEGORY(__first));
}

template <class _BidirectionalIter, class _OutputIter>
__STL_INLINE_LOOP
_OutputIter reverse_copy(_BidirectionalIter __first,
                            _BidirectionalIter __last,
                            _OutputIter __result) {
  __stl_debug_check(__check_range(__first, __last));
  while (__first != __last) {
    --__last;
    *__result = *__last;
    ++__result;
  }
  return __result;
}

// rotate and rotate_copy, and their auxiliary functions

template <class _EuclideanRingElement>
__STL_INLINE_LOOP
_EuclideanRingElement __gcd(_EuclideanRingElement __m,
                            _EuclideanRingElement __n)
{
  while (__n != 0) {
    _EuclideanRingElement __t = __m % __n;
    __m = __n;
    __n = __t;
  }
  return __m;
}

template <class _ForwardIter, class _Distance>
_ForwardIter __rotate(_ForwardIter __first,
                      _ForwardIter __middle,
                      _ForwardIter __last,
                      _Distance*,
                      forward_iterator_tag);

template <class _BidirectionalIter, class _Distance>
_BidirectionalIter __rotate(_BidirectionalIter __first,
                            _BidirectionalIter __middle,
                            _BidirectionalIter __last,
                            _Distance*,
                            bidirectional_iterator_tag);

template <class _RandomAccessIter, class _Distance, class _Tp>
_RandomAccessIter __rotate(_RandomAccessIter __first,
                           _RandomAccessIter __middle,
                           _RandomAccessIter __last,
                           _Distance *, _Tp *);


template <class _RandomAccessIter, class _Distance>
inline _RandomAccessIter __rotate(_RandomAccessIter __first,
                           _RandomAccessIter __middle,
                           _RandomAccessIter __last,
                           _Distance * __dis,
			   random_access_iterator_tag) {
  return __rotate(__first, __middle, __last,
		  __dis, __VALUE_TYPE(__first));
}

template <class _ForwardIter>
inline _ForwardIter rotate(_ForwardIter __first, _ForwardIter __middle,
                           _ForwardIter __last) {
  __stl_debug_check(__check_range(__first, __middle));
  __stl_debug_check(__check_range(__middle, __last));
  return __rotate(__first, __middle, __last,
                  __DISTANCE_TYPE(__first),
                  __ITERATOR_CATEGORY(__first));
}

template <class _ForwardIter, class _OutputIter>
inline _OutputIter rotate_copy(_ForwardIter __first, _ForwardIter __middle,
			       _ForwardIter __last, _OutputIter __result) {
  return copy(__first, __middle, copy(__middle, __last, __result));
}


// random_shuffle

template <class _RandomAccessIter>
void random_shuffle(_RandomAccessIter __first,
		    _RandomAccessIter __last);

template <class _RandomAccessIter, class _RandomNumberGenerator>
void random_shuffle(_RandomAccessIter __first, _RandomAccessIter __last,
                    _RandomNumberGenerator& __rand);

// random_sample and random_sample_n (extensions, not part of the standard).

template <class _ForwardIter, class _OutputIter, class _Distance>
_OutputIter random_sample_n(_ForwardIter __first, _ForwardIter __last,
                            _OutputIter __out, const _Distance __n);

template <class _ForwardIter, class _OutputIter, class _Distance,
          class _RandomNumberGenerator>
_OutputIter random_sample_n(_ForwardIter __first, _ForwardIter __last,
                            _OutputIter __out, const _Distance __n,
                            _RandomNumberGenerator& __rand);

template <class _InputIter, class _RandomAccessIter, class _Distance>
_RandomAccessIter __random_sample(_InputIter __first, _InputIter __last,
                                  _RandomAccessIter __out,
                                  const _Distance __n);

template <class _InputIter, class _RandomAccessIter,
          class _RandomNumberGenerator, class _Distance>
_RandomAccessIter __random_sample(_InputIter __first, _InputIter __last,
                                  _RandomAccessIter __out,
                                  _RandomNumberGenerator& __rand,
                                  const _Distance __n);

template <class _InputIter, class _RandomAccessIter>
inline _RandomAccessIter
random_sample(_InputIter __first, _InputIter __last,
              _RandomAccessIter __out_first, _RandomAccessIter __out_last) 
{
  return __random_sample(__first, __last,
                         __out_first, __out_last - __out_first);
}


template <class _InputIter, class _RandomAccessIter, 
          class _RandomNumberGenerator>
inline _RandomAccessIter
random_sample(_InputIter __first, _InputIter __last,
              _RandomAccessIter __out_first, _RandomAccessIter __out_last,
              _RandomNumberGenerator& __rand) 
{
  return __random_sample(__first, __last,
                         __out_first, __rand,
                         __out_last - __out_first);
}

// partition, stable_partition, and their auxiliary functions

template <class _ForwardIter, class _Predicate>
_ForwardIter __partition(_ForwardIter __first,
		         _ForwardIter __last,
			 _Predicate   __pred,
			 forward_iterator_tag);

template <class _BidirectionalIter, class _Predicate>
_BidirectionalIter __partition(_BidirectionalIter __first,
                               _BidirectionalIter __last,
			       _Predicate __pred,
			       bidirectional_iterator_tag);


# if defined (__STL_NONTEMPL_BASE_MATCH_BUG)
template <class _BidirectionalIter, class _Predicate>
inline
_BidirectionalIter __partition(_BidirectionalIter __first,
                               _BidirectionalIter __last,
			       _Predicate __pred,
			       random_access_iterator_tag) {
  return __partition(__first, __last, __pred, bidirectional_iterator_tag());
}
# endif

template <class _ForwardIter, class _Predicate>
inline _ForwardIter partition(_ForwardIter __first,
   			      _ForwardIter __last,
			      _Predicate   __pred) {
  __stl_debug_check(__check_range(__first, __last));
  return __partition(__first, __last, __pred, __ITERATOR_CATEGORY(__first));
}

template <class _ForwardIter, class _Predicate, class _Distance>
_ForwardIter __inplace_stable_partition(_ForwardIter __first,
                                        _ForwardIter __last,
                                        _Predicate __pred, _Distance __len);


template <class _ForwardIter, class _Pointer, class _Predicate, 
          class _Distance>
_ForwardIter __stable_partition_adaptive(_ForwardIter __first,
                                         _ForwardIter __last,
                                         _Predicate __pred, _Distance __len,
                                         _Pointer __buffer,
                                         _Distance __buffer_size);


template <class _ForwardIter, class _Predicate, class _Tp, class _Distance>
inline _ForwardIter
__stable_partition_aux(_ForwardIter __first, _ForwardIter __last, 
                       _Predicate __pred, _Tp*, _Distance*)
{
  _Temporary_buffer<_ForwardIter, _Tp> __buf(__first, __last);
  return (__buf.size() > 0) ?
    __stable_partition_adaptive(__first, __last, __pred,
				_Distance(__buf.requested_size()),
				__buf.begin(), __buf.size())  :
    __inplace_stable_partition(__first, __last, __pred, 
			       _Distance(__buf.requested_size()));
}

template <class _ForwardIter, class _Predicate>
inline _ForwardIter stable_partition(_ForwardIter __first,
                                     _ForwardIter __last, 
                                     _Predicate __pred) {
  __stl_debug_check(__check_range(__first, __last));
  if (__first == __last)
    return __first;
  else
    return __stable_partition_aux(__first, __last, __pred,
                                  __VALUE_TYPE(__first),
                                  __DISTANCE_TYPE(__first));
}

template <class _RandomAccessIter, class _Tp>
_RandomAccessIter __unguarded_partition(_RandomAccessIter __first, 
                                        _RandomAccessIter __last, 
                                        _Tp __pivot); 

template <class _RandomAccessIter, class _Tp, class _Compare>
_RandomAccessIter __unguarded_partition(_RandomAccessIter __first, 
                                        _RandomAccessIter __last, 
                                        _Tp __pivot, _Compare __comp);

// sort() and its auxiliary functions. 

template <class _RandomAccessIter>
void __final_insertion_sort(_RandomAccessIter __first, 
                            _RandomAccessIter __last);

template <class _RandomAccessIter, class _Compare>
void __final_insertion_sort(_RandomAccessIter __first, 
                            _RandomAccessIter __last, _Compare __comp);


template <class _Size>
inline _Size __lg(_Size __n) {
  _Size __k;
  for (__k = 0; __n != 1; __n >>= 1) ++__k;
  return __k;
}

template <class _RandomAccessIter, class _Tp, class _Size>
void __introsort_loop(_RandomAccessIter __first,
                      _RandomAccessIter __last, _Tp*,
                      _Size __depth_limit);

template <class _RandomAccessIter, class _Tp, class _Size, class _Compare>
void __introsort_loop(_RandomAccessIter __first,
                      _RandomAccessIter __last, _Tp*,
                      _Size __depth_limit, _Compare __comp);


template <class _RandomAccessIter>
inline void sort(_RandomAccessIter __first, _RandomAccessIter __last) {
  __stl_debug_check(__check_range(__first, __last));
  if (__first != __last) {
    __introsort_loop(__first, __last,
                     __VALUE_TYPE(__first),
                     __lg(__last - __first) * 2);
    __final_insertion_sort(__first, __last);
  }
}

template <class _RandomAccessIter, class _Compare>
inline void sort(_RandomAccessIter __first, _RandomAccessIter __last,
                 _Compare __comp) {
  if (__first != __last) {
    __introsort_loop(__first, __last,
                     __VALUE_TYPE(__first),
                     __lg(__last - __first) * 2,
                     __comp);
    __final_insertion_sort(__first, __last, __comp);
  }
}

// stable_sort() and its auxiliary functions.
template <class _RandomAccessIter>
void stable_sort(_RandomAccessIter __first,
		 _RandomAccessIter __last);

template <class _RandomAccessIter, class _Compare>
void stable_sort(_RandomAccessIter __first,
		 _RandomAccessIter __last, _Compare __comp);

// partial_sort, partial_sort_copy, and auxiliary functions.

template <class _RandomAccessIter, class _Tp>
void __partial_sort(_RandomAccessIter __first, _RandomAccessIter __middle,
                    _RandomAccessIter __last, _Tp*);

template <class _RandomAccessIter>
inline void partial_sort(_RandomAccessIter __first,
                         _RandomAccessIter __middle,
                         _RandomAccessIter __last) {
  __stl_debug_check(__check_range(__first, __middle));
  __stl_debug_check(__check_range(__middle, __last));
  __partial_sort(__first, __middle, __last, __VALUE_TYPE(__first));
}

template <class _RandomAccessIter, class _Tp, class _Compare>
void __partial_sort(_RandomAccessIter __first, _RandomAccessIter __middle,
                    _RandomAccessIter __last, _Tp*, _Compare __comp);

template <class _RandomAccessIter, class _Compare>
inline void partial_sort(_RandomAccessIter __first,
                         _RandomAccessIter __middle,
                         _RandomAccessIter __last, _Compare __comp) {
  __stl_debug_check(__check_range(__first, __middle));
  __stl_debug_check(__check_range(__middle, __last));
  __partial_sort(__first, __middle, __last, __VALUE_TYPE(__first), __comp);
}

template <class _InputIter, class _RandomAccessIter, class _Distance,
          class _Tp>
_RandomAccessIter __partial_sort_copy(_InputIter __first,
                                         _InputIter __last,
                                         _RandomAccessIter __result_first,
                                         _RandomAccessIter __result_last, 
                                         _Distance*, _Tp*);

template <class _InputIter, class _RandomAccessIter>
inline _RandomAccessIter
partial_sort_copy(_InputIter __first, _InputIter __last,
                  _RandomAccessIter __result_first,
                  _RandomAccessIter __result_last) {
  __stl_debug_check(__check_range(__first, __last));
  __stl_debug_check(__check_range(__result_first, __result_last));
  return __partial_sort_copy(__first, __last, __result_first, __result_last, 
                             __DISTANCE_TYPE(__result_first),
                             __VALUE_TYPE(__first));
}

template <class _InputIter, class _RandomAccessIter, class _Compare,
          class _Distance, class _Tp>
_RandomAccessIter __partial_sort_copy(_InputIter __first,
                                         _InputIter __last,
                                         _RandomAccessIter __result_first,
                                         _RandomAccessIter __result_last,
                                         _Compare __comp, _Distance*, _Tp*);

template <class _InputIter, class _RandomAccessIter, class _Compare>
inline _RandomAccessIter
partial_sort_copy(_InputIter __first, _InputIter __last,
                  _RandomAccessIter __result_first,
                  _RandomAccessIter __result_last, _Compare __comp) {
  __stl_debug_check(__check_range(__first, __last));
  __stl_debug_check(__check_range(__result_first, __result_last));
  return __partial_sort_copy(__first, __last, __result_first, __result_last,
                             __comp,
                             __DISTANCE_TYPE(__result_first),
                             __VALUE_TYPE(__first));
}

// nth_element() and its auxiliary functions.  

template <class _RandomAccessIter, class _Tp>
void __nth_element(_RandomAccessIter __first, _RandomAccessIter __nth,
                   _RandomAccessIter __last, _Tp*);

template <class _RandomAccessIter>
inline void nth_element(_RandomAccessIter __first, _RandomAccessIter __nth,
                        _RandomAccessIter __last) {
  __stl_debug_check(__check_range(__first, __nth));
  __stl_debug_check(__check_range(__nth, __last));
  __nth_element(__first, __nth, __last, __VALUE_TYPE(__first));
}

template <class _RandomAccessIter, class _Tp, class _Compare>
void __nth_element(_RandomAccessIter __first, _RandomAccessIter __nth,
                   _RandomAccessIter __last, _Tp*, _Compare __comp);

template <class _RandomAccessIter, class _Compare>
inline void nth_element(_RandomAccessIter __first, _RandomAccessIter __nth,
                 _RandomAccessIter __last, _Compare __comp) {
  __stl_debug_check(__check_range(__first, __nth));
  __stl_debug_check(__check_range(__nth, __last));
  __nth_element(__first, __nth, __last, __VALUE_TYPE(__first), __comp);
}


// Binary search (lower_bound, upper_bound, equal_range, binary_search).

template <class _ForwardIter, class _Tp, class _Distance>
_ForwardIter __lower_bound(_ForwardIter __first, _ForwardIter __last,
                           const _Tp& __val, _Distance*);


template <class _ForwardIter, class _Tp>
inline _ForwardIter lower_bound(_ForwardIter __first, _ForwardIter __last,
                                   const _Tp& __val) {
  __stl_debug_check(__check_range(__first, __last));
  return __lower_bound(__first, __last, __val,
                       __DISTANCE_TYPE(__first));
}

template <class _ForwardIter, class _Tp, class _Compare, class _Distance>
_ForwardIter __lower_bound(_ForwardIter __first, _ForwardIter __last,
                              const _Tp& __val, _Compare __comp, _Distance*);


template <class _ForwardIter, class _Tp, class _Compare>
inline _ForwardIter lower_bound(_ForwardIter __first, _ForwardIter __last,
                                const _Tp& __val, _Compare __comp) {
  __stl_debug_check(__check_range(__first, __last));
  return __lower_bound(__first, __last, __val, __comp,
                       __DISTANCE_TYPE(__first));
}

template <class _ForwardIter, class _Tp, class _Distance>
_ForwardIter __upper_bound(_ForwardIter __first, _ForwardIter __last,
                           const _Tp& __val, _Distance*);

template <class _ForwardIter, class _Tp>
inline _ForwardIter upper_bound(_ForwardIter __first, _ForwardIter __last,
                                const _Tp& __val) {
  __stl_debug_check(__check_range(__first, __last));
  return __upper_bound(__first, __last, __val,
                       __DISTANCE_TYPE(__first));
}

template <class _ForwardIter, class _Tp, class _Compare, class _Distance>
_ForwardIter __upper_bound(_ForwardIter __first, _ForwardIter __last,
                           const _Tp& __val, _Compare __comp, _Distance*);

template <class _ForwardIter, class _Tp, class _Compare>
inline _ForwardIter upper_bound(_ForwardIter __first, _ForwardIter __last,
                                const _Tp& __val, _Compare __comp) {
  __stl_debug_check(__check_range(__first, __last));
  return __upper_bound(__first, __last, __val, __comp,
                       __DISTANCE_TYPE(__first));
}

template <class _ForwardIter, class _Tp, class _Distance>
pair<_ForwardIter, _ForwardIter>
__equal_range(_ForwardIter __first, _ForwardIter __last, const _Tp& __val,
              _Distance*);

template <class _ForwardIter, class _Tp>
inline pair<_ForwardIter, _ForwardIter>
equal_range(_ForwardIter __first, _ForwardIter __last, const _Tp& __val) {
  __stl_debug_check(__check_range(__first, __last));
  return __equal_range(__first, __last, __val,
                       __DISTANCE_TYPE(__first));
}

template <class _ForwardIter, class _Tp, class _Compare, class _Distance>
pair<_ForwardIter, _ForwardIter>
__equal_range(_ForwardIter __first, _ForwardIter __last, const _Tp& __val,
              _Compare __comp, _Distance*);

template <class _ForwardIter, class _Tp, class _Compare>
inline pair<_ForwardIter, _ForwardIter>
equal_range(_ForwardIter __first, _ForwardIter __last, const _Tp& __val,
            _Compare __comp) {
  __stl_debug_check(__check_range(__first, __last));
  return __equal_range(__first, __last, __val, __comp,
                       __DISTANCE_TYPE(__first));
} 

template <class _ForwardIter, class _Tp>
inline bool binary_search(_ForwardIter __first, _ForwardIter __last,
                   const _Tp& __val) {
  _ForwardIter __i = lower_bound(__first, __last, __val);
  return __i != __last && !(__val < *__i);
}

template <class _ForwardIter, class _Tp, class _Compare>
inline bool binary_search(_ForwardIter __first, _ForwardIter __last,
                   const _Tp& __val,
                   _Compare __comp) {
  _ForwardIter __i = lower_bound(__first, __last, __val, __comp);
  return __i != __last && !__comp(__val, *__i);
}

// merge, with and without an explicitly supplied comparison function.

template <class _InputIter1, class _InputIter2, class _OutputIter>
_OutputIter merge(_InputIter1 __first1, _InputIter1 __last1,
                  _InputIter2 __first2, _InputIter2 __last2,
                  _OutputIter __result);
 
template <class _InputIter1, class _InputIter2, class _OutputIter,
          class _Compare>
_OutputIter merge(_InputIter1 __first1, _InputIter1 __last1,
                  _InputIter2 __first2, _InputIter2 __last2,
                  _OutputIter __result, _Compare __comp);


// inplace_merge and its auxiliary functions. 


template <class _BidirectionalIter>
void inplace_merge(_BidirectionalIter __first,
		   _BidirectionalIter __middle,
		   _BidirectionalIter __last) ;

template <class _BidirectionalIter, class _Compare>
void inplace_merge(_BidirectionalIter __first,
		   _BidirectionalIter __middle,
		   _BidirectionalIter __last, _Compare __comp);

// Set algorithms: includes, set_union, set_intersection, set_difference,
// set_symmetric_difference.  All of these algorithms have the precondition
// that their input ranges are sorted and the postcondition that their output
// ranges are sorted.

template <class _InputIter1, class _InputIter2>
bool includes(_InputIter1 __first1, _InputIter1 __last1,
              _InputIter2 __first2, _InputIter2 __last2);

template <class _InputIter1, class _InputIter2, class _Compare>
bool includes(_InputIter1 __first1, _InputIter1 __last1,
              _InputIter2 __first2, _InputIter2 __last2, _Compare __comp);
 
template <class _InputIter1, class _InputIter2, class _OutputIter>
_OutputIter set_union(_InputIter1 __first1, _InputIter1 __last1,
                      _InputIter2 __first2, _InputIter2 __last2,
                      _OutputIter __result);

template <class _InputIter1, class _InputIter2, class _OutputIter,
          class _Compare>
_OutputIter set_union(_InputIter1 __first1, _InputIter1 __last1,
                      _InputIter2 __first2, _InputIter2 __last2,
                      _OutputIter __result, _Compare __comp);

template <class _InputIter1, class _InputIter2, class _OutputIter>
_OutputIter set_intersection(_InputIter1 __first1, _InputIter1 __last1,
                             _InputIter2 __first2, _InputIter2 __last2,
                             _OutputIter __result);

template <class _InputIter1, class _InputIter2, class _OutputIter,
          class _Compare>
_OutputIter set_intersection(_InputIter1 __first1, _InputIter1 __last1,
                             _InputIter2 __first2, _InputIter2 __last2,
                             _OutputIter __result, _Compare __comp);



template <class _InputIter1, class _InputIter2, class _OutputIter>
_OutputIter set_difference(_InputIter1 __first1, _InputIter1 __last1,
                           _InputIter2 __first2, _InputIter2 __last2,
                           _OutputIter __result);

template <class _InputIter1, class _InputIter2, class _OutputIter, 
          class _Compare>
_OutputIter set_difference(_InputIter1 __first1, _InputIter1 __last1,
                           _InputIter2 __first2, _InputIter2 __last2, 
                           _OutputIter __result, _Compare __comp);

template <class _InputIter1, class _InputIter2, class _OutputIter>
_OutputIter 
set_symmetric_difference(_InputIter1 __first1, _InputIter1 __last1,
                         _InputIter2 __first2, _InputIter2 __last2,
                         _OutputIter __result);


template <class _InputIter1, class _InputIter2, class _OutputIter,
          class _Compare>
_OutputIter 
set_symmetric_difference(_InputIter1 __first1, _InputIter1 __last1,
                         _InputIter2 __first2, _InputIter2 __last2,
                         _OutputIter __result,
                         _Compare __comp);


// min_element and max_element, with and without an explicitly supplied
// comparison function.

template <class _ForwardIter>
_ForwardIter max_element(_ForwardIter __first, _ForwardIter __last);
template <class _ForwardIter, class _Compare>
_ForwardIter max_element(_ForwardIter __first, _ForwardIter __last,
                            _Compare __comp);

template <class _ForwardIter>
_ForwardIter min_element(_ForwardIter __first, _ForwardIter __last);

template <class _ForwardIter, class _Compare>
_ForwardIter min_element(_ForwardIter __first, _ForwardIter __last,
                            _Compare __comp);

// next_permutation and prev_permutation, with and without an explicitly 
// supplied comparison function.

template <class _BidirectionalIter>
bool next_permutation(_BidirectionalIter __first, _BidirectionalIter __last);

template <class _BidirectionalIter, class _Compare>
bool next_permutation(_BidirectionalIter __first, _BidirectionalIter __last,
                      _Compare __comp);


template <class _BidirectionalIter>
bool prev_permutation(_BidirectionalIter __first, _BidirectionalIter __last);


template <class _BidirectionalIter, class _Compare>
bool prev_permutation(_BidirectionalIter __first, _BidirectionalIter __last,
                      _Compare __comp);


// find_first_of, with and without an explicitly supplied comparison function.

template <class _InputIter, class _ForwardIter>
_InputIter find_first_of(_InputIter __first1, _InputIter __last1,
                         _ForwardIter __first2, _ForwardIter __last2);

template <class _InputIter, class _ForwardIter, class _BinaryPredicate>
_InputIter find_first_of(_InputIter __first1, _InputIter __last1,
                         _ForwardIter __first2, _ForwardIter __last2,
                         _BinaryPredicate __comp);


// find_end, with and without an explicitly supplied comparison function.
// Search [first2, last2) as a subsequence in [first1, last1), and return
// the *last* possible match.  Note that find_end for bidirectional iterators
// is much faster than for forward iterators.


// Dispatching functions for find_end.

template <class _ForwardIter1, class _ForwardIter2>
_ForwardIter1 
find_end(_ForwardIter1 __first1, _ForwardIter1 __last1, 
         _ForwardIter2 __first2, _ForwardIter2 __last2);

template <class _ForwardIter1, class _ForwardIter2, 
          class _BinaryPredicate>
_ForwardIter1 
find_end(_ForwardIter1 __first1, _ForwardIter1 __last1, 
         _ForwardIter2 __first2, _ForwardIter2 __last2,
         _BinaryPredicate __comp);

// is_heap, a predicate testing whether or not a range is
// a heap.  This function is an extension, not part of the C++
// standard.

template <class _RandomAccessIter>
bool is_heap(_RandomAccessIter __first, _RandomAccessIter __last);

template <class _RandomAccessIter, class _StrictWeakOrdering>
bool is_heap(_RandomAccessIter __first, _RandomAccessIter __last,
	     _StrictWeakOrdering __comp);

// is_sorted, a predicated testing whether a range is sorted in
// nondescending order.  This is an extension, not part of the C++
// standard.

template <class _ForwardIter>
bool is_sorted(_ForwardIter __first, _ForwardIter __last);

template <class _ForwardIter, class _StrictWeakOrdering>
bool is_sorted(_ForwardIter __first, _ForwardIter __last,
               _StrictWeakOrdering __comp);

#if defined(__sgi) && !defined(__GNUC__) && (_MIPS_SIM != _MIPS_SIM_ABI32)
#pragma reset woff 1209
#endif

__STL_END_NAMESPACE

# if !defined (__STL_LINK_TIME_INSTANTIATION)
#  include <stl_algo.c>
# endif

#endif /* __SGI_STL_INTERNAL_ALGO_H */

// Local Variables:
// mode:C++
// End:

