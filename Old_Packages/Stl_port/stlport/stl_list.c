/*
 *
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
#ifndef __STL_LIST_C
#define __STL_LIST_C

#if defined(__sgi) && !defined(__GNUC__) && (_MIPS_SIM != _MIPS_SIM_ABI32)
#pragma set woff 1174
#pragma set woff 1375
#endif

# undef list
# define  list  __WORKAROUND_RENAME(list)

__STL_BEGIN_NAMESPACE

template <class _Dummy>
void 
_List_global<_Dummy>::_Transfer(_List_node_base* __position, 
				_List_node_base* __first, _List_node_base* __last) {
  if (__position != __last) {
    // Remove [first, last) from its old position.
    ((_Node*) (__last->_M_prev))->_M_next = __position;
    ((_Node*) (__first->_M_prev))->_M_next    = __last;
    ((_Node*) (__position->_M_prev))->_M_next = __first; 
    
    // Splice [first, last) into its new position.
    _Node* __tmp = (_Node*) (__position->_M_prev);
    __position->_M_prev = __last->_M_prev;
    __last->_M_prev      = __first->_M_prev; 
    __first->_M_prev    = __tmp;
  }
}

template <class _Tp, class _Alloc>
void 
_List_base<_Tp,_Alloc>::clear() 
{
  __stl_debug_do(_Invalidate_all());
  _List_node<_Tp>* __cur = (_List_node<_Tp>*) _M_node._M_data->_M_next;
  while (__cur != _M_node._M_data) {
    _List_node<_Tp>* __tmp = __cur;
    __cur = (_List_node<_Tp>*) __cur->_M_next;
    destroy(&__tmp->_M_data);
    _M_node.deallocate(__tmp, 1);
  }
  _M_node._M_data->_M_next = _M_node._M_data;
  _M_node._M_data->_M_prev = _M_node._M_data;
}

# if defined (__STL_NESTED_TYPE_PARAM_BUG) 
#  define __iterator__   _List_iterator<_Tp, _Nonconst_traits<_Tp> >
#  define iterator       __iterator__
#  define const_iterator _List_iterator<_Tp, _Const_traits<_Tp> >
#  define size_type      size_t
# else
#  define __iterator__ __STL_TYPENAME_ON_RETURN_TYPE list<_Tp,_Alloc>::iterator
# endif

# if defined (__STL_MEMBER_TEMPLATES) && ! defined (__STL_INLINE_MEMBER_TEMPLATES)

template <class _Tp, class _Alloc> template <class _Predicate>
void list<_Tp, _Alloc>::remove_if(_Predicate __pred) 
{
  iterator __first = begin();
  iterator __last = end();
  while (__first != __last) {
    iterator __next = __first;
    ++__next;
    if (__pred(*__first)) erase(__first);
    __first = __next;
  }
}

template <class _Tp, class _Alloc>  template <class _BinaryPredicate>
void list<_Tp, _Alloc>::unique(_BinaryPredicate __binary_pred) {
  iterator __first = begin();
  iterator __last = end();
  if (__first == __last) return;
  iterator __next = __first;
  while (++__next != __last) {
    if (__binary_pred(*__first, *__next))
      erase(__next);
    else
      __first = __next;
    __next = __first;
  }
}

template <class _Tp, class _Alloc>  template <class _StrictWeakOrdering>
void list<_Tp, _Alloc>::merge(list<_Tp, _Alloc>& __x,
			      _StrictWeakOrdering __comp) {
  iterator __first1 = begin();
  iterator __last1 = end();
  iterator __first2 = __x.begin();
  iterator __last2 = __x.end();
  while (__first1 != __last1 && __first2 != __last2)
    if (__comp(*__first2, *__first1)) {
      iterator __next = __first2;
      _List_global_inst::_Transfer(__first1._M_node, __first2._M_node, (++__next)._M_node);
      __first2 = __next;
    }
    else
      ++__first1;
  if (__first2 != __last2) _List_global_inst::_Transfer(__last1._M_node, __first2._M_node, __last2._M_node);
}

template <class _Tp, class _Alloc> template <class _StrictWeakOrdering>
void list<_Tp, _Alloc>::sort(_StrictWeakOrdering __comp) {
  // Do nothing if the list has length 0 or 1.
  if (_M_node._M_data->_M_next != _M_node._M_data &&
      ((_Node*) (_M_node._M_data->_M_next))->_M_next != _M_node._M_data) {
    list<_Tp, _Alloc> __carry;
    list<_Tp, _Alloc> __counter[64];
    int __fill = 0;
    while (!empty()) {
      __carry.splice(__carry.begin(), *this, begin());
      int __i = 0;
      while(__i < __fill && !__counter[__i].empty()) {
	__counter[__i].merge(__carry, __comp);
	__carry.swap(__counter[__i++]);
      }
      __carry.swap(__counter[__i]);         
      if (__i == __fill) ++__fill;
    } 

    for (int __i = 1; __i < __fill; ++__i) 
      __counter[__i].merge(__counter[__i-1], __comp);
    swap(__counter[__fill-1]);
  }
  __stl_debug_do(_Invalidate_all());
}
# endif

#ifndef __STL_MEMBER_TEMPLATES

template <class _Tp, class _Alloc>
void 
list<_Tp, _Alloc>::insert(iterator __position, 
                          const _Tp* __first, const _Tp* __last)
{
  for ( ; __first != __last; ++__first)
    insert(__position, *__first);
}

template <class _Tp, class _Alloc>
void 
list<_Tp, _Alloc>::insert(iterator __position,
                         const_iterator __first, const_iterator __last)
{
  for ( ; __first != __last; ++__first)
    insert(__position, *__first);
}

#endif

template <class _Tp, class _Alloc>
void 
list<_Tp, _Alloc>::_M_fill_insert(iterator __position, size_type __n, const _Tp& __x)
{
  for ( ; __n > 0; --__n)
    insert(__position, __x);
}

template <class _Tp, class _Alloc>
__iterator__ list<_Tp, _Alloc>::erase(iterator __first, 
				      iterator __last)
{
  while (__first != __last)
    erase(__first++);
  return __last;
}

template <class _Tp, class _Alloc>
void list<_Tp, _Alloc>::resize(size_type __new_size, const _Tp& __x)
{
  iterator __i = begin();
  size_type __len = 0;
  for ( ; __i != end() && __len < __new_size; ++__i, ++__len)
    ;
  if (__len == __new_size)
    erase(__i, end());
  else                          // __i == end()
    insert(end(), __new_size - __len, __x);
}

template <class _Tp, class _Alloc>
list<_Tp, _Alloc>& list<_Tp, _Alloc>::operator=(const list<_Tp, _Alloc>& __x)
{
  if (this != &__x) {
    iterator __first1 = begin();
    iterator __last1 = end();
    const_iterator __first2 = __x.begin();
    const_iterator __last2 = __x.end();
    while (__first1 != __last1 && __first2 != __last2) 
      *__first1++ = *__first2++;
    if (__first2 == __last2)
      erase(__first1, __last1);
    else
      insert(__last1, __first2, __last2);
  }
  __stl_debug_do(_Invalidate_all());
  return *this;
}

template <class _Tp, class _Alloc>
void list<_Tp, _Alloc>::_M_fill_assign(size_type __n, const _Tp& __val) {
  iterator __i = begin();
  for ( ; __i != end() && __n > 0; ++__i, --__n)
    *__i = __val;
  if (__n > 0)
    insert(end(), __n, __val);
  else
    erase(__i, end());
}

template <class _Tp, class _Alloc>
void list<_Tp, _Alloc>::remove(const _Tp& __value)
{
  iterator __first = begin();
  iterator __last = end();
  while (__first != __last) {
    iterator __next = __first;
    ++__next;
    if (*__first == __value) 
      erase(__first);
    __first = __next;
  }
}

template <class _Tp, class _Alloc>
void list<_Tp, _Alloc>::unique()
{
  iterator __first = begin();
  iterator __last = end();
  if (__first == __last) return;
  iterator __next = __first;
  while (++__next != __last) {
    if (*__first == *__next)
      erase(__next);
    else
      __first = __next;
    __next = __first;
  }
}

template <class _Tp, class _Alloc>
void list<_Tp, _Alloc>::merge(list<_Tp, _Alloc>& __x)
{
  iterator __first1 = begin();
  iterator __last1 = end();
  iterator __first2 = __x.begin();
  iterator __last2 = __x.end();
  while (__first1 != __last1 && __first2 != __last2)
    if (*__first2 < *__first1) {
      iterator __next = __first2;
      _List_global_inst::_Transfer(__first1._M_node, __first2._M_node, (++__next)._M_node);
      __first2 = __next;
    }
    else
      ++__first1;
  if (__first2 != __last2) _List_global_inst::_Transfer(__last1._M_node, __first2._M_node, __last2._M_node);
  __stl_debug_do(__x._Invalidate_all());
}

template <class _Tp, class _Alloc>
void list<_Tp, _Alloc>::reverse() 
{
  // Do nothing if the list has length 0 or 1.
  if (_M_node._M_data->_M_next != _M_node._M_data &&
      ((_Node*) (_M_node._M_data->_M_next))->_M_next != _M_node._M_data) {
    _Node* __first = (_Node*)((_Node*)_M_node._M_data->_M_next)->_M_next;
    //    __first = __first->_M_next;
    while (__first != _M_node._M_data) {
      _Node* __old = __first;
      __first = (_Node*)__first->_M_next;
      _List_global_inst::_Transfer(( _Node*)_M_node._M_data->_M_next, __old, __first);
    }
  }
  __stl_debug_do(_Invalidate_all());
}    

template <class _Tp, class _Alloc>
void list<_Tp, _Alloc>::sort()
{
  // Do nothing if the list has length 0 or 1.
  if (_M_node._M_data->_M_next != _M_node._M_data &&
      ((_Node*) (_M_node._M_data->_M_next))->_M_next != _M_node._M_data) {
    list<_Tp, _Alloc> __carry;
    list<_Tp, _Alloc> __counter[64];
    int __fill = 0;
    while (!empty()) {
      __carry.splice(__carry.begin(), *this, begin());
      int __i = 0;
      while(__i < __fill && !__counter[__i].empty()) {
        __counter[__i].merge(__carry);
        __carry.swap(__counter[__i++]);
      }
      __carry.swap(__counter[__i]);         
      if (__i == __fill) ++__fill;
    } 

    for (int __i = 1; __i < __fill; ++__i)
      __counter[__i].merge(__counter[__i-1]);
    swap(__counter[__fill-1]);
  }
  __stl_debug_do(_Invalidate_all());
}

# undef list

# undef  iterator
# undef  const_iterator
# undef  size_type
# undef __iterator__

__STL_END_NAMESPACE

#endif /*  __STL_LIST_C */

// Local Variables:
// mode:C++
// End:
