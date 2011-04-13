/*
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
#ifndef __STL_SLIST_C
#define __STL_SLIST_C

#if defined(__sgi) && !defined(__GNUC__) && (_MIPS_SIM != _MIPS_SIM_ABI32)
#pragma set woff 1174
#pragma set woff 1375
#endif

# if defined (__STL_DEBUG)
#  define _Make_iterator(__l) iterator(&_M_iter_list, (__l))
#  define _Make_const_iterator(__l) const_iterator(&_M_iter_list, (__l))
# else
#  define _Make_iterator(__l) iterator((__l))
#  define _Make_const_iterator(__l) const_iterator((__l))
# endif
#  define  slist  __WORKAROUND_RENAME(slist)

__STL_BEGIN_NAMESPACE

template <class _Tp, class _Alloc> 
_Slist_node_base*
_Slist_base<_Tp,_Alloc>::_M_erase_after(_Slist_node_base* __before_first,
                                        _Slist_node_base* __last_node) {
  _Slist_node<_Tp>* __cur = (_Slist_node<_Tp>*) (__before_first->_M_next);
  while (__cur != __last_node) {
    _Slist_node<_Tp>* __tmp = __cur;
    __cur = (_Slist_node<_Tp>*) __cur->_M_next;
    destroy(&__tmp->_M_data);
    _M_head.deallocate(__tmp,1);
  }
  __before_first->_M_next = __last_node;
  return __last_node;
}

#if defined (__STL_MEMBER_TEMPLATES) && ! defined (__STL_INLINE_MEMBER_TEMPLATES)
template <class _Tp, class _Alloc> template <class _InputIter>
void
slist<_Tp,_Alloc>::_M_assign_dispatch(_InputIter __first, _InputIter __last,
		   __false_type)
{
  _Node_base* __prev = &this->_M_head._M_data;
  _Node* __node = (_Node*) this->_M_head._M_data._M_next;
  while (__node != 0 && __first != __last) {
    __node->_M_data = *__first;
    __prev = __node;
    __node = (_Node*) __node->_M_next;
    ++__first;
  }
  if (__first != __last)
    _M_insert_after_range(__prev, __first, __last);
  else
    this->_M_erase_after(__prev, 0);
}

template <class _Tp, class _Alloc>   template <class _Predicate>
void slist<_Tp,_Alloc>::remove_if(_Predicate __pred) {
  _Node_base* __cur = &this->_M_head._M_data;
  while (__cur->_M_next) {
    if (__pred(((_Node*) __cur->_M_next)->_M_data))
      this->_M_erase_after(__cur);
    else
      __cur = __cur->_M_next;
  }
}

template <class _Tp, class _Alloc>   template <class _BinaryPredicate> 
void slist<_Tp,_Alloc>::unique(_BinaryPredicate __pred) {
  _Node* __cur = (_Node*) this->_M_head._M_data._M_next;
  if (__cur) {
    while (__cur->_M_next) {
      if (__pred(((_Node*)__cur)->_M_data, 
		 ((_Node*)(__cur->_M_next))->_M_data))
	this->_M_erase_after(__cur);
      else
	__cur = (_Node*) __cur->_M_next;
    }
  }
}

template <class _Tp, class _Alloc>   template <class _StrictWeakOrdering>
void slist<_Tp,_Alloc>::merge(slist<_Tp,_Alloc>& __x,
			      _StrictWeakOrdering __comp) {
  _Node_base* __n1 = &this->_M_head._M_data;
  while (__n1->_M_next && __x._M_head._M_data._M_next) {
    if (__comp(((_Node*) __x._M_head._M_data._M_next)->_M_data,
	       ((_Node*)       __n1->_M_next)->_M_data))
      _Sl_global_inst::__splice_after(__n1, &__x._M_head._M_data, __x._M_head._M_data._M_next);
    __n1 = __n1->_M_next;
  }
  if (__x._M_head._M_data._M_next) {
    __n1->_M_next = __x._M_head._M_data._M_next;
    __x._M_head._M_data._M_next = 0;
  }
}

template <class _Tp, class _Alloc>   template <class _StrictWeakOrdering> 
void slist<_Tp,_Alloc>::sort(_StrictWeakOrdering __comp) {
  if (this->_M_head._M_data._M_next && this->_M_head._M_data._M_next->_M_next) {
    slist __carry;
    slist __counter[64];
    int __fill = 0;
    while (!empty()) {
      _Sl_global_inst::__splice_after(&__carry._M_head._M_data, &this->_M_head._M_data, this->_M_head._M_data._M_next);
      int __i = 0;
      while (__i < __fill && !__counter[__i].empty()) {
	__counter[__i].merge(__carry, __comp);
	__carry.swap(__counter[__i]);
	++__i;
      }
      __carry.swap(__counter[__i]);
      if (__i == __fill)
	++__fill;
    }
    
    for (int __i = 1; __i < __fill; ++__i)
	__counter[__i].merge(__counter[__i-1], __comp);
    this->swap(__counter[__fill-1]);
  }
}

#endif  /* __STL_MEMBER_TEMPLATES */

# if defined (__STL_NESTED_TYPE_PARAM_BUG) 
#  define __size_type__      size_t
# else
#  define __size_type__      size_type
# endif

template <class _Tp, class _Alloc>
slist<_Tp,_Alloc>& slist<_Tp,_Alloc>::operator=(const slist<_Tp,_Alloc>& __x)
{
  if (&__x != this) {
    _Node_base* __p1 = &this->_M_head._M_data;
    _Node* __n1 = (_Node*) this->_M_head._M_data._M_next;
    const _Node* __n2 = (const _Node*) __x._M_head._M_data._M_next;
    while (__n1 && __n2) {
      __n1->_M_data = __n2->_M_data;
      __p1 = __n1;
      __n1 = (_Node*) __n1->_M_next;
      __n2 = (const _Node*) __n2->_M_next;
    }
    if (__n2 == 0)
      this->_M_erase_after(__p1, 0);
    else
      _M_insert_after_range(__p1, _Make_const_iterator((_Node*)__n2), 
                                  _Make_const_iterator(0));
  }
  return *this;
}

template <class _Tp, class _Alloc>
void slist<_Tp, _Alloc>::_M_fill_assign(__size_type__ __n, const _Tp& __val) {
  _Node_base* __prev = &this->_M_head._M_data;
  _Node* __node = (_Node*) this->_M_head._M_data._M_next;
  for ( ; __node != 0 && __n > 0 ; --__n) {
    __node->_M_data = __val;
    __prev = __node;
    __node = (_Node*) __node->_M_next;
  }
  if (__n > 0)
    _M_insert_after_fill(__prev, __n, __val);
  else
    this->_M_erase_after(__prev, 0);
}


template <class _Tp, class _Alloc>
void slist<_Tp,_Alloc>::resize(__size_type__ __len, const _Tp& __x)
{
  _Node_base* __cur = &this->_M_head._M_data;
  while (__cur->_M_next != 0 && __len > 0) {
    --__len;
    __cur = __cur->_M_next;
  }
  if (__cur->_M_next) 
    this->_M_erase_after(__cur, 0);
  else
    _M_insert_after_fill(__cur, __len, __x);
}

template <class _Tp, class _Alloc>
void slist<_Tp,_Alloc>::remove(const _Tp& __val)
{
  _Node_base* __cur = &this->_M_head._M_data;
  while (__cur && __cur->_M_next) {
    if (((_Node*) __cur->_M_next)->_M_data == __val)
      this->_M_erase_after(__cur);
    else
      __cur = __cur->_M_next;
  }
}

template <class _Tp, class _Alloc> 
void slist<_Tp,_Alloc>::unique()
{
  _Node_base* __cur = this->_M_head._M_data._M_next;
  if (__cur) {
    while (__cur->_M_next) {
      if (((_Node*)__cur)->_M_data == 
          ((_Node*)(__cur->_M_next))->_M_data)
        this->_M_erase_after(__cur);
      else
        __cur = __cur->_M_next;
    }
  }
}

template <class _Tp, class _Alloc>
void slist<_Tp,_Alloc>::merge(slist<_Tp,_Alloc>& __x)
{
  _Node_base* __n1 = &this->_M_head._M_data;
  while (__n1->_M_next && __x._M_head._M_data._M_next) {
    if (((_Node*) __x._M_head._M_data._M_next)->_M_data < 
        ((_Node*)       __n1->_M_next)->_M_data) 
      _Sl_global_inst::__splice_after(__n1, &__x._M_head._M_data, __x._M_head._M_data._M_next);
    __n1 = __n1->_M_next;
  }
  if (__x._M_head._M_data._M_next) {
    __n1->_M_next = __x._M_head._M_data._M_next;
    __x._M_head._M_data._M_next = 0;
  }
  __stl_debug_do(__x._Invalidate_all());
}

template <class _Tp, class _Alloc>
void slist<_Tp,_Alloc>::sort()
{
  if (this->_M_head._M_data._M_next && this->_M_head._M_data._M_next->_M_next) {
    _Self __carry;
    _Self __counter[64];
    int __fill = 0;
    while (!empty()) {
      _Sl_global_inst::__splice_after(&__carry._M_head._M_data, &this->_M_head._M_data, this->_M_head._M_data._M_next);
      int __i = 0;
      while (__i < __fill && !__counter[__i].empty()) {
        __counter[__i].merge(__carry);
        __carry.swap(__counter[__i]);
        ++__i;
      }
      __carry.swap(__counter[__i]);
      if (__i == __fill)
        ++__fill;
    }

    for (int __i = 1; __i < __fill; ++__i)
      __counter[__i].merge(__counter[__i-1]);
    this->swap(__counter[__fill-1]);
  }
}

# undef _Make_iterator
# undef _Make_const_iterator
# undef slist
# undef __size_type__

__STL_END_NAMESPACE

#endif /*  __STL_SLIST_C */

// Local Variables:
// mode:C++
// End:
