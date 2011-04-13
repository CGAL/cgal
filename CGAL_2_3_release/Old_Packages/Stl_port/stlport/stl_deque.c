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
#ifndef __STL_DEQUE_C
#define __STL_DEQUE_C

#if defined(__sgi) && !defined(__GNUC__) && (_MIPS_SIM != _MIPS_SIM_ABI32)
#pragma set woff 1174
#pragma set woff 1375
#endif

# undef deque
# if defined ( __STL_NO_DEFAULT_NON_TYPE_PARAM )
#  define deque __deque
# else
#  define deque __WORKAROUND_RENAME(deque)
# endif

__STL_BEGIN_NAMESPACE

# ifdef __STL_DEBUG
// this hack is horrible, but, given no Alloc parameter
// for _Deque_iterator, we are not able to restore full deque structure anyways

template <class _Tp>
struct _Deq_iter_guts : public __owned_link {
  _Tp* _M_cur;
  _Tp* _M_first;
  _Tp* _M_last;
  _Tp** _M_node; 
  bool _M_unsafe;
};

// We do not send __ptr as _Tp*, as only void* is guaranteed to hold any pointer
template <class _Tp>
bool __Deq_dereferenceable(const void* __ptr, _Tp*) {
  typedef _Deq_iter_guts<_Tp> _Guts;
  const _Guts * __guts = (const _Guts*)(void*)__ptr;
  __stl_verbose_return(__guts->_Valid(), _StlMsg_INVALID_ITERATOR);
  
  if (__guts->_M_unsafe) return true;
  
  const _Guts* __start = (const _Guts*)(__guts->_Owner()->_Owner());
  const _Guts* __finish = __start+1;
  
  __stl_verbose_return(
		       ((__guts->_M_node == __finish->_M_node) ? /* *__guts < *__finish */
			(__guts->_M_cur < __finish->_M_cur) : (__guts->_M_node < __finish->_M_node)) &&
		       
		       !((__guts->_M_node == __start->_M_node) ?     /* ! (*__guts < *__start) */
			 (__guts->_M_cur < __start->_M_cur) : (__guts->_M_node < __start->_M_node)),
		       _StlMsg_NOT_DEREFERENCEABLE); 
  return true;  
}

template <class _Tp>
bool __Deq_nonsingular(const void* __ptr, _Tp*) {
  typedef _Deq_iter_guts<_Tp> _Guts;
  const _Guts * __guts = (const _Guts*)(void*)__ptr;
  __stl_verbose_return(__guts->_Valid(), _StlMsg_INVALID_ITERATOR);
  
  if (__guts->_M_unsafe) return true;
  
  const _Guts* __start = (const _Guts*)(__guts->_Owner()->_Owner());
  const _Guts* __finish = __start+1;
  
  __stl_verbose_return(
		       !(((__guts->_M_node == __finish->_M_node) ? /* *__guts > *__finish */
			(__guts->_M_cur > __finish->_M_cur) : (__guts->_M_node > __finish->_M_node)) ||
		       
		       ((__guts->_M_node == __start->_M_node) ?     /* (*__guts < *__start) */
			 (__guts->_M_cur < __start->_M_cur) : (__guts->_M_node < __start->_M_node))),
		       _StlMsg_SINGULAR_ITERATOR); 
  return true;  
}

#  endif

// Non-inline member functions from _Deque_base.

template <class _Tp, class _Alloc, size_t __bufsiz>
_Deque_base<_Tp,_Alloc,__bufsiz>::~_Deque_base() {
  if (_M_map._M_data) {
    _M_destroy_nodes(_M_start._M_node, _M_finish._M_node + 1);
    _M_map.deallocate(_M_map._M_data, _M_map_size._M_data);
  }
  // should be done here instead of ~deque to ensure 
  // no detach is ever possible
  __stl_debug_do(_M_start._Invalidate());
  __stl_debug_do(_M_finish._Invalidate());
}

template <class _Tp, class _Alloc, size_t __bufsiz>
void
_Deque_base<_Tp,_Alloc,__bufsiz>::_M_initialize_map(size_t __num_elements)
{
  size_t __num_nodes = 
    __num_elements / __deque_buf_size(__bufsiz, sizeof(_Tp)) + 1;

  _M_map_size._M_data = max((size_t) _S_initial_map_size, __num_nodes + 2);
  _M_map._M_data = _M_map.allocate(_M_map_size._M_data);

  _Tp** __nstart = _M_map._M_data + (_M_map_size._M_data - __num_nodes) / 2;
  _Tp** __nfinish = __nstart + __num_nodes;
    
  __STL_TRY {
    _M_create_nodes(__nstart, __nfinish);
  }
  __STL_UNWIND((_M_map.deallocate(_M_map._M_data, _M_map_size._M_data), 
                _M_map._M_data = 0, _M_map_size._M_data = 0));
  _M_start._M_set_node(__nstart);
  _M_finish._M_set_node(__nfinish - 1);
  _M_start._M_cur = _M_start._M_first;
  _M_finish._M_cur = _M_finish._M_first +
               __num_elements % __deque_buf_size(__bufsiz, sizeof(_Tp));
}

template <class _Tp, class _Alloc, size_t __bufsiz>
void
_Deque_base<_Tp,_Alloc,__bufsiz>::_M_create_nodes(_Tp** __nstart,
                                                  _Tp** __nfinish)
{
  _Tp** __cur;
  __STL_TRY {
    for (__cur = __nstart; __cur < __nfinish; ++__cur)
      *__cur = _M_map_size.allocate(__buf_traits::_buf_size);
  }
  __STL_UNWIND(_M_destroy_nodes(__nstart, __cur));
}

template <class _Tp, class _Alloc, size_t __bufsiz>
void 
_Deque_base<_Tp,_Alloc,__bufsiz>::_M_destroy_nodes(_Tp** __nstart,
                                                   _Tp** __nfinish)
{
  for (_Tp** __n = __nstart; __n < __nfinish; ++__n)
    _M_map_size.deallocate(*__n, __buf_traits::_buf_size);
}



// Non-inline member functions

# if defined ( __STL_NESTED_TYPE_PARAM_BUG )
// qualified references 
#   define __iterator__           _Deque_iterator<_Tp, _Nonconst_traits<_Tp>, _Buf_size_traits<_Tp, __bufsiz> >
#   define const_iterator         _Deque_iterator<_Tp, _Const_traits<_Tp>, _Buf_size_traits<_Tp, __bufsiz> > 
#   define iterator               __iterator__
#   define size_type              size_t
#   define value_type             _Tp
# else
#  define __iterator__           __STL_TYPENAME_ON_RETURN_TYPE deque<_Tp, _Alloc, __bufsiz>::iterator
# endif

template <class _Tp, class _Alloc, size_t __bufsiz>
deque<_Tp, _Alloc, __bufsiz>&  
deque<_Tp, _Alloc, __bufsiz>::operator= (const deque<_Tp, _Alloc, __bufsiz>& __x) {
  const size_type __len = size();
  if (&__x != this) {
    if (__len >= __x.size())
      erase(copy(__x.begin(), __x.end(), _M_start), _M_finish);
    else {
      const_iterator __mid = __x.begin() + difference_type(__len);
      copy(__x.begin(), __mid, _M_start);
      insert(_M_finish, __mid, __x.end());
    }
  }
  __stl_debug_do(_Invalidate_all());
  return *this;
}        

template <class _Tp, class _Alloc, size_t __bufsiz>
void 
deque<_Tp, _Alloc, __bufsiz>::_M_fill_insert(iterator __pos,
					     size_type __n, const value_type& __x)
{
  __stl_debug_check(__check_if_owner(&_M_iter_list, __pos));
  if (__pos._M_cur == _M_start._M_cur) {
    iterator __new_start = _M_reserve_elements_at_front(__n);
    __STL_TRY {
      uninitialized_fill(__new_start, _M_start, __x);
    }
    __STL_UNWIND(_M_destroy_nodes(__new_start._M_node, _M_start._M_node));
    _M_start = __new_start;
    __stl_debug_do(_M_orphan_start());
  }
  else if (__pos._M_cur == _M_finish._M_cur) {
    iterator __new_finish = _M_reserve_elements_at_back(__n);
    __STL_TRY {
      uninitialized_fill(_M_finish, __new_finish, __x);
    }
    __STL_UNWIND(_M_destroy_nodes(_M_finish._M_node+1, __new_finish._M_node+1));
    _M_finish = __new_finish;
    __stl_debug_do(_M_orphan_finish());
  }
  else 
    _M_insert_aux(__pos, __n, __x);
}

#ifndef __STL_MEMBER_TEMPLATES  

template <class _Tp, class _Alloc, size_t __bufsiz>
void deque<_Tp, _Alloc, __bufsiz>::insert(iterator __pos,
                                           const value_type* __first,
                                           const value_type* __last) {
  __stl_debug_check(__check_if_owner(&_M_iter_list, __pos));
  size_type __n = __last - __first;
  if (__pos._M_cur == _M_start._M_cur) {
    iterator __new_start = _M_reserve_elements_at_front(__n);
    __STL_TRY {
      uninitialized_copy(__first, __last, __new_start);
    }
    __STL_UNWIND(_M_destroy_nodes(__new_start._M_node, _M_start._M_node));
    _M_start = __new_start;
    __stl_debug_do(_M_orphan_start());
  }
  else if (__pos._M_cur == _M_finish._M_cur) {
    iterator __new_finish = _M_reserve_elements_at_back(__n);
    __STL_TRY {
      uninitialized_copy(__first, __last, _M_finish);
    }
    __STL_UNWIND(_M_destroy_nodes(_M_finish._M_node + 1, 
                                  __new_finish._M_node + 1));
    _M_finish = __new_finish;
    __stl_debug_do(_M_orphan_finish());
  }
  else
    _M_insert_aux(__pos, __first, __last, __n);
}

template <class _Tp, class _Alloc, size_t __bufsiz>
void deque<_Tp,_Alloc,__bufsiz>::insert(iterator __pos,
                                         const_iterator __first,
                                         const_iterator __last)
{
  __stl_debug_check(__check_if_owner(&_M_iter_list, __pos));
  size_type __n = __last - __first;
  if (__pos._M_cur == _M_start._M_cur) {
    iterator __new_start = _M_reserve_elements_at_front(__n);
    __STL_TRY {
      uninitialized_copy(__first, __last, __new_start);
    }
    __STL_UNWIND(_M_destroy_nodes(__new_start._M_node, _M_start._M_node));
    _M_start = __new_start;
    __stl_debug_do(_M_orphan_start());
  }
  else if (__pos._M_cur == _M_finish._M_cur) {
    iterator __new_finish = _M_reserve_elements_at_back(__n);
    __STL_TRY {
      uninitialized_copy(__first, __last, _M_finish);
    }
    __STL_UNWIND(_M_destroy_nodes(_M_finish._M_node + 1,__new_finish._M_node + 1));
    _M_finish = __new_finish;
    __stl_debug_do(_M_orphan_finish());
  }
  else
    _M_insert_aux(__pos, __first, __last, __n);
}

#endif /* __STL_MEMBER_TEMPLATES */

template <class _Tp, class _Alloc, size_t __bufsiz>
__iterator__ 
deque<_Tp,_Alloc,__bufsiz>::erase(iterator __first, iterator __last)
{
  __stl_debug_check(__check_if_owner(&_M_iter_list, __first) && __check_range(__first,__last));
  if (__first == _M_start && __last == _M_finish) {
    clear();
    return _M_finish;
  }
  else {
    difference_type __n = __last - __first;
    difference_type __elems_before = __first - _M_start;
    if (__elems_before < difference_type(size() - __n) / 2) {
      copy_backward(_M_start, __first, __last);
      iterator __new_start = _M_start + __n;
      destroy(_M_start, __new_start);
      _M_destroy_nodes(__new_start._M_node, _M_start._M_node);
      _M_start = __new_start;
      __stl_debug_do(_M_orphan_start());
    }
    else {
      copy(__last, _M_finish, __first);
      iterator __new_finish = _M_finish - __n;
      destroy(__new_finish, _M_finish);
      _M_destroy_nodes(__new_finish._M_node + 1, _M_finish._M_node + 1);
      _M_finish = __new_finish;
      __stl_debug_do(_M_orphan_finish());
    }
    return _M_start + __elems_before;
  }
}

template <class _Tp, class _Alloc, size_t __bufsiz>
void deque<_Tp,_Alloc,__bufsiz>::clear()
{
  __stl_debug_do(_Invalidate_all());
  for (_Map_pointer __node = _M_start._M_node + 1;
       __node < _M_finish._M_node;
       ++__node) {
    destroy(*__node, *__node + __buf_traits::_buf_size);
    _M_map_size.deallocate(*__node, __buf_traits::_buf_size);
  }

  if (_M_start._M_node != _M_finish._M_node) {
    destroy(_M_start._M_cur, _M_start._M_last);
    destroy(_M_finish._M_first, _M_finish._M_cur);
    _M_map_size.deallocate(_M_finish._M_first, __buf_traits::_buf_size);
  }
  else
    destroy(_M_start._M_cur, _M_finish._M_cur);

  _M_finish = _M_start;
  __stl_debug_do(_M_orphan_finish());
}

// Precondition: _M_start and _M_finish have already been initialized,
// but none of the deque's elements have yet been constructed.
template <class _Tp, class _Alloc, size_t __bufsiz>
void 
deque<_Tp,_Alloc,__bufsiz>::_M_fill_initialize(const value_type& __value) {
  _Map_pointer __cur;
  __stl_debug_do(_M_iter_list._Safe_init(&_M_start));
  __stl_debug_do(_Init_bounds());
  __STL_TRY {
    for (__cur = _M_start._M_node; __cur < _M_finish._M_node; ++__cur)
      uninitialized_fill(*__cur, *__cur + __buf_traits::_buf_size, __value);
    uninitialized_fill(_M_finish._M_first, _M_finish._M_cur, __value);
  }
# ifdef __STL_DEBUG
  __STL_UNWIND(destroy(_M_start, iterator(&_M_iter_list,*__cur, __cur)));
# else
  __STL_UNWIND(destroy(_M_start, iterator(*__cur, __cur)));
# endif
}


// Called only if _M_finish._M_cur == _M_finish._M_last - 1.
template <class _Tp, class _Alloc, size_t __bufsiz>
void
deque<_Tp,_Alloc,__bufsiz>::_M_push_back_aux(const value_type& __t)
{
  value_type __t_copy = __t;
  _M_reserve_map_at_back();
  *(_M_finish._M_node + 1) = _M_map_size.allocate(__buf_traits::_buf_size);
  __STL_TRY {
    construct(_M_finish._M_cur, __t_copy);
    _M_finish._M_set_node(_M_finish._M_node + 1);
    _M_finish._M_cur = _M_finish._M_first;
  }
  __STL_UNWIND(_M_map_size.deallocate(*(_M_finish._M_node + 1), 
				      __buf_traits::_buf_size));
}

// Called only if _M_finish._M_cur == _M_finish._M_last - 1.
template <class _Tp, class _Alloc, size_t __bufsiz>
void
deque<_Tp,_Alloc,__bufsiz>::_M_push_back_aux()
{
  _M_reserve_map_at_back();
  *(_M_finish._M_node + 1) = _M_map_size.allocate(__buf_traits::_buf_size);
  __STL_TRY {
    construct(_M_finish._M_cur);
    _M_finish._M_set_node(_M_finish._M_node + 1);
    _M_finish._M_cur = _M_finish._M_first;
  }
  __STL_UNWIND(_M_map_size.deallocate(*(_M_finish._M_node + 1), 
				      __buf_traits::_buf_size));
}

// Called only if _M_start._M_cur == _M_start._M_first.
template <class _Tp, class _Alloc, size_t __bufsiz>
void 
deque<_Tp,_Alloc,__bufsiz>::_M_push_front_aux(const value_type& __t)
{
  value_type __t_copy = __t;
  _M_reserve_map_at_front();
  *(_M_start._M_node - 1) = _M_map_size.allocate(__buf_traits::_buf_size);
  __STL_TRY {
    _M_start._M_set_node(_M_start._M_node - 1);
    _M_start._M_cur = _M_start._M_last - 1;
    construct(_M_start._M_cur, __t_copy);
  }
  __STL_UNWIND((++_M_start, 
		_M_map_size.deallocate(*(_M_start._M_node - 1), __buf_traits::_buf_size)));
} 

// Called only if _M_start._M_cur == _M_start._M_first.
template <class _Tp, class _Alloc, size_t __bufsiz>
void 
deque<_Tp,_Alloc,__bufsiz>::_M_push_front_aux()
{
  _M_reserve_map_at_front();
  *(_M_start._M_node - 1) = _M_map_size.allocate(__buf_traits::_buf_size);
  __STL_TRY {
    _M_start._M_set_node(_M_start._M_node - 1);
    _M_start._M_cur = _M_start._M_last - 1;
    construct(_M_start._M_cur);
  }
  __STL_UNWIND((++_M_start, _M_map_size.deallocate(*(_M_start._M_node - 1), 
						   __buf_traits::_buf_size )));
} 

// Called only if _M_finish._M_cur == _M_finish._M_first.
template <class _Tp, class _Alloc, size_t __bufsiz>
void 
deque<_Tp,_Alloc,__bufsiz>::_M_pop_back_aux()
{
  _M_map_size.deallocate(_M_finish._M_first, __buf_traits::_buf_size);
  _M_finish._M_set_node(_M_finish._M_node - 1);
  _M_finish._M_cur = _M_finish._M_last - 1;
  destroy(_M_finish._M_cur);
}

// Called only if _M_start._M_cur == _M_start._M_last - 1.  Note that 
// if the deque has at least one element (a precondition for this member 
// function), and if _M_start._M_cur == _M_start._M_last, then the deque 
// must have at least two nodes.
template <class _Tp, class _Alloc, size_t __bufsiz>
void 
deque<_Tp,_Alloc,__bufsiz>::_M_pop_front_aux()
{
  destroy(_M_start._M_cur);
  _M_map_size.deallocate(_M_start._M_first, __buf_traits::_buf_size);
  _M_start._M_set_node(_M_start._M_node + 1);
  _M_start._M_cur = _M_start._M_first;
}      



template <class _Tp, class _Alloc, size_t __bufsiz>
__iterator__
deque<_Tp,_Alloc,__bufsiz>::_M_insert_aux_prepare(__iterator__ __pos) {
  difference_type __index = __pos - _M_start;
  if (__index < difference_type(size() / 2)) {
    push_front(front());
    iterator __front1 = _M_start;
    ++__front1;
    iterator __front2 = __front1;
    ++__front2;
    __pos = _M_start + __index;
    iterator __pos1 = __pos;
    ++__pos1;
    copy(__front2, __pos1, __front1);
  }
  else {
    push_back(back());
    iterator __back1 = _M_finish;
    --__back1;
    iterator __back2 = __back1;
    --__back2;
    __pos = _M_start + __index;
    copy_backward(__pos, __back2, __back1);
  }
  return __pos;
}

template <class _Tp, class _Alloc, size_t __bufsiz>
__iterator__
deque<_Tp,_Alloc,__bufsiz>::_M_insert_aux(__iterator__ __pos,
                                           const value_type& __x) {
  value_type __x_copy = __x;
  __pos = _M_insert_aux_prepare(__pos);
  *__pos = __x_copy;
  return __pos;
}

template <class _Tp, class _Alloc, size_t __bufsiz>
__iterator__
deque<_Tp,_Alloc,__bufsiz>::_M_insert_aux(__iterator__ __pos)
{
  __pos = _M_insert_aux_prepare(__pos);
  *__pos = value_type();
  return __pos;
}

template <class _Tp, class _Alloc, size_t __bufsiz>
void
deque<_Tp,_Alloc,__bufsiz>::_M_insert_aux(iterator __pos,
                                           size_type __n,
                                           const value_type& __x)
{
  const difference_type __elems_before = __pos - _M_start;
  size_type __length = size();
  value_type __x_copy = __x;
  if (__elems_before < difference_type(__length / 2)) {
    iterator __new_start = _M_reserve_elements_at_front(__n);
    iterator __old_start = _M_start;
    __pos = _M_start + __elems_before;
    __STL_TRY {
      if (__elems_before >= difference_type(__n)) {
        iterator __start_n = _M_start + difference_type(__n);
        uninitialized_copy(_M_start, __start_n, __new_start);
        _M_start = __new_start;
	__stl_debug_do(_M_orphan_start());
        copy(__start_n, __pos, __old_start);
        fill(__pos - difference_type(__n), __pos, __x_copy);
      }
      else {
        __uninitialized_copy_fill(_M_start, __pos, __new_start, 
	                          _M_start, __x_copy);
        _M_start = __new_start;
	__stl_debug_do(_M_orphan_start());
        fill(__old_start, __pos, __x_copy);
      }
    }
    __STL_UNWIND(_M_destroy_nodes(__new_start._M_node, _M_start._M_node));
  }
  else {
    iterator __new_finish = _M_reserve_elements_at_back(__n);
    iterator __old_finish = _M_finish;
    const difference_type __elems_after = 
      difference_type(__length) - __elems_before;
    __pos = _M_finish - __elems_after;
    __STL_TRY {
      if (__elems_after > difference_type(__n)) {
        iterator __finish_n = _M_finish - difference_type(__n);
        uninitialized_copy(__finish_n, _M_finish, _M_finish);
        _M_finish = __new_finish;
	__stl_debug_do(_M_orphan_finish());
        copy_backward(__pos, __finish_n, __old_finish);
        fill(__pos, __pos + difference_type(__n), __x_copy);
      }
      else {
        __uninitialized_fill_copy(_M_finish, __pos + difference_type(__n),
                                  __x_copy, __pos, _M_finish);
        _M_finish = __new_finish;
	__stl_debug_do(_M_orphan_finish());	
        fill(__pos, __old_finish, __x_copy);
      }
    }
    __STL_UNWIND(_M_destroy_nodes(_M_finish._M_node + 1, __new_finish._M_node + 1));
  }
  __stl_debug_do(_Invalidate_all());        
}

#ifndef __STL_MEMBER_TEMPLATES 
template <class _Tp, class _Alloc, size_t __bufsiz>
void 
deque<_Tp,_Alloc,__bufsiz>::_M_insert_aux(iterator __pos,
                                           const value_type* __first,
                                           const value_type* __last,
                                           size_type __n)
{

  const difference_type __elemsbefore = __pos - _M_start;
  size_type __length = size();
  if (__elemsbefore < difference_type(__length / 2)) {
    iterator __new_start = _M_reserve_elements_at_front(__n);
    iterator __old_start = _M_start;
    __pos = _M_start + __elemsbefore;
    __STL_TRY {
      if (__elemsbefore >= difference_type(__n)) {
        iterator __start_n = _M_start + difference_type(__n);
        uninitialized_copy(_M_start, __start_n, __new_start);
        _M_start = __new_start;
	__stl_debug_do(_M_orphan_start());
        copy(__start_n, __pos, __old_start);
        copy(__first, __last, __pos - difference_type(__n));
      }
      else {
        const value_type* __mid = 
	  __first + (difference_type(__n) - __elemsbefore);
        __uninitialized_copy_copy(_M_start, __pos, __first, __mid,
                                  __new_start);
        _M_start = __new_start;
	__stl_debug_do(_M_orphan_start());
        copy(__mid, __last, __old_start);
      }
    }
    __STL_UNWIND(_M_destroy_nodes(__new_start._M_node, _M_start._M_node));
  }
  else {
    iterator __new_finish = _M_reserve_elements_at_back(__n);
    iterator __old_finish = _M_finish;
    const difference_type __elemsafter = 
      difference_type(__length) - __elemsbefore;
    __pos = _M_finish - __elemsafter;
    __STL_TRY {
      if (__elemsafter > difference_type(__n)) {
        iterator __finish_n = _M_finish - difference_type(__n);
        uninitialized_copy(__finish_n, _M_finish, _M_finish);
        _M_finish = __new_finish;
	__stl_debug_do(_M_orphan_finish());
        copy_backward(__pos, __finish_n, __old_finish);
        copy(__first, __last, __pos);
      }
      else {
        const value_type* __mid = __first + __elemsafter;
        __uninitialized_copy_copy(__mid, __last, __pos, _M_finish, _M_finish);
        _M_finish = __new_finish;
	__stl_debug_do(_M_orphan_finish());
        copy(__first, __mid, __pos);
      }
    }
    __STL_UNWIND(_M_destroy_nodes(_M_finish._M_node + 1, __new_finish._M_node + 1));
  }
  __stl_debug_do(_Invalidate_all());        
}

template <class _Tp, class _Alloc, size_t __bufsiz>
void
deque<_Tp,_Alloc,__bufsiz>::_M_insert_aux(iterator __pos,
                                           const_iterator __first,
                                           const_iterator __last,
                                           size_type __n)
{
  const difference_type __elemsbefore = __pos - _M_start;
  size_type __length = size();
  if (__elemsbefore < difference_type(__length / 2)) {
    iterator __new_start = _M_reserve_elements_at_front(__n);
    iterator __old_start = _M_start;
    __pos = _M_start + __elemsbefore;
    __STL_TRY {
      if (__elemsbefore >= difference_type(__n)) {
        iterator __start_n = _M_start + __n;
        uninitialized_copy(_M_start, __start_n, __new_start);
        _M_start = __new_start;
	__stl_debug_do(_M_orphan_start());
        copy(__start_n, __pos, __old_start);
        copy(__first, __last, __pos - difference_type(__n));
      }
      else {
        const_iterator __mid = __first + (__n - __elemsbefore);
        __uninitialized_copy_copy(_M_start, __pos, __first, __mid,
                                  __new_start);
        _M_start = __new_start;
	__stl_debug_do(_M_orphan_start());
        copy(__mid, __last, __old_start);
      }
    }
    __STL_UNWIND(_M_destroy_nodes(__new_start._M_node, _M_start._M_node));
  }
  else {
    iterator __new_finish = _M_reserve_elements_at_back(__n);
    iterator __old_finish = _M_finish;
    const difference_type __elemsafter = __length - __elemsbefore;
    __pos = _M_finish - __elemsafter;
    __STL_TRY {
      if (__elemsafter > difference_type(__n)) {
        iterator __finish_n = _M_finish - difference_type(__n);
        uninitialized_copy(__finish_n, _M_finish, _M_finish);
        _M_finish = __new_finish;
	__stl_debug_do(_M_orphan_finish());
        copy_backward(__pos, __finish_n, __old_finish);
        copy(__first, __last, __pos);
      }
      else {
        const_iterator __mid = __first + __elemsafter;
        __uninitialized_copy_copy(__mid, __last, __pos, _M_finish, _M_finish);
        _M_finish = __new_finish;
	__stl_debug_do(_M_orphan_finish());
        copy(__first, __mid, __pos);
      }
    }
    __STL_UNWIND(_M_destroy_nodes(_M_finish._M_node + 1, __new_finish._M_node + 1));
  }
  __stl_debug_do(_Invalidate_all());        
}

#endif /* __STL_MEMBER_TEMPLATES */

template <class _Tp, class _Alloc, size_t __bufsiz>
void 
deque<_Tp,_Alloc,__bufsiz>::_M_new_elements_at_front(size_type __new_elems)
{
  size_type __new_nodes
      = (__new_elems + __buf_traits::_buf_size - 1) / __buf_traits::_buf_size;
  _M_reserve_map_at_front(__new_nodes);
  size_type __i;
  __STL_TRY {
    for (__i = 1; __i <= __new_nodes; ++__i)
      *(_M_start._M_node - __i) = _M_map_size.allocate(__buf_traits::_buf_size);
  }
#       ifdef __STL_USE_EXCEPTIONS
  catch(...) {
    for (size_type __j = 1; __j < __i; ++__j)
      _M_map_size.deallocate(*(_M_start._M_node - __j), __buf_traits::_buf_size);
    throw;
  }
#       endif /* __STL_USE_EXCEPTIONS */
}

template <class _Tp, class _Alloc, size_t __bufsiz>
void 
deque<_Tp,_Alloc,__bufsiz>::_M_new_elements_at_back(size_type __new_elems)
{
  size_type __new_nodes
      = (__new_elems + __buf_traits::_buf_size - 1) / __buf_traits::_buf_size;
  _M_reserve_map_at_back(__new_nodes);
  size_type __i;
  __STL_TRY {
    for (__i = 1; __i <= __new_nodes; ++__i)
      *(_M_finish._M_node + __i) = _M_map_size.allocate(__buf_traits::_buf_size);
  }
#       ifdef __STL_USE_EXCEPTIONS
  catch(...) {
    for (size_type __j = 1; __j < __i; ++__j)
      _M_map_size.deallocate(*(_M_finish._M_node + __j), __buf_traits::_buf_size);
    throw;
  }
#       endif /* __STL_USE_EXCEPTIONS */
}

template <class _Tp, class _Alloc, size_t __bufsiz>
void 
deque<_Tp,_Alloc,__bufsiz>::_M_reallocate_map(size_type __nodes_to_add,
                                              bool __add_at_front)
{
  size_type __old_num_nodes = _M_finish._M_node - _M_start._M_node + 1;
  size_type __new_num_nodes = __old_num_nodes + __nodes_to_add;

  _Map_pointer __new_nstart;
  if (_M_map_size._M_data > 2 * __new_num_nodes) {
    __new_nstart = _M_map._M_data + (_M_map_size._M_data - __new_num_nodes) / 2 
                     + (__add_at_front ? __nodes_to_add : 0);
    if (__new_nstart < _M_start._M_node)
      copy(_M_start._M_node, _M_finish._M_node + 1, __new_nstart);
    else
      copy_backward(_M_start._M_node, _M_finish._M_node + 1, 
                    __new_nstart + __old_num_nodes);
  }
  else {
    size_type __new_map_size = 
      _M_map_size._M_data + max((size_t)_M_map_size._M_data, __nodes_to_add) + 2;

    _Map_pointer __new_map = _M_map.allocate(__new_map_size);
    __new_nstart = __new_map + (__new_map_size - __new_num_nodes) / 2
                         + (__add_at_front ? __nodes_to_add : 0);
    copy(_M_start._M_node, _M_finish._M_node + 1, __new_nstart);
    _M_map.deallocate(_M_map._M_data, _M_map_size._M_data);

    _M_map._M_data = __new_map;
    _M_map_size._M_data = __new_map_size;
  }

  _M_start._M_set_node(__new_nstart);
  _M_finish._M_set_node(__new_nstart + __old_num_nodes - 1);
  __stl_debug_do(_Invalidate_all());        
}


#if defined (__STL_MEMBER_TEMPLATES) && ! defined (__STL_INLINE_MEMBER_TEMPLATES )
template <class _Tp, class _Alloc, size_t __bufsize>
template <class _InputIterator>
void deque<_Tp,_Alloc,__bufsize>::_M_range_initialize(_InputIterator __first,
						      _InputIterator __last,
						      input_iterator_tag) {
    _M_initialize_map(0);
    __stl_debug_do(_M_iter_list._Safe_init(&_M_start));
    __stl_debug_do(_Init_bounds());
    __STL_TRY {
      for ( ; __first != __last; ++__first)
	push_back(*__first);
    }
    __STL_UNWIND(clear());
  }
  
template <class _Tp, class _Alloc, size_t __bufsize>
template <class _ForwardIterator>
void  deque<_Tp,_Alloc,__bufsize>::_M_range_initialize(_ForwardIterator __first,
						       _ForwardIterator __last,
						       forward_iterator_tag) {
   size_type __n = 0;
   distance(__first, __last, __n);
   __stl_debug_do(_M_iter_list._Safe_init(&_M_start));
   __stl_debug_do(_Init_bounds());
   _M_initialize_map(__n);
   _Map_pointer __cur_node;
   __STL_TRY {
    for (__cur_node = _M_start._M_node; 
         __cur_node < _M_finish._M_node; 
	 ++__cur_node) {
      _ForwardIterator __mid = __first;
      advance(__mid, __buf_traits::_buf_size);
      uninitialized_copy(__first, __mid, *__cur_node);
      __first = __mid;
    }
    uninitialized_copy(__first, __last, _M_finish._M_first);
   }
# ifdef __STL_DEBUG
  __STL_UNWIND(destroy(_M_start, iterator(&_M_iter_list, *__cur_node, __cur_node)));
# else
  __STL_UNWIND(destroy(_M_start, iterator(*__cur_node, __cur_node)));
# endif
 }

template <class _Tp, class _Alloc, size_t __bufsize>
template <class _ForwardIterator>
void  deque<_Tp,_Alloc,__bufsize>::insert(iterator __pos,
					  _ForwardIterator __first,
					  _ForwardIterator __last,
					  forward_iterator_tag)
 {
  __stl_debug_check(__check_if_owner(&_M_iter_list, __pos));
  size_type __n = 0;
  distance(__first, __last, __n);
  if (__pos._M_cur == _M_start._M_cur) {
    iterator __new_start = _M_reserve_elements_at_front(__n);
    __STL_TRY {
      uninitialized_copy(__first, __last, __new_start);
      _M_start = __new_start;
      __stl_debug_do(_M_orphan_start());
    }
    __STL_UNWIND(_M_destroy_nodes(__new_start._M_node, _M_start._M_node));
  }
  else if (__pos._M_cur == _M_finish._M_cur) {
    iterator __new_finish = _M_reserve_elements_at_back(__n);
    __STL_TRY {
      uninitialized_copy(__first, __last, _M_finish);
      _M_finish = __new_finish;
      __stl_debug_do(_M_orphan_finish());
    }
    __STL_UNWIND(_M_destroy_nodes(_M_finish._M_node + 1, __new_finish._M_node + 1));
  }
  else
    _M_insert_aux(__pos, __first, __last, __n);
}

template <class _Tp, class _Alloc, size_t __bufsize>
template <class _ForwardIterator>
void deque<_Tp,_Alloc,__bufsize>::_M_insert_aux(iterator __pos,
						_ForwardIterator __first,
						_ForwardIterator __last,
						size_type __n)
{
    
    const difference_type __elemsbefore = __pos - _M_start;
    size_type __length = size();
    if (__elemsbefore < difference_type(__length / 2)) {
      iterator __new_start = _M_reserve_elements_at_front(__n);
      iterator __old_start = _M_start;
      __pos = _M_start + __elemsbefore;
      __STL_TRY {
	if (__elemsbefore >= difference_type(__n)) {
	  iterator __start_n = _M_start + difference_type(__n); 
	  uninitialized_copy(_M_start, __start_n, __new_start);
	  _M_start = __new_start;
	  __stl_debug_do(_M_orphan_start());
	  copy(__start_n, __pos, __old_start);
	  copy(__first, __last, __pos - difference_type(__n));
	}
	else {
	  _ForwardIterator __mid = __first;
	  advance(__mid, difference_type(__n) - __elemsbefore);
	  __uninitialized_copy_copy(_M_start, __pos, __first, __mid,
				    __new_start);
	  _M_start = __new_start;
	  __stl_debug_do(_M_orphan_start());
	  copy(__mid, __last, __old_start);
	}
      }
      __STL_UNWIND(_M_destroy_nodes(__new_start._M_node, _M_start._M_node));
    }
    else {
      iterator __new_finish = _M_reserve_elements_at_back(__n);
      iterator __old_finish = _M_finish;
      const difference_type __elemsafter = 
	difference_type(__length) - __elemsbefore;
      __pos = _M_finish - __elemsafter;
      __STL_TRY {
      if (__elemsafter > difference_type(__n)) {
        iterator __finish_n = _M_finish - difference_type(__n);
        uninitialized_copy(__finish_n, _M_finish, _M_finish);
        _M_finish = __new_finish;
	__stl_debug_do(_M_orphan_finish());
        copy_backward(__pos, __finish_n, __old_finish);
        copy(__first, __last, __pos);
      }
      else {
        _ForwardIterator __mid = __first;
        advance(__mid, __elemsafter);
        __uninitialized_copy_copy(__mid, __last, __pos, _M_finish, _M_finish);
        _M_finish = __new_finish;
	__stl_debug_do(_M_orphan_finish());
        copy(__first, __mid, __pos);
      }
      }
      __STL_UNWIND(_M_destroy_nodes(_M_finish._M_node + 1, __new_finish._M_node + 1));
    }
    __stl_debug_do(_Invalidate_all());        
  }
# endif /* __STL_MEMBER_TEMPLATES */

__STL_END_NAMESPACE

# undef deque
# undef __iterator__
# undef iterator
# undef const_iterator
# undef size_type
# undef value_type

#if defined(__sgi) && !defined(__GNUC__) && (_MIPS_SIM != _MIPS_SIM_ABI32)
#pragma reset woff 1174
#pragma reset woff 1375
#endif

#endif /*  __STL_DEQUE_C */

// Local Variables:
// mode:C++
// End:
