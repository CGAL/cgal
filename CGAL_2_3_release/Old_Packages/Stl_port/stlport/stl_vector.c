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
#ifndef __STL_VECTOR_C
#define __STL_VECTOR_C

#if defined(__sgi) && !defined(__GNUC__) && (_MIPS_SIM != _MIPS_SIM_ABI32)
#pragma set woff 1174
#pragma set woff 1375
#endif

#  undef  vector
#  define vector __original_vector
//#  define vector __WORKAROUND_RENAME(vector)

# if defined ( __STL_NESTED_TYPE_PARAM_BUG )
#  define __pointer__             _Tp*
#  define __const_pointer__       const _Tp*
#  define __size_type__           size_t
#  define __difference_type__     ptrdiff_t
# else
#  define __pointer__         pointer
#  define __const_pointer__   const_pointer
#  define __size_type__       size_type
#  define __difference_type__ difference_type
# endif

# undef _Make_ptr
# if defined (__STL_DEBUG)
#  define _Make_iterator(__i) iterator(&_M_iter_list, __i)
#  define _Make_const_iterator(__i) const_iterator(&_M_iter_list, __i)
#  define _Make_ptr(__i)   __i._M_iterator
#  define __iterator__       _Vec_iter<_Tp, _Nonconst_traits<_Tp> >
#  define __const_iterator__ _Vec_iter<_Tp, _Const_traits<_Tp> >
# else
#  define __iterator__       __pointer__
#  define __const_iterator__ __const_pointer__  
#  define _Make_iterator(__i) __i
#  define _Make_const_iterator(__i) __i
#  define _Make_ptr(__i)   __i
# endif

__STL_BEGIN_NAMESPACE

template <class _Tp, class _Alloc>
void 
vector<_Tp, _Alloc>::reserve(__size_type__ __n) {
  if (capacity() < __n) {
    const size_type __old_size = size();
    pointer __tmp;
    if (_M_start) {
      __tmp = _M_allocate_and_copy(__n, _M_start, _M_finish);
      __STLPORT_STD::destroy(_M_start, _M_finish);
      _M_end_of_storage.deallocate(_M_start, _M_end_of_storage._M_data - _M_start);
     } else {
      __tmp = _M_end_of_storage.allocate(__n);
    }
    _M_start = __tmp;
    _M_finish = __tmp + __old_size;
    _M_end_of_storage._M_data = _M_start + __n;
  }
}

#if defined (__STL_MEMBER_TEMPLATES) && ! defined (__STL_INLINE_MEMBER_TEMPLATES)

template <class _Tp, class _Alloc>  template <class _ForwardIter>
void 
vector<_Tp, _Alloc>::_M_assign_aux(_ForwardIter __first, _ForwardIter __last,
				   forward_iterator_tag) {
  size_type __len = 0;
  distance(__first, __last, __len);
    
  if (__len > capacity()) {
    iterator __tmp = _M_allocate_and_copy(__len, __first, __last);
    __STLPORT_STD::destroy(_M_start, _M_finish);
    _M_end_of_storage.deallocate(_M_start, _M_end_of_storage._M_data - _M_start);
    _M_start = __tmp;
    _M_end_of_storage._M_data = _M_finish = _M_start + __len;
  }
  else if (size() >= __len) {
    iterator __new_finish = copy(__first, __last, _M_start);
    __STLPORT_STD::destroy(__new_finish, _M_finish);
    _M_finish = __new_finish;
  }
  else {
    _ForwardIter __mid = __first;
    advance(__mid, size());
    copy(__first, __mid, _M_start);
    _M_finish = uninitialized_copy(__mid, __last, _M_finish);
  }
}

#endif /* __STL_MEMBER_TEMPLATES */

#ifndef __STL_INLINE_MEMBER_TEMPLATES
#ifndef __STL_MEMBER_TEMPLATES
template <class _Tp, class _Alloc>
void 
vector<_Tp, _Alloc>::insert(
# if defined ( __STL_DEBUG )
			    __iterator__ __pos, 
			    __const_pointer__ __first, 
			    __const_pointer__ __last
# else
			    __iterator__ __position, 
			    __const_iterator__ __first, 
			    __const_iterator__ __last
# endif
			    )
#else /* MEMBER_TEMPLATES */

template <class _Tp, class _Alloc> template <class _ForwardIterator>
void  vector<_Tp, _Alloc>::_M_range_insert(
# ifdef __STL_DEBUG
					   __iterator__ __pos,
# else
					   __iterator__ __position,
# endif
					   _ForwardIterator __first,
					   _ForwardIterator __last,
					   forward_iterator_tag)
# endif /* MEMBER_TEMPLATES */

{

#ifdef __STL_DEBUG
  __stl_debug_check(__check_if_owner(&_M_iter_list, __pos));
  pointer __position(_Make_ptr(__pos));
#endif
  if (__first != __last) {
    __stl_debug_check(__check_range(__first,__last));
    size_type __n = 0;
    distance(__first, __last, __n);
    if (size_type(_M_end_of_storage._M_data - _M_finish) >= __n) {
      const size_type __elems_after = _M_finish - __position;
      pointer __old_finish = _M_finish;
      if (__elems_after > __n) {
	uninitialized_copy(_M_finish - __n, _M_finish, _M_finish);
	_M_finish += __n;
	copy_backward(__position, __old_finish - __n, __old_finish);
	copy(__first, __last, __position);
      }
      else {
# ifdef __STL_MEMBER_TEMPLATES
	_ForwardIterator __mid = __first;
	advance(__mid, __elems_after);
# else
	__const_pointer__ __mid = __first + __elems_after;
# endif
	uninitialized_copy(__mid, __last, _M_finish);
	_M_finish += __n - __elems_after;
	uninitialized_copy(__position, __old_finish, _M_finish);
	_M_finish += __elems_after;
	copy(__first, __mid, __position);
      }
      __stl_debug_do(__invalidate_range(&_M_iter_list, __pos, end()));
    }
    else {
      const size_type __old_size = size();
      const size_type __len = __old_size + max(__old_size, __n);
      pointer __new_start = _M_end_of_storage.allocate(__len);
      pointer __new_finish = __new_start;
      __STL_TRY {
	__new_finish = uninitialized_copy(_M_start, __position, __new_start);
	__new_finish = uninitialized_copy(__first, __last, __new_finish);
	__new_finish
	  = uninitialized_copy(__position, _M_finish, __new_finish);
      }
      __STL_UNWIND((__STLPORT_STD::destroy(__new_start,__new_finish), 
		    _M_end_of_storage.deallocate(__new_start,__len)));
      __STLPORT_STD::destroy(_M_start, _M_finish);
      _M_end_of_storage.deallocate(_M_start, _M_end_of_storage._M_data - _M_start);
      _M_start = __new_start;
      _M_finish = __new_finish;
      _M_end_of_storage._M_data = __new_start + __len;
    }
  }
}

#endif /* __STL_INLINE_MEMBER_TEMPLATES */


template <class _Tp, class _Alloc>
vector<_Tp,_Alloc>& 
vector<_Tp,_Alloc>::operator=(const vector<_Tp, _Alloc>& __x)
{
  if (&__x != this) {
    const size_type __xlen = __x.size();
    __stl_debug_do(_M_iter_list._Invalidate_all());
    if (__xlen > capacity()) {
      pointer __tmp = _M_allocate_and_copy(__xlen, (const_pointer)__x._M_start+0, (const_pointer)__x._M_finish+0);
      __STLPORT_STD::destroy(_M_start, _M_finish);
      _M_end_of_storage.deallocate(_M_start, _M_end_of_storage._M_data - _M_start);
      _M_start = __tmp;
      _M_end_of_storage._M_data = _M_start + __xlen;
    }
    else if (size() >= __xlen) {
      pointer __i = copy((const_pointer)__x._M_start+0, (const_pointer)__x._M_finish+0, _M_start);
      __STLPORT_STD::destroy(__i, _M_finish);
    }
    else {
      copy((const_pointer)__x._M_start, (const_pointer)__x._M_start + size(), _M_start);
      uninitialized_copy((const_pointer)__x._M_start + size(), (const_pointer)__x._M_finish+0, _M_finish);
    }
    _M_finish = _M_start + __xlen;
  }
  return *this;
}

template <class _Tp, class _Alloc>
void vector<_Tp, _Alloc>::_M_fill_assign(size_t __n, const _Tp& __val) {
  if (__n > capacity()) {
    vector<_Tp, _Alloc> __tmp(__n, __val, get_allocator());
    __tmp.swap(*this);
  }
  else if (__n > size()) {
    fill(begin(), end(), __val);
    _M_finish = uninitialized_fill_n(_M_finish, __n - size(), __val);
  }
  else
    erase(fill_n(begin(), __n, __val), end());
}

template <class _Tp, class _Alloc>
void 
vector<_Tp, _Alloc>::_M_insert_overflow(_Tp* __position, const _Tp& __x, size_type __fill_len)
{
  const size_type __old_size = size();
  const size_type __len = __old_size + max(__old_size, __fill_len);

  pointer __new_start = _M_end_of_storage.allocate(__len);
  pointer __new_finish = __new_start;
  __STL_TRY {
    __new_finish = uninitialized_copy(_M_start, __position, __new_start);
    // handle insertion
    if (__fill_len==1) {
      __STLPORT_STD::construct(__new_finish, __x);
      ++__new_finish;
    }
    else
      __new_finish = uninitialized_fill_n(__new_finish, __fill_len, __x);
    // copy remainder
    __new_finish = uninitialized_copy(__position, _M_finish, __new_finish);
  }
  __STL_UNWIND((__STLPORT_STD::destroy(__new_start,__new_finish), 
		_M_end_of_storage.deallocate(__new_start,__len)));
  __STLPORT_STD::destroy(_M_start, _M_finish);
  _M_end_of_storage.deallocate(_M_start, _M_end_of_storage._M_data - _M_start);
  _M_start = __new_start;
  _M_finish = __new_finish;
  _M_end_of_storage._M_data = __new_start + __len;
}

template <class _Tp, class _Alloc>
void 
vector<_Tp, _Alloc>::_M_fill_insert(
# if defined ( __STL_DEBUG )
				    __iterator__ __pos, 
# else
				    __iterator__ __position, 
# endif
				    __size_type__ __n, const _Tp& __x) {
# if defined ( __STL_DEBUG )
  __stl_debug_check(__check_if_owner(&_M_iter_list, __pos));
  pointer __position=_Make_ptr(__pos);
# endif
  if (__n != 0) {
    if (size_type(_M_end_of_storage._M_data - _M_finish) >= __n) {
      _Tp __x_copy = __x;
      const size_type __elems_after = _M_finish - __position;
      pointer __old_finish = _M_finish;
      if (__elems_after > __n) {
	uninitialized_copy(_M_finish - __n, _M_finish, _M_finish);
	_M_finish += __n;
	copy_backward(__position, __old_finish - __n, __old_finish);
	fill(__position, __position + __n, __x_copy);
      }
      else {
	uninitialized_fill_n(_M_finish, __n - __elems_after, __x_copy);
	_M_finish += __n - __elems_after;
	uninitialized_copy(__position, __old_finish, _M_finish);
	_M_finish += __elems_after;
	fill(__position, __old_finish, __x_copy);
      }
      __stl_debug_do(__invalidate_range(&_M_iter_list, __pos, end()));
    }
    else 
      _M_insert_overflow(__position, __x, __n);
  }
}


__STL_END_NAMESPACE

# undef vector

# undef __pointer__
# undef __const_pointer__
# undef __size_type__
# undef __difference_type__
# undef _Make_iterator
# undef _Make_const_iterator
# undef _Make_ptr
# undef __iterator__
# undef __const_iterator__

#endif /*  __STL_VECTOR_C */

      // Local Variables:
	// mode:C++
	// End:
