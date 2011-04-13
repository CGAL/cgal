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
#ifndef __STL_STRING_C
#define __STL_STRING_C


#ifndef __SGI_STDEXCEPT
#  include <stdexcept>      
#endif

# if !defined (__STL_LINK_TIME_INSTANTIATION)
#  include <stl_string_fwd.c>
# endif

# if defined (__STL_USE_NEW_IOSTREAMS) && ! defined (__STLPORT_NEW_IOSTREAMS) && \
       !defined (__STL_MSVC) && !defined (__STL_USE_MSIPL) && !defined (__BORLANDC__)
#  include <locale>
# endif

# if defined (__STL_USE_OWN_NAMESPACE) || !defined (__STL_USE_NATIVE_STRING)

#if defined(__sgi) && !defined(__GNUC__) && (_MIPS_SIM != _MIPS_SIM_ABI32)
#pragma set woff 1174
#pragma set woff 1375
#endif

# undef _Make_ptr
# if defined (__STL_DEBUG)
#  define _Make_ptr(__i)   __i._M_iterator
# else
#  define _Make_ptr(__i)   __i
# endif

# if defined (__STL_NESTED_TYPE_PARAM_BUG)
#  define __size_type__ size_t
#  define size_type size_t
#  if defined (__STL_DEBUG)
#   define __iterator__       _Vec_iter<_CharT, _Nonconst_traits<_CharT> >
#  else
#   define __iterator__  _CharT*
#  endif
#   define iterator      __iterator__
# else
#  define __size_type__ __STL_TYPENAME_ON_RETURN_TYPE basic_string<_CharT,_Traits,_Alloc>::size_type
#  define __iterator__  __STL_TYPENAME_ON_RETURN_TYPE basic_string<_CharT,_Traits,_Alloc>::iterator
# endif

__STL_BEGIN_NAMESPACE

#if defined (__STL_MEMBER_TEMPLATES) && ! defined (__STL_INLINE_MEMBER_TEMPLATES)
template <class _CharT, class _Traits, class _Alloc>  
template <class _ForwardIter>
basic_string<_CharT, _Traits, _Alloc>& 
basic_string<_CharT, _Traits, _Alloc>::append(_ForwardIter __first, _ForwardIter __last, 
					      forward_iterator_tag)
{
  __stl_debug_do(__check_range(__first, __last));
  if (__first != __last) {
    const size_type __old_size = size();
    difference_type __n = 0;
    distance(__first, __last, __n);
    if (__STATIC_CAST(size_type,__n) > max_size() || 
	__old_size > max_size() - __STATIC_CAST(size_type,__n))
      _M_throw_length_error();
    if (__old_size + __n > capacity()) {
      const size_type __len = __old_size +
	max(__old_size, __STATIC_CAST(size_type,__n)) + 1;
      pointer __new_start = _M_end_of_storage.allocate(__len);
      pointer __new_finish = __new_start;
      __STL_TRY {
	__new_finish = uninitialized_copy(_M_start, _M_finish, __new_start);
	__new_finish = uninitialized_copy(__first, __last, __new_finish);
	_M_construct_null(__new_finish);
      }
      __STL_UNWIND((destroy(__new_start,__new_finish),
		    _M_end_of_storage.deallocate(__new_start,__len)));
      destroy(_M_start, _M_finish + 1);
      _M_deallocate_block();
      _M_start = __new_start;
      _M_finish = __new_finish;
      _M_end_of_storage._M_data = __new_start + __len; 
    }
    else {
      _ForwardIter __f1 = __first;
      ++__f1;
      uninitialized_copy(__f1, __last, _M_finish + 1);
      __STL_TRY {
	_M_construct_null(_M_finish + __n);
      }
      __STL_UNWIND(destroy(_M_finish + 1, _M_finish + __n));
      _Traits::assign(*_M_finish, *__first);
      _M_finish += __n;
    }
  }
  return *this;  
}

template <class _CharT, class _Traits, class _Alloc>  template <class _ForwardIter>
void 
basic_string<_CharT, _Traits, _Alloc>::insert(iterator __position, 
					      _ForwardIter __first, _ForwardIter __last, 
					      forward_iterator_tag)  {
  __stl_debug_do(__check_range(__first,__last));
  if (__first != __last) {
    difference_type __n = 0;
    distance(__first, __last, __n);
    if (_M_end_of_storage._M_data - _M_finish >= __n + 1) {
      const difference_type __elems_after = _M_finish - _Make_ptr(__position);
      pointer __old_finish = _M_finish;
      if (__elems_after >= __n) {
	uninitialized_copy((_M_finish - __n) + 1, _M_finish + 1,
			   _M_finish + 1);
	_M_finish += __n;
	_Traits::move(_Make_ptr(__position) + __n,
		      _Make_ptr(__position), (__elems_after - __n) + 1);
	_M_copy(__first, __last, _Make_ptr(__position));
      }
      else {
	_ForwardIter __mid = __first;
	advance(__mid, __elems_after + 1);
	uninitialized_copy(__mid, __last, _M_finish + 1);
	_M_finish += __n - __elems_after;
	__STL_TRY {
	  uninitialized_copy(_Make_ptr(__position), __old_finish + 1, _M_finish);
	  _M_finish += __elems_after;
	}
	__STL_UNWIND((destroy(__old_finish + 1, _M_finish), 
		      _M_finish = __old_finish));
	_M_copy(__first, __mid, _Make_ptr(__position));
      }
    }
    else {
      const size_type __old_size = size();        
      const size_type __len
	= __old_size + max(__old_size, __STATIC_CAST(size_type,__n)) + 1;
      pointer __new_start = _M_end_of_storage.allocate(__len);
      pointer __new_finish = __new_start;
      __STL_TRY {
	__new_finish = uninitialized_copy(_M_start, _Make_ptr(__position), __new_start);
	__new_finish = uninitialized_copy(__first, __last, __new_finish);
	__new_finish
	  = uninitialized_copy(_Make_ptr(__position), _M_finish, __new_finish);
	_M_construct_null(__new_finish);
      }
      __STL_UNWIND((destroy(__new_start,__new_finish),
		    _M_end_of_storage.deallocate(__new_start,__len)));
      destroy(_M_start, _M_finish + 1);
      _M_deallocate_block();
      _M_start = __new_start;
      _M_finish = __new_finish;
      _M_end_of_storage._M_data = __new_start + __len; 
    }
  }
}

# endif /* __STL_INLINE_MEMBER_TEMPLATES */


// ------------------------------------------------------------
// Non-inline declarations.

# if (__STL_STATIC_TEMPLATE_DATA > 0)
template <class _CharT, class _Traits, class _Alloc> 
const __size_type__
basic_string<_CharT,_Traits,_Alloc>::npos = (size_t)-1;
# else /* __STL_STATIC_TEMPLATE_DATA */
template class basic_string<char,char_traits<char>,allocator<char> >;
const size_t string::npos __STL_WEAK = (size_t)-1;
#  ifdef __STL_HAS_WCHAR_T
template class basic_string<wchar_t,char_traits<wchar_t>,allocator<wchar_t> >;
const size_t wstring::npos __STL_WEAK = (size_t)-1;
#  endif
# endif /* __STL_STATIC_TEMPLATE_DATA */

// _String_base methods
template <class _Tp, class _Alloc> 
void _String_base<_Tp,_Alloc>::_M_throw_length_error() const {
  __STL_THROW(length_error(string("basic_string")));
}

template <class _Tp, class _Alloc> 
void _String_base<_Tp, _Alloc>::_M_throw_out_of_range() const {
  __STL_THROW(out_of_range(string("basic_string")));
}

// Change the string's capacity so that it is large enough to hold
//  at least __res_arg elements, plus the terminating _CharT().  Note that,
//  if __res_arg < capacity(), this member function may actually decrease
//  the string's capacity.
template <class _CharT, class _Traits, class _Alloc> 
void basic_string<_CharT,_Traits,_Alloc>::reserve(__size_type__ __res_arg) {
  if (__res_arg > max_size())
    _M_throw_length_error();

  size_type __n = max(__res_arg, size()) + 1;
  pointer __new_start = _M_end_of_storage.allocate(__n);
  pointer __new_finish = __new_start;

  __STL_TRY {
    __new_finish = uninitialized_copy(_M_start, _M_finish, __new_start);
    _M_construct_null(__new_finish);
  }
  __STL_UNWIND((destroy(__new_start, __new_finish), 
                _M_end_of_storage.deallocate(__new_start, __n)));

  destroy(_M_start, _M_finish + 1);
  _M_deallocate_block();
  _M_start = __new_start;
  _M_finish = __new_finish;
  _M_end_of_storage._M_data = __new_start + __n;
}

template <class _CharT, class _Traits, class _Alloc> 
basic_string<_CharT,_Traits,_Alloc>& 
basic_string<_CharT,_Traits,_Alloc>::append(__size_type__ __n, _CharT __c) {
  if (__n > max_size() || size() > max_size() - __n)
    _M_throw_length_error();
  if (size() + __n > capacity())
    reserve(size() + max(size(), __n));
  if (__n > 0) {
    uninitialized_fill_n(_M_finish + 1, __n - 1, __c);
    __STL_TRY {
      _M_construct_null(_M_finish + __n);
    }
    __STL_UNWIND(destroy(_M_finish + 1, _M_finish + __n));
    _Traits::assign(*_M_finish, __c);
    _M_finish += __n;
  }
  return *this;
}

#ifndef __STL_MEMBER_TEMPLATES

template <class _CharT, class _Traits, class _Alloc> 
basic_string<_CharT, _Traits, _Alloc>& 
basic_string<_CharT, _Traits, _Alloc>::append(const _CharT* __first,
					      const _CharT* __last)
{
  if (__first != __last) {
    const size_type __old_size = size();
    ptrdiff_t __n = __last - __first;
    if ((size_type)__n > max_size() || __old_size > max_size() - __n)
      _M_throw_length_error();
    if (__old_size + __n > capacity()) {
      const size_type __len = __old_size + max(__old_size, (size_t) __n) + 1;
      pointer __new_start = _M_end_of_storage.allocate(__len);
      pointer __new_finish = __new_start;
      __STL_TRY {
        __new_finish = uninitialized_copy(_M_start, _M_finish, __new_start);
        __new_finish = uninitialized_copy(__first, __last, __new_finish);
        _M_construct_null(__new_finish);
      }
      __STL_UNWIND((destroy(__new_start,__new_finish),
                    _M_end_of_storage.deallocate(__new_start,__len)));
      destroy(_M_start, _M_finish + 1);
      _M_deallocate_block();
      _M_start = __new_start;
      _M_finish = __new_finish;
      _M_end_of_storage._M_data = __new_start + __len; 
    }
    else {
      const _CharT* __f1 = __first;
      ++__f1;
      uninitialized_copy(__f1, __last, _M_finish + 1);
      __STL_TRY {
        _M_construct_null(_M_finish + __n);
      }
      __STL_UNWIND(destroy(_M_finish + 1, _M_finish + __n));
      _Traits::assign(*_M_finish, *__first);
      _M_finish += __n;
    }
  }
  return *this;  
}

#endif /* __STL_MEMBER_TEMPLATES */

template <class _CharT, class _Traits, class _Alloc> 
basic_string<_CharT,_Traits,_Alloc>& 
basic_string<_CharT,_Traits,_Alloc>::assign(__size_type__ __n, _CharT __c) {
  if (__n <= size()) {
    _Traits::assign(_M_start, __n, __c);
    erase(begin() + __n, end());
  }
  else {
    _Traits::assign(_M_start, size(), __c);
    append(__n - size(), __c);
  }
  return *this;
}

template <class _CharT, class _Traits, class _Alloc> 
basic_string<_CharT,_Traits,_Alloc>& 
basic_string<_CharT,_Traits,_Alloc>::assign(const _CharT* __f, 
                                            const _CharT* __l)
{
  __stl_debug_do(__check_range(__f, __l));
  ptrdiff_t __n = __l - __f;
  if (__STATIC_CAST(size_type,__n) <= size()) {
    _Traits::copy(_M_start, __f, __n);
    erase(begin() + __n, end());
  }
  else {
    _Traits::copy(_M_start, __f, size());
    append(__f + size(), __l);
  }
  return *this;
}

template <class _CharT, class _Traits, class _Alloc>
_CharT* 
basic_string<_CharT,_Traits,_Alloc>
  ::_M_insert_aux(_CharT* __p,
                  _CharT __c)
{
  pointer __new_pos = __p;
  if (_M_finish + 1 < _M_end_of_storage._M_data) {
    _M_construct_null(_M_finish + 1);
    _Traits::move(__p + 1, __p, _M_finish - __p);
    _Traits::assign(*__p, __c);
    ++_M_finish;
  }
  else {
    const size_type __old_len = size();
    const size_type __len = __old_len +
                            max(__old_len, __STATIC_CAST(size_type,1)) + 1;
    pointer __new_start = _M_end_of_storage.allocate(__len);
    pointer __new_finish = __new_start;
    __STL_TRY {
      __new_pos = uninitialized_copy(_M_start, __p, __new_start);
      construct(__new_pos, __c);
      __new_finish = __new_pos + 1;
      __new_finish = uninitialized_copy(__p, _M_finish, __new_finish);
      _M_construct_null(__new_finish);
    }
    __STL_UNWIND((destroy(__new_start,__new_finish), 
                  _M_end_of_storage.deallocate(__new_start,__len)));
    destroy(_M_start, _M_finish + 1);
    _M_deallocate_block();
    _M_start = __new_start;
    _M_finish = __new_finish;
    _M_end_of_storage._M_data = __new_start + __len;
  }
  return __new_pos;
}

template <class _CharT, class _Traits, class _Alloc>
void basic_string<_CharT,_Traits,_Alloc>::insert(__iterator__ __position,
           size_t __n, _CharT __c)
{
  __stl_debug_do(__check_if_owner(&_M_iter_list,__position));
  if (__n != 0) {
    if (size_type(_M_end_of_storage._M_data - _M_finish) >= __n + 1) {
      const size_type __elems_after = _M_finish - _Make_ptr(__position);
      pointer __old_finish = _M_finish;
      if (__elems_after >= __n) {
        uninitialized_copy((_M_finish - __n) + 1, _M_finish + 1,
                           _M_finish + 1);
        _M_finish += __n;
        _Traits::move(_Make_ptr(__position) + __n,
                      _Make_ptr(__position), (__elems_after - __n) + 1);
        _Traits::assign(_Make_ptr(__position), __n, __c);
      }
      else {
        uninitialized_fill_n(_M_finish + 1, __n - __elems_after - 1, __c);
        _M_finish += __n - __elems_after;
        __STL_TRY {
          uninitialized_copy(_Make_ptr(__position), __old_finish + 1, _M_finish);
          _M_finish += __elems_after;
        }
        __STL_UNWIND((destroy(__old_finish + 1, _M_finish), 
                      _M_finish = __old_finish));
        _Traits::assign(_Make_ptr(__position), __elems_after + 1, __c);
      }
    }
    else {
      const size_type __old_size = size();        
      const size_type __len = __old_size + max(__old_size, __n) + 1;
      pointer __new_start = _M_end_of_storage.allocate(__len);
      pointer __new_finish = __new_start;
      __STL_TRY {
        __new_finish = uninitialized_copy(_M_start, _Make_ptr(__position), __new_start);
        __new_finish = uninitialized_fill_n(__new_finish, __n, __c);
        __new_finish = uninitialized_copy(_Make_ptr(__position), _M_finish,
                                          __new_finish);
        _M_construct_null(__new_finish);
      }
      __STL_UNWIND((destroy(__new_start,__new_finish),
                    _M_end_of_storage.deallocate(__new_start,__len)));
      destroy(_M_start, _M_finish + 1);
      _M_deallocate_block();
      _M_start = __new_start;
      _M_finish = __new_finish;
      _M_end_of_storage._M_data = __new_start + __len;    
    }
  }
}

#ifndef __STL_MEMBER_TEMPLATES

template <class _CharT, class _Traits, class _Alloc>
void 
basic_string<_CharT,_Traits,_Alloc>::insert(__iterator__ __position,
                                            const _CharT* __first, 
                                            const _CharT* __last)
{
  __stl_debug_do(__check_if_owner(&_M_iter_list,__position) && 
		 __check_range(__first, __last));
  if (__first != __last) {
    const ptrdiff_t __n = __last - __first;
    if (_M_end_of_storage._M_data - _M_finish >= __n + 1) {
      const ptrdiff_t __elems_after = _M_finish - _Make_ptr(__position);
      pointer __old_finish = _M_finish;
      if (__elems_after >= __n) {
        uninitialized_copy((_M_finish - __n) + 1, _M_finish + 1,
                           _M_finish + 1);
        _M_finish += __n;
        _Traits::move(_Make_ptr(__position) + __n,
                      _Make_ptr(__position), (__elems_after - __n) + 1);
        _M_copy(__first, __last, _Make_ptr(__position));
      }
      else {
        const _CharT* __mid = __first;
        advance(__mid, __elems_after + 1);
        uninitialized_copy(__mid, __last, _M_finish + 1);
        _M_finish += __n - __elems_after;
        __STL_TRY {
          uninitialized_copy(_Make_ptr(__position), __old_finish + 1, _M_finish);
          _M_finish += __elems_after;
        }
        __STL_UNWIND((destroy(__old_finish + 1, _M_finish), 
                      _M_finish = __old_finish));
        _M_copy(__first, __mid, _Make_ptr(__position));
      }
    }
    else {
      size_type __old_size = size();        
      size_type __len
        = __old_size + max(__old_size, __STATIC_CAST(size_type,__n)) + 1;
      pointer __new_start = _M_end_of_storage.allocate(__len);
      pointer __new_finish = __new_start;
      __STL_TRY {
        __new_finish = uninitialized_copy(_M_start, _Make_ptr(__position), __new_start);
        __new_finish = uninitialized_copy(__first, __last, __new_finish);
        __new_finish
          = uninitialized_copy(_Make_ptr(__position), _M_finish, __new_finish);
        _M_construct_null(__new_finish);
      }
      __STL_UNWIND((destroy(__new_start,__new_finish),
                    _M_end_of_storage.deallocate(__new_start,__len)));
      destroy(_M_start, _M_finish + 1);
      _M_deallocate_block();
      _M_start = __new_start;
      _M_finish = __new_finish;
      _M_end_of_storage._M_data = __new_start + __len; 
    }
  }
}

#endif /* __STL_MEMBER_TEMPLATES */

template <class _CharT, class _Traits, class _Alloc>
basic_string<_CharT,_Traits,_Alloc>&
basic_string<_CharT,_Traits,_Alloc>
  ::replace(iterator __first, iterator __last, __size_type__ __n, _CharT __c)
{
  __stl_debug_do(__check_if_owner(&_M_iter_list,__first) 
		 && __check_range(__first, __last));
  const size_type __len = __STATIC_CAST(size_type,(__last - __first));
  if (__len >= __n) {
    _Traits::assign(_Make_ptr(__first), __n, __c);
    erase(__first + __n, __last);
  }
  else {
    _Traits::assign(_Make_ptr(__first), __len, __c);
    insert(__last, __n - __len, __c);
  }
  return *this;
}

#ifndef __STL_MEMBER_TEMPLATES


template <class _CharT, class _Traits, class _Alloc>
basic_string<_CharT,_Traits,_Alloc>&
basic_string<_CharT,_Traits,_Alloc>
  ::replace(iterator __first, iterator __last,
            const _CharT* __f, const _CharT* __l)
{
  __stl_debug_do(__check_if_owner(&_M_iter_list,__first) && 
		 __check_range(__first, __last) &&
		 __check_range(__f, __l));
  const ptrdiff_t         __n = __l - __f;
  const difference_type __len = __last - __first;
  if (__len >= __n) {
    _M_copy(__f, __l, _Make_ptr(__first));
    erase(__first + __n, __last);
  }
  else {
    const _CharT* __m = __f + __len;
    _M_copy(__f, __m, _Make_ptr(__first));
    insert(__last, __m, __l);
  }
  return *this;
}

#endif /* __STL_MEMBER_TEMPLATES */

template <class _CharT, class _Traits, class _Alloc>
__size_type__
basic_string<_CharT,_Traits,_Alloc>
  ::find(const _CharT* __s, size_type __pos, size_type __n) const 
{
  if (__pos + __n >= size())
    return npos;
  else {
    const const_pointer __result =
      search((const _CharT*)_M_start + __pos, (const _CharT*)_M_finish, 
             __s, __s + __n, _Eq_traits<_Traits>());
    return __result != _M_finish ? __result - _M_start : npos;
  }
}

template <class _CharT, class _Traits, class _Alloc>
__size_type__
basic_string<_CharT,_Traits,_Alloc>
  ::find(_CharT __c, size_type __pos) const 
{
  if (__pos >= size())
    return npos;
  else {
    const const_pointer __result =
      find_if((const _CharT*)_M_start + __pos, (const _CharT*)_M_finish,
              bind2nd(_Eq_traits<_Traits>(), __c));
    return __result != _M_finish ? __result - _M_start : npos;
  }
}    

template <class _CharT, class _Traits, class _Alloc>
__size_type__
basic_string<_CharT,_Traits,_Alloc>
  ::rfind(const _CharT* __s, size_type __pos, size_type __n) const 
{
  const size_t __len = size();

  if (__n > __len)
    return npos;
  else if (__n == 0)
    return min(__len, __pos);
  else {
    const_pointer __last = _M_start + min(__len - __n, __pos) + __n;
    const_pointer __result = find_end((const_pointer)_M_start, __last,
				      __s, __s + __n,
				      _Eq_traits<_Traits>());
    return __result != __last ? __result - _M_start : npos;
  }
}

template <class _CharT, class _Traits, class _Alloc>
__size_type__
basic_string<_CharT,_Traits,_Alloc>
  ::rfind(_CharT __c, size_type __pos) const 
{
  const size_type __len = size();

  if (__len < 1)
    return npos;
  else {
    const const_iterator __last = begin() + min(__len - 1, __pos) + 1;
    const_reverse_iterator __rresult =
      find_if(const_reverse_iterator(__last), rend(),
              bind2nd(_Eq_traits<_Traits>(), __c));
    return __rresult != rend() ? (__rresult.base() - 1) - begin() : npos;
  }
}

template <class _CharT, class _Traits, class _Alloc>
__size_type__
basic_string<_CharT,_Traits,_Alloc>
  ::find_first_of(const _CharT* __s, size_type __pos, size_type __n) const
{
  if (__pos >= size())
    return npos;
  else {
    const_iterator __result = __STLPORT_STD::find_first_of(begin() + __pos, end(),
                                                   __s, __s + __n,
                                                   _Eq_traits<_Traits>());
    return __result != end() ? __result - begin() : npos;
  }
}


template <class _CharT, class _Traits, class _Alloc>
__size_type__
basic_string<_CharT,_Traits,_Alloc>
  ::find_last_of(const _CharT* __s, size_type __pos, size_type __n) const
{
  const size_type __len = size();

  if (__len < 1)
    return npos;
  else {
    const const_iterator __last = begin() + min(__len - 1, __pos) + 1;
    const const_reverse_iterator __rresult =
      __STLPORT_STD::find_first_of(const_reverse_iterator(__last), rend(),
                           __s, __s + __n,
                           _Eq_traits<_Traits>());
    return __rresult != rend() ? (__rresult.base() - 1) - begin() : npos;
  }
}


template <class _CharT, class _Traits, class _Alloc>
__size_type__
basic_string<_CharT,_Traits,_Alloc>
  ::find_first_not_of(const _CharT* __s, size_type __pos, size_type __n) const
{
  typedef typename _Traits::char_type _CharType;
  if (__pos > size())
    return npos;
  else {
    const_pointer __result = find_if((const _CharT*)_M_start + __pos, 
				      (const _CharT*)_M_finish,
                                _Not_within_traits<_Traits>((const _CharType*)__s, 
							    (const _CharType*)__s + __n));
    return __result != _M_finish ? __result - _M_start : npos;
  }
}

template <class _CharT, class _Traits, class _Alloc>
__size_type__
basic_string<_CharT,_Traits,_Alloc>
  ::find_first_not_of(_CharT __c, size_type __pos) const
{
  if (__pos > size())
    return npos;
  else {
    const_pointer __result = find_if((const _CharT*)_M_start + __pos, (const _CharT*)_M_finish,
				     not1(bind2nd(_Eq_traits<_Traits>(), __c)));
    return __result != _M_finish ? __result - _M_start : npos;
  }
}    

template <class _CharT, class _Traits, class _Alloc>
__size_type__
basic_string<_CharT,_Traits,_Alloc>
  ::find_last_not_of(const _CharT* __s, size_type __pos, size_type __n) const 
{
  typedef typename _Traits::char_type _CharType;
  const size_type __len = size();

  if (__len < 1)
    return npos;
  else {
    const_iterator __last = begin() + min(__len - 1, __pos) + 1;
    const_reverse_iterator __rlast = const_reverse_iterator(__last);
    const_reverse_iterator __rresult =
      find_if(__rlast, rend(),
              _Not_within_traits<_Traits>((const _CharType*)__s, 
					  (const _CharType*)__s + __n));
    return __rresult != rend() ? (__rresult.base() - 1) - begin() : npos;
  }
}

template <class _CharT, class _Traits, class _Alloc>
__size_type__
basic_string<_CharT, _Traits, _Alloc>
  ::find_last_not_of(_CharT __c, size_type __pos) const 
{
  const size_type __len = size();

  if (__len < 1)
    return npos;
  else {
    const_iterator __last = begin() + min(__len - 1, __pos) + 1;
    const_reverse_iterator __rlast = const_reverse_iterator(__last);
    const_reverse_iterator __rresult =
      find_if(__rlast, rend(),
              not1(bind2nd(_Eq_traits<_Traits>(), __c)));
    return __rresult != rend() ? (__rresult.base() - 1) - begin() : npos;
  }
}


#if defined (__STL_USE_NEW_IOSTREAMS)

template <class _CharT, class _Traits>
inline bool
__sgi_string_fill(basic_ostream<_CharT, _Traits>& __os,
                  basic_streambuf<_CharT, _Traits>* __buf,
                  size_t __n)
{
  _CharT __f = __os.fill();
  size_t __i;
  bool __ok = true;

  for (__i = 0; __i < __n; ++__i)
    __ok == __ok && !_Traits::eq_int_type(__buf->sputc(__f), _Traits::eof());
  return __ok;
}

# if !(defined(__STLPORT_NEW_IOSTREAMS) || (defined (__BORLANDC__) && ! defined (__STL_USE_OWN_NAMESPACE)))

template <class _CharT, class _Traits, class _Alloc>
basic_ostream<_CharT, _Traits>&
operator<<(basic_ostream<_CharT, _Traits>& __os, 
           const basic_string<_CharT,_Traits,_Alloc>& __s)
{
  typename __STL_VENDOR_STD::basic_ostream<_CharT, _Traits>::sentry __sentry(__os);
  bool __ok = false;

  if (__sentry) {
    __ok = true;
    size_t __n = __s.size();
    size_t __pad_len = 0;
    const bool __left = (__os.flags() & ios::left) != 0;
    const size_t __w = __os.width(0);
    basic_streambuf<_CharT, _Traits>* __buf = __os.rdbuf();

    if (__w > 0) {
      __n = min(__w, __n);
      __pad_len = __w - __n;
    }
    
    if (!__left)
      __ok = __sgi_string_fill(__os, __buf, __pad_len);    

    __ok = __ok && ((size_t)__buf->sputn(__s.data(), __n) == __n);

    if (__left)
      __ok = __ok && __sgi_string_fill(__os, __buf, __pad_len);
  }

  if (!__ok)
    __os.setstate(ios_base::failbit);

  return __os;
}

template <class _CharT, class _Traits, class _Alloc>
basic_istream<_CharT, _Traits>& 
operator>>(basic_istream<_CharT, _Traits>& __is,
           basic_string<_CharT,_Traits, _Alloc>& __s)
{
# ifndef __STL_HAS_NO_NAMESPACES
using namespace __STL_VENDOR_STD;
# endif
  typename basic_istream<_CharT, _Traits>::sentry __sentry(__is);

  if (__sentry) {
    basic_streambuf<_CharT, _Traits>* __buf = __is.rdbuf();
    typedef ctype<_CharT> _C_type;
    const locale& __loc = __is.getloc();
#ifdef __STLPORT_NEW_IOSTREAMS
    const _C_type& _Ctype = use_facet<_C_type>(__loc);
#else
# if defined (__STL_MSVC) && (__STL_MSVC <= 1200 )
    const _C_type& _Ctype = use_facet(__loc , ( _C_type * ) 0, true);
# elif defined (__SUNPRO_CC)
    const _C_type& _Ctype = use_facet(__loc , ( _C_type * ) 0);
# else
    const _C_type& _Ctype = use_facet<_C_type>(__loc);
# endif
#endif
    __s.clear();
    size_t __n = __is.width(0);
    if (__n == 0)
      __n = __STATIC_CAST(size_t,-1);
    else
      __s.reserve(__n);
    

    while (__n-- > 0) {
      typename _Traits::int_type __c1 = __buf->sbumpc();
      if (_Traits::eq_int_type(__c1, _Traits::eof())) {
        __is.setstate(ios_base::eofbit);
        break;
      }
      else {
        _CharT __c = _Traits::to_char_type(__c1);

        if (_Ctype.is(_C_type::space, __c)) {
          if (_Traits::eq_int_type(__buf->sputbackc(__c), _Traits::eof()))
            __is.setstate(ios_base::failbit);
          break;
        }
        else
          __s.push_back(__c);
      }
    }
    
    // If we have read no characters, then set failbit.
    if (__s.size() == 0)
      __is.setstate(ios_base::failbit);
  }
  else
    __is.setstate(ios_base::failbit);

  return __is;
}

template <class _CharT, class _Traits, class _Alloc>    
basic_istream<_CharT, _Traits>& 
getline(basic_istream<_CharT, _Traits>& __is,
        basic_string<_CharT,_Traits,_Alloc>& __s,
        _CharT __delim)
{
  size_t __nread = 0;
  typename __STL_VENDOR_STD::basic_istream<_CharT, _Traits>::sentry __sentry(__is, true);
  if (__sentry) {
    basic_streambuf<_CharT, _Traits>* __buf = __is.rdbuf();
    __s.clear();

    while (__nread < __s.max_size()) {
      int __c1 = __buf->sbumpc();
      if (_Traits::eq_int_type(__c1, _Traits::eof())) {
        __is.setstate(ios_base::eofbit);
        break;
      }
      else {
        ++__nread;
        _CharT __c = _Traits::to_char_type(__c1);
        if (!_Traits::eq(__c, __delim)) 
          __s.push_back(__c);
        else
          break;              // Character is extracted but not appended.
      }
    }
  }
  if (__nread == 0 || __nread >= __s.max_size())
    __is.setstate(ios_base::failbit);

  return __is;
}
# endif /* __BORLANDC */

#elif ! defined ( __STL_USE_NO_IOSTREAMS )

inline void 
__sgi_string_fill(ostream& __os, streambuf* __buf, size_t __n)
{
  char __f = __os.fill();
  size_t __i;

  for (__i = 0; __i < __n; ++__i) __buf->sputc(__f);
}

template <class _CharT, class _Traits, class _Alloc>
ostream& operator<<(ostream& __os, 
                    const basic_string<_CharT,_Traits,_Alloc>& __s)
{
  streambuf* __buf = __os.rdbuf();
  if (__buf) {
    size_t __n = __s.size();
    size_t __pad_len = 0;
    const bool __left = (__os.flags() & ios::left) !=0;
    const size_t __w = __os.width();

    if (__w > 0) {
      __n = min(__w, __n);
      __pad_len = __w - __n;
    }
    
    if (!__left)
      __sgi_string_fill(__os, __buf, __pad_len);
  
    const size_t __nwritten = __buf->sputn(__s.data(), __n);

    if (__left)
      __sgi_string_fill(__os, __buf, __pad_len);

    if (__nwritten != __n)
      __os.clear(__os.rdstate() | ios::failbit);

    __os.width(0);
  }
  else
    __os.clear(__os.rdstate() | ios::badbit);

  return __os;
}

template <class _CharT, class _Traits, class _Alloc>
istream& operator>>(istream& __is, basic_string<_CharT,_Traits,_Alloc>& __s)
{
  if (!__is)
    return __is;

  streambuf* __buf = __is.rdbuf();
  if (__buf) {
#ifdef __USLC__
/* Jochen Schlick '1999  - operator >> modified. Work-around to get the 
 *                         output buffer flushed (necessary when using 
 *                         "cout" (without endl or flushing) followed by
 *                         "cin >>" ...)
 */
    if (__is.flags() & ios::skipws) {
      _CharT __c;
      do 
         __is.get(__c);
      while (__is && isspace(__c));
      if (__is)
         __is.putback(__c);
    }
#else

    if (__is.flags() & ios::skipws) {
      _CharT __c;
      do {
        int __c1 = __buf->sbumpc();
        if (__c1 == EOF) {
          __is.clear(__is.rdstate() | ios::eofbit | ios::failbit);
          break;
        }
        else
          __c = _Traits::to_char_type(__c1);
      }
      while (isspace((unsigned char) __c));

      if (__buf->sputbackc(__c) == EOF)
        __is.clear(__is.rdstate() | ios::failbit);
    }
# endif
    // If we arrive at end of file (or fail for some other reason) while
    // still discarding whitespace, then we don't try to read the string.
    if (__is) {
      __s.clear();

      size_t __n = __is.width();
      if (__n == 0)
        __n = __STATIC_CAST(size_t,-1);
      else
        __s.reserve(__n);

      while (__n-- > 0) {
        int __c1 = __buf->sbumpc();
        if (__c1 == EOF) {
          __is.clear(__is.rdstate() | ios::eofbit);
          break;
        }
        else {
          _CharT __c = _Traits::to_char_type(__c1);

          if (isspace((unsigned char) __c)) {
            if (__buf->sputbackc(__c) == EOF)
              __is.clear(__is.rdstate() | ios::failbit);
            break;
          }
          else
            __s.push_back(__c);
        }
      }
    
      // If we have read no characters, then set failbit.
      if (__s.size() == 0)
        __is.clear(__is.rdstate() | ios::failbit);
    }

    __is.width(0);
  }
  else                          // We have no streambuf.
    __is.clear(__is.rdstate() | ios::badbit);

  return __is;
}

template <class _CharT, class _Traits, class _Alloc>    
istream& getline(istream& __is,
                 basic_string<_CharT,_Traits,_Alloc>& __s,
                 _CharT __delim)
{
  streambuf* __buf = __is.rdbuf();
  if (__buf) {
    size_t __nread = 0;
    if (__is) {
      __s.clear();

      while (__nread < __s.max_size()) {
        int __c1 = __buf->sbumpc();
        if (__c1 == EOF) {
          __is.clear(__is.rdstate() | ios::eofbit);
          break;
        }
        else {
          ++__nread;
          _CharT __c = _Traits::to_char_type(__c1);
          if (!_Traits::eq(__c, __delim)) 
            __s.push_back(__c);
          else
            break;              // Character is extracted but not appended.
        }
      }
    }

    if (__nread == 0 || __nread >= __s.max_size())
      __is.clear(__is.rdstate() | ios::failbit);
  }
  else
    __is.clear(__is.rdstate() | ios::badbit);

  return __is;
}

# endif /* __STL_NEW_IOSTREAMS */

template <class _CharT, class _Traits, class _Alloc>
void _S_string_copy(const basic_string<_CharT,_Traits,_Alloc>& __s,
                    _CharT* __buf,
                    size_t __n)
{
  if (__n > 0) {
    __n = min(__n - 1, __s.size());
    copy(__s.begin(), __s.begin() + __n, __buf);
    __buf[__n] = _CharT();
  }
}

__STL_END_NAMESPACE

# undef __size_type__
# undef size_type
# undef _Make_ptr
# undef __iterator__
# undef iterator


#if defined(__sgi) && !defined(__GNUC__) && (_MIPS_SIM != _MIPS_SIM_ABI32)
#pragma reset woff 1174
#pragma reset woff 1375
#endif

# endif /* NATIVE */

#endif /*  __STL_STRING_C */

// Local Variables:
// mode:C++
// End:
