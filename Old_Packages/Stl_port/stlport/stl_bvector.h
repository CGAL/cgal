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

#ifndef __SGI_STL_INTERNAL_BVECTOR_H
#define __SGI_STL_INTERNAL_BVECTOR_H

#ifndef __SGI_STL_INTERNAL_VECTOR_H
#  include <stl_vector.h>
#endif

__STL_BEGIN_NAMESPACE 

#define __WORD_BIT (int(CHAR_BIT*sizeof(unsigned int)))

#if defined(__sgi) && !defined(__GNUC__) && (_MIPS_SIM != _MIPS_SIM_ABI32)
#pragma set woff 1174
#pragma set woff 1375
#endif

struct _Bit_reference {
  unsigned int* _M_p;
  unsigned int _M_mask;
  _Bit_reference(unsigned int* __x, unsigned int __y) 
    : _M_p(__x), _M_mask(__y) {}

public:
  _Bit_reference() : _M_p(0), _M_mask(0) {}
  operator bool() const { return !(!(*_M_p & _M_mask)); }
  _Bit_reference& operator=(bool __x)
  {
    if (__x)  *_M_p |= _M_mask;
    else      *_M_p &= ~_M_mask;
    return *this;
  }
  _Bit_reference& operator=(const _Bit_reference& __x) 
    { return *this = bool(__x); }
  bool operator==(const _Bit_reference& __x) const
    { return bool(*this) == bool(__x); }
  bool operator<(const _Bit_reference& __x) const {
    return !bool(*this) && bool(__x);
  }
  void flip() { *_M_p ^= _M_mask; }
};

__STL_END_NAMESPACE

# if defined (__SGI_STL_NO_ARROW_OPERATOR) && ! defined (__STL_NO_PROXY_ARROW_OPERATOR)

__STL_TEMPLATE_NULL struct __arrow_op_dispatch<_Bit_reference, _Bit_reference*> {
  __arrow_op_dispatch(_Bit_reference) {}
  __arrow_op_dummy operator ->() const { return __arrow_op_dummy(); }
};

__STL_TEMPLATE_NULL struct __arrow_op_dispatch<bool, const bool*> {
  __arrow_op_dispatch(bool) {}
  __arrow_op_dummy operator ->() const { return __arrow_op_dummy(); }
};

# endif

__STL_BEGIN_NAMESPACE

inline void swap(_Bit_reference __x, _Bit_reference __y)
{
  bool __tmp = __x;
  __x = __y;
  __y = __tmp;
}

struct _Bit_iterator_base;

# ifdef __STL_DEBUG
template <class _Dummy> 
struct _Bv_global {
  static bool _Dereferenceable(const _Bit_iterator_base& __x);
  static bool _Nonsingular(const _Bit_iterator_base& __x );
};
typedef _Bv_global<bool> _Bv_global_inst;
# endif

struct _Bit_iterator_base
# if defined (__STL_DEBUG)
 : public __owned_link
# endif
{
  typedef ptrdiff_t difference_type;

  unsigned int* _M_p;
  unsigned int _M_offset;

  void _M_bump_up() {
    if (_M_offset++ == __WORD_BIT - 1) {
      _M_offset = 0;
      ++_M_p;
    }
    __stl_debug_check(_Bv_global_inst::_Nonsingular(*this));  
  }

  void _M_bump_down() {
    if (_M_offset-- == 0) {
      _M_offset = __WORD_BIT - 1;
      --_M_p;
    }
    __stl_debug_check(_Bv_global_inst::_Nonsingular(*this));
  }

# if defined ( __STL_DEBUG )
  bool _M_unsafe;
  bool _Overrun_ok() const { return _M_unsafe; } 
  void _Set_overrun(bool __val) { _M_unsafe = __val; } 
  _Bit_iterator_base() : __owned_link(0), _M_p(0),_M_offset(0), _M_unsafe(false) {}
  _Bit_iterator_base(const __owned_list* __root,unsigned int* __x, 
		 unsigned int __y, bool __over = false) : 
    __owned_link(__root), _M_p(__x), _M_offset(__y), _M_unsafe(__over) {}
  /*
  // these are trivial copy constructor and assignment. egcs chokes
  // on compiler-generated ones.
  _Bit_iterator_base(const _Bit_iterator_base& __it) : __owned_link(__it){
    _M_p = __it._M_p;
    _M_offset = __it._M_offset;
    _M_unsafe = __it._M_unsafe;
  }
 _Bit_iterator_base& operator=(const _Bit_iterator_base& __it) {
   __owned_link::operator=(__it);
   _M_p = __it._M_p;
   _M_offset = __it._M_offset;
   _M_unsafe = __it._M_unsafe;
   return *this;
 }
 */
# else
  _Bit_iterator_base() : _M_p(0), _M_offset(0) {}
  _Bit_iterator_base(unsigned int* __x, unsigned int __y) : _M_p(__x), _M_offset(__y) {}
# endif

  void _M_advance (difference_type __i) {
    difference_type __n = __i + _M_offset;
    _M_p += __n / __WORD_BIT;
    __n = __n % __WORD_BIT;
    if (__n < 0) {
      _M_offset = (unsigned int) __n + __WORD_BIT;
      --_M_p;
    } else
      _M_offset = (unsigned int) __n;
    __stl_debug_check(_Bv_global_inst::_Nonsingular(*this));
  }

  difference_type _M_subtract(const _Bit_iterator_base& __x) const {
    __stl_debug_check(__check_same_owner(*this,__x));    
    return __WORD_BIT * (_M_p - __x._M_p) + _M_offset - __x._M_offset;
  }
};

inline bool operator==(const _Bit_iterator_base& __x, const _Bit_iterator_base& __y) {
  __stl_debug_check(__check_same_owner_or_null(__x, __y));    
  return __y._M_p == __x._M_p && __y._M_offset == __x._M_offset;
}
inline bool operator!=(const _Bit_iterator_base& __x, const _Bit_iterator_base& __y) {
  __stl_debug_check(__check_same_owner_or_null(__x, __y));    
  return __y._M_p != __x._M_p || __y._M_offset != __x._M_offset;
}

inline bool operator<(const _Bit_iterator_base& __x, const _Bit_iterator_base& __y) {
  __stl_debug_check(__check_same_owner(__y,__x));    
  return __x._M_p < __y._M_p || (__x._M_p == __y._M_p && __x._M_offset < __y._M_offset);
}

inline bool operator>(const _Bit_iterator_base& __x, const _Bit_iterator_base& __y)  { 
  return __STLPORT_STD::operator <(__y , __x); 
}
inline bool operator<=(const _Bit_iterator_base& __x, const _Bit_iterator_base& __y) { 
  return !(__y < __x); 
}
inline bool operator>=(const _Bit_iterator_base& __x, const _Bit_iterator_base& __y) { 
  return !(__x < __y); 
}

# ifdef __STL_DEBUG
template <class _Dummy>
bool _Bv_global<_Dummy>::_Dereferenceable(const _Bit_iterator_base& __x) {
  __stl_verbose_return(__x._Valid(), _StlMsg_INVALID_ITERATOR);
  // hack : assumes _M_start is followed by _M_finish 
  _Bit_iterator_base* __start = (_Bit_iterator_base*)(__x._Owner()->_Owner());
  if (__x._Overrun_ok()) return true;
  __stl_verbose_return(__x >= *__start && __x < *(__start+1),
		       _StlMsg_NOT_DEREFERENCEABLE);    
  return true;
}

template <class _Dummy>
bool _Bv_global<_Dummy>::_Nonsingular(const _Bit_iterator_base& __x ) {
    __stl_verbose_return(__x._Valid(), _StlMsg_INVALID_ITERATOR);
    // hack : assumes _M_start is followed by _M_finish 
    _Bit_iterator_base* __start = (_Bit_iterator_base*)(__x._Owner()->_Owner());
    if (__x._Overrun_ok()) return true;
    __stl_verbose_return(__x >= *__start && __x <= *(__start+1),
			 _StlMsg_SINGULAR_ITERATOR);    
    return true;
}
# endif /* __STL_DEBUG */

template <class _Ref, class _Ptr>
struct _Bit_iter : public _Bit_iterator_base
{
  typedef _Ref  reference;
  typedef _Ptr  pointer;
  typedef _Bit_iter<_Bit_reference, _Bit_reference*> iterator;
  typedef _Bit_iter<bool, const bool*> const_iterator;
  typedef _Bit_iter<_Ref, _Ptr> _Self;
  typedef random_access_iterator_tag iterator_category;
  typedef bool value_type;
  typedef ptrdiff_t difference_type;
  typedef size_t size_type;

# if defined ( __STL_DEBUG )
  _Bit_iter(const __owned_list* __root,unsigned int* __x, 
	    unsigned int __y, bool __over = false) : _Bit_iterator_base(__root, __x, __y, __over) {}
# else
  _Bit_iter(unsigned int* __x, unsigned int __y) : _Bit_iterator_base(__x, __y) {}
# endif
  _Bit_iter() {}
  _Bit_iter(const _Bit_iter<_Bit_reference, _Bit_reference*>& __x): 
    _Bit_iterator_base((const _Bit_iterator_base&)__x) {}
  reference operator*() const { 
    __stl_debug_check(_Bv_global_inst::_Dereferenceable(*this));    
    return _Bit_reference(_M_p, 1UL << _M_offset); 
  }
  _Self& operator++() {
    _M_bump_up();
    return *this;
  }
  _Self operator++(int) {
    _Self __tmp = *this;
    _M_bump_up();
    return __tmp;
  }
  _Self& operator--() {
    _M_bump_down();
    return *this;
  }
  _Self operator--(int) {
    _Self __tmp = *this;
    _M_bump_down();
    return __tmp;
  }
  _Self& operator+=(difference_type __i) {
    _M_advance(__i);
    return *this;
  }
  _Self& operator-=(difference_type __i) {
    *this += -__i;
    return *this;
  }
  _Self operator+(difference_type __i) const {
    _Self __tmp = *this;
    return __tmp += __i;
  }
  _Self operator-(difference_type __i) const {
    _Self __tmp = *this;
    return __tmp -= __i;
  }
  difference_type operator-(const _Self& __x) const {
    return _M_subtract(__x);
  }
  reference operator[](difference_type __i) { return *(*this + __i); }
};

template <class _Ref, class _Ptr>
inline _Bit_iter<_Ref,_Ptr> operator+(ptrdiff_t __n, const _Bit_iter<_Ref, _Ptr>& __x) {
   return __x + __n;
}

# ifndef __STL_CLASS_PARTIAL_SPECIALIZATION
inline random_access_iterator_tag 
iterator_category(const _Bit_iterator_base&) {return random_access_iterator_tag();}
inline ptrdiff_t* 
distance_type(const _Bit_iterator_base&) {return (ptrdiff_t*)0;}
inline bool* value_type(const _Bit_iter<_Bit_reference, _Bit_reference*>&) {return (bool*)0;}
inline bool* value_type(const _Bit_iter<bool, const bool*>&) {return (bool*)0;}
# endif

typedef _Bit_iter<bool, const bool*> _Bit_const_iterator;
typedef _Bit_iter<_Bit_reference, _Bit_reference*> _Bit_iterator;



// Bit-vector base class, which encapsulates the difference between
//  old SGI-style allocators and standard-conforming allocators.


template <class _Alloc>
class _Bvector_base
{
public:
  typedef typename _Alloc_traits<bool, _Alloc>::allocator_type allocator_type;
  typedef unsigned int __chunk_type;
  typedef typename _Alloc_traits<__chunk_type, 
          _Alloc>::allocator_type __chunk_allocator_type;
  allocator_type get_allocator() const { 
    return __STL_CONVERT_ALLOCATOR((const __chunk_allocator_type&)_M_end_of_storage, bool); 
  }
  _Bvector_base(const allocator_type& __a)
    : _M_start(), _M_finish(), _M_end_of_storage(__STL_CONVERT_ALLOCATOR(__a, __chunk_type),
						 (__chunk_type*)0) {
      __stl_debug_do(_M_iter_list._Safe_init(&_M_start));        
      __stl_debug_do(_Init_bounds());
  }
  ~_Bvector_base() { _M_deallocate();
  __stl_debug_do(_M_start._Invalidate()); 
  __stl_debug_do(_M_finish._Invalidate());
  }

protected:

  unsigned int* _M_bit_alloc(size_t __n) 
    { return _M_end_of_storage.allocate((__n + __WORD_BIT - 1)/__WORD_BIT); }
  void _M_deallocate() {
    if (_M_start._M_p)
      _M_end_of_storage.deallocate(_M_start._M_p,
				   _M_end_of_storage._M_data - _M_start._M_p);
    __stl_debug_do(_Invalidate_all());
  }

  _Bit_iterator _M_start;
  _Bit_iterator _M_finish;
  _STL_alloc_proxy<__chunk_type*, __chunk_type, __chunk_allocator_type> _M_end_of_storage;  

# if defined (__STL_DEBUG)
  __owned_list _M_iter_list;
  void _Init_bounds() {
    _M_iter_list._Orphan(_M_start);
    _M_iter_list._Orphan(_M_finish);
  }
  void _Invalidate_iterator(const _Bit_iterator& __it) { 
    __invalidate_iterator(&_M_iter_list,__it);
  }
  void _Invalidate_all() {
    _M_iter_list._Invalidate_all();
  }
# endif
};


// The next few lines are confusing.  What we're doing is declaring a
//  partial specialization of vector<T, Alloc> if we have the necessary
//  compiler support.  Otherwise, we define a class bit_vector which uses
//  the default allocator. 

#if defined(__STL_CLASS_PARTIAL_SPECIALIZATION) \
 && !defined(__STL_PARTIAL_SPECIALIZATION_BUG) && ! defined(__STL_NO_BOOL) && \
  ! defined(CGAL_LIMITED_ITERATOR_TRAITS_SUPPORT)
# define __SGI_STL_VECBOOL_TEMPLATE
# define __BVEC_TMPL_HEADER template <class _Alloc>
#else
# undef __SGI_STL_VECBOOL_TEMPLATE
# ifdef __STL_NO_BOOL
#  define __BVEC_TMPL_HEADER
# else
#  define __BVEC_TMPL_HEADER __STL_TEMPLATE_NULL
# endif
# define _Alloc __STL_DEFAULT_ALLOCATOR(bool)
#endif

#ifdef __STL_NO_BOOL
# define __BVECTOR_QUALIFIED bit_vector
# define __BVECTOR           bit_vector
#else
# ifdef __SGI_STL_VECBOOL_TEMPLATE
#  define __BVECTOR_QUALIFIED vector<bool, _Alloc>
# else
#  define __BVECTOR_QUALIFIED vector<bool>
# endif
#ifdef __STL_PARTIAL_SPEC_NEEDS_TEMPLATE_ARGS
# define __BVECTOR __BVECTOR_QUALIFIED
#else
# define __BVECTOR vector
#endif
#endif


__BVEC_TMPL_HEADER
class __BVECTOR_QUALIFIED : public _Bvector_base<_Alloc >
{
  typedef _Bvector_base<_Alloc > _Base;
  typedef __BVECTOR_QUALIFIED _Self;
public:
  typedef bool value_type;
  typedef size_t size_type;
  typedef ptrdiff_t difference_type; 
  typedef _Bit_reference reference;
  typedef bool const_reference;
  typedef _Bit_reference* pointer;
  typedef const bool* const_pointer;

  typedef _Bit_iterator                iterator;
  typedef _Bit_const_iterator          const_iterator;

#if (defined ( __STL_CLASS_PARTIAL_SPECIALIZATION ) && \
  ! defined (__STL_PARTIAL_SPECIALIZATION_BUG)) || \
  defined(CGAL_LIMITED_ITERATOR_TRAITS_SUPPORT)
  typedef __STLPORT_STD::reverse_iterator<const_iterator> const_reverse_iterator;
  typedef __STLPORT_STD::reverse_iterator<iterator> reverse_iterator;
#else /* __STL_CLASS_PARTIAL_SPECIALIZATION */
#error " wrong reverse_iterator! "
# if defined (__STL_MSVC50_COMPATIBILITY)
  typedef __STLPORT_STD::reverse_iterator<const_iterator, value_type, const_reference, 
  const_pointer, difference_type> const_reverse_iterator;
  typedef __STLPORT_STD::reverse_iterator<iterator, value_type, reference, reference*, 
  difference_type> reverse_iterator;
# else
#error " wrong reverse_iterator! "
  typedef __STLPORT_STD::reverse_iterator<const_iterator, value_type, const_reference, 
                                  difference_type> const_reverse_iterator;
  typedef __STLPORT_STD::reverse_iterator<iterator, value_type, reference, difference_type>
          reverse_iterator;
# endif
#endif /* __STL_CLASS_PARTIAL_SPECIALIZATION */

# ifdef __SGI_STL_VECBOOL_TEMPLATE
  typedef typename _Base::allocator_type allocator_type;
# else
  typedef _Base::allocator_type allocator_type;
# endif

#if defined( __STL_HAS_NAMESPACES )
  __STL_USING_BASE_MEMBER _Bvector_base<_Alloc >::get_allocator;
protected:
  __STL_USING_BASE_MEMBER _Bvector_base<_Alloc >::_M_bit_alloc;
  __STL_USING_BASE_MEMBER _Bvector_base<_Alloc >::_M_deallocate;
  __STL_USING_BASE_MEMBER _Bvector_base<_Alloc >::_M_start;
  __STL_USING_BASE_MEMBER _Bvector_base<_Alloc >::_M_finish;
  __STL_USING_BASE_MEMBER _Bvector_base<_Alloc >::_M_end_of_storage;
# ifdef __STL_DEBUG
  __STL_USING_BASE_MEMBER _Bvector_base<_Alloc >::_M_iter_list;
# endif
#endif /* __STL_HAS_NAMESPACES */

protected:
# if defined (__STL_DEBUG)
# define _Make_iterator(__x,__y,__z) iterator(&_M_iter_list,__x,__y,__z)
# define _Make_const_iterator(__x,__y,__z) const_iterator(&_M_iter_list,__x,__y,__z)
# else
# define _Make_iterator(__x,__y,__z) iterator(__x,__y)
# define _Make_const_iterator(__x,__y,__z) const_iterator(__x,__y)
# endif

  void _M_initialize(size_type __n) {
    unsigned int* __q = _M_bit_alloc(__n);
    _M_end_of_storage._M_data = __q + (__n + __WORD_BIT - 1)/__WORD_BIT;
    __stl_debug_do(_M_iter_list._Safe_init(&_M_start));
    _M_start = _Make_iterator(__q, 0, false);
    __stl_debug_do(_M_start._Set_overrun(true));
    _M_finish = _M_start + difference_type(__n);
    __stl_debug_do(_Init_bounds());
  }
  void _M_insert_aux(iterator __position, bool __x) {
    __stl_debug_do(_M_finish._Set_overrun(true));
    if (_M_finish._M_p != _M_end_of_storage._M_data) {
      copy_backward(__position, _M_finish, _M_finish + 1);
      *__position = __x;
      ++_M_finish;
    }
    else {
      size_type __len = size() ? 2 * size() : __WORD_BIT;
      unsigned int* __q = _M_bit_alloc(__len);
      iterator __i = copy(begin(), __position, _Make_iterator(__q, 0, true));
      *__i++ = __x;
      _M_finish = copy(__position, end(), __i);
      __stl_debug_do(_M_iter_list._Orphan(_M_finish));
      _M_deallocate();
      _M_end_of_storage._M_data = __q + (__len + __WORD_BIT - 1)/__WORD_BIT;
      _M_start = _Make_iterator(__q, 0, false);
      __stl_debug_do(_M_iter_list._Orphan(_M_start));
    }
    __stl_debug_do(_M_finish._Set_overrun(false));
  }

#ifdef __STL_MEMBER_TEMPLATES
  template <class _InputIterator>
  void _M_initialize_range(_InputIterator __first, _InputIterator __last,
			input_iterator_tag) {
    __stl_debug_check(__check_range(__first,__last));
    __stl_debug_do(_Safe_init(&_M_start));
    _M_start = iterator();
    _M_finish = iterator();
    __stl_debug_do(_Init_bounds());
    _M_end_of_storage._M_data = 0;
# if defined (__STL_DEBUG) && defined (__STL_USE_EXCEPTIONS)
    try {
# endif
      for ( ; __first != __last; ++__first) 
        push_back(*__first);
# if defined (__STL_DEBUG) && defined (__STL_USE_EXCEPTIONS)
    }
    catch (...) {
      __stl_debug_do(_M_start._Invalidate());
      __stl_debug_do(_M_finish._Invalidate());
      __stl_debug_do(_M_iter_list._Invalidate());        
    }
# endif
  }

  template <class _ForwardIterator>
  void _M_initialize_range(_ForwardIterator __first, _ForwardIterator __last,
                           forward_iterator_tag) {
    size_type __n = 0;
    distance(__first, __last, __n);
    _M_initialize(__n);
    copy(__first, __last, _M_start);
  }

  template <class _InputIterator>
  void _M_insert_range(iterator __pos,
                       _InputIterator __first, _InputIterator __last,
                       input_iterator_tag) {
    for ( ; __first != __last; ++__first) {
      __pos = insert(__pos, *__first);
      ++__pos;
    }
  }

  template <class _ForwardIterator>
  void _M_insert_range(iterator __position,
                       _ForwardIterator __first, _ForwardIterator __last,
                       forward_iterator_tag) {
    __stl_debug_check(__check_range(__first,__last));
    if (__first != __last) {
      size_type __n = 0;
      distance(__first, __last, __n);
      __stl_debug_do(_M_finish._Set_overrun(true));
      __stl_debug_do(__position._Set_overrun(true));
      if (capacity() - size() >= __n) {
        copy_backward(__position, end(), _M_finish + difference_type(__n));
        copy(__first, __last, __position);
        _M_finish += difference_type(__n);
      }
      else {
        size_type __len = size() + max(size(), __n);
        unsigned int* __q = _M_bit_alloc(__len);
        iterator __i = copy(begin(), __position, _Make_iterator(__q, 0, true));
        __i = copy(__first, __last, __i);
        _M_finish = copy(__position, end(), __i);
	__stl_debug_do(_M_iter_list._Orphan(_M_finish));
        _M_deallocate();
        _M_end_of_storage._M_data = __q + (__len + __WORD_BIT - 1)/__WORD_BIT;
        _M_start = _Make_iterator(__q, 0, false);
	__stl_debug_do(_M_iter_list._Orphan(_M_start));
      }
    }
    __stl_debug_do(_M_finish._Set_overrun(false));
  }      

#endif /* __STL_MEMBER_TEMPLATES */

public:
  iterator begin() { return _M_start; }
  const_iterator begin() const { return *(const const_iterator*)&_M_start; }
  iterator end() { return _M_finish; }
  const_iterator end() const { return *(const const_iterator*)&_M_finish; }

  reverse_iterator rbegin() { return reverse_iterator(end()); }
  const_reverse_iterator rbegin() const { 
    return const_reverse_iterator(end()); 
  }
  reverse_iterator rend() { return reverse_iterator(begin()); }
  const_reverse_iterator rend() const { 
    return const_reverse_iterator(begin()); 
  }

  size_type size() const { return size_type(end() - begin()); }
  size_type max_size() const { return size_type(-1); }
  size_type capacity() const {
    return size_type(_Make_const_iterator(_M_end_of_storage._M_data, 0, false) - begin());
  }
  bool empty() const { return begin() == end(); }
  reference operator[](size_type __n) 
  { return *(begin() + difference_type(__n)); }
  const_reference operator[](size_type __n) const 
  { return *(begin() + difference_type(__n)); }

#ifdef __STL_THROW_RANGE_ERRORS
  void _M_range_check(size_type __n) const {
    if (__n >= this->size())
      __stl_throw_range_error("vector<bool>");
  }

  reference at(size_type __n)
    { _M_range_check(__n); return (*this)[__n]; }
  const_reference at(size_type __n) const
    { _M_range_check(__n); return (*this)[__n]; }
#endif /* __STL_THROW_RANGE_ERRORS */


  explicit __BVECTOR(const allocator_type& __a = 
		     __STL_ALLOC_INSTANCE(allocator_type))
    : _Bvector_base<_Alloc >(__a) {}

  __BVECTOR(size_type __n, bool __value,
            const allocator_type& __a = 
	    __STL_ALLOC_INSTANCE(allocator_type))
    : _Bvector_base<_Alloc >(__a)
  {
    _M_initialize(__n);
    fill(_M_start._M_p, (__chunk_type*)_M_end_of_storage._M_data, __value ? ~0 : 0);
  }

  explicit __BVECTOR(size_type __n)
    : _Bvector_base<_Alloc >(allocator_type())
  {
    _M_initialize(__n);
    fill(_M_start._M_p, (__chunk_type*)_M_end_of_storage._M_data, 0);
  }

  __BVECTOR(const _Self& __x) : _Bvector_base<_Alloc >(__x.get_allocator()) {
    _M_initialize(__x.size());
    copy(__x.begin(), __x.end(), _M_start);
  }

#if defined (__STL_MEMBER_TEMPLATES) && !(defined (__GNUC__) && __GNUC_MINOR__ < 90) 
  // Check whether it's an integral type.  If so, it's not an iterator.
  template <class _InputIterator>
  __BVECTOR(_InputIterator __first, _InputIterator __last,
            const allocator_type& __a = __STL_ALLOC_INSTANCE(allocator_type))
    : _Base(__a)
  {
    typedef typename _Is_integer<_InputIterator>::_Integral _Integral;
    _M_initialize_dispatch(__first, __last, _Integral());
  }
    
  template <class _Integer>
  void _M_initialize_dispatch(_Integer __n, _Integer __x, __true_type) {
    _M_initialize(__n);
    fill(_M_start._M_p, _M_end_of_storage._M_data, __x ? ~0 : 0);
  }
    
  template <class _InputIterator>
  void _M_initialize_dispatch(_InputIterator __first, _InputIterator __last,
                              __false_type) {
    _M_initialize_range(__first, __last, __ITERATOR_CATEGORY(__first));
  }
#else /* __STL_MEMBER_TEMPLATES */
  __BVECTOR(const_iterator __first, const_iterator __last,
            const allocator_type& __a = __STL_ALLOC_INSTANCE(allocator_type))
    : _Bvector_base<_Alloc >(__a)
  {
    size_type __n = 0;
    distance(__first, __last, __n);
    _M_initialize(__n);
    copy(__first, __last, _M_start);
  }
  __BVECTOR(const bool* __first, const bool* __last,
            const allocator_type& __a = __STL_ALLOC_INSTANCE(allocator_type))
    : _Bvector_base<_Alloc >(__a)
  {
    size_type __n = 0;
    distance(__first, __last, __n);
    _M_initialize(__n);
    copy(__first, __last, _M_start);
  }
#endif /* __STL_MEMBER_TEMPLATES */

  ~__BVECTOR() { }

  __BVECTOR_QUALIFIED& operator=(const __BVECTOR_QUALIFIED& __x) {
    if (&__x == this) return *this;
    __stl_debug_do(_Invalidate_all());
    if (__x.size() > capacity()) {
      _M_deallocate();
      _M_initialize(__x.size());
    }
    copy(__x.begin(), __x.end(), begin());
    _M_finish = begin() + difference_type(__x.size());
    __stl_debug_do(_M_iter_list._Orphan(_M_finish));
    return *this;
  }

  // assign(), a generalized assignment member function.  Two
  // versions: one that takes a count, and one that takes a range.
  // The range version is a member template, so we dispatch on whether
  // or not the type is an integer.

  void _M_fill_assign(size_t __n, bool __x) {
    if (__n > size()) {
      fill(_M_start._M_p, (__chunk_type*)_M_end_of_storage._M_data, __x ? ~0 : 0);
      insert(end(), __n - size(), __x);
    }
    else {
      erase(begin() + __n, end());
      fill(_M_start._M_p, (__chunk_type*)_M_end_of_storage._M_data, __x ? ~0 : 0);
    }
  }
  void assign(size_t __n, bool __x) { _M_fill_assign(__n, __x); }

#ifdef __STL_MEMBER_TEMPLATES

  template <class _InputIterator>
  void assign(_InputIterator __first, _InputIterator __last) {
    typedef typename _Is_integer<_InputIterator>::_Integral _Integral;
    _M_assign_dispatch(__first, __last, _Integral());
  }

  template <class _Integer>
  void _M_assign_dispatch(_Integer __n, _Integer __val, __true_type)
    { _M_fill_assign((size_t) __n, (bool) __val); }

  template <class _InputIter>
  void _M_assign_dispatch(_InputIter __first, _InputIter __last, __false_type)
    { _M_assign_aux(__first, __last, __ITERATOR_CATEGORY(__first)); }

  template <class _InputIterator>
  void _M_assign_aux(_InputIterator __first, _InputIterator __last,
                     input_iterator_tag) {
    iterator __cur = begin();
    for ( ; __first != __last && __cur != end(); ++__cur, ++__first)
      *__cur = *__first;
    if (__first == __last)
      erase(__cur, end());
    else
      insert(end(), __first, __last);
  }

  template <class _ForwardIterator>
  void _M_assign_aux(_ForwardIterator __first, _ForwardIterator __last,
                     forward_iterator_tag) {
    size_type __len = 0;
    distance(__first, __last, __len);
    if (__len < size())
      erase(copy(__first, __last, begin()), end());
    else {
      _ForwardIterator __mid = __first;
      advance(__mid, size());
      copy(__first, __mid, begin());
      insert(end(), __mid, __last);
    }
  }    

#endif /* __STL_MEMBER_TEMPLATES */

  void reserve(size_type __n) {
    if (capacity() < __n) {
      unsigned int* __q = _M_bit_alloc(__n);
# ifdef __STL_DEBUG
      _Bit_iterator __z(&_M_iter_list, __q, 0, true);
# else
      _Bit_iterator __z(__q, 0);
# endif
      _M_finish = copy(begin(), end(), __z);
      __stl_debug_do(_M_finish._Set_overrun(false));
      __stl_debug_do(_M_iter_list._Orphan(_M_finish));
      _M_deallocate();
      _M_start = _Make_iterator(__q, 0, false);
      __stl_debug_do(_M_iter_list._Orphan(_M_start));
      _M_end_of_storage._M_data = __q + (__n + __WORD_BIT - 1)/__WORD_BIT;
    }
  }

  reference front() { return *begin(); }
  const_reference front() const { return *begin(); }
  reference back() { return *(end() - 1); }
  const_reference back() const { return *(end() - 1); }
  void push_back(bool __x) {
    if (_M_finish._M_p != _M_end_of_storage._M_data) {
      __stl_debug_do(_M_finish._Set_overrun(true));
      *_M_finish++ = __x;
      __stl_debug_do(_M_finish._Set_overrun(false));
    }
    else
      _M_insert_aux(end(), __x);
  }
  void swap(__BVECTOR_QUALIFIED& __x) {
    __STLPORT_STD::swap(_M_start, __x._M_start);
    __STLPORT_STD::swap(_M_finish, __x._M_finish);
    __STLPORT_STD::swap(_M_end_of_storage._M_data, __x._M_end_of_storage._M_data);
    __stl_debug_do(_M_iter_list._Swap_owners(__x._M_iter_list));
  }
  iterator insert(iterator __position, bool __x = bool()) {
    __stl_debug_check(__check_if_owner(&_M_iter_list,__position));
    difference_type __n = __position - begin();
    if (_M_finish._M_p != _M_end_of_storage._M_data && __position == end()) {
      __stl_debug_do(_M_finish._Set_overrun(true));
      *_M_finish++ = __x;
      __stl_debug_do(_M_finish._Set_overrun(false));
    }
    else
      _M_insert_aux(__position, __x);
    return begin() + __n;
  }

#if defined ( __STL_MEMBER_TEMPLATES ) && !(defined (__GNUC__) && (__GNUC_MINOR__ < 90)) 

  template <class _Integer>
  void _M_insert_dispatch(iterator __pos, _Integer __n, _Integer __x,
                          __true_type) {
    _M_fill_insert(__pos, (size_type) __n, (bool) __x);
  }

  template <class _InputIterator>
  void _M_insert_dispatch(iterator __pos,
                          _InputIterator __first, _InputIterator __last,
                          __false_type) {
    _M_insert_range(__pos, __first, __last, __ITERATOR_CATEGORY(__first));
  }

  // Check whether it's an integral type.  If so, it's not an iterator.
  template <class _InputIterator>
  void insert(iterator __position,
              _InputIterator __first, _InputIterator __last) {
    typedef typename _Is_integer<_InputIterator>::_Integral _Is_Integral;
    _M_insert_dispatch(__position, __first, __last, _Is_Integral());
  }
#else /* __STL_MEMBER_TEMPLATES */
  void insert(iterator __position,
              const_iterator __first, const_iterator __last) {
    __stl_debug_check(__check_if_owner(&_M_iter_list,__position));
    if (__first == __last) return;
    size_type __n = 0;
    distance(__first, __last, __n);
    __stl_debug_do(_M_finish._Set_overrun(true));
    __stl_debug_do(__position._Set_overrun(true));
    if (capacity() - size() >= __n) {
      copy_backward(__position, end(), _M_finish + __n);
      copy(__first, __last, __position);
      _M_finish += __n;
    }
    else {
      size_type __len = size() + max(size(), __n);
      unsigned int* __q = _M_bit_alloc(__len);
      iterator __i = copy(begin(), __position, _Make_iterator(__q, 0, true));
      __i = copy(__first, __last, __i);
      _M_finish = copy(__position, end(), __i);
      __stl_debug_do(_M_iter_list._Orphan(_M_finish));
      _M_deallocate();
      _M_end_of_storage._M_data = __q + (__len + __WORD_BIT - 1)/__WORD_BIT;
      _M_start = _Make_iterator(__q, 0, false);
      __stl_debug_do(_M_iter_list._Orphan(_M_start));
    }
    __stl_debug_do(_M_finish._Set_overrun(false));
  }

  void insert(iterator __position, const bool* __first, const bool* __last) {
    __stl_debug_check(__check_if_owner(&_M_iter_list,__position));
    if (__first == __last) return;
    size_type __n = 0;
    distance(__first, __last, __n);
    __stl_debug_do(_M_finish._Set_overrun(true));
    __stl_debug_do(__position._Set_overrun(true));
    if (capacity() - size() >= __n) {
      copy_backward(__position, end(), _M_finish + __n);
      copy(__first, __last, __position);
      _M_finish += __n;
    }
    else {
      size_type __len = size() + max(size(), __n);
      unsigned int* __q = _M_bit_alloc(__len);
      iterator __i = copy(begin(), __position, _Make_iterator(__q, 0, true));
      __i = copy(__first, __last, __i);
      _M_finish = copy(__position, end(), __i);
      __stl_debug_do(_M_iter_list._Orphan(_M_finish));
      _M_deallocate();
      _M_end_of_storage._M_data = __q + (__len + __WORD_BIT - 1)/__WORD_BIT;
      _M_start = _Make_iterator(__q, 0, false);
      __stl_debug_do(_M_iter_list._Orphan(_M_start));
    }
    __stl_debug_do(_M_finish._Set_overrun(false));
  }
#endif /* __STL_MEMBER_TEMPLATES */
  
  void _M_fill_insert(iterator __position, size_type __n, bool __x) {
    __stl_debug_check(__check_if_owner(&_M_iter_list,__position));
    if (__n == 0) return;
    __stl_debug_do(_M_finish._Set_overrun(true));
    __stl_debug_do(__position._Set_overrun(true));
    if (capacity() - size() >= __n) {
      copy_backward(__position, end(), _M_finish + difference_type(__n));
      fill(__position, __position + difference_type(__n), __x);
      _M_finish += difference_type(__n);
    }
    else {
      size_type __len = size() + max(size(), __n);
      unsigned int* __q = _M_bit_alloc(__len);
      iterator __i = copy(begin(), __position, _Make_iterator(__q, 0, true));
      fill_n(__i, __n, __x);
      _M_finish = copy(__position, end(), __i + difference_type(__n));
      __stl_debug_do(_M_iter_list._Orphan(_M_finish));
      _M_deallocate();
      _M_end_of_storage._M_data = __q + (__len + __WORD_BIT - 1)/__WORD_BIT;
      _M_start = _Make_iterator(__q, 0, false);
      __stl_debug_do(_M_iter_list._Orphan(_M_start));
    }
    __stl_debug_do(_M_finish._Set_overrun(false));
  }

  void insert(iterator __position, size_type __n, bool __x) {
    _M_fill_insert(__position, __n, __x);
  }

  void pop_back() { 
        __stl_verbose_assert(!empty(), _StlMsg_EMPTY_CONTAINER);
	--_M_finish; 
  }
  iterator erase(iterator __position) {
    __stl_debug_check(__check_if_owner(&_M_iter_list,__position));
    __stl_verbose_assert(__position!= end(),_StlMsg_ERASE_PAST_THE_END);
    if (__position + 1 != end())
      copy(__position + 1, end(), __position);
      --_M_finish;
    return __position;
  }
  iterator erase(iterator __first, iterator __last) {
    __stl_debug_check(__check_if_owner(&_M_iter_list,__first));
    _M_finish = copy(__last, end(), __first);
    __stl_debug_do(_M_iter_list._Orphan(_M_finish));
    return __first;
  }
  void resize(size_type __new_size, bool __x = bool()) {
    if (__new_size < size()) 
      erase(begin() + difference_type(__new_size), end());
    else
      insert(end(), __new_size - size(), __x);
  }
  void flip() {
    for (unsigned int* __p = _M_start._M_p; __p != _M_end_of_storage._M_data; ++__p)
      *__p = ~*__p;
  }

  void clear() { erase(begin(), end()); }
#ifdef __BORLANDC__
  inline bool 
  operator==(const __BVECTOR_QUALIFIED& __y) const
    {
      return (size() == __y.size() && 
	      equal(begin(), end(), __y.begin()));
    }
  
  inline bool 
  operator<(const __BVECTOR_QUALIFIED& __y) const
    {
      return lexicographical_compare(begin(), end(), 
				     __y.begin(), __y.end());
    }
#endif
};

#if !defined (__STL_NO_BOOL)

// This typedef is non-standard.  It is provided for backward compatibility.
# ifdef __SGI_STL_VECBOOL_TEMPLATE
typedef vector<bool, allocator<bool> > bit_vector;
# else
typedef vector<bool> bit_vector;
# endif
#endif

__BVEC_TMPL_HEADER
inline void swap(__BVECTOR_QUALIFIED& __x, __BVECTOR_QUALIFIED& __y) {
  __x.swap(__y);
}

#ifndef __BORLANDC__
__BVEC_TMPL_HEADER
inline bool 
operator==(const __BVECTOR_QUALIFIED& __x, const __BVECTOR_QUALIFIED& __y)
{
  return (__x.size() == __y.size() && 
          equal(__x.begin(), __x.end(), __y.begin()));
}


__BVEC_TMPL_HEADER
inline bool 
operator<(const __BVECTOR_QUALIFIED& __x, const __BVECTOR_QUALIFIED& __y)
{
  return lexicographical_compare(__x.begin(), __x.end(), 
                                 __y.begin(), __y.end());
}
# endif /* __BORLANDC__ */

# ifdef __STL_USE_SEPARATE_RELOPS_NAMESPACE

__BVEC_TMPL_HEADER
inline bool 
operator!=(const __BVECTOR_QUALIFIED& __x, const __BVECTOR_QUALIFIED& __y)
{
  return !(__x == __y);
}

__BVEC_TMPL_HEADER
inline bool operator>(const __BVECTOR_QUALIFIED& __x, const __BVECTOR_QUALIFIED& __y)
{
  return __y < __x;
}

__BVEC_TMPL_HEADER
inline bool operator<=(const __BVECTOR_QUALIFIED& __x, const __BVECTOR_QUALIFIED& __y)
{
  return !(__y < __x);
}

__BVEC_TMPL_HEADER
inline bool operator>=(const __BVECTOR_QUALIFIED& __x, const __BVECTOR_QUALIFIED& __y)
{
  return !(__x < __y);
}
#endif

#undef _Alloc
#undef __SGI_STL_VECBOOL_TEMPLATE
#undef __BVECTOR
#undef __BVECTOR_QUALIFIED
#undef __BVEC_TMPL_HEADER


#if defined(__sgi) && !defined(__GNUC__) && (_MIPS_SIM != _MIPS_SIM_ABI32)
#pragma reset woff 1174
#pragma reset woff 1375
#endif

# undef __WORD_BIT
# undef _Make_iterator
# undef _Make_const_iterator

__STL_END_NAMESPACE

#endif /* __SGI_STL_INTERNAL_BVECTOR_H */

// Local Variables:
// mode:C++
// End:
