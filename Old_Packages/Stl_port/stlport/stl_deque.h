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

#ifndef __SGI_STL_INTERNAL_DEQUE_H
#define __SGI_STL_INTERNAL_DEQUE_H

# ifndef __SGI_STL_INTERNAL_ALGOBASE_H
#  include <stl_algobase.h>
# endif

# ifndef __SGI_STL_INTERNAL_ALLOC_H
#  include <stl_alloc.h>
# endif

# ifndef __SGI_STL_INTERNAL_ITERATOR_H
#  include <stl_iterator.h>
# endif

# ifndef __SGI_STL_INTERNAL_UNINITIALIZED_H
#  include <stl_uninitialized.h>
# endif

# ifndef __STL_RANGE_ERRORS_H
#  include <stl_range_errors.h>
# endif

/* Class invariants:
 *  For any nonsingular iterator i:
 *    i.node is the address of an element in the map array.  The
 *      contents of i.node is a pointer to the beginning of a node.
 *    i.first == *(i.node) 
 *    i.last  == i.first + node_size
 *    i.cur is a pointer in the range [i.first, i.last).  NOTE:
 *      the implication of this is that i.cur is always a dereferenceable
 *      pointer, even if i is a past-the-end iterator.
 *  Start and Finish are always nonsingular iterators.  NOTE: this means
 *    that an empty deque must have one node, and that a deque
 *    with N elements, where N is the buffer size, must have two nodes.
 *  For every node other than start.node and finish.node, every element
 *    in the node is an initialized object.  If start.node == finish.node,
 *    then [start.cur, finish.cur) are initialized objects, and
 *    the elements outside that range are uninitialized storage.  Otherwise,
 *    [start.cur, start.last) and [finish.first, finish.cur) are initialized
 *    objects, and [start.first, start.cur) and [finish.cur, finish.last)
 *    are uninitialized storage.
 *  [map, map + map_size) is a valid, non-empty range.  
 *  [start.node, finish.node] is a valid range contained within 
 *    [map, map + map_size).  
 *  A pointer in the range [map, map + map_size) points to an allocated node
 *    if and only if the pointer is in the range [start.node, finish.node].
 */


/*
 * In previous versions of deque, node_size was fixed by the 
 * implementation.  In this version, however, users can select
 * the node size.  Deque has three template parameters; the third,
 * a number of type size_t, is the number of elements per node.
 * If the third template parameter is 0 (which is the default), 
 * then deque will use a default node size.
 *
 * The only reason for using an alternate node size is if your application
 * requires a different performance tradeoff than the default.  If,
 * for example, your program contains many deques each of which contains
 * only a few elements, then you might want to save memory (possibly
 * by sacrificing some speed) by using smaller nodes.
 *
 * Unfortunately, some compilers have trouble with non-type template 
 * parameters; stl_config.h defines __STL_NON_TYPE_TMPL_PARAM_BUG if
 * that is the case.  If your compiler is one of them, then you will
 * not be able to use alternate node sizes; you will have to use the
 * default value.
 */

# undef __deque__
# undef deque
# if defined ( __STL_NO_DEFAULT_NON_TYPE_PARAM )
#  define deque __deque
#  define __deque__ __deque
# else
#  define __deque__ __FULL_NAME(deque)
#  define deque __WORKAROUND_RENAME(deque)
# endif

__STL_BEGIN_NAMESPACE

#if defined(__sgi) && !defined(__GNUC__) && (_MIPS_SIM != _MIPS_SIM_ABI32)
#pragma set woff 1174
#pragma set woff 1375
#endif

# if defined ( __STL_USE_ABBREVS )
#  define __deque_iterator         _dQ__It
#  define _Buf_traits              _dQ__BTr
#  define _Deque_iterator          _Dq__It
# endif

// Note: this function is simply a kludge to work around several compilers'
//  bugs in handling constant expressions.
inline size_t
__deque_buf_size(size_t __n, size_t size)
{
  return __n != 0 ? __n : (size < 512 ? size_t(512 / size) : size_t(1));
}

// this helper class is needed to pass deque not-type parameter
// to its iterators. 

template <class _Tp, size_t _BufSize>
struct _Buf_size_traits {
public:
  enum _Constants { 
    _blocksize = 512, 
    _buf_size = (_BufSize != 0 ? _BufSize : (sizeof(_Tp) < (size_t)_blocksize ?
					     ( (size_t)_blocksize / sizeof(_Tp)) : size_t(1))), 
    _byte_buf_size = sizeof(_Tp)*_buf_size 
  };
  static size_t buffer_size() { return (size_t)_buf_size; } 
};

# ifdef __STL_DEBUG
template <class _Tp>
bool __Deq_dereferenceable(const void* __ptr, _Tp*);
template <class _Tp>
bool __Deq_nonsingular(const void* __ptr, _Tp*);
# endif

template <class _Tp, class __bufsiz>
# if defined ( __STL_DEBUG )
struct _Deque_iterator_base : public __owned_link {
# else
struct _Deque_iterator_base {
# endif

  enum { __buffer_size = __bufsiz::_buf_size } ;

  typedef random_access_iterator_tag iterator_category;

  typedef _Tp value_type;
  typedef size_t size_type;
  typedef ptrdiff_t difference_type;

  typedef value_type** _Map_pointer;

  typedef _Deque_iterator_base< _Tp, __bufsiz> _Self;

  value_type* _M_cur;
  value_type* _M_first;
  value_type* _M_last;
  _Map_pointer _M_node;

# if defined ( __STL_DEBUG )
  bool _M_unsafe;

  _Deque_iterator_base(const __owned_list* __root, value_type* __x, _Map_pointer __y) 
    : __owned_link(__root),_M_cur(__x), _M_first(*__y), _M_last(*__y + __buffer_size), 
      _M_node(__y), _M_unsafe(false) {}
  _Deque_iterator_base() : __owned_link(0), _M_cur(0), _M_first(0), _M_last(0), 
          _M_node(0), _M_unsafe(false) {}
# else
  _Deque_iterator_base(value_type* __x, _Map_pointer __y) 
    : _M_cur(__x), _M_first(*__y),
      _M_last(*__y + __buffer_size), _M_node(__y) {}
  _Deque_iterator_base() : _M_cur(0), _M_first(0), _M_last(0), _M_node(0) {}
# endif

  difference_type _M_subtract(const _Self& __x) const {
    __stl_debug_check(__check_same_owner(*this,__x));    
    return difference_type(__buffer_size) * (_M_node - __x._M_node - 1) +
      (_M_cur - _M_first) + (__x._M_last - __x._M_cur);
  }

  void _M_increment() {
    //    ++_M_cur;
    if (++_M_cur == _M_last) {
      _M_set_node(_M_node + 1);
      _M_cur = _M_first;
    }
    __stl_debug_check(__Deq_nonsingular(this,(value_type*)0));
  }

  void _M_decrement() {
    if (_M_cur == _M_first) {
      _M_set_node(_M_node - 1);
      _M_cur = _M_last;
    }
    --_M_cur;
    __stl_debug_check(__Deq_nonsingular(this,(value_type*)0));
  }

  void _M_advance(difference_type __n)
  {
    difference_type __offset = __n + (_M_cur - _M_first);
    if (__offset >= 0 && __offset < difference_type(__buffer_size))
      _M_cur += __n;
    else {
      difference_type __node_offset =
        __offset > 0 ? __offset / __buffer_size
                   : -difference_type((-__offset - 1) / __buffer_size) - 1;
      _M_set_node(_M_node + __node_offset);
      _M_cur = _M_first + 
        (__offset - __node_offset * difference_type(__buffer_size));
    }
    __stl_debug_check(__Deq_nonsingular(this,(_Tp*)0));
  }

  void _M_set_node(_Map_pointer __new_node) {
    _M_last = (_M_first = *(_M_node = __new_node)) + difference_type(__buffer_size);
  }
};



template <class _Tp, class _Traits, class __bufsiz>
struct _Deque_iterator : public _Deque_iterator_base< _Tp, __bufsiz> {

  typedef random_access_iterator_tag iterator_category;
  typedef _Tp value_type;
  typedef typename _Traits::reference  reference;
  typedef typename _Traits::pointer    pointer;
  typedef size_t size_type;
  typedef ptrdiff_t difference_type;
  typedef value_type** _Map_pointer;

  typedef _Deque_iterator_base< _Tp, __bufsiz> _Base;
  typedef _Deque_iterator<_Tp, _Traits, __bufsiz> _Self;
  typedef _Deque_iterator<_Tp, _Nonconst_traits<_Tp>,__bufsiz> _Nonconst_self;
  typedef _Deque_iterator<_Tp, _Const_traits<_Tp>,__bufsiz> _Const_self;
# ifdef __STL_HAS_NAMESPACES
  __STL_USING_BASE_MEMBER _Base::_M_cur;
  __STL_USING_BASE_MEMBER _Base::_M_first;
  __STL_USING_BASE_MEMBER _Base::_M_last;
  __STL_USING_BASE_MEMBER _Base::_M_node;
# endif

# if defined ( __STL_DEBUG )
#  ifdef __STL_HAS_NAMESPACES
    __STL_USING_BASE_MEMBER _Base::_M_unsafe;
#  endif
  _Deque_iterator(const __owned_list* __root, value_type* __x, _Map_pointer __y):
    _Deque_iterator_base<value_type, __bufsiz>(__root,__x,__y) {}
# else /* __STL_DEBUG */
  _Deque_iterator(value_type* __x, _Map_pointer __y) :
    _Deque_iterator_base<value_type, __bufsiz>(__x,__y) {}
# endif
  _Deque_iterator() {}
  _Deque_iterator(const _Nonconst_self& __x) : 
    _Deque_iterator_base<value_type, __bufsiz>(__x) {}

  reference operator*() const { 
      __stl_debug_check(__Deq_dereferenceable((const void*)this, (value_type*)0));    
      return *_M_cur; 
  }

  __STL_DEFINE_ARROW_OPERATOR

  difference_type operator-(const _Self& __x) const { return _M_subtract(__x); }

  _Self& operator++() { _M_increment(); return *this; }
  _Self operator++(int)  {
    _Self __tmp = *this;
    ++*this;
    return __tmp;
  }

  _Self& operator--() { _M_decrement(); return *this; }
  _Self operator--(int) {
    _Self __tmp = *this;
    --*this;
    return __tmp;
  }

  _Self& operator+=(difference_type __n) { _M_advance(__n); return *this; }
  _Self operator+(difference_type __n) const
  {
    _Self __tmp = *this;
    return __tmp += __n;
  }

  _Self& operator-=(difference_type __n) { return *this += -__n; }
  _Self operator-(difference_type __n) const {
    _Self __tmp = *this;
    return __tmp -= __n;
  }

  reference operator[](difference_type __n) const { return *(*this + __n); }
};

template <class _Tp, class _Traits, class __bufsiz>
inline _Deque_iterator<_Tp, _Traits, __bufsiz>
operator+(ptrdiff_t __n, const _Deque_iterator<_Tp, _Traits, __bufsiz>& __x)
{
   return __x + __n;
}


#ifdef __STL_USE_SEPARATE_RELOPS_NAMESPACE

template <class _Tp, class __bufsiz>
inline bool 
operator==(const _Deque_iterator_base<_Tp,__bufsiz>& __x,
	   const _Deque_iterator_base<_Tp,__bufsiz>& __y) { 
    __stl_debug_check(__check_same_owner_or_null(__x, __y));    
    return __x._M_cur == __y._M_cur; 
}

template <class _Tp, class __bufsiz>
inline bool 
operator < (const _Deque_iterator_base<_Tp,__bufsiz>& __x,
	    const _Deque_iterator_base<_Tp,__bufsiz>& __y) { 
  __stl_debug_check(__check_same_owner(__x, __y));    
  return (__x._M_node == __y._M_node) ? 
    (__x._M_cur < __y._M_cur) : (__x._M_node < __y._M_node);
}

template <class _Tp, class __bufsiz>
inline bool 
operator!=(const _Deque_iterator_base<_Tp,__bufsiz>& __x,
	   const _Deque_iterator_base<_Tp,__bufsiz>& __y) { 
    __stl_debug_check(__check_same_owner_or_null(__x, __y));    
    return __x._M_cur != __y._M_cur; 
}
template <class _Tp,  class __bufsiz>
inline bool 
operator>(const _Deque_iterator_base<_Tp,__bufsiz>& __x,
	  const _Deque_iterator_base<_Tp,__bufsiz>& __y) { 
    return __y < __x;
}
template <class _Tp, class __bufsiz>
inline bool operator>=(const _Deque_iterator_base<_Tp,__bufsiz>& __x,
		       const _Deque_iterator_base<_Tp,__bufsiz>& __y) { 
    return !(__x < __y);
}
template <class _Tp, class __bufsiz>
inline bool operator<=(const _Deque_iterator_base<_Tp,__bufsiz>& __x,
		       const _Deque_iterator_base<_Tp,__bufsiz>& __y) { 
    return !(__y < __x);
}
# else

template <class _Tp, class _Traits1, class _Traits2, class __bufsiz>
inline bool 
operator==(const _Deque_iterator<_Tp, _Traits1,__bufsiz>& __x,
	   const _Deque_iterator<_Tp, _Traits2,__bufsiz>& __y) { 
    __stl_debug_check(__check_same_owner_or_null(__x, __y));    
    return __x._M_cur == __y._M_cur; 
}

template <class _Tp, class _Traits1, class _Traits2, class __bufsiz>
inline bool 
operator < (const _Deque_iterator<_Tp, _Traits1,__bufsiz>& __x,
	    const _Deque_iterator<_Tp, _Traits2,__bufsiz>& __y) { 
  __stl_debug_check(__check_same_owner(__x, __y));    
  return (__x._M_node == __y._M_node) ? 
    (__x._M_cur < __y._M_cur) : (__x._M_node < __y._M_node);
}

template <class _Tp, class __bufsiz>
inline bool 
operator!=(const _Deque_iterator<_Tp, _Nonconst_traits<_Tp>,__bufsiz>& __x,
	   const _Deque_iterator<_Tp, _Const_traits<_Tp>,__bufsiz>& __y) { 
    __stl_debug_check(__check_same_owner_or_null(__x, __y));    
    return __x._M_cur != __y._M_cur; 
}
template <class _Tp,  class __bufsiz>
inline bool 
operator>(const _Deque_iterator<_Tp, _Nonconst_traits<_Tp>,__bufsiz>& __x,
	  const _Deque_iterator<_Tp, _Const_traits<_Tp>, __bufsiz>& __y) { 
    return __y < __x;
}
template <class _Tp, class __bufsiz>
inline bool operator>=(const _Deque_iterator<_Tp, _Nonconst_traits<_Tp>,__bufsiz>& __x,
		       const _Deque_iterator<_Tp, _Const_traits<_Tp>, __bufsiz>& __y) { 
    return !(__x < __y);
}
template <class _Tp, class __bufsiz>
inline bool operator<=(const _Deque_iterator<_Tp, _Nonconst_traits<_Tp>,__bufsiz>& __x,
		       const _Deque_iterator<_Tp, _Const_traits<_Tp>, __bufsiz>& __y) { 
    return !(__y < __x);
}
# endif

#ifndef __STL_CLASS_PARTIAL_SPECIALIZATION

template <class _Tp, class _Traits, class __bufsiz>
inline _Tp*
value_type(const _Deque_iterator<_Tp, _Traits,__bufsiz>&) {
  return (_Tp*)0;
}

template <class _Tp, class _Traits, class __bufsiz>
inline random_access_iterator_tag
iterator_category(const _Deque_iterator<_Tp, _Traits,__bufsiz>&) {
  return random_access_iterator_tag();
}

template <class _Tp, class _Traits, class __bufsiz>
inline ptrdiff_t*
distance_type(const _Deque_iterator<_Tp, _Traits,__bufsiz>&) {
  return 0;
}

#endif /* __STL_CLASS_PARTIAL_SPECIALIZATION */


// Deque base class.  It has two purposes.  First, its constructor
//  and destructor allocate (but don't initialize) storage.  This makes
//  exception safety easier.  Second, the base class encapsulates all of
//  the differences between SGI-style allocators and standard-conforming
//  allocators.

template <class _Tp, class _Alloc, size_t __bufsiz>
class _Deque_base {
public:
  typedef _Tp value_type;
  typedef typename _Alloc_traits<_Tp,_Alloc>::allocator_type  allocator_type;
  typedef typename _Alloc_traits<_Tp*, _Alloc>::allocator_type _Map_alloc_type;

  typedef _Buf_size_traits <_Tp, __bufsiz > __buf_traits;

  typedef _Deque_iterator<_Tp, _Nonconst_traits<_Tp>,__buf_traits> iterator;
  typedef _Deque_iterator<_Tp, _Const_traits<_Tp>, __buf_traits>   const_iterator;

  _Deque_base(const allocator_type& __a, size_t __num_elements)
    : _M_start(), _M_finish(), _M_map(__STL_CONVERT_ALLOCATOR(__a, _Tp*), (_Tp**)0),
      _M_map_size(__a, (size_t)0) {
	_M_initialize_map(__num_elements);
  }
  _Deque_base(const allocator_type& __a)
    : _M_start(), _M_finish(), _M_map(__STL_CONVERT_ALLOCATOR(__a, _Tp*), (_Tp**)0), 
      _M_map_size(__a, (size_t)0) {
  }
  ~_Deque_base();    
  allocator_type get_allocator() const { return _M_map_size; }

protected:
  void _M_initialize_map(size_t);
  void _M_create_nodes(_Tp** __nstart, _Tp** __nfinish);
  void _M_destroy_nodes(_Tp** __nstart, _Tp** __nfinish);
  enum { _S_initial_map_size = 8 };

protected:
  iterator _M_start;
  iterator _M_finish;
  _STL_alloc_proxy<value_type**, value_type*, _Map_alloc_type>  _M_map;
  _STL_alloc_proxy<size_t, value_type,  allocator_type>   _M_map_size;  
};


// See __deque_buf_size().  The only reason that the default value is 0
//  is as a workaround for bugs in the way that some compilers handle
//  constant expressions.
# if defined ( __STL_NO_DEFAULT_NON_TYPE_PARAM )
template <class _Tp, class _Alloc, size_t __bufsiz> 
# else
template <class _Tp, __STL_DEFAULT_ALLOCATOR_SELECT(_Tp),
  __DFL_NON_TYPE_PARAM(size_t, __bufsiz, 0)>
# endif 
class deque : protected _Deque_base<_Tp, _Alloc, __bufsiz> {
  typedef _Deque_base<_Tp, _Alloc, __bufsiz> _Base;
  typedef deque<_Tp, _Alloc, __bufsiz> _Self;
public:                         // Basic types
  typedef _Tp value_type;
  typedef value_type* pointer;
  typedef const value_type* const_pointer;
  typedef value_type& reference;
  typedef const value_type& const_reference;
  typedef size_t size_type;
  typedef ptrdiff_t difference_type;

  typedef typename _Base::allocator_type allocator_type;

public:                         // Iterators
  typedef typename _Base::iterator       iterator;
  typedef typename _Base::const_iterator const_iterator;

#if defined ( __STL_CLASS_PARTIAL_SPECIALIZATION ) && \
! defined (__STL_PARTIAL_SPECIALIZATION_BUG) && \
! defined (CGAL_LIMITED_ITERATOR_TRAITS_SUPPORT)
  typedef __STLPORT_STD::reverse_iterator<const_iterator> const_reverse_iterator;
  typedef __STLPORT_STD::reverse_iterator<iterator> reverse_iterator;
#else /* __STL_CLASS_PARTIAL_SPECIALIZATION */
# if defined (__STL_MSVC50_COMPATIBILITY)
    typedef __STLPORT_STD::reverse_iterator<const_iterator, value_type, const_reference, 
    const value_type*, difference_type>  const_reverse_iterator;
    typedef __STLPORT_STD::reverse_iterator<iterator, value_type, reference, pointer, 
        difference_type> reverse_iterator; 
# else
  typedef __STLPORT_STD::reverse_iterator<const_iterator, value_type, const_reference, 
                           difference_type>  
          const_reverse_iterator;
  typedef __STLPORT_STD::reverse_iterator<iterator, value_type, reference, difference_type>
          reverse_iterator; 
# endif /* __STL_MSVC50_COMPATIBILITY */
#endif /* __STL_CLASS_PARTIAL_SPECIALIZATION */

protected:                      // Internal typedefs
  typedef pointer* _Map_pointer;

  typedef _Buf_size_traits < _Tp, __bufsiz > __buf_traits;

protected:
#if defined( __STL_HAS_NAMESPACES )
  __STL_USING_BASE_MEMBER _Deque_base<_Tp, _Alloc, __bufsiz>::_M_initialize_map;
  __STL_USING_BASE_MEMBER _Deque_base<_Tp, _Alloc, __bufsiz>::_M_create_nodes;
# ifndef __STL_DEBUG
  __STL_USING_BASE_MEMBER _Deque_base<_Tp, _Alloc, __bufsiz>::_M_destroy_nodes;
# endif
  __STL_USING_BASE_MEMBER _Deque_base<_Tp, _Alloc, __bufsiz>::_M_map;
  __STL_USING_BASE_MEMBER _Deque_base<_Tp, _Alloc, __bufsiz>::_M_map_size;
  __STL_USING_BASE_MEMBER _Deque_base<_Tp, _Alloc, __bufsiz>::_M_start;
  __STL_USING_BASE_MEMBER _Deque_base<_Tp, _Alloc, __bufsiz>::_M_finish;
#endif /* __STL_HAS_NAMESPACES */
# if defined (__STL_DEBUG)
protected:
  __owned_list _M_iter_list;
  void _Init_bounds() {
    _M_orphan_start();
    _M_orphan_finish();
  }
  void _Invalidate_iterator(const iterator& __it) { 
    __invalidate_iterator(&_M_iter_list,__it);
  }
  void _Invalidate_all() {
      _M_iter_list._Invalidate_all();
  }
  void _M_destroy_nodes(_Tp** __nstart, _Tp** __nfinish) {
    _Base::_M_destroy_nodes(__nstart, __nfinish);
    // to simplify unwind handling
      _M_finish._M_unsafe= false;
      _M_start._M_unsafe = false;
  }
  void _M_orphan_start() {
    _M_iter_list._Orphan(_M_start);
    _M_start._M_unsafe = false;
  }
  void _M_orphan_finish() {
    _M_iter_list._Orphan(_M_finish);
    _M_finish._M_unsafe= false;
  }
# endif
public:                         // Basic accessors
  iterator begin() { return _M_start; }
  iterator end() { return _M_finish; }
  const_iterator begin() const { return *(const const_iterator*)&_M_start; }
  const_iterator end() const { return *(const const_iterator*)&_M_finish; }

  reverse_iterator rbegin() { return reverse_iterator(_M_finish); }
  reverse_iterator rend() { return reverse_iterator(_M_start); }
  const_reverse_iterator rbegin() const 
    { return const_reverse_iterator(_M_finish); }
  const_reverse_iterator rend() const 
    { return const_reverse_iterator(_M_start); }

  reference operator[](size_type __n)
    { return _M_start[difference_type(__n)]; }
  const_reference operator[](size_type __n) const 
    { return _M_start[difference_type(__n)]; }

#ifdef __STL_THROW_RANGE_ERRORS
  void _M_range_check(size_type __n) const {
    if (__n >= this->size())
      __stl_throw_range_error("deque");
  }

  reference at(size_type __n)
    { _M_range_check(__n); return (*this)[__n]; }
  const_reference at(size_type __n) const
    { _M_range_check(__n); return (*this)[__n]; }
#endif /* __STL_THROW_RANGE_ERRORS */

  reference front() { return *_M_start; }
  reference back() {
    iterator __tmp = _M_finish;
    --__tmp;
    return *__tmp;
  }
  const_reference front() const { return *_M_start; }
  const_reference back() const {
    const_iterator __tmp = _M_finish;
    --__tmp;
    return *__tmp;
  }

  size_type size() const { return _M_finish - _M_start; }
  size_type max_size() const { return size_type(-1); }
  bool empty() const { return _M_finish == _M_start; }

public:                         // Constructor, destructor.
  explicit deque(const allocator_type& __a = __STL_ALLOC_INSTANCE(allocator_type)) 
    : _Deque_base<_Tp, _Alloc, __bufsiz>(__a, 0) {
    __stl_debug_do(_M_iter_list._Safe_init(&_M_start));
    __stl_debug_do(_Init_bounds());
}
  deque(const _Self& __x) : 
    _Deque_base<_Tp, _Alloc, __bufsiz>(__x.get_allocator(), __x.size()) 
    { 
      __stl_debug_do(_M_iter_list._Safe_init(&_M_start));
      __stl_debug_do(_Init_bounds());
      uninitialized_copy(__x.begin(), __x.end(), _M_start); 
    }

  deque(size_type __n, const value_type& __value,
        const allocator_type& __a = __STL_ALLOC_INSTANCE(allocator_type)) : 
    _Deque_base<_Tp, _Alloc, __bufsiz>(__a, __n)
    { _M_fill_initialize(__value); }
  // int,long variants may be needed 
  explicit deque(size_type __n) : _Deque_base<_Tp, _Alloc, __bufsiz>(allocator_type(), __n)
    { _M_fill_initialize(value_type()); }

#ifdef __STL_MEMBER_TEMPLATES

  // Check whether it's an integral type.  If so, it's not an iterator.
  template <class _InputIterator>
  deque(_InputIterator __first, _InputIterator __last,
        const allocator_type& __a = __STL_ALLOC_INSTANCE(allocator_type)) : 
    _Deque_base<_Tp, _Alloc, __bufsiz>(__a) {
    typedef typename _Is_integer<_InputIterator>::_Integral _Integral;
    _M_initialize_dispatch(__first, __last, _Integral());
  }

  template <class _Integer>
  void _M_initialize_dispatch(_Integer __n, _Integer __x, __true_type) {
    _M_initialize_map(__n);
    _M_fill_initialize(__x);
  }

  template <class _InputIter>
  void _M_initialize_dispatch(_InputIter __first, _InputIter __last,
                              __false_type) {
    _M_range_initialize(__first, __last, __ITERATOR_CATEGORY(__first));
  }

#else /* __STL_MEMBER_TEMPLATES */

  deque(const value_type* __first, const value_type* __last,
        const allocator_type& __a = __STL_ALLOC_INSTANCE(allocator_type)) 
    : _Deque_base<_Tp, _Alloc, __bufsiz>(__a, __last - __first)
    { 
      __stl_debug_do(_M_iter_list._Safe_init(&_M_start));
      __stl_debug_do(_Init_bounds());
      uninitialized_copy(__first, __last, _M_start); 
    }
  deque(const_iterator __first, const_iterator __last,
        const allocator_type& __a = __STL_ALLOC_INSTANCE(allocator_type)) 
    : _Deque_base<_Tp, _Alloc, __bufsiz>(__a, __last - __first)
    { 
      __stl_debug_do(_M_iter_list._Safe_init(&_M_start));
      __stl_debug_do(_Init_bounds());
      uninitialized_copy(__first, __last, _M_start); 
    }

#endif /* __STL_MEMBER_TEMPLATES */

  ~deque() { 
    destroy(_M_start, _M_finish); 
  }

  _Self& operator= (const _Self& __x);

  void swap(_Self& __x) {
    __STLPORT_STD::swap(_M_start, __x._M_start);
    __STLPORT_STD::swap(_M_finish, __x._M_finish);
    __STLPORT_STD::swap(_M_map, __x._M_map);
    __STLPORT_STD::swap(_M_map_size, __x._M_map_size);
    __stl_debug_do(_M_iter_list._Swap_owners(__x._M_iter_list));
  }

public: 
  // assign(), a generalized assignment member function.  Two
  // versions: one that takes a count, and one that takes a range.
  // The range version is a member template, so we dispatch on whether
  // or not the type is an integer.

  void _M_fill_assign(size_type __n, const _Tp& __val) {
    if (__n > size()) {
      fill(begin(), end(), __val);
      insert(end(), __n - size(), __val);
    }
    else {
      erase(begin() + __n, end());
      fill(begin(), end(), __val);
    }
  }

  void assign(size_type __n, const _Tp& __val) {
    _M_fill_assign(__n, __val);
  }

#ifdef __STL_MEMBER_TEMPLATES

  template <class _InputIterator>
  void assign(_InputIterator __first, _InputIterator __last) {
    typedef typename _Is_integer<_InputIterator>::_Integral _Integral;
    _M_assign_dispatch(__first, __last, _Integral());
  }

private:                        // helper functions for assign() 

  template <class _Integer>
  void _M_assign_dispatch(_Integer __n, _Integer __val, __true_type)
    { _M_fill_assign((size_type) __n, (_Tp) __val); }

  template <class _InputIterator>
  void _M_assign_dispatch(_InputIterator __first, _InputIterator __last,
                          __false_type) {
    _M_assign_aux(__first, __last, __ITERATOR_CATEGORY(__first));
  }

  template <class _InputIter>
  void _M_assign_aux(_InputIter __first, _InputIter __last, input_iterator_tag) {
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
    if (__len > size()) {
      _ForwardIterator __mid = __first;
      advance(__mid, size());
      copy(__first, __mid, begin());
      insert(end(), __mid, __last);
    }
    else
      erase(copy(__first, __last, begin()), end());
  }

#endif /* __STL_MEMBER_TEMPLATES */

public:                         // push_* and pop_*
  
  void push_back(const value_type& __t) {
    if (_M_finish._M_cur != _M_finish._M_last - 1) {
      construct(_M_finish._M_cur, __t);
      ++_M_finish._M_cur;
    }
    else
      _M_push_back_aux(__t);
    __stl_debug_do(_Invalidate_all());
  }

  void push_back() {
    if (_M_finish._M_cur != _M_finish._M_last - 1) {
      construct(_M_finish._M_cur);
      ++_M_finish._M_cur;
    }
    else
      _M_push_back_aux();
    __stl_debug_do(_Invalidate_all());
  }

  void push_front(const value_type& __t) {
    if (_M_start._M_cur != _M_start._M_first) {
      construct(_M_start._M_cur - 1, __t);
      --_M_start._M_cur;
    }
    else
      _M_push_front_aux(__t);
    __stl_debug_do(_Invalidate_all());
  }

  void push_front() {
    if (_M_start._M_cur != _M_start._M_first) {
      construct(_M_start._M_cur - 1);
      --_M_start._M_cur;
    }
    else
      _M_push_front_aux();
    __stl_debug_do(_Invalidate_all());
  }


  void pop_back() {
    __stl_debug_do(_Invalidate_iterator(_M_finish));
    if (_M_finish._M_cur != _M_finish._M_first) {
      --_M_finish._M_cur;
      destroy(_M_finish._M_cur);
    }
    else
      _M_pop_back_aux();
  }

  void pop_front() {
    __stl_debug_do(_Invalidate_iterator(_M_start));        
    if (_M_start._M_cur != _M_start._M_last - 1) {
      destroy(_M_start._M_cur);
      ++_M_start._M_cur;
    }
    else 
      _M_pop_front_aux();
  }

public:                         // Insert

  iterator insert(iterator __position, const value_type& __x) {
    __stl_debug_check(__check_if_owner(&_M_iter_list, __position));
    if (__position._M_cur == _M_start._M_cur) {
      push_front(__x);
      return _M_start;
    }
    else if (__position._M_cur == _M_finish._M_cur) {
      push_back(__x);
      iterator __tmp = _M_finish;
      --__tmp;
      return __tmp;
    }
    else {
      return _M_insert_aux(__position, __x);
    }
  }

  iterator insert(iterator __position)
    { return insert(__position, value_type()); }

  void insert(iterator __pos, size_type __n, const value_type& __x) {
    _M_fill_insert(__pos, __n, __x);
  }

  void _M_fill_insert(iterator __pos, size_type __n, const value_type& __x);

#ifdef __STL_MEMBER_TEMPLATES  

  // Check whether it's an integral type.  If so, it's not an iterator.
  template <class _InputIterator>
  void insert(iterator __pos, _InputIterator __first, _InputIterator __last) {
    typedef typename _Is_integer<_InputIterator>::_Integral _Integral;
    _M_insert_dispatch(__pos, __first, __last, _Integral());
  }

  template <class _Integer>
  void _M_insert_dispatch(iterator __pos, _Integer __n, _Integer __x,
                          __true_type) {
    _M_fill_insert(__pos, (size_type) __n, (value_type) __x);
  }

  template <class _InputIterator>
  void _M_insert_dispatch(iterator __pos,
                          _InputIterator __first, _InputIterator __last,
                          __false_type) {
    insert(__pos, __first, __last, __ITERATOR_CATEGORY(__first));
  }

#else /* __STL_MEMBER_TEMPLATES */

  void insert(iterator __pos,
              const value_type* __first, const value_type* __last);
  void insert(iterator __pos,
              const_iterator __first, const_iterator __last);

#endif /* __STL_MEMBER_TEMPLATES */

  void resize(size_type __new_size, const value_type& __x) {
    const size_type __len = size();
    if (__new_size < __len) 
      erase(_M_start + __new_size, _M_finish);
    else
      insert(_M_finish, __new_size - __len, __x);
  }

  void resize(size_type new_size) { resize(new_size, value_type()); }

public:                         // Erase
  iterator erase(iterator __pos) {
    __stl_debug_check(__check_if_owner(&_M_iter_list, __pos) && (__pos != end()));
    iterator __next = __pos;
    ++__next;
    difference_type __index = __pos - _M_start;
    if (__index < difference_type(size() >> 1)) {
      copy_backward(_M_start, __pos, __next);
      pop_front();
    }
    else {
      copy(__next, _M_finish, __pos);
      pop_back();
    }
    return _M_start + __index;
  }

  iterator erase(iterator __first, iterator __last);
  void clear(); 

protected:                        // Internal construction/destruction

  void _M_fill_initialize(const value_type& __value);

#ifdef __STL_MEMBER_TEMPLATES 

  template <class _InputIterator>
  void _M_range_initialize(_InputIterator __first,
			   _InputIterator __last,
			   input_iterator_tag)
#ifndef __STL_INLINE_MEMBER_TEMPLATES
; 
#else
 {
    _M_initialize_map(0);
    __stl_debug_do(_M_iter_list._Safe_init(&_M_start));
    __stl_debug_do(_Init_bounds());
    __STL_TRY {
      for ( ; __first != __last; ++__first)
	push_back(*__first);
    }
    __STL_UNWIND(clear());
  }
# endif /* __STL_OUTLINE_MEMBER_TEMPLATES */
  
 template <class _ForwardIterator>
 void  _M_range_initialize(_ForwardIterator __first,
			   _ForwardIterator __last,
			   forward_iterator_tag) 
# ifndef __STL_INLINE_MEMBER_TEMPLATES
;
# else
{
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
# endif /* __STL_INLINE_MEMBER_TEMPLATES */
#endif /* __STL_MEMBER_TEMPLATES */

protected:                        // Internal push_* and pop_*

  void _M_push_back_aux(const value_type&);
  void _M_push_back_aux();
  void _M_push_front_aux(const value_type&);
  void _M_push_front_aux();
  void _M_pop_back_aux();
  void _M_pop_front_aux();

protected:                        // Internal insert functions

#ifdef __STL_MEMBER_TEMPLATES

template <class _InputIterator>
void 
insert(iterator __pos,
       _InputIterator __first,
       _InputIterator __last,
       input_iterator_tag)
{
  __stl_debug_check(__check_if_owner(&_M_iter_list, __pos) &&
                    __check_range(__first,__last));    
  copy(__first, __last, inserter(*this, __pos));
}

template <class _ForwardIterator>
void  insert(iterator __pos,
	     _ForwardIterator __first,
	     _ForwardIterator __last,
	     forward_iterator_tag)
# ifndef __STL_INLINE_MEMBER_TEMPLATES
;
# else
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
# endif /* __STL_INLINE_MEMBER_TEMPLATES */
#endif /* __STL_MEMBER_TEMPLATES */

  iterator _M_insert_aux(iterator __pos, const value_type& __x);
  iterator _M_insert_aux(iterator __pos);
  iterator _M_insert_aux_prepare(iterator __pos);

  void _M_insert_aux(iterator __pos, size_type __n, const value_type& __x);

#ifdef __STL_MEMBER_TEMPLATES  
  template <class _ForwardIterator>
  void _M_insert_aux(iterator __pos,
					    _ForwardIterator __first,
					    _ForwardIterator __last,
					    size_type __n)
# ifndef __STL_INLINE_MEMBER_TEMPLATES
;
# else
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
# endif /* __STL_INLINE_MEMBER_TEMPLATES */
#else /* __STL_MEMBER_TEMPLATES */
  
  void _M_insert_aux(iterator __pos,
                     const value_type* __first, const value_type* __last,
                     size_type __n);

  void _M_insert_aux(iterator __pos, 
                     const_iterator __first, const_iterator __last,
                     size_type __n);
 
#endif /* __STL_MEMBER_TEMPLATES */

  iterator _M_reserve_elements_at_front(size_type __n) {
    size_type __vacancies = _M_start._M_cur - _M_start._M_first;
    if (__n > __vacancies) 
      _M_new_elements_at_front(__n - __vacancies);
    __stl_debug_do(_M_start._M_unsafe=true);
    return _M_start - difference_type(__n);
  }

  iterator _M_reserve_elements_at_back(size_type __n) {
    size_type __vacancies = (_M_finish._M_last - _M_finish._M_cur) - 1;
    if (__n > __vacancies)
      _M_new_elements_at_back(__n - __vacancies);
    __stl_debug_do(_M_finish._M_unsafe=true);
    return _M_finish + difference_type(__n);
  }

  void _M_new_elements_at_front(size_type __new_elements);
  void _M_new_elements_at_back(size_type __new_elements);

protected:                      // Allocation of _M_map and nodes

  // Makes sure the _M_map has space for new nodes.  Does not actually
  //  add the nodes.  Can invalidate _M_map pointers.  (And consequently, 
  //  deque iterators.)

  void _M_reserve_map_at_back (size_type __nodes_to_add = 1) {
    if (__nodes_to_add + 1 > _M_map_size._M_data - (_M_finish._M_node - _M_map._M_data))
      _M_reallocate_map(__nodes_to_add, false);
  }

  void _M_reserve_map_at_front (size_type __nodes_to_add = 1) {
    if (__nodes_to_add > size_type(_M_start._M_node - _M_map._M_data))
      _M_reallocate_map(__nodes_to_add, true);
  }

  void _M_reallocate_map(size_type __nodes_to_add, bool __add_at_front);
 
#ifdef __STL_NON_TYPE_TMPL_PARAM_BUG
public:
  bool operator==(const deque<_Tp,_Alloc,__bufsiz>& __x) const {
    return size() == __x.size() && equal(begin(), end(), __x.begin());
  }
  bool operator!=(const deque<_Tp,_Alloc,__bufsiz>& __x) const {
    return size() != __x.size() || !equal(begin(), end(), __x.begin());
  }
  bool operator<(const deque<_Tp,_Alloc,__bufsiz>& __x) const {
    return lexicographical_compare(begin(), end(), __x.begin(), __x.end());
  }
  bool operator>(const deque<_Tp,_Alloc,__bufsiz>& __x) const {
    return __x < *this;
  }
  bool operator<=(const deque<_Tp,_Alloc,__bufsiz>& __x) const {
    return !(__x < *this);
  }
  bool operator>=(const deque<_Tp,_Alloc,__bufsiz>& __x) const {
    return !(*this < __x);
  }
# endif /* __STL_NON_TYPE_TMPL_PARAM_BUG */

};

// Nonmember functions.

#ifndef __STL_NON_TYPE_TMPL_PARAM_BUG

template <class _Tp, class _Alloc, size_t __bufsiz>
inline bool operator==(const deque<_Tp, _Alloc, __bufsiz>& __x,
		       const deque<_Tp, _Alloc, __bufsiz>& __y)
{
  return __x.size() == __y.size() &&
  equal(__x.begin(), __x.end(), __y.begin());
}

template <class _Tp, class _Alloc, size_t __bufsiz>
inline bool operator<(const deque<_Tp, _Alloc, __bufsiz>& __x,
		      const deque<_Tp, _Alloc, __bufsiz>& __y)
{
  return lexicographical_compare(__x.begin(), __x.end(), 
                                 __y.begin(), __y.end());
}

#if defined(__STL_USE_SEPARATE_RELOPS_NAMESPACE)

template <class _Tp, class _Alloc, size_t __bufsiz>
inline bool operator>(const deque<_Tp, _Alloc, __bufsiz>& __x,
		      const deque<_Tp, _Alloc, __bufsiz>& __y)
{
  return __y < __x; 
}

template <class _Tp, class _Alloc, size_t __bufsiz>
inline bool operator>=(const deque<_Tp, _Alloc, __bufsiz>& __x,
		       const deque<_Tp, _Alloc, __bufsiz>& __y)
{
  return !(__x < __y); 
}

template <class _Tp, class _Alloc, size_t __bufsiz>
inline bool operator<=(const deque<_Tp, _Alloc, __bufsiz>& __x,
		       const deque<_Tp, _Alloc, __bufsiz>& __y)
{
 return !(__y < __x); 
}
# endif /* __STL_FUNCTION_TMPL_PARTIAL_ORDER */
#endif /* __STL_NON_TYPE_TMPL_PARAM_BUG */


// do a cleanup
# undef deque

# if defined (__STL_USE_WRAPPER_FOR_ALLOC_PARAM) || \
     defined (__STL_NO_DEFAULT_NON_TYPE_PARAM)
// provide a "default" deque adaptor
# if (defined (__STL_DEFAULT_TYPE_PARAM) && defined (__STL_USE_SGI_ALLOCATORS)) || \
      ! defined (__STL_LIMITED_DEFAULT_TEMPLATES)
template <class _Tp, __STL_DEFAULT_ALLOCATOR_SELECT(_Tp) >
class deque : public __deque__<_Tp,_Alloc,size_t(0)> {
 public:
#   define _DEQUE_SUPER __deque__<_Tp,_Alloc,size_t(0)>
  typedef deque<_Tp, _Alloc> _Self;
# else
template <class _Tp>
class deque : public __deque__<_Tp, __STL_DEFAULT_ALLOCATOR(_Tp), size_t(0)> {
public:
#   define _DEQUE_SUPER __deque__<_Tp, __STL_DEFAULT_ALLOCATOR(_Tp), size_t(0)>
  typedef deque<_Tp> _Self;
# endif
  typedef _DEQUE_SUPER _Super;
    __IMPORT_WITH_REVERSE_ITERATORS(_Super)
    __IMPORT_SUPER_COPY_ASSIGNMENT(deque, _Self, _DEQUE_SUPER)
    deque() : _DEQUE_SUPER() { }
    deque(size_type __n, const _Tp& __value) : _DEQUE_SUPER(__n, __value) { }
    explicit deque(size_type __n) : _DEQUE_SUPER(__n) { }
    deque(const _Tp* __first, const _Tp* __last) : _DEQUE_SUPER(__first, __last) { }
    deque(const_iterator __first, const_iterator __last) : _DEQUE_SUPER(__first, __last) { }
    ~deque() { }
};

#  if defined (__STL_BASE_MATCH_BUG)
template <class _Tp>
inline bool 
operator==(const deque<_Tp>& __x, const deque<_Tp>& __y) {
    return __x.size() == __y.size() && equal(__x.begin(), __x.end(), __y.begin());
}

template <class _Tp>
inline bool 
operator<(const deque<_Tp>& __x, const deque<_Tp>& __y) {
    return lexicographical_compare(__x.begin(), __x.end(), __y.begin(), __y.end());
}
#  endif /* BASE_MATCH_BUG */

# undef _DEQUE_SUPER
# endif /* __STL_DEFAULT_TYPE_PARAM */

#if defined(__STL_FUNCTION_TMPL_PARTIAL_ORDER) && \
    !defined(__STL_NON_TYPE_TMPL_PARAM_BUG)

template <class _Tp, class _Alloc, size_t __bufsiz>
inline void 
swap(__deque__<_Tp,_Alloc,__bufsiz>& __x, __deque__<_Tp,_Alloc,__bufsiz>& __y)
{
  __x.swap(__y);
}

#endif

#if defined(__sgi) && !defined(__GNUC__) && (_MIPS_SIM != _MIPS_SIM_ABI32)
#pragma reset woff 1174
#pragma reset woff 1375
#endif
          
__STL_END_NAMESPACE 

# if !defined (__STL_LINK_TIME_INSTANTIATION)
#  include <stl_deque.c>
# endif
  
#endif /* __SGI_STL_INTERNAL_DEQUE_H */

// Local Variables:
// mode:C++
// End:
