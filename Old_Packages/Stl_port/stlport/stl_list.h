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

#ifndef __SGI_STL_INTERNAL_LIST_H
#define __SGI_STL_INTERNAL_LIST_H

# ifndef __SGI_STL_INTERNAL_ALGOBASE_H
#  include <stl_algobase.h>
# endif

# ifndef __SGI_STL_INTERNAL_ALLOC_H
#  include <stl_alloc.h>
# endif

# ifndef __SGI_STL_INTERNAL_ITERATOR_H
#  include <stl_iterator.h>
# endif

# ifndef __SGI_STL_INTERNAL_CONSTRUCT_H
#  include <stl_construct.h>
# endif

__STL_BEGIN_NAMESPACE

#if defined(__sgi) && !defined(__GNUC__) && (_MIPS_SIM != _MIPS_SIM_ABI32)
#pragma set woff 1174
#pragma set woff 1375
#endif

#  define  list  __WORKAROUND_RENAME(list)

# if defined ( __STL_USE_ABBREVS )
#  define __list_iterator         _L__It
# endif

struct _List_node_base {
  void* _M_next;
  void* _M_prev;
};

template <class _Dummy>
struct _List_global {
  typedef _List_node_base _Node;
  static void _Transfer(_List_node_base* __position, 
			_List_node_base* __first, _List_node_base* __last);
};

typedef _List_global<bool> _List_global_inst;

template <class _Tp>
struct _List_node : public _List_node_base {
  _Tp _M_data;
  __TRIVIAL_STUFF(_List_node)
};

template<class _Tp, class _Traits>
# if defined ( __STL_DEBUG )
struct _List_iterator : public __owned_link {
# else
struct _List_iterator {
# endif
  typedef _Tp value_type;
  typedef typename _Traits::pointer    pointer;
  typedef typename _Traits::reference  reference;

  typedef _List_iterator<_Tp, _Nonconst_traits<_Tp> > iterator;
  typedef _List_iterator<_Tp, _Const_traits<_Tp> >    const_iterator;
  typedef _List_iterator<_Tp, _Traits>                       _Self;

  typedef bidirectional_iterator_tag iterator_category;
  typedef _List_node<_Tp> _Node;
  typedef size_t size_type;
  typedef ptrdiff_t difference_type;

  _Node* _M_node;

  //  operator const const_iterator& () const { return *(const const_iterator*)this; } 

# if defined ( __STL_DEBUG )
  _Node* _Owner_node() const {
      const __owned_list* __ptr = _Owner();
      return __ptr ? (_Node*)__ptr->_Owner() : (_Node*)0; 
  }
  _List_iterator(const __owned_list* __root, _Node* __x) : 
    __owned_link(__root), _M_node(__x) {}
  _List_iterator() : __owned_link(0) {}
  _List_iterator(const iterator& __x) : __owned_link(__x), _M_node(__x._M_node) {}
# else
  _List_iterator(_Node* __x) : _M_node(__x) {}
  _List_iterator() {}
  _List_iterator(const iterator& __x) : _M_node(__x._M_node) {}
# endif

  reference operator*() const { 
            __stl_verbose_assert(_Valid() && _M_node!=_Owner_node(), 
				 _StlMsg_NOT_DEREFERENCEABLE); 
            return (*_M_node)._M_data; 
  }

  __STL_DEFINE_ARROW_OPERATOR

  _Self& operator++() { 
    __stl_verbose_assert(_M_node!=_Owner_node(), _StlMsg_INVALID_ADVANCE); 
    _M_node = (_Node*)(_M_node->_M_next);
    return *this;
  }
  _Self operator++(int) { 
    _Self __tmp = *this;
    ++*this;
    return __tmp;
  }
  _Self& operator--() { 
    _M_node = (_Node*)(_M_node->_M_prev);
    __stl_verbose_assert(_M_node!=_Owner_node(), _StlMsg_INVALID_ADVANCE); 
    return *this;
  }
  _Self operator--(int) { 
    _Self __tmp = *this;
    --*this;
    return __tmp;
  }
};


template<class _Tp, class _Traits, class _Traits1>
inline  bool operator==(const _List_iterator<_Tp, _Traits>& __x,
			const _List_iterator<_Tp, _Traits1>& __y ) { 
  __stl_debug_check(__check_same_owner_or_null(__x,__y));                         
  return __x._M_node == __y._M_node; 
}

#ifdef __STL_USE_SEPARATE_RELOPS_NAMESPACE
template<class _Tp, class _Traits, class _Traits1>
inline  bool operator!=(const _List_iterator<_Tp, _Traits>& __x,
			const _List_iterator<_Tp, _Traits1>& __y ) { 
    __stl_debug_check(__check_same_owner_or_null(__x, __y));                         
    return __x._M_node != __y._M_node; 
}
#else

template<class _Tp>
inline  bool operator!=(const _List_iterator<_Tp, _Nonconst_traits<_Tp> >& __x,
			const _List_iterator<_Tp, _Const_traits<_Tp> >& __y ) { 
    __stl_debug_check(__check_same_owner_or_null(__x, __y));                         
    return __x._M_node != __y._M_node; 
}

# ifdef __SUNPRO_CC
template<class _Tp>
inline  bool operator!=(const _List_iterator<_Tp, _Const_traits<_Tp> >& __x,
			const _List_iterator<_Tp, _Nonconst_traits<_Tp> >& __y ) { 
    __stl_debug_check(__check_same_owner_or_null(__x, __y));                         
    return __x._M_node != __y._M_node; 
}
# endif

#endif

#ifndef __STL_CLASS_PARTIAL_SPECIALIZATION

template <class _Tp, class _Traits>
inline _Tp*
value_type(const _List_iterator<_Tp, _Traits>&) { return 0; }

template <class _Tp, class _Traits>
inline bidirectional_iterator_tag
iterator_category(const _List_iterator<_Tp, _Traits>&) { return bidirectional_iterator_tag();}

template <class _Tp, class _Traits>
inline ptrdiff_t* distance_type(const _List_iterator<_Tp, _Traits>&) { return 0; }

#endif /* __STL_CLASS_PARTIAL_SPECIALIZATION */


// Base class that encapsulates details of allocators and helps 
// to simplify EH

template <class _Tp, class _Alloc>
class _List_base 
{
protected:
  typedef _List_node<_Tp> _Node;
  typedef typename _Alloc_traits<_Node, _Alloc>::allocator_type
           _Node_allocator_type;
public:
  typedef typename _Alloc_traits<_Tp, _Alloc>::allocator_type
          allocator_type;

  allocator_type get_allocator() const { 
    return __STL_CONVERT_ALLOCATOR((const _Node_allocator_type&)_M_node, _Tp);
  }

  _List_base(const allocator_type& __a) : _M_node(__STL_CONVERT_ALLOCATOR(__a, _Node), (_Node*)0) {
    _M_node._M_data = _M_node.allocate(1);
    _M_node._M_data->_M_next = _M_node._M_data;
    _M_node._M_data->_M_prev = _M_node._M_data;
    __stl_debug_do(_M_iter_list._Safe_init(_M_node._M_data));
  }
  ~_List_base() {
    clear();
    _M_node.deallocate(_M_node._M_data, 1);
  }

  void clear();

public:
  _STL_alloc_proxy<_Node*, _Node, _Node_allocator_type>  _M_node;
# if defined (__STL_DEBUG)
protected:
    mutable __owned_list _M_iter_list;
    void _Invalidate_all() { _M_iter_list._Invalidate_all();}
# endif
};

template <class _Tp, __STL_DEFAULT_ALLOCATOR_SELECT(_Tp) >
class list : protected _List_base<_Tp, _Alloc> {
  typedef _List_base<_Tp, _Alloc> _Base;
  typedef list<_Tp, _Alloc> _Self;
protected:
  typedef void* _Void_pointer;
public:      
  typedef _Tp value_type;
  typedef value_type* pointer;
  typedef const value_type* const_pointer;
  typedef value_type& reference;
  typedef const value_type& const_reference;
  typedef _List_node<_Tp> _Node;
  typedef size_t size_type;
  typedef ptrdiff_t difference_type;
  typedef typename _Base::allocator_type allocator_type;

public:
  typedef _List_iterator<_Tp, _Nonconst_traits<_Tp> > iterator;
  typedef _List_iterator<_Tp, _Const_traits<_Tp> >    const_iterator;

#if defined ( __STL_CLASS_PARTIAL_SPECIALIZATION ) && \
  ! defined (__STL_PARTIAL_SPECIALIZATION_BUG) && \
  ! defined(CGAL_LIMITED_ITERATOR_TRAITS_SUPPORT)
    typedef __STLPORT_STD::reverse_iterator<const_iterator> const_reverse_iterator;
    typedef __STLPORT_STD::reverse_iterator<iterator> reverse_iterator;
#else /* __STL_CLASS_PARTIAL_SPECIALIZATION */
# if defined (__STL_MSVC50_COMPATIBILITY)
    typedef reverse_bidirectional_iterator<const_iterator, value_type,
                                           const_reference, const value_type*, difference_type>
            const_reverse_iterator;
    typedef reverse_bidirectional_iterator<iterator, value_type, reference,
                                           pointer, difference_type>
            reverse_iterator; 
# else
  typedef reverse_bidirectional_iterator<const_iterator,value_type,
                                         const_reference,difference_type>
          const_reverse_iterator;
  typedef reverse_bidirectional_iterator<iterator,value_type,reference,
                                         difference_type>
          reverse_iterator; 
# endif /* __STL_MSVC50_COMPATIBILITY */
#endif /* __STL_CLASS_PARTIAL_SPECIALIZATION */

  // protected:
#if defined( __STL_HAS_NAMESPACES )
  __STL_USING_BASE_MEMBER _List_base<_Tp, _Alloc>::_M_node;
#endif /* __STL_HAS_NAMESPACES */
public:
  __STL_USING_BASE_MEMBER _List_base<_Tp, _Alloc>::get_allocator;
  __STL_USING_BASE_MEMBER _List_base<_Tp, _Alloc>::clear;

protected:
  _Node* _M_create_node(const _Tp& __x)
  {
    _Node* __p = _M_node.allocate(1);
    __STL_TRY {
      construct(&__p->_M_data, __x);
    }
    __STL_UNWIND(_M_node.deallocate(__p, 1));
    return __p;
  }

  _Node* _M_create_node()
  {
    _Node* __p = _M_node.allocate(1);
    __STL_TRY {
      construct(&__p->_M_data);
    }
    __STL_UNWIND(_M_node.deallocate(__p, 1));
    return __p;
  }

public:
  explicit list(const allocator_type& __a = __STL_ALLOC_INSTANCE(allocator_type)) :
    _List_base<_Tp, _Alloc>(__a) {}

# if defined (__STL_DEBUG)
  iterator begin()             { return iterator(&_M_iter_list, (_Node*)(_M_node._M_data->_M_next)); }
  const_iterator begin() const { return const_iterator(&_M_iter_list, (_Node*)(_M_node._M_data->_M_next)); }

  iterator end()             { return iterator(&_M_iter_list, _M_node._M_data); }
  const_iterator end() const { return const_iterator(&_M_iter_list, _M_node._M_data); }
# else
  iterator begin()             { return iterator((_Node*)(_M_node._M_data->_M_next)); }
  const_iterator begin() const { return const_iterator((_Node*)(_M_node._M_data->_M_next)); }

  iterator end()             { return _M_node._M_data; }
  const_iterator end() const { return _M_node._M_data; }
# endif

  reverse_iterator rbegin() 
    { return reverse_iterator(end()); }
  const_reverse_iterator rbegin() const 
    { return const_reverse_iterator(end()); }

  reverse_iterator rend()
    { return reverse_iterator(begin()); }
  const_reverse_iterator rend() const
    { return const_reverse_iterator(begin()); }

  bool empty() const { return _M_node._M_data->_M_next == _M_node._M_data; }
  size_type size() const {
    size_type __result = 0;
    distance(begin(), end(), __result);
    return __result;
  }
  size_type max_size() const { return size_type(-1); }

  reference front() { return *begin(); }
  const_reference front() const { return *begin(); }
  reference back() { return *(--end()); }
  const_reference back() const { return *(--end()); }

  void swap(list<_Tp, _Alloc>& __x) {
    __stl_debug_do(_M_iter_list._Swap_owners(__x._M_iter_list, true));
    __STLPORT_STD::swap(_M_node, __x._M_node); 
  }

  iterator insert(iterator __position, const _Tp& __x) {
    __stl_debug_check(__check_if_owner(&_M_iter_list,__position));
    _Node* __tmp = _M_create_node(__x);
    __tmp->_M_next = __position._M_node;
    __tmp->_M_prev = __position._M_node->_M_prev;
    ((_Node*) (__position._M_node->_M_prev))->_M_next = __tmp;
    __position._M_node->_M_prev = __tmp;
#  if defined ( __STL_DEBUG )
      return iterator(&_M_iter_list,__tmp);
#  else
    return __tmp;
#  endif
  }

  iterator insert(iterator __position) { return insert(__position, _Tp()); }
#ifdef __STL_MEMBER_TEMPLATES
  // Check whether it's an integral type.  If so, it's not an iterator.

  template<class _Integer>
  void _M_insert_dispatch(iterator __pos, _Integer __n, _Integer __x,
                          __true_type) {
    _M_fill_insert(__pos, (size_type) __n, (_Tp) __x);
  }

  template <class _InputIter>
  void 
  _M_insert_dispatch(iterator __position,
		     _InputIter __first, _InputIter __last,
		     __false_type) {
    for ( ; __first != __last; ++__first)
      insert(__position, *__first);
  }

  template <class _InputIterator>
  void insert(iterator __pos, _InputIterator __first, _InputIterator __last) {
    typedef typename _Is_integer<_InputIterator>::_Integral _Integral;
    _M_insert_dispatch(__pos, __first, __last, _Integral());
  }

#else /* __STL_MEMBER_TEMPLATES */
  void insert(iterator __position, const _Tp* __first, const _Tp* __last);
  void insert(iterator __position,
              const_iterator __first, const_iterator __last);
#endif /* __STL_MEMBER_TEMPLATES */
  void insert(iterator __pos, size_type __n, const _Tp& __x)
  { _M_fill_insert(__pos, __n, __x); }
  void _M_fill_insert(iterator __pos, size_type __n, const _Tp& __x);
 
  void push_front(const _Tp& __x) { insert(begin(), __x); }
  void push_front() {insert(begin());}
  void push_back(const _Tp& __x) { insert(end(), __x); }
  void push_back() {insert(end());}

  iterator erase(iterator __position) {
    __stl_debug_check(__check_if_owner(&_M_iter_list,__position));
    __stl_verbose_assert(__position._M_node!=_M_node._M_data, _StlMsg_ERASE_PAST_THE_END);
    _Node* __next_node = (_Node*) (__position._M_node->_M_next);
    _Node* __prev_node = (_Node*) (__position._M_node->_M_prev);
    __prev_node->_M_next = __next_node;
    __next_node->_M_prev = __prev_node;
    __stl_debug_do(__invalidate_iterator(&_M_iter_list, __position));
    destroy(&__position._M_node->_M_data);
    _M_node.deallocate(__position._M_node, 1);
#  if defined ( __STL_DEBUG )
      return iterator(&_M_iter_list,__next_node);
#  else
      return iterator(__next_node);
#  endif
    }
  iterator erase(iterator __first, iterator __last);
  //  void clear() { _Base::clear(); }

  void resize(size_type __new_size, const _Tp& __x);
  void resize(size_type __new_size) { resize(__new_size, _Tp()); }

  void pop_front() { erase(begin()); }
  void pop_back() { 
    iterator __tmp = end();
    erase(--__tmp);
  }
  list(size_type __n, const _Tp& __value,
       const allocator_type& __a = __STL_ALLOC_INSTANCE(allocator_type))
    : _List_base<_Tp, _Alloc>(__a)
    { insert(begin(), __n, __value); }
  explicit list(size_type __n)
    : _List_base<_Tp, _Alloc>(allocator_type())
    { insert(begin(), __n, _Tp()); }

#ifdef __STL_MEMBER_TEMPLATES

  // We don't need any dispatching tricks here, because insert does all of
  // that anyway.  
  template <class _InputIterator>
  list(_InputIterator __first, _InputIterator __last,
       const allocator_type& __a = __STL_ALLOC_INSTANCE(allocator_type))
    : _List_base<_Tp, _Alloc>(__a)
    { insert(begin(), __first, __last); }

#else /* __STL_MEMBER_TEMPLATES */

  list(const _Tp* __first, const _Tp* __last,
       const allocator_type& __a = __STL_ALLOC_INSTANCE(allocator_type))
    : _List_base<_Tp, _Alloc>(__a)
    { insert(begin(), __first, __last); }
  list(const_iterator __first, const_iterator __last,
       const allocator_type& __a)
    : _List_base<_Tp, _Alloc>(__a)
    { insert(begin(), __first, __last); }
  list(const_iterator __first, const_iterator __last)
    : _List_base<_Tp, _Alloc>(allocator_type())
    { insert(begin(), __first, __last); }

#endif /* __STL_MEMBER_TEMPLATES */
  list(const list<_Tp, _Alloc>& __x) : _List_base<_Tp, _Alloc>(__x.get_allocator())
    { insert(begin(), __x.begin(), __x.end()); }

  ~list() { }

  list<_Tp, _Alloc>& operator=(const list<_Tp, _Alloc>& __x);

public:
  // assign(), a generalized assignment member function.  Two
  // versions: one that takes a count, and one that takes a range.
  // The range version is a member template, so we dispatch on whether
  // or not the type is an integer.

  void assign(size_type __n, const _Tp& __val) { _M_fill_assign(__n, __val); }

  void _M_fill_assign(size_type __n, const _Tp& __val);

#ifdef __STL_MEMBER_TEMPLATES

  template <class _InputIterator>
  void assign(_InputIterator __first, _InputIterator __last) {
    typedef typename _Is_integer<_InputIterator>::_Integral _Integral;
    _M_assign_dispatch(__first, __last, _Integral());
  }

  template <class _Integer>
  void _M_assign_dispatch(_Integer __n, _Integer __val, __true_type)
    { assign((size_type) __n, (_Tp) __val); }

  template <class _InputIterator>
  void _M_assign_dispatch(_InputIterator __first2, _InputIterator __last2,
                          __false_type) {
    iterator __first1 = begin();
    iterator __last1 = end();
    for ( ; __first1 != __last1 && __first2 != __last2; ++__first1, ++__first2)
      *__first1 = *__first2;
    if (__first2 == __last2)
      erase(__first1, __last1);
    else
      insert(__last1, __first2, __last2);
  }

#endif /* __STL_MEMBER_TEMPLATES */

protected:
public:
  void splice(iterator __position, _Self& __x) {
    __stl_verbose_assert(&__x!=this, _StlMsg_INVALID_ARGUMENT);
    __stl_debug_check(__check_if_owner(&_M_iter_list,__position));
    if (!__x.empty()) 
      _List_global_inst::_Transfer(__position._M_node, __x.begin()._M_node, __x.end()._M_node);
    __stl_debug_do(__x._Invalidate_all());
  }
# if defined ( __STL_DEBUG )
  void splice(iterator __position, _Self& __x, iterator __i) {
# else
  void splice(iterator __position, _Self&, iterator __i) {
# endif
    __stl_debug_check(__check_if_owner(&_M_iter_list,__position) &&
		      __check_if_owner(&__x._M_iter_list ,__i));
    __stl_verbose_assert(__i._M_node!=__i._Owner_node(), _StlMsg_NOT_DEREFERENCEABLE); 
    iterator __j = __i;
    ++__j;
#  if defined ( __STL_DEBUG )
    if (( &__x == this ) && (__position == __i || __position == __j)) return;
#  else
    if (__position == __i || __position == __j) return;
#  endif
    _List_global_inst::_Transfer(__position._M_node, __i._M_node, __j._M_node);
    __stl_debug_do(__invalidate_iterator(&__x._M_iter_list, __i));
  }
#  if defined ( __STL_DEBUG )
  void splice(iterator __position, _Self& __x, iterator __first, iterator __last) {
    __stl_debug_check(__check_if_owner(&_M_iter_list, __position));
    __stl_verbose_assert(__first._Owner()==&__x._M_iter_list && __last._Owner()==&__x._M_iter_list, 
			 _StlMsg_NOT_OWNER);
#  else
  void splice(iterator __position, _Self&, iterator __first, iterator __last) {
#  endif
    if (__first != __last) 
      _List_global_inst::_Transfer(__position._M_node, __first._M_node, __last._M_node);
  }
  void remove(const _Tp& __value);
  void unique();
  void merge(_Self& __x);
  void reverse();
  void sort();

#ifdef __STL_MEMBER_TEMPLATES
# ifndef __STL_INLINE_MEMBER_TEMPLATES
  template <class _Predicate> void remove_if(_Predicate);
  template <class _BinaryPredicate> void unique(_BinaryPredicate);
  template <class _StrictWeakOrdering> void merge(list<_Tp, _Alloc>&, _StrictWeakOrdering);
  template <class _StrictWeakOrdering> void sort(_StrictWeakOrdering);
# else
  template <class _Predicate>
    void remove_if(_Predicate __pred) 
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
  template <class _BinaryPredicate>
    void unique(_BinaryPredicate __binary_pred) {
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

  template <class _StrictWeakOrdering>
    void merge(list<_Tp, _Alloc>& __x,
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

  template <class _StrictWeakOrdering>
    void sort(_StrictWeakOrdering __comp) {
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
# endif /* __STL_INLINE_MEMBER_TEMPLATES */
#endif /* __STL_MEMBER_TEMPLATES */

};

template <class _Tp, class _Alloc>
__STL_INLINE_LOOP bool operator==(const list<_Tp,_Alloc>& __x,
                       const list<_Tp,_Alloc>& __y)
{
  typedef typename list<_Tp,_Alloc>::const_iterator const_iterator;
  const_iterator __end1 = __x.end();
  const_iterator __end2 = __y.end();

  const_iterator __i1 = __x.begin();
  const_iterator __i2 = __y.begin();
  while (__i1 != __end1 && __i2 != __end2 && *__i1 == *__i2) {
    ++__i1;
    ++__i2;
  }
  return __i1 == __end1 && __i2 == __end2;
}

template <class _Tp, class _Alloc>
inline bool operator<(const list<_Tp,_Alloc>& __x,
                      const list<_Tp,_Alloc>& __y)
{
  return lexicographical_compare(__x.begin(), __x.end(),
                                 __y.begin(), __y.end());
}

#ifdef __STL_USE_SEPARATE_RELOPS_NAMESPACE

template <class _Tp, class _Alloc>
inline void 
swap(list<_Tp, _Alloc>& __x, list<_Tp, _Alloc>& __y)
{
  __x.swap(__y);
}

#endif /* __STL_USE_SEPARATE_RELOPS_NAMESPACE */

// do a cleanup
# undef list
# define __list__ __FULL_NAME(list)

# if defined (__STL_USE_WRAPPER_FOR_ALLOC_PARAM)
// provide a "default" list adaptor
template <class _Tp>
class list : public __list__<_Tp, __STL_DEFAULT_ALLOCATOR(_Tp) >
{
#   define __LIST_SUPER __list__<_Tp, __STL_DEFAULT_ALLOCATOR(_Tp) >
public:
    typedef __LIST_SUPER _Super;
    __IMPORT_WITH_REVERSE_ITERATORS(_Super)
    __IMPORT_SUPER_COPY_ASSIGNMENT(list, list<_Tp>, __LIST_SUPER)
    list() { }
    explicit list(size_type __n, const _Tp& __value) : __LIST_SUPER(__n, __value) { }
    explicit list(size_type __n) :  __LIST_SUPER(__n) { } 
    list(const _Tp* __first, const _Tp* __last) : __LIST_SUPER(__first, __last) { } 
    list(const_iterator __first, const_iterator __last) : __LIST_SUPER(__first, __last) { }
# undef __LIST_SUPER
};

#  if defined (__STL_BASE_MATCH_BUG)
template <class _Tp>
inline bool operator==(const list<_Tp>& __x, const list<_Tp>& __y) {
    typedef typename list<_Tp>::_Super _Super;
    return operator == ((const _Super&)__x,(const _Super&)__y);
}

template <class _Tp>
inline bool operator<(const list<_Tp>& __x, const list<_Tp>& __y) {
  return lexicographical_compare(__x.begin(), __x.end(),
                                 __y.begin(), __y.end());
}
#  endif
# endif /*  WRAPPER */

__STL_END_NAMESPACE 

#if defined(__sgi) && !defined(__GNUC__) && (_MIPS_SIM != _MIPS_SIM_ABI32)
#pragma reset woff 1174
#pragma reset woff 1375
#endif


# if !defined (__STL_LINK_TIME_INSTANTIATION)
#  include <stl_list.c>
# endif

#endif /* __SGI_STL_INTERNAL_LIST_H */

// Local Variables:
// mode:C++
// End:
