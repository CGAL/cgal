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

#ifndef __SGI_STL_INTERNAL_SET_H
#define __SGI_STL_INTERNAL_SET_H

#if defined(__sgi) && !defined(__GNUC__) && (_MIPS_SIM != _MIPS_SIM_ABI32)
#pragma set woff 1174
#pragma set woff 1375
#endif
#define set __WORKAROUND_RENAME(set)

__STL_BEGIN_NAMESPACE

template <class _Key, __DFL_TMPL_PARAM(_Compare,less<_Key>), 
                     __STL_DEFAULT_ALLOCATOR_SELECT(_Key) >
class set {
public:
// typedefs:
  typedef _Key     key_type;
  typedef _Key     value_type;
  typedef _Compare key_compare;
  typedef _Compare value_compare;
private:
  typedef _Rb_tree<key_type, value_type, 
                  _Identity<value_type>, key_compare, _Alloc> _Rep_type;
public:
  typedef typename _Rep_type::const_pointer pointer;
  typedef typename _Rep_type::const_pointer const_pointer;
  typedef typename _Rep_type::const_reference reference;
  typedef typename _Rep_type::const_reference const_reference;
  typedef typename _Rep_type::const_iterator const_iterator;
  typedef const_iterator iterator;
  typedef typename _Rep_type::const_reverse_iterator reverse_iterator;
  typedef typename _Rep_type::const_reverse_iterator const_reverse_iterator;
  typedef typename _Rep_type::size_type size_type;
  typedef typename _Rep_type::difference_type difference_type;
  typedef typename _Rep_type::allocator_type allocator_type;

private:
  _Rep_type _M_t;  // red-black tree representing set
public:

  // allocation/deallocation

  set() : _M_t(_Compare(), allocator_type()) {}
  explicit set(const _Compare& __comp,
	       const allocator_type& __a = __STL_ALLOC_INSTANCE(allocator_type))
    : _M_t(__comp, __a) {}

#ifdef __STL_MEMBER_TEMPLATES
  template <class _InputIterator>
  set(_InputIterator __first, _InputIterator __last)
    : _M_t(_Compare(), allocator_type())
    { _M_t.insert_unique(__first, __last); }

  template <class _InputIterator>
  set(_InputIterator __first, _InputIterator __last, const _Compare& __comp,
      const allocator_type& __a = __STL_ALLOC_INSTANCE(allocator_type))
    : _M_t(__comp, __a) { _M_t.insert_unique(__first, __last); }
#else
  set(const value_type* __first, const value_type* __last) 
    : _M_t(_Compare(), allocator_type()) 
    { _M_t.insert_unique(__first, __last); }

  set(const value_type* __first, 
      const value_type* __last, const _Compare& __comp,
      const allocator_type& __a = __STL_ALLOC_INSTANCE(allocator_type))
    : _M_t(__comp, __a) { _M_t.insert_unique(__first, __last); }

  set(const_iterator __first, const_iterator __last)
    : _M_t(_Compare(), allocator_type()) 
    { _M_t.insert_unique(__first, __last); }

  set(const_iterator __first, const_iterator __last, const _Compare& __comp,
      const allocator_type& __a = __STL_ALLOC_INSTANCE(allocator_type))
    : _M_t(__comp, __a) { _M_t.insert_unique(__first, __last); }
#endif /* __STL_MEMBER_TEMPLATES */

  set(const set<_Key,_Compare,_Alloc>& __x) : _M_t(__x._M_t) {}
  set<_Key,_Compare,_Alloc>& operator=(const set<_Key, _Compare, _Alloc>& __x)
  { 
    _M_t = __x._M_t; 
    return *this;
  }

  // accessors:

  key_compare key_comp() const { return _M_t.key_comp(); }
  value_compare value_comp() const { return _M_t.key_comp(); }
  allocator_type get_allocator() const { return _M_t.get_allocator(); }

  iterator begin() const { return _M_t.begin(); }
  iterator end() const { return _M_t.end(); }
  reverse_iterator rbegin() const { return _M_t.rbegin(); } 
  reverse_iterator rend() const { return _M_t.rend(); }
  bool empty() const { return _M_t.empty(); }
  size_type size() const { return _M_t.size(); }
  size_type max_size() const { return _M_t.max_size(); }
  void swap(set<_Key,_Compare,_Alloc>& __x) { _M_t.swap(__x._M_t); }

  // insert/erase
  pair<iterator,bool> insert(const value_type& __x) { 
    typedef typename _Rep_type::iterator _Rep_iterator;
    pair<_Rep_iterator, bool> __p = _M_t.insert_unique(__x); 
    return pair<iterator, bool>((const iterator&)__p.first, __p.second);
  }
  iterator insert(iterator __position, const value_type& __x) {
    typedef typename _Rep_type::iterator _Rep_iterator;
    return _M_t.insert_unique((_Rep_iterator&)__position, __x);
  }
#ifdef __STL_MEMBER_TEMPLATES
  template <class _InputIterator>
  void insert(_InputIterator __first, _InputIterator __last) {
    _M_t.insert_unique(__first, __last);
  }
#else
  void insert(const_iterator __first, const_iterator __last) {
    _M_t.insert_unique(__first, __last);
  }
  void insert(const value_type* __first, const value_type* __last) {
    _M_t.insert_unique(__first, __last);
  }
#endif /* __STL_MEMBER_TEMPLATES */
  void erase(iterator __position) { 
    typedef typename _Rep_type::iterator _Rep_iterator;
    _M_t.erase((_Rep_iterator&)__position); 
  }
  size_type erase(const key_type& __x) { 
    return _M_t.erase(__x); 
  }
  void erase(iterator __first, iterator __last) { 
    typedef typename _Rep_type::iterator _Rep_iterator;
    _M_t.erase((_Rep_iterator&)__first, (_Rep_iterator&)__last); 
  }
  void clear() { _M_t.clear(); }

  // set operations:

  iterator find(const key_type& __x) const { return _M_t.find(__x); }
  size_type count(const key_type& __x) const { return _M_t.count(__x); }
  iterator lower_bound(const key_type& __x) const {
    return _M_t.lower_bound(__x);
  }
  iterator upper_bound(const key_type& __x) const {
    return _M_t.upper_bound(__x); 
  }
  pair<iterator,iterator> equal_range(const key_type& __x) const {
    return _M_t.equal_range(__x);
  }
};

template <class _Key, class _Compare, class _Alloc>
inline bool operator==(const set<_Key,_Compare,_Alloc>& __x, 
                       const set<_Key,_Compare,_Alloc>& __y) {
  return __x.size() == __y.size() &&
         equal(__x.begin(), __x.end(), __y.begin());
}

template <class _Key, class _Compare, class _Alloc>
inline bool operator<(const set<_Key,_Compare,_Alloc>& __x, 
                      const set<_Key,_Compare,_Alloc>& __y) {
  return lexicographical_compare(__x.begin(), __x.end(), 
                                 __y.begin(), __y.end());
}

#ifdef __STL_USE_SEPARATE_RELOPS_NAMESPACE

template <class _Key, class _Compare, class _Alloc>
inline bool operator!=(const set<_Key,_Compare,_Alloc>& __x, 
                       const set<_Key,_Compare,_Alloc>& __y) {
  return !(__x == __y);
}

template <class _Key, class _Compare, class _Alloc>
inline bool operator>(const set<_Key,_Compare,_Alloc>& __x, 
                      const set<_Key,_Compare,_Alloc>& __y) {
  return __y < __x;
}

template <class _Key, class _Compare, class _Alloc>
inline bool operator<=(const set<_Key,_Compare,_Alloc>& __x, 
                       const set<_Key,_Compare,_Alloc>& __y) {
  return !(__y < __x);
}

template <class _Key, class _Compare, class _Alloc>
inline bool operator>=(const set<_Key,_Compare,_Alloc>& __x, 
                       const set<_Key,_Compare,_Alloc>& __y) {
  return !(__x < __y);
}
#endif /* __STL_USE_SEPARATE_RELOPS_NAMESPACE */

#ifdef __STL_FUNCTION_TMPL_PARTIAL_ORDER

template <class _Key, class _Compare, class _Alloc>
inline void swap(set<_Key,_Compare,_Alloc>& __x, 
                 set<_Key,_Compare,_Alloc>& __y) {
  __x.swap(__y);
}

#endif /* __STL_FUNCTION_TMPL_PARTIAL_ORDER */

// do a cleanup
# undef set
// provide a way to access full funclionality 
# define __set__  __FULL_NAME(set)

# ifdef __STL_USE_WRAPPER_FOR_ALLOC_PARAM

#  if defined (__STL_MINIMUM_DEFAULT_TEMPLATE_PARAMS)
#   define __SET_TEMPLATE_HEADER  template <class _Key>
#   define __SET_ARGUMENTS        _Key
#   define _Compare less<_Key>
#  else
#   define __SET_TEMPLATE_HEADER  template <class _Key, class _Compare >
#   define __SET_ARGUMENTS        _Key, _Compare
#  endif

#   define __SET_SUPER  __set< _Key, _Compare, __STL_DEFAULT_ALLOCATOR(_Key) >

// provide a "default" set adaptor
__SET_TEMPLATE_HEADER
class set : public __SET_SUPER
{
  typedef set< __SET_ARGUMENTS > _Self;
public:
    typedef __SET_SUPER _Super;
    __IMPORT_WITH_REVERSE_ITERATORS(_Super)
    // copy & assignment from super
    __IMPORT_SUPER_COPY_ASSIGNMENT(set,_Self,__SET_SUPER)
    // specific constructors
    explicit set() : __SET_SUPER(_Compare()) {}
    explicit set(const _Compare& __comp) : __SET_SUPER(__comp) {}
    set(const value_type* __first, const value_type* __last) : 
        __SET_SUPER(__first, __last, _Compare()) { }
    set(const value_type* __first, const value_type* __last, 
        const _Compare& __comp) : __SET_SUPER(__first, __last, __comp) { }
    set(const_iterator __first, const_iterator __last) : 
        __SET_SUPER(__first, __last, _Compare()) { }
    set(const_iterator __first, const_iterator __last, 
        const _Compare& __comp) : __SET_SUPER(__first, __last, __comp) { }
};

#  if defined (__STL_BASE_MATCH_BUG)
__SET_TEMPLATE_HEADER 
inline bool operator==(const set< __SET_ARGUMENTS >& __x, 
                       const set< __SET_ARGUMENTS >& __y) {
  typedef __SET_SUPER _Super;
  return operator==((const _Super&)__x,(const _Super&)__y);
}

__SET_TEMPLATE_HEADER 
inline bool operator<(const set< __SET_ARGUMENTS >& __x, 
                      const set< __SET_ARGUMENTS >& __y) {
  typedef __SET_SUPER _Super;
  return operator < ((const _Super&)__x , (const _Super&)__y);
}
#  endif

# undef __SET_TEMPLATE_HEADER
# undef __SET_ARGUMENTS
# undef __SET_SUPER 
# undef _Compare

# endif /*  WRAPPER */

#if defined(__sgi) && !defined(__GNUC__) && (_MIPS_SIM != _MIPS_SIM_ABI32)
#pragma reset woff 1174
#pragma reset woff 1375
#endif

__STL_END_NAMESPACE

#endif /* __SGI_STL_INTERNAL_SET_H */

// Local Variables:
// mode:C++
// End:
