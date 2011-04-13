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

#ifndef __SGI_STL_INTERNAL_MULTIMAP_H
#define __SGI_STL_INTERNAL_MULTIMAP_H

__STL_BEGIN_NAMESPACE

#if defined(__sgi) && !defined(__GNUC__) && (_MIPS_SIM != _MIPS_SIM_ABI32)
#pragma set woff 1174
#pragma set woff 1375
#endif

#define multimap __WORKAROUND_RENAME(multimap)

template <class _Key, class _Tp, __DFL_TMPL_PARAM(_Compare, less<_Key> ), 
          __STL_DEFAULT_PAIR_ALLOCATOR_SELECT(const _Key, _Tp) >
class multimap {
public:

// typedefs:

  typedef _Key                  key_type;
  typedef _Tp                   data_type;
  typedef _Tp                   mapped_type;
  typedef pair<const _Key, _Tp> value_type;
  typedef _Compare              key_compare;

  class value_compare : public binary_function<value_type, value_type, bool> {
  friend class multimap<_Key,_Tp,_Compare,_Alloc>;
  protected:
    _Compare _M_comp;
    value_compare(_Compare __c) : _M_comp(__c) {}
  public:
    bool operator()(const value_type& __x, const value_type& __y) const {
      return _M_comp(__x.first, __y.first);
    }
  };

private:
# ifdef __STL_MULTI_CONST_TEMPLATE_ARG_BUG
  typedef rb_tree<key_type, value_type, 
                  __select1st_hint<value_type, _Key>, key_compare, _Alloc> _Rep_type;
# else
  typedef rb_tree<key_type, value_type, 
                  select1st<value_type>, key_compare, _Alloc> _Rep_type;
# endif
  _Rep_type _M_t;  // red-black tree representing multimap
public:
  typedef typename _Rep_type::pointer pointer;
  typedef typename _Rep_type::const_pointer const_pointer;
  typedef typename _Rep_type::reference reference;
  typedef typename _Rep_type::const_reference const_reference;
  typedef typename _Rep_type::iterator iterator;
  typedef typename _Rep_type::const_iterator const_iterator; 
  typedef typename _Rep_type::reverse_iterator reverse_iterator;
  typedef typename _Rep_type::const_reverse_iterator const_reverse_iterator;
  typedef typename _Rep_type::size_type size_type;
  typedef typename _Rep_type::difference_type difference_type;
  typedef typename _Rep_type::allocator_type allocator_type;

// allocation/deallocation

  multimap() : _M_t(_Compare(), allocator_type()) { }
  explicit multimap(const _Compare& __comp,
                    const allocator_type& __a = __STL_ALLOC_INSTANCE(allocator_type))
    : _M_t(__comp, __a) { }

#ifdef __STL_MEMBER_TEMPLATES  
  template <class _InputIterator>
  multimap(_InputIterator __first, _InputIterator __last)
    : _M_t(_Compare(), allocator_type())
    { _M_t.insert_equal(__first, __last); }

  template <class _InputIterator>
  multimap(_InputIterator __first, _InputIterator __last,
           const _Compare& __comp,
           const allocator_type& __a = __STL_ALLOC_INSTANCE(allocator_type))
    : _M_t(__comp, __a) { _M_t.insert_equal(__first, __last); }
#else
  multimap(const value_type* __first, const value_type* __last)
    : _M_t(_Compare(), allocator_type())
    { _M_t.insert_equal(__first, __last); }
  multimap(const value_type* __first, const value_type* __last,
           const _Compare& __comp,
           const allocator_type& __a = __STL_ALLOC_INSTANCE(allocator_type))
    : _M_t(__comp, __a) { _M_t.insert_equal(__first, __last); }

  multimap(const_iterator __first, const_iterator __last)
    : _M_t(_Compare(), allocator_type())
    { _M_t.insert_equal(__first, __last); }
  multimap(const_iterator __first, const_iterator __last,
           const _Compare& __comp,
           const allocator_type& __a = __STL_ALLOC_INSTANCE(allocator_type))
    : _M_t(__comp, __a) { _M_t.insert_equal(__first, __last); }
#endif /* __STL_MEMBER_TEMPLATES */

  multimap(const multimap<_Key,_Tp,_Compare,_Alloc>& __x) : _M_t(__x._M_t) { }
  multimap<_Key,_Tp,_Compare,_Alloc>&
  operator=(const multimap<_Key,_Tp,_Compare,_Alloc>& __x) {
    _M_t = __x._M_t;
    return *this; 
  }

  // accessors:

  key_compare key_comp() const { return _M_t.key_comp(); }
  value_compare value_comp() const { return value_compare(_M_t.key_comp()); }
  allocator_type get_allocator() const { return _M_t.get_allocator(); }

  iterator begin() { return _M_t.begin(); }
  const_iterator begin() const { return _M_t.begin(); }
  iterator end() { return _M_t.end(); }
  const_iterator end() const { return _M_t.end(); }
  reverse_iterator rbegin() { return _M_t.rbegin(); }
  const_reverse_iterator rbegin() const { return _M_t.rbegin(); }
  reverse_iterator rend() { return _M_t.rend(); }
  const_reverse_iterator rend() const { return _M_t.rend(); }
  bool empty() const { return _M_t.empty(); }
  size_type size() const { return _M_t.size(); }
  size_type max_size() const { return _M_t.max_size(); }
  void swap(multimap<_Key,_Tp,_Compare,_Alloc>& __x) { _M_t.swap(__x._M_t); }

  // insert/erase

  iterator insert(const value_type& __x) { return _M_t.insert_equal(__x); }
  iterator insert(iterator __position, const value_type& __x) {
    return _M_t.insert_equal(__position, __x);
  }
#ifdef __STL_MEMBER_TEMPLATES  
  template <class _InputIterator>
  void insert(_InputIterator __first, _InputIterator __last) {
    _M_t.insert_equal(__first, __last);
  }
#else
  void insert(const value_type* __first, const value_type* __last) {
    _M_t.insert_equal(__first, __last);
  }
  void insert(const_iterator __first, const_iterator __last) {
    _M_t.insert_equal(__first, __last);
  }
#endif /* __STL_MEMBER_TEMPLATES */
  void erase(iterator __position) { _M_t.erase(__position); }
  size_type erase(const key_type& __x) { return _M_t.erase(__x); }
  void erase(iterator __first, iterator __last)
    { _M_t.erase(__first, __last); }
  void clear() { _M_t.clear(); }

  // multimap operations:

  iterator find(const key_type& __x) { return _M_t.find(__x); }
  const_iterator find(const key_type& __x) const { return _M_t.find(__x); }
  size_type count(const key_type& __x) const { return _M_t.count(__x); }
  iterator lower_bound(const key_type& __x) {return _M_t.lower_bound(__x); }
  const_iterator lower_bound(const key_type& __x) const {
    return _M_t.lower_bound(__x); 
  }
  iterator upper_bound(const key_type& __x) {return _M_t.upper_bound(__x); }
  const_iterator upper_bound(const key_type& __x) const {
    return _M_t.upper_bound(__x); 
  }
   pair<iterator,iterator> equal_range(const key_type& __x) {
    return _M_t.equal_range(__x);
  }
  pair<const_iterator,const_iterator> equal_range(const key_type& __x) const {
    return _M_t.equal_range(__x);
  }
};

template <class _Key, class _Tp, class _Compare, class _Alloc>
inline bool operator==(const multimap<_Key,_Tp,_Compare,_Alloc>& __x, 
                       const multimap<_Key,_Tp,_Compare,_Alloc>& __y) {
  return __x.size() == __y.size() &&
         equal(__x.begin(), __x.end(), __y.begin());
}

template <class _Key, class _Tp, class _Compare, class _Alloc>
inline bool operator<(const multimap<_Key,_Tp,_Compare,_Alloc>& __x, 
                      const multimap<_Key,_Tp,_Compare,_Alloc>& __y) {
  return lexicographical_compare(__x.begin(), __x.end(), 
                                 __y.begin(), __y.end());
}

#ifdef __STL_USE_SEPARATE_RELOPS_NAMESPACE

template <class _Key, class _Tp, class _Compare, class _Alloc>
inline bool operator!=(const multimap<_Key,_Tp,_Compare,_Alloc>& __x, 
                       const multimap<_Key,_Tp,_Compare,_Alloc>& __y) {
  return !(__x == __y);
}

template <class _Key, class _Tp, class _Compare, class _Alloc>
inline bool operator>(const multimap<_Key,_Tp,_Compare,_Alloc>& __x, 
                      const multimap<_Key,_Tp,_Compare,_Alloc>& __y) {
  return __y < __x;
}

template <class _Key, class _Tp, class _Compare, class _Alloc>
inline bool operator<=(const multimap<_Key,_Tp,_Compare,_Alloc>& __x, 
                       const multimap<_Key,_Tp,_Compare,_Alloc>& __y) {
  return !(__y < __x);
}

template <class _Key, class _Tp, class _Compare, class _Alloc>
inline bool operator>=(const multimap<_Key,_Tp,_Compare,_Alloc>& __x, 
                       const multimap<_Key,_Tp,_Compare,_Alloc>& __y) {
  return !(__x < __y);
}

#endif /* __STL_USE_SEPARATE_RELOPS_NAMESPACE */

#ifdef __STL_FUNCTION_TMPL_PARTIAL_ORDER
template <class _Key, class _Tp, class _Compare, class _Alloc>
inline void swap(multimap<_Key,_Tp,_Compare,_Alloc>& __x, 
                 multimap<_Key,_Tp,_Compare,_Alloc>& __y) {
  __x.swap(__y);
}
#endif /* __STL_FUNCTION_TMPL_PARTIAL_ORDER */

// do a cleanup
#  undef multimap
// provide a way to access full funclionality 
#  define __multimap__  __FULL_NAME(multimap)

# ifdef __STL_USE_WRAPPER_FOR_ALLOC_PARAM

#  if defined (__STL_MINIMUM_DEFAULT_TEMPLATE_PARAMS)
#   define __MMAP_TEMPLATE_HEADER  template <class _Key, class _Tp>
#   define __MMAP_ARGUMENTS        _Key, _Tp
#   define _Compare less<_Key>
#  else
#   define __MMAP_TEMPLATE_HEADER  template <class _Key, class _Tp, class _Compare >
#   define __MMAP_ARGUMENTS        _Key, _Tp, _Compare
#  endif

#   define __MMAP_SUPER  __multimap< _Key, _Tp, _Compare, __STL_DEFAULT_PAIR_ALLOCATOR(const _Key, _Tp) >

// provide a "default" multimap adaptor
__MMAP_TEMPLATE_HEADER
class multimap : public __MMAP_SUPER
{
  typedef multimap< __MMAP_ARGUMENTS > _Self;
public:
    typedef __MMAP_SUPER  _Super;
    __IMPORT_WITH_REVERSE_ITERATORS(_Super)
    // copy & assignment from super
    __IMPORT_SUPER_COPY_ASSIGNMENT(multimap, _Self, __MMAP_SUPER)
    multimap() : __MMAP_SUPER(_Compare()) {}
    explicit multimap(const _Compare& __comp) : __MMAP_SUPER(__comp) {}
    multimap(const typename _Super::value_type* __first, 
	     const typename _Super::value_type* __last) : 
        __MMAP_SUPER(__first, __last, _Compare()) { }
    multimap(const typename _Super::value_type* __first,
	     const typename _Super::value_type* __last, 
	     const _Compare& __comp) : __MMAP_SUPER(__first, __last, __comp) { }
    multimap(typename _Super::const_iterator __first, 
	     typename _Super::const_iterator __last) : 
      __MMAP_SUPER(__first, __last, _Compare()) { }
    multimap(typename _Super::const_iterator __first, 
	     typename _Super::const_iterator __last, 
	     const _Compare& __comp) : __MMAP_SUPER(__first, __last, __comp) { }
};

#  if defined (__STL_BASE_MATCH_BUG)
__MMAP_TEMPLATE_HEADER
inline bool operator==(const multimap< __MMAP_ARGUMENTS >& __x, 
                       const multimap< __MMAP_ARGUMENTS >& __y) {
  typedef __MMAP_SUPER  _Super;
  return (const _Super&)__x == (const _Super&)__y;
}

__MMAP_TEMPLATE_HEADER
inline bool operator<(const multimap< __MMAP_ARGUMENTS >& __x, 
                      const multimap< __MMAP_ARGUMENTS >& __y) {
  typedef __MMAP_SUPER  _Super;
  return (const _Super&)__x < (const _Super&)__y;
}
#  endif

# undef __MMAP_TEMPLATE_HEADER
# undef __MMAP_ARGUMENTS
# undef __MMAP_SUPER
# undef _Compare

# endif /*  WRAPPER */

#if defined(__sgi) && !defined(__GNUC__) && (_MIPS_SIM != _MIPS_SIM_ABI32)
#pragma reset woff 1174
#pragma reset woff 1375
#endif

__STL_END_NAMESPACE

#endif /* __SGI_STL_INTERNAL_MULTIMAP_H */

// Local Variables:
// mode:C++
// End:
