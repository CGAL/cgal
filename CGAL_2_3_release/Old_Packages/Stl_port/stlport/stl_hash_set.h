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

#ifndef __SGI_STL_INTERNAL_HASH_SET_H
#define __SGI_STL_INTERNAL_HASH_SET_H

__STL_BEGIN_NAMESPACE

#if defined(__sgi) && !defined(__GNUC__) && (_MIPS_SIM != _MIPS_SIM_ABI32)
#pragma set woff 1174
#pragma set woff 1375
#endif

# define  hash_set      __WORKAROUND_RENAME(hash_set)
# define  hash_multiset __WORKAROUND_RENAME(hash_multiset)

template <class _Value, __DFL_TMPL_PARAM(_HashFcn,hash<_Value>),
          __DFL_TMPL_PARAM(_EqualKey,equal_to<_Value>),
          __STL_DEFAULT_ALLOCATOR_SELECT(_Value) >
class hash_set
{
private:
  typedef hashtable<_Value, _Value, _HashFcn, _Identity<_Value>, 
                    _EqualKey, _Alloc> _Ht;
  typedef hash_set<_Value, _HashFcn, _EqualKey, _Alloc> _Self;
  typedef typename _Ht::iterator _ht_iterator;
public:
  typedef typename _Ht::key_type key_type;
  typedef typename _Ht::value_type value_type;
  typedef typename _Ht::hasher hasher;
  typedef typename _Ht::key_equal key_equal;

  typedef typename _Ht::size_type size_type;
  typedef typename _Ht::difference_type difference_type;
  typedef typename _Ht::const_pointer pointer;
  typedef typename _Ht::const_pointer const_pointer;
  typedef typename _Ht::const_reference reference;
  typedef typename _Ht::const_reference const_reference;

  // SunPro bug
  typedef typename _Ht::const_iterator const_iterator;
  typedef const_iterator iterator;

  typedef typename _Ht::allocator_type allocator_type;

  hasher hash_funct() const { return _M_ht.hash_funct(); }
  key_equal key_eq() const { return _M_ht.key_eq(); }
  allocator_type get_allocator() const { return _M_ht.get_allocator(); }

private:
  _Ht _M_ht;

public:
  hash_set()
    : _M_ht(100, hasher(), key_equal(), allocator_type()) {}
  explicit hash_set(size_type __n)
    : _M_ht(__n, hasher(), key_equal(), allocator_type()) {}
  hash_set(size_type __n, const hasher& __hf)
    : _M_ht(__n, __hf, key_equal(), allocator_type()) {}
  hash_set(size_type __n, const hasher& __hf, const key_equal& __eql,
           const allocator_type& __a = __STL_ALLOC_INSTANCE(allocator_type))
    : _M_ht(__n, __hf, __eql, __a) {}

#ifdef __STL_MEMBER_TEMPLATES
  template <class _InputIterator>
  hash_set(_InputIterator __f, _InputIterator __l)
    : _M_ht(100, hasher(), key_equal(), allocator_type())
    { _M_ht.insert_unique(__f, __l); }
  template <class _InputIterator>
  hash_set(_InputIterator __f, _InputIterator __l, size_type __n)
    : _M_ht(__n, hasher(), key_equal(), allocator_type())
    { _M_ht.insert_unique(__f, __l); }
  template <class _InputIterator>
  hash_set(_InputIterator __f, _InputIterator __l, size_type __n,
           const hasher& __hf)
    : _M_ht(__n, __hf, key_equal(), allocator_type())
    { _M_ht.insert_unique(__f, __l); }
  template <class _InputIterator>
  hash_set(_InputIterator __f, _InputIterator __l, size_type __n,
           const hasher& __hf, const key_equal& __eql,
           const allocator_type& __a = __STL_ALLOC_INSTANCE(allocator_type))
    : _M_ht(__n, __hf, __eql, __a)
    { _M_ht.insert_unique(__f, __l); }
#else

  hash_set(const value_type* __f, const value_type* __l)
    : _M_ht(100, hasher(), key_equal(), allocator_type())
    { _M_ht.insert_unique(__f, __l); }
  hash_set(const value_type* __f, const value_type* __l, size_type __n)
    : _M_ht(__n, hasher(), key_equal(), allocator_type())
    { _M_ht.insert_unique(__f, __l); }
  hash_set(const value_type* __f, const value_type* __l, size_type __n,
           const hasher& __hf)
    : _M_ht(__n, __hf, key_equal(), allocator_type())
    { _M_ht.insert_unique(__f, __l); }
  hash_set(const value_type* __f, const value_type* __l, size_type __n,
           const hasher& __hf, const key_equal& __eql,
           const allocator_type& __a = __STL_ALLOC_INSTANCE(allocator_type))
    : _M_ht(__n, __hf, __eql, __a)
    { _M_ht.insert_unique(__f, __l); }

  hash_set(const_iterator __f, const_iterator __l)
    : _M_ht(100, hasher(), key_equal(), allocator_type())
    { _M_ht.insert_unique(__f, __l); }
  hash_set(const_iterator __f, const_iterator __l, size_type __n)
    : _M_ht(__n, hasher(), key_equal(), allocator_type())
    { _M_ht.insert_unique(__f, __l); }
  hash_set(const_iterator __f, const_iterator __l, size_type __n,
           const hasher& __hf)
    : _M_ht(__n, __hf, key_equal(), allocator_type())
    { _M_ht.insert_unique(__f, __l); }
  hash_set(const_iterator __f, const_iterator __l, size_type __n,
           const hasher& __hf, const key_equal& __eql,
           const allocator_type& __a = __STL_ALLOC_INSTANCE(allocator_type))
    : _M_ht(__n, __hf, __eql, __a)
    { _M_ht.insert_unique(__f, __l); }
#endif /*__STL_MEMBER_TEMPLATES */

public:
  size_type size() const { return _M_ht.size(); }
  size_type max_size() const { return _M_ht.max_size(); }
  bool empty() const { return _M_ht.empty(); }
  void swap(_Self& __hs) { _M_ht.swap(__hs._M_ht); }

  iterator begin() const { return _M_ht.begin(); }
  iterator end() const { return _M_ht.end(); }

public:
  pair<iterator, bool> insert(const value_type& __obj)
    {
      pair<_ht_iterator, bool> __p = _M_ht.insert_unique(__obj);
      return pair<iterator,bool>((const iterator&)__p.first, __p.second);
    }
#ifdef __STL_MEMBER_TEMPLATES
  template <class _InputIterator>
  void insert(_InputIterator __f, _InputIterator __l) 
    { _M_ht.insert_unique(__f,__l); }
#else
  void insert(const value_type* __f, const value_type* __l) {
    _M_ht.insert_unique(__f,__l);
  }
  void insert(const_iterator __f, const_iterator __l) 
    {_M_ht.insert_unique(__f, __l); }
#endif /*__STL_MEMBER_TEMPLATES */
  pair<iterator, bool> insert_noresize(const value_type& __obj)
  {
    pair<_ht_iterator, bool> __p = 
      _M_ht.insert_unique_noresize(__obj);
    return pair<iterator, bool>(__p.first, __p.second);
  }

  iterator find(const key_type& __key) const { return _M_ht.find(__key); }

  size_type count(const key_type& __key) const { return _M_ht.count(__key); }
  
  pair<iterator, iterator> equal_range(const key_type& __key) const
    { return _M_ht.equal_range(__key); }

  size_type erase(const key_type& __key) {return _M_ht.erase(__key); }
  void erase(iterator __it) { _M_ht.erase(__it); }
  void erase(iterator __f, iterator __l) { _M_ht.erase(__f, __l); }
  void clear() { _M_ht.clear(); }

public:
  void resize(size_type __hint) { _M_ht.resize(__hint); }
  size_type bucket_count() const { return _M_ht.bucket_count(); }
  size_type max_bucket_count() const { return _M_ht.max_bucket_count(); }
  size_type elems_in_bucket(size_type __n) const
    { return _M_ht.elems_in_bucket(__n); }

  static bool _M_equal (const _Self& __x, const _Self& __y) {
    return _Ht::_M_equal(__x._M_ht,__y._M_ht);
  }

};


template <class _Value, class _HashFcn, class _EqualKey, class _Alloc>
inline bool 
operator==(const hash_set<_Value,_HashFcn,_EqualKey,_Alloc>& __hs1,
           const hash_set<_Value,_HashFcn,_EqualKey,_Alloc>& __hs2)
{
  return hash_set<_Value,_HashFcn,_EqualKey,_Alloc>::_M_equal(__hs1, __hs2);
}

#ifdef __STL_USE_SEPARATE_RELOPS_NAMESPACE

template <class _Value, class _HashFcn, class _EqualKey, class _Alloc>
inline bool 
operator!=(const hash_set<_Value,_HashFcn,_EqualKey,_Alloc>& __hs1,
           const hash_set<_Value,_HashFcn,_EqualKey,_Alloc>& __hs2) {
  return !(__hs1 == __hs2);
}

template <class _Val, class _HashFcn, class _EqualKey, class _Alloc>
inline void 
swap(hash_set<_Val,_HashFcn,_EqualKey,_Alloc>& __hs1,
     hash_set<_Val,_HashFcn,_EqualKey,_Alloc>& __hs2)
{
  __hs1.swap(__hs2);
}

#endif /* __STL_USE_SEPARATE_RELOPS_NAMESPACE */


template <class _Value, __DFL_TMPL_PARAM(_HashFcn,hash<_Value>),
          __DFL_TMPL_PARAM(_EqualKey,equal_to<_Value>),
          __STL_DEFAULT_ALLOCATOR_SELECT(_Value) >
class hash_multiset
{
private:
  typedef hashtable<_Value, _Value, _HashFcn, _Identity<_Value>, 
                    _EqualKey, _Alloc> _Ht;
  typedef hash_multiset<_Value, _HashFcn, _EqualKey, _Alloc> _Self;

public:
  typedef typename _Ht::key_type key_type;
  typedef typename _Ht::value_type value_type;
  typedef typename _Ht::hasher hasher;
  typedef typename _Ht::key_equal key_equal;

  typedef typename _Ht::size_type size_type;
  typedef typename _Ht::difference_type difference_type;
  typedef typename _Ht::const_pointer pointer;
  typedef typename _Ht::const_pointer const_pointer;
  typedef typename _Ht::const_reference reference;
  typedef typename _Ht::const_reference const_reference;

  typedef typename _Ht::const_iterator const_iterator;
  // SunPro bug
  typedef const_iterator iterator;

  typedef typename _Ht::allocator_type allocator_type;

  hasher hash_funct() const { return _M_ht.hash_funct(); }
  key_equal key_eq() const { return _M_ht.key_eq(); }
  allocator_type get_allocator() const { return _M_ht.get_allocator(); }

private:
  _Ht _M_ht;

public:
  hash_multiset()
    : _M_ht(100, hasher(), key_equal(), allocator_type()) {}
  explicit hash_multiset(size_type __n)
    : _M_ht(__n, hasher(), key_equal(), allocator_type()) {}
  hash_multiset(size_type __n, const hasher& __hf)
    : _M_ht(__n, __hf, key_equal(), allocator_type()) {}
  hash_multiset(size_type __n, const hasher& __hf, const key_equal& __eql,
                const allocator_type& __a = __STL_ALLOC_INSTANCE(allocator_type))
    : _M_ht(__n, __hf, __eql, __a) {}

#ifdef __STL_MEMBER_TEMPLATES
  template <class _InputIterator>
  hash_multiset(_InputIterator __f, _InputIterator __l)
    : _M_ht(100, hasher(), key_equal(), allocator_type())
    { _M_ht.insert_equal(__f, __l); }
  template <class _InputIterator>
  hash_multiset(_InputIterator __f, _InputIterator __l, size_type __n)
    : _M_ht(__n, hasher(), key_equal(), allocator_type())
    { _M_ht.insert_equal(__f, __l); }
  template <class _InputIterator>
  hash_multiset(_InputIterator __f, _InputIterator __l, size_type __n,
                const hasher& __hf)
    : _M_ht(__n, __hf, key_equal(), allocator_type())
    { _M_ht.insert_equal(__f, __l); }
  template <class _InputIterator>
  hash_multiset(_InputIterator __f, _InputIterator __l, size_type __n,
                const hasher& __hf, const key_equal& __eql,
                const allocator_type& __a = __STL_ALLOC_INSTANCE(allocator_type))
    : _M_ht(__n, __hf, __eql, __a)
    { _M_ht.insert_equal(__f, __l); }
#else

  hash_multiset(const value_type* __f, const value_type* __l)
    : _M_ht(100, hasher(), key_equal(), allocator_type())
    { _M_ht.insert_equal(__f, __l); }
  hash_multiset(const value_type* __f, const value_type* __l, size_type __n)
    : _M_ht(__n, hasher(), key_equal(), allocator_type())
    { _M_ht.insert_equal(__f, __l); }
  hash_multiset(const value_type* __f, const value_type* __l, size_type __n,
                const hasher& __hf)
    : _M_ht(__n, __hf, key_equal(), allocator_type())
    { _M_ht.insert_equal(__f, __l); }
  hash_multiset(const value_type* __f, const value_type* __l, size_type __n,
                const hasher& __hf, const key_equal& __eql,
                const allocator_type& __a = __STL_ALLOC_INSTANCE(allocator_type))
    : _M_ht(__n, __hf, __eql, __a)
    { _M_ht.insert_equal(__f, __l); }

  hash_multiset(const_iterator __f, const_iterator __l)
    : _M_ht(100, hasher(), key_equal(), allocator_type())
    { _M_ht.insert_equal(__f, __l); }
  hash_multiset(const_iterator __f, const_iterator __l, size_type __n)
    : _M_ht(__n, hasher(), key_equal(), allocator_type())
    { _M_ht.insert_equal(__f, __l); }
  hash_multiset(const_iterator __f, const_iterator __l, size_type __n,
                const hasher& __hf)
    : _M_ht(__n, __hf, key_equal(), allocator_type())
    { _M_ht.insert_equal(__f, __l); }
  hash_multiset(const_iterator __f, const_iterator __l, size_type __n,
                const hasher& __hf, const key_equal& __eql,
                const allocator_type& __a = __STL_ALLOC_INSTANCE(allocator_type))
    : _M_ht(__n, __hf, __eql, __a)
    { _M_ht.insert_equal(__f, __l); }
#endif /*__STL_MEMBER_TEMPLATES */

public:
  size_type size() const { return _M_ht.size(); }
  size_type max_size() const { return _M_ht.max_size(); }
  bool empty() const { return _M_ht.empty(); }
  void swap(_Self& hs) { _M_ht.swap(hs._M_ht); }

  iterator begin() const { return _M_ht.begin(); }
  iterator end() const { return _M_ht.end(); }

public:
  iterator insert(const value_type& __obj)
    { return _M_ht.insert_equal(__obj); }
#ifdef __STL_MEMBER_TEMPLATES
  template <class _InputIterator>
  void insert(_InputIterator __f, _InputIterator __l) 
    { _M_ht.insert_equal(__f,__l); }
#else
  void insert(const value_type* __f, const value_type* __l) {
    _M_ht.insert_equal(__f,__l);
  }
  void insert(const_iterator __f, const_iterator __l) 
    { _M_ht.insert_equal(__f, __l); }
#endif /*__STL_MEMBER_TEMPLATES */
  iterator insert_noresize(const value_type& __obj)
    { return _M_ht.insert_equal_noresize(__obj); }    

  iterator find(const key_type& __key) const { return _M_ht.find(__key); }

  size_type count(const key_type& __key) const { return _M_ht.count(__key); }
  
  pair<iterator, iterator> equal_range(const key_type& __key) const
    { return _M_ht.equal_range(__key); }

  size_type erase(const key_type& __key) {return _M_ht.erase(__key); }
  void erase(iterator __it) { _M_ht.erase(__it); }
  void erase(iterator __f, iterator __l) { _M_ht.erase(__f, __l); }
  void clear() { _M_ht.clear(); }

public:
  void resize(size_type __hint) { _M_ht.resize(__hint); }
  size_type bucket_count() const { return _M_ht.bucket_count(); }
  size_type max_bucket_count() const { return _M_ht.max_bucket_count(); }
  size_type elems_in_bucket(size_type __n) const
    { return _M_ht.elems_in_bucket(__n); }
  static bool _M_equal (const _Self& __x, const _Self& __y) {
    return _Ht::_M_equal(__x._M_ht,__y._M_ht);
  }
};

template <class _Value, class _HashFcn, class _EqualKey, class _Alloc>
inline bool 
operator==(const hash_multiset<_Value,_HashFcn,_EqualKey,_Alloc>& __hs1,
           const hash_multiset<_Value,_HashFcn,_EqualKey,_Alloc>& __hs2)
{
  return hash_multiset<_Value,_HashFcn,_EqualKey,_Alloc>::_M_equal(__hs1, __hs2);
}

#ifdef __STL_USE_SEPARATE_RELOPS_NAMESPACE

template <class _Value, class _HashFcn, class _EqualKey, class _Alloc>
inline bool 
operator!=(const hash_multiset<_Value,_HashFcn,_EqualKey,_Alloc>& __hs1,
           const hash_multiset<_Value,_HashFcn,_EqualKey,_Alloc>& __hs2) {
  return !(__hs1 == __hs2);
}

template <class _Val, class _HashFcn, class _EqualKey, class _Alloc>
inline void 
swap(hash_multiset<_Val,_HashFcn,_EqualKey,_Alloc>& __hs1,
     hash_multiset<_Val,_HashFcn,_EqualKey,_Alloc>& __hs2) {
  __hs1.swap(__hs2);
}

#endif /* __STL_USE_SEPARATE_RELOPS_NAMESPACE */

// Specialization of insert_iterator so that it will work for hash_set
// and hash_multiset.

#if defined( __STL_CLASS_PARTIAL_SPECIALIZATION) && \
!defined(CGAL_LIMITED_ITERATOR_TRAITS_SUPPORT)

template <class _Value, class _HashFcn, class _EqualKey, class _Alloc>
class insert_iterator<hash_set<_Value, _HashFcn, _EqualKey, _Alloc> > {
protected:
  typedef hash_set<_Value, _HashFcn, _EqualKey, _Alloc> _Container;
  _Container* container;
public:
  typedef _Container          container_type;
  typedef output_iterator_tag iterator_category;
  typedef void                value_type;
  typedef void                difference_type;
  typedef void                pointer;
  typedef void                reference;

  insert_iterator(_Container& __x) : container(&__x) {}
  insert_iterator(_Container& __x, typename _Container::iterator)
    : container(&__x) {}
  insert_iterator<_Container>&
  operator=(const typename _Container::value_type& __value) { 
    container->insert(__value);
    return *this;
  }
  insert_iterator<_Container>& operator*() { return *this; }
  insert_iterator<_Container>& operator++() { return *this; }
  insert_iterator<_Container>& operator++(int) { return *this; }
};

template <class _Value, class _HashFcn, class _EqualKey, class _Alloc>
class insert_iterator<hash_multiset<_Value, _HashFcn, _EqualKey, _Alloc> > {
protected:
  typedef hash_multiset<_Value, _HashFcn, _EqualKey, _Alloc> _Container;
  _Container* container;
  typename _Container::iterator iter;
public:
  typedef _Container          container_type;
  typedef output_iterator_tag iterator_category;
  typedef void                value_type;
  typedef void                difference_type;
  typedef void                pointer;
  typedef void                reference;

  insert_iterator(_Container& __x) : container(&__x) {}
  insert_iterator(_Container& __x, typename _Container::iterator)
    : container(&__x) {}
  insert_iterator<_Container>&
  operator=(const typename _Container::value_type& __value) { 
    container->insert(__value);
    return *this;
  }
  insert_iterator<_Container>& operator*() { return *this; }
  insert_iterator<_Container>& operator++() { return *this; }
  insert_iterator<_Container>& operator++(int) { return *this; }
};

#endif /* __STL_CLASS_PARTIAL_SPECIALIZATION */

// do a cleanup
#  undef hash_set
#  undef hash_multiset
// provide a uniform way to access full funclionality 
#  define __hash_set__       __FULL_NAME(hash_set)
#  define __hash_multiset__  __FULL_NAME(hash_multiset)

# if defined ( __STL_USE_WRAPPER_FOR_ALLOC_PARAM )

#  if defined (__STL_MINIMUM_DEFAULT_TEMPLATE_PARAMS)
#   define __HS_TEMPLATE_HEADER  template <class _Value>
#   define __HS_ARGUMENTS        _Value
#   define __HS_BASE_ARGUMENTS   _Value, hash<_Value>, equal_to<_Value>, __STL_DEFAULT_ALLOCATOR(_Value)
#  else
#   define __HS_TEMPLATE_HEADER  template <class _Value, class _HashFcn, class _EqualKey >
#   define __HS_ARGUMENTS        _Value, _HashFcn, _EqualKey
#   define __HS_BASE_ARGUMENTS   _Value, _HashFcn, _EqualKey, __STL_DEFAULT_ALLOCATOR(_Value)
#  endif

# define __HS_SUPER  __hash_set< __HS_BASE_ARGUMENTS >
# define __HMS_SUPER __hash_multiset< __HS_BASE_ARGUMENTS >

// provide a "default" hash_set adaptor
__HS_TEMPLATE_HEADER
class hash_set : public __HS_SUPER 
{
  typedef hash_set< __HS_ARGUMENTS > _Self;
public:
  typedef  __HS_SUPER _Super;
  __IMPORT_WITH_ITERATORS(_Super)
  typedef typename _Super::key_type key_type;
  typedef typename _Super::hasher hasher;
  typedef typename _Super::key_equal key_equal;
  hash_set() {}
  hash_set(size_type n) : __HS_SUPER(n) {}
  hash_set(size_type n, const hasher& hf) : __HS_SUPER(n, hf) {}
  hash_set(size_type n, const hasher& hf, const key_equal& eql): __HS_SUPER(n, hf, eql) {}

  hash_set(const value_type* f, const value_type* l) : __HS_SUPER(f,l) {}
  hash_set(const value_type* f, const value_type* l, size_type n): __HS_SUPER(f,l,n) {}
  hash_set(const value_type* f, const value_type* l, size_type n,
           const hasher& hf) : __HS_SUPER(f,l,n,hf) {}
  hash_set(const value_type* f, const value_type* l, size_type n,
           const hasher& hf, const key_equal& eql) : __HS_SUPER(f,l,n,hf, eql) {}

  hash_set(const_iterator f, const_iterator l) : __HS_SUPER(f,l) { }
  hash_set(const_iterator f, const_iterator l, size_type n) : __HS_SUPER(f,l,n) { }
  hash_set(const_iterator f, const_iterator l, size_type n,
           const hasher& hf) : __HS_SUPER(f, l, n, hf) { }
  hash_set(const_iterator f, const_iterator l, size_type n,
           const hasher& hf, const key_equal& eql) : __HS_SUPER(f, l, n, hf, eql) { }
# if defined (__STL_BASE_MATCH_BUG)
    friend inline bool operator== __STL_NULL_TMPL_ARGS (const _Self& hs1, const _Self& hs2);
# endif
};

# if defined (__STL_BASE_MATCH_BUG)
__HS_TEMPLATE_HEADER
inline bool operator==(const hash_set< __HS_ARGUMENTS >& hs1, 
                       const hash_set< __HS_ARGUMENTS >& hs2)
{
    typedef __HS_SUPER _Super;
    return (const _Super&)hs1 == (const _Super&)hs2; 
}
# endif

// provide a "default" hash_multiset adaptor
__HS_TEMPLATE_HEADER
class hash_multiset : public __HMS_SUPER
{
  typedef hash_multiset< __HS_ARGUMENTS > _Self;
public:
  typedef __HMS_SUPER _Super;
  __IMPORT_WITH_ITERATORS(_Super)
  typedef typename _Super::key_type key_type;
  typedef typename _Super::hasher hasher;
  typedef typename _Super::key_equal key_equal;

  hash_multiset() {}
  hash_multiset(size_type __n) : __HMS_SUPER(__n) {}
  hash_multiset(size_type __n, const hasher& __hf) : __HMS_SUPER(__n, __hf) {}
  hash_multiset(size_type __n, const hasher& __hf, const key_equal& __eql): __HMS_SUPER(__n, __hf, __eql) {}

  hash_multiset(const value_type* __f, const value_type* __l) : __HMS_SUPER(__f,__l) {}
  hash_multiset(const value_type* __f, const value_type* __l, size_type __n): __HMS_SUPER(__f,__l,__n) {}
  hash_multiset(const value_type* __f, const value_type* __l, size_type __n,
           const hasher& __hf) : __HMS_SUPER(__f,__l,__n,__hf) {}
  hash_multiset(const value_type* __f, const value_type* __l, size_type __n,
           const hasher& __hf, const key_equal& __eql) : __HMS_SUPER(__f,__l,__n,__hf, __eql) {}

  hash_multiset(const_iterator __f, const_iterator __l) : __HMS_SUPER(__f,__l) { }
  hash_multiset(const_iterator __f, const_iterator __l, size_type __n) : __HMS_SUPER(__f,__l,__n) { }
  hash_multiset(const_iterator __f, const_iterator __l, size_type __n,
           const hasher& __hf) : __HMS_SUPER(__f, __l, __n, __hf) { }
  hash_multiset(const_iterator __f, const_iterator __l, size_type __n,
           const hasher& __hf, const key_equal& __eql) : __HMS_SUPER(__f, __l, __n, __hf, __eql) { }
# if defined (__STL_BASE_MATCH_BUG)
  friend inline bool operator== __STL_NULL_TMPL_ARGS (const _Self& __hs1, const _Self& __hs2);
# endif
};

# if defined (__STL_BASE_MATCH_BUG)
__HS_TEMPLATE_HEADER
inline bool operator==(const hash_multiset< __HS_ARGUMENTS >& __hs1, 
                       const hash_multiset< __HS_ARGUMENTS >& __hs2)
{
    typedef __HMS_SUPER  __s;
    return __STLPORT_STD::operator==((const __s&)__hs1,(const __s&)__hs2);
}
# endif

# undef __HS_SUPER
# undef __HMS_SUPER
# undef __HS_ARGUMENTS
# undef __HS_BASE_ARGUMENTS
# undef __HS_TEMPLATE_HEADER

# endif /*  WRAPPER */

#if defined(__sgi) && !defined(__GNUC__) && (_MIPS_SIM != _MIPS_SIM_ABI32)
#pragma reset woff 1174
#pragma reset woff 1375
#endif

__STL_END_NAMESPACE

#endif /* __SGI_STL_INTERNAL_HASH_SET_H */

// Local Variables:
// mode:C++
// End:
