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

#ifndef __SGI_STL_INTERNAL_HASHTABLE_H
#define __SGI_STL_INTERNAL_HASHTABLE_H

# ifndef __SGI_STL_INTERNAL_VECTOR_H
#  include <stl_vector.h>
# endif

# ifndef __SGI_STL_INTERNAL_ITERATOR_H
#  include <stl_iterator.h>
# endif

# ifndef __SGI_STL_INTERNAL_FUNCTION_H
#  include <stl_function.h>
# endif

# ifndef __SGI_STL_INTERNAL_ALGO_H
#  include <stl_algo.h>
# endif

# ifndef __SGI_STL_HASH_FUN_H
#  include <stl_hash_fun.h>
# endif

// Hashtable class, used to implement the hashed associative containers
// hash_set, hash_map, hash_multiset, and hash_multimap.


__STL_BEGIN_NAMESPACE

# if defined ( __STL_USE_ABBREVS )
#  define _Hashtable_iterator         _hT__It
#  define _Hashtable_const_iterator   _hT__cIt
#  define _Hashtable_node             _hT__N
#  define _Hashtable_base             _hT__B
#  define hashtable                   _h__T
#  define _Ht_iterator _Ht_It
# endif


template <class _Val>
struct _Hashtable_node
{
  typedef _Hashtable_node<_Val> _Self;
  _Self* _M_next;
  _Val _M_val;
  __TRIVIAL_STUFF(_Hashtable_node)
};  

// some compilers require the names of template parameters to be the same
template <class _Val, class _Key, class _HF,
          class _ExK, class _EqK, class _All>
class hashtable;

template <class _Val, class _Key, class _HF,
          class _ExK, class _EqK, class _All>
struct _Hashtable_iterator
# if defined ( __STL_DEBUG )
    : public __owned_link 
# endif
{
  typedef hashtable<_Val,_Key,_HF,_ExK,_EqK,_All>
          _Hashtable;
  typedef _Hashtable_node<_Val> _Node;

  typedef forward_iterator_tag iterator_category;
  typedef _Val value_type;
  typedef ptrdiff_t difference_type;
  typedef size_t size_type;
  typedef _Val& reference;
  typedef _Val* pointer;

  _Node* _M_cur;
  _Hashtable* _M_ht;

# if defined ( __STL_DEBUG )
  _Hashtable_iterator(_Node* __n, _Hashtable* __tab) : 
      __owned_link(__tab), _M_cur(__n), _M_ht(__tab) {}
  _Hashtable_iterator() : __owned_link(0) {}
# else
  _Hashtable_iterator(_Node* __n, _Hashtable* __tab) 
    : _M_cur(__n), _M_ht(__tab) {}
  _Hashtable_iterator() {}
# endif

  _Node* _M_skip_to_next();
};


template <class _Val, class _Traits, class _Key, class _HF,
          class _ExK, class _EqK, class _All>
struct _Ht_iterator : public _Hashtable_iterator< _Val, _Key,_HF, _ExK,_EqK,_All>
{
  
  typedef _Hashtable_iterator<_Val,_Key,_HF,_ExK,_EqK,_All> _Base;

  typedef _Ht_iterator<_Val, _Nonconst_traits<_Val>,_Key,_HF,_ExK,_EqK,_All> iterator;
  typedef _Ht_iterator<_Val, _Const_traits<_Val>,_Key,_HF,_ExK,_EqK,_All> const_iterator;
  typedef _Ht_iterator<_Val, _Traits,_Key,_HF,_ExK,_EqK,_All> _Self;

  typedef hashtable<_Val,_Key,_HF,_ExK,_EqK,_All> _Hashtable;
  typedef _Hashtable_node<_Val> _Node;

  typedef _Val value_type;
  typedef forward_iterator_tag iterator_category;
  typedef ptrdiff_t difference_type;
  typedef size_t size_type;
  typedef typename _Traits::reference reference;
  typedef typename _Traits::pointer   pointer;

# ifdef __STL_HAS_NAMESPACES
  __STL_USING_BASE_MEMBER _Hashtable_iterator<_Val,_Key,_HF,_ExK,_EqK,_All>::_M_skip_to_next;  
  __STL_USING_BASE_MEMBER _Hashtable_iterator<_Val,_Key,_HF,_ExK,_EqK,_All>::_M_cur;  
  __STL_USING_BASE_MEMBER _Hashtable_iterator<_Val,_Key,_HF,_ExK,_EqK,_All>::_M_ht;  
# endif

  _Ht_iterator(const _Node* __n, const _Hashtable* __tab) :
    _Hashtable_iterator<_Val,_Key,_HF,_ExK,_EqK,_All>((_Node*)__n, (_Hashtable*)__tab) {}
  _Ht_iterator() {}
  _Ht_iterator(const iterator& __it) : _Hashtable_iterator<_Val,_Key,_HF,_ExK,_EqK,_All>(__it) {}

  reference operator*() const { 
      __stl_verbose_assert(_Valid() && _M_cur!=0,_StlMsg_NOT_DEREFERENCEABLE);
      return _M_cur->_M_val; 
  }
  __STL_DEFINE_ARROW_OPERATOR

  _Self& operator++() {
    _Node* __n = _M_cur->_M_next;
    _M_cur =  (__n !=0 ? __n : _M_skip_to_next());
    return *this;
  }
  inline  _Self operator++(int) {
     _Self __tmp = *this;
    ++*this;
    return __tmp;
  }
};

template <class _Val, class _Traits, class _Traits1, class _Key, class _HF,
          class _ExK, class _EqK, class _All>
inline bool 
operator==(const _Ht_iterator<_Val, _Traits,_Key,_HF,_ExK,_EqK,_All>& __x, 
	   const _Ht_iterator<_Val, _Traits1,_Key,_HF,_ExK,_EqK,_All>& __y) { 
  __stl_debug_check(__check_same_owner_or_null(__x,__y));
  return __x._M_cur == __y._M_cur; 
}

#ifdef __STL_USE_SEPARATE_RELOPS_NAMESPACE
template <class _Val, class _Key, class _HF,
          class _ExK, class _EqK, class _All>
inline bool 
operator!=(const _Hashtable_iterator<_Val,_Key,_HF,_ExK,_EqK,_All>& __x, 
	   const _Hashtable_iterator<_Val,_Key,_HF,_ExK,_EqK,_All>& __y) { 
  __stl_debug_check(__check_same_owner_or_null(__x,__y));
  return __x._M_cur != __y._M_cur; 
}
#else
template <class _Val, class _Key, class _HF,
          class _ExK, class _EqK, class _All>
inline bool 
operator!=(const _Ht_iterator<_Val, _Nonconst_traits<_Val>,_Key,_HF,_ExK,_EqK,_All>& __x, 
	   const _Ht_iterator<_Val, _Const_traits<_Val>,_Key,_HF,_ExK,_EqK,_All>& __y) { 
  __stl_debug_check(__check_same_owner_or_null(__x,__y));
  return __x._M_cur != __y._M_cur; 
}
#endif

#ifndef __STL_CLASS_PARTIAL_SPECIALIZATION

template <class _Val, class _Traits, class _Key, class _HF, class _ExK, class _EqK, class _All>
inline _Val*
value_type(const _Ht_iterator<_Val, _Traits,_Key,_HF,_ExK,_EqK,_All>&)
{
  return (_Val*) 0;
}

template <class _Val, class _Traits, class _Key, class _HF, class _ExK, class _EqK, class _All>
inline forward_iterator_tag
iterator_category(const _Ht_iterator<_Val, _Traits,_Key,_HF,_ExK,_EqK,_All>&)
{
  return forward_iterator_tag();
}

template <class _Val, class _Traits, class _Key, class _HF, class _ExK, class _EqK, 
          class _All>
inline ptrdiff_t*
distance_type(const _Ht_iterator<_Val,_Traits,_Key,_HF,_ExK,_EqK,_All>&)
{
  return (ptrdiff_t*) 0;
}

#endif /* __STL_CLASS_PARTIAL_SPECIALIZATION */

// Note: assumes long is at least 32 bits.
# define __stl_num_primes  28

# if ( __STL_STATIC_TEMPLATE_DATA > 0 ) /* && ! defined (__GNUC__) */
#  define __stl_prime_list _Stl_prime<const unsigned long>::_M_list
   template <class _Tp>
   struct _Stl_prime {
   public:
       static _Tp _M_list[__stl_num_primes];
   };
#  else
#  if ( __STL_WEAK_ATTRIBUTE > 0 )
      extern const unsigned long __stl_prime_list[__stl_num_primes] __attribute__((weak));
#  else
      // give up
      static const unsigned long __stl_prime_list[__stl_num_primes];
#  endif /* __STL_WEAK_ATTRIBUTE */
#endif /* __STL_STATIC_TEMPLATE_DATA */


// Hashtables handle allocators a bit differently than other containers
//  do.  If we're using standard-conforming allocators, then a hashtable
//  unconditionally has a member variable to hold its allocator, even if
//  it so happens that all instances of the allocator type are identical.
// This is because, for hashtables, this extra storage is negligible.  
//  Additionally, a base class wouldn't serve any other purposes; it 
//  wouldn't, for example, simplify the exception-handling code.
template <class _Val, class _Key, class _HF,
          class _ExK, class _EqK, class _All>
# if defined ( __STL_DEBUG )
class hashtable : public __owned_list {
# else
class hashtable {
# endif
  typedef hashtable<_Val, _Key, _HF, _ExK, _EqK, _All> _Self;
public:
  typedef _Key key_type;
  typedef _Val value_type;
  typedef _HF hasher;
  typedef _EqK key_equal;

  typedef size_t            size_type;
  typedef ptrdiff_t         difference_type;
  typedef value_type*       pointer;
  typedef const value_type* const_pointer;
  typedef value_type&       reference;
  typedef const value_type& const_reference;

  hasher hash_funct() const { return _M_hash; }
  key_equal key_eq() const { return _M_equals; }

private:
  typedef _Hashtable_node<_Val> _Node;

# define __HASH_ALLOC_PARAM   allocator_type

private:
  typedef typename _Alloc_traits<_Node, _All>::allocator_type _M_node_allocator_type;
  typedef typename _Alloc_traits<void*, _All>::allocator_type _M_node_ptr_allocator_type;
//  typedef __vector__<void*, _M_node_ptr_allocator_type> _BucketVector;
  typedef __original_vector<void*, _M_node_ptr_allocator_type> _BucketVector;
public:
  typedef typename _Alloc_traits<_Val,_All>::allocator_type allocator_type;
  allocator_type get_allocator() const { 
    return __STL_CONVERT_ALLOCATOR((const _M_node_allocator_type&)_M_num_elements, _Val); 
  }
private:
  hasher                _M_hash;
  key_equal             _M_equals;
  _ExK                  _M_get_key;
  _BucketVector         _M_buckets;
  _STL_alloc_proxy<size_type, _Node, _M_node_allocator_type>  _M_num_elements;
  const _Node* _M_get_bucket(size_t __n) const { return (_Node*)_M_buckets[__n]; }

public:
  typedef _Const_traits<_Val> __const_val_traits;
  typedef _Nonconst_traits<_Val> __nonconst_val_traits;
  typedef _Ht_iterator<_Val, __const_val_traits,_Key,_HF,_ExK,_EqK, _All> const_iterator;
  typedef _Ht_iterator<_Val, __nonconst_val_traits,_Key,_HF,_ExK,_EqK,_All> iterator;
  friend struct _Hashtable_iterator<_Val,_Key,_HF,_ExK,_EqK,_All>;
  friend struct _Ht_iterator<_Val, _Nonconst_traits<_Val>,_Key,_HF,_ExK,_EqK,_All>;
  friend struct _Ht_iterator<_Val, _Const_traits<_Val>,_Key,_HF,_ExK,_EqK, _All>;

public:
  hashtable(size_type __n,
            const _HF&  __hf,
            const _EqK& __eql,
            const _ExK& __ext,
            const allocator_type& __a = __STL_ALLOC_INSTANCE(allocator_type))
    :
      _M_hash(__hf),
      _M_equals(__eql),
      _M_get_key(__ext),
      _M_buckets(__STL_CONVERT_ALLOCATOR(__a,void*)),
      _M_num_elements(__STL_CONVERT_ALLOCATOR(__a,_Node), (size_type)0)
  {
    _M_initialize_buckets(__n);
  }

  hashtable(size_type __n,
            const _HF&    __hf,
            const _EqK&   __eql,
            const allocator_type& __a = __STL_ALLOC_INSTANCE(allocator_type))
    :
      _M_hash(__hf),
      _M_equals(__eql),
      _M_get_key(_ExK()),
      _M_buckets(__STL_CONVERT_ALLOCATOR(__a,void*)),
      _M_num_elements(__STL_CONVERT_ALLOCATOR(__a,_Node), (size_type)0)
  {
    _M_initialize_buckets(__n);
  }

  hashtable(const _Self& __ht)
    :
      _M_hash(__ht._M_hash),
      _M_equals(__ht._M_equals),
      _M_get_key(__ht._M_get_key),
      _M_buckets(__STL_CONVERT_ALLOCATOR(__ht.get_allocator(),void*)),
      _M_num_elements((const _M_node_allocator_type&)__ht._M_num_elements, (size_type)0)
  {
    __stl_debug_do(_Safe_init(this));
    _M_copy_from(__ht);
  }

  _Self& operator= (const _Self& __ht)
  {
    if (&__ht != this) {
      clear();
      _M_hash = __ht._M_hash;
      _M_equals = __ht._M_equals;
      _M_get_key = __ht._M_get_key;
      _M_copy_from(__ht);
    }
    return *this;
  }

  ~hashtable() { clear(); }

  size_type size() const { return _M_num_elements._M_data; }
  size_type max_size() const { return size_type(-1); }
  bool empty() const { return size() == 0; }

  void swap(_Self& __ht)
  {
    __STLPORT_STD::swap(_M_hash, __ht._M_hash);
    __STLPORT_STD::swap(_M_equals, __ht._M_equals);
    __STLPORT_STD::swap(_M_get_key, __ht._M_get_key);
    _M_buckets.swap(__ht._M_buckets);
    __STLPORT_STD::swap(_M_num_elements, __ht._M_num_elements);
    __stl_debug_do(_Swap_owners(__ht));
  }

  iterator begin()
  { 
    for (size_type __n = 0; __n < _M_buckets.size(); ++__n)
      if (_M_buckets[__n])
        return iterator((_Node*)_M_buckets[__n], this);
    return end();
  }

  iterator end() { return iterator((_Node*)0, this); }

  const_iterator begin() const
  {
    for (size_type __n = 0; __n < _M_buckets.size(); ++__n)
      if (_M_buckets[__n])
        return const_iterator((_Node*)_M_buckets[__n], this);
    return end();
  }

  const_iterator end() const { return const_iterator((_Node*)0, this); }

  static bool _M_equal (const hashtable<_Val, _Key, _HF, _ExK, _EqK, _All>&,
			const hashtable<_Val, _Key, _HF, _ExK, _EqK, _All>&);

public:

  size_type bucket_count() const { return _M_buckets.size(); }

  size_type max_bucket_count() const
    { return __stl_prime_list[(int)__stl_num_primes - 1]; } 

  size_type elems_in_bucket(size_type __bucket) const
  {
    size_type __result = 0;
    for (_Node* __cur = (_Node*)_M_buckets[__bucket]; __cur; __cur = __cur->_M_next)
      __result += 1;
    return __result;
  }

  pair<iterator, bool> insert_unique(const value_type& __obj)
  {
    resize(_M_num_elements._M_data + 1);
    return insert_unique_noresize(__obj);
  }

  iterator insert_equal(const value_type& __obj)
  {
    resize(_M_num_elements._M_data + 1);
    return insert_equal_noresize(__obj);
  }

  pair<iterator, bool> insert_unique_noresize(const value_type& __obj);
  iterator insert_equal_noresize(const value_type& __obj);
 
#ifdef __STL_MEMBER_TEMPLATES
  template <class _InputIterator>
  void insert_unique(_InputIterator __f, _InputIterator __l)
  {
    insert_unique(__f, __l, __ITERATOR_CATEGORY(__f));
  }

  template <class _InputIterator>
  void insert_equal(_InputIterator __f, _InputIterator __l)
  {
    insert_equal(__f, __l, __ITERATOR_CATEGORY(__f));
  }

  template <class _InputIterator>
  void insert_unique(_InputIterator __f, _InputIterator __l,
                     input_iterator_tag)
  {
    for ( ; __f != __l; ++__f)
      insert_unique(*__f);
  }

  template <class _InputIterator>
  void insert_equal(_InputIterator __f, _InputIterator __l,
                    input_iterator_tag)
  {
    for ( ; __f != __l; ++__f)
      insert_equal(*__f);
  }

  template <class _ForwardIterator>
  void insert_unique(_ForwardIterator __f, _ForwardIterator __l,
                     forward_iterator_tag)
  {
    size_type __n = 0;
    distance(__f, __l, __n);
    resize(_M_num_elements._M_data + __n);
    for ( ; __n > 0; --__n, ++__f)
      insert_unique_noresize(*__f);
  }

  template <class _ForwardIterator>
  void insert_equal(_ForwardIterator __f, _ForwardIterator __l,
                    forward_iterator_tag)
  {
    size_type __n = 0;
    distance(__f, __l, __n);
    resize(_M_num_elements._M_data + __n);
    for ( ; __n > 0; --__n, ++__f)
      insert_equal_noresize(*__f);
  }

#else /* __STL_MEMBER_TEMPLATES */
  void insert_unique(const value_type* __f, const value_type* __l)
  {
    size_type __n = __l - __f;
    resize(_M_num_elements._M_data + __n);
    for ( ; __n > 0; --__n, ++__f)
      insert_unique_noresize(*__f);
  }

  void insert_equal(const value_type* __f, const value_type* __l)
  {
    size_type __n = __l - __f;
    resize(_M_num_elements._M_data + __n);
    for ( ; __n > 0; --__n, ++__f)
      insert_equal_noresize(*__f);
  }

  void insert_unique(const_iterator __f, const_iterator __l)
  {
    size_type __n = 0;
    distance(__f, __l, __n);
    resize(_M_num_elements._M_data + __n);
    for ( ; __n > 0; --__n, ++__f)
      insert_unique_noresize(*__f);
  }

  void insert_equal(const_iterator __f, const_iterator __l)
  {
    size_type __n = 0;
    distance(__f, __l, __n);
    resize(_M_num_elements._M_data + __n);
    for ( ; __n > 0; --__n, ++__f)
      insert_equal_noresize(*__f);
  }
#endif /*__STL_MEMBER_TEMPLATES */

  reference find_or_insert(const value_type& __obj);

  iterator find(const key_type& __key) 
  {
    size_type __n = _M_bkt_num_key(__key);
    _Node* __first;
    for ( __first = (_Node*)_M_buckets[__n];
          __first && !_M_equals(_M_get_key(__first->_M_val), __key);
          __first = __first->_M_next)
      {}
    return iterator(__first, this);
  } 

  const_iterator find(const key_type& __key) const
  {
    size_type __n = _M_bkt_num_key(__key);
    const _Node* __first;
    for ( __first = (_Node*)_M_buckets[__n];
          __first && !_M_equals(_M_get_key(__first->_M_val), __key);
          __first = __first->_M_next)
      {}
    return const_iterator(__first, this);
  } 

  size_type count(const key_type& __key) const
  {
    const size_type __n = _M_bkt_num_key(__key);
    size_type __result = 0;

    for (const _Node* __cur = (_Node*)_M_buckets[__n]; __cur; __cur = __cur->_M_next)
      if (_M_equals(_M_get_key(__cur->_M_val), __key))
        ++__result;
    return __result;
  }

  pair<iterator, iterator> 
  equal_range(const key_type& __key);

  pair<const_iterator, const_iterator> 
  equal_range(const key_type& __key) const;

  size_type erase(const key_type& __key);
  //   void erase(const iterator& __it); `
  void erase(const const_iterator& __it) ;

  //  void erase(const const_iterator& __first, const const_iterator __last) {
  //     erase((const iterator&)__first, (const iterator&)__last);
  //  }
  void erase(const_iterator __first, const_iterator __last);
  void resize(size_type __num_elements_hint);
  void clear();
# if defined (__STL_DEBUG)
    void _Invalidate_node(_Node* __it) { 
      __invalidate_iterator((__owned_list*)this,iterator(__it, this));
    }
# endif

private:

  size_type _M_next_size(size_type __n) const
    { 
      const size_type* __first = (const size_type*)__stl_prime_list;
      const size_type* __last =  (const size_type*)__stl_prime_list + (int)__stl_num_primes;
      const size_type* pos = lower_bound(__first, __last, __n);
      return (pos == __last ? *(__last - 1) : *pos);
    }

  void _M_initialize_buckets(size_type __n)
  {
    const size_type __n_buckets = _M_next_size(__n);
    _M_buckets.reserve(__n_buckets);
    _M_buckets.insert(_M_buckets.end(), __n_buckets, (void*) 0);
    _M_num_elements._M_data = 0;
    __stl_debug_do(_Safe_init(this));
  }

  size_type _M_bkt_num_key(const key_type& __key) const
  {
    return _M_bkt_num_key(__key, _M_buckets.size());
  }

  size_type _M_bkt_num(const value_type& __obj) const
  {
    return _M_bkt_num_key(_M_get_key(__obj));
  }

  size_type _M_bkt_num_key(const key_type& __key, size_t __n) const
  {
    return _M_hash(__key) % __n;
  }

  size_type _M_bkt_num(const value_type& __obj, size_t __n) const
  {
    return _M_bkt_num_key(_M_get_key(__obj), __n);
  }

  _Node* _M_new_node(const value_type& __obj)
  {
    _Node* __n = _M_num_elements.allocate(1);
    __n->_M_next = 0;
    __STL_TRY {
      construct(&__n->_M_val, __obj);
      //      return __n;
    }
    __STL_UNWIND(_M_num_elements.deallocate(__n, 1));
    return __n;
  }
  
  void _M_delete_node(_Node* __n)
  {
    __stl_debug_do(_Invalidate_node(__n));
    destroy(&__n->_M_val);
    _M_num_elements.deallocate(__n, 1);
  }

  void _M_erase_bucket(const size_type __n, _Node* __first, _Node* __last);
  void _M_erase_bucket(const size_type __n, _Node* __last);

  void _M_copy_from(const _Self& __ht);
};

template <class _Val, class _Key, class _HF, class _ExK, class _EqK, class _All>
inline bool operator==(const hashtable<_Val,_Key,_HF,_ExK,_EqK,_All>& __ht1,
                       const hashtable<_Val,_Key,_HF,_ExK,_EqK,_All>& __ht2)
{
  return hashtable<_Val,_Key,_HF,_ExK,_EqK,_All>::_M_equal( __ht1, __ht2 );
}

#ifdef __STL_USE_SEPARATE_RELOPS_NAMESPACE

template <class _Val, class _Key, class _HF, class _Ex, class _Eq, class _All>
inline bool operator!=(const hashtable<_Val,_Key,_HF,_Ex,_Eq,_All>& __ht1,
                       const hashtable<_Val,_Key,_HF,_Ex,_Eq,_All>& __ht2) {
  return !(__ht1 == __ht2);
}

template <class _Val, class _Key, class _HF, class _ExK, class _EqK, 
          class _All>
inline void swap(hashtable<_Val, _Key, _HF, _ExK, _EqK, _All>& __ht1,
                 hashtable<_Val, _Key, _HF, _ExK, _EqK, _All>& __ht2) {
  __ht1.swap(__ht2);
}

#endif /* __STL_USE_SEPARATE_RELOPS_NAMESPACE */

__STL_END_NAMESPACE

# undef __stl_prime_list

# if !defined (__STL_LINK_TIME_INSTANTIATION)
#  include <stl_hashtable.c>
# endif

#endif /* __SGI_STL_INTERNAL_HASHTABLE_H */

// Local Variables:
// mode:C++
// End:


