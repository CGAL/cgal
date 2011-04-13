/*
 * Copyright (c) 1997-1999
 * Silicon Graphics Computer Systems, Inc.
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

#ifndef __SGI_STL_STRING_H
#define __SGI_STL_STRING_H

# ifndef __STL_CONFIG_H
#  include <stl_config.h>
# endif

#if !defined (__STLPORT_DEBUG_H) && (defined  (__STL_DEBUG) || defined (__STL_ASSERTIONS))
# include <stldebug.h>
#endif

# if defined (__STL_DEBUG)
#  define _Make_iterator(__i) iterator(&_M_iter_list, __i)
#  define _Make_const_iterator(__i) const_iterator(&_M_iter_list, __i)
#  define _Make_ptr(__i)   (__i)._M_iterator
# else
#  define _Make_iterator(__i) __i
#  define _Make_const_iterator(__i) __i
#  define _Make_ptr(__i)   __i
# endif

# ifndef __STLPORT_CCTYPE
#  include <cctype>
# endif

#ifndef __SGI_STL_STRING_FWD_H
#  include <stl_string_fwd.h>
#endif

#ifndef __SGI_STL_INTERNAL_FUNCTION_H
# include <stl_function.h>
#endif

# include <stl_ctraits_fns.h>

//#ifndef __SGI_STDEXCEPT
//# include <stdexcept>      
//#endif

#ifndef __SGI_STL_MEMORY
# include <memory>
#endif

#ifndef __SGI_STL_INTERNAL_ALGO_H
# include <stl_algo.h>
#endif

# ifndef __STLPORT_IOSFWD
#  include <iosfwd>
# endif

#if defined (__STL_DEBUG) && ! defined (__STLPORT_VEC_ITERATOR_H)
// string uses the same debug iterator as vector
#  include <stl_vec_iterator.h>
#endif

#if defined( __MWERKS__ ) && ! defined (__STL_USE_OWN_NAMESPACE)

// MSL implementation classes expect to see the definition of streampos
// when this header is included. We expect this to be fixed in later MSL
// implementations
# if !defined( __MSL_CPP__ ) || __MSL_CPP__ < 0x4105
#  include <msl_string.h>
# endif

#endif // __MWERKS__

// Standard C++ string class.  This class has performance
// characteristics very much like vector<>, meaning, for example, that
// it does not perform reference-count or copy-on-write, and that
// concatenation of two strings is an O(N) operation. 

// There are three reasons why basic_string is not identical to
// vector.  First, basic_string always stores a null character at the
// end; this makes it possible for c_str to be a fast operation.
// Second, the C++ standard requires basic_string to copy elements
// using char_traits<>::assign, char_traits<>::copy, and
// char_traits<>::move.  This means that all of vector<>'s low-level
// operations must be rewritten.  Third, basic_string<> has a lot of
// extra functions in its interface that are convenient but, strictly
// speaking, redundant.

// Additionally, the C++ standard imposes a major restriction: according
// to the standard, the character type _CharT must be a POD type.  This
// implementation weakens that restriction, and allows _CharT to be a
// a user-defined non-POD type.  However, _CharT must still have a
// default constructor.

__STL_BEGIN_NAMESPACE

#if defined(__sgi) && !defined(__GNUC__) && (_MIPS_SIM != _MIPS_SIM_ABI32)
#pragma set woff 1174
#pragma set woff 1375
#endif

// A helper class to use a char_traits as a function object.

template <class _Traits>
struct _Not_within_traits
  : public unary_function<typename _Traits::char_type, bool>
{
  typedef typename _Traits::char_type _CharT;
  const _CharT* _M_first;
  const _CharT* _M_last;

  _Not_within_traits(const typename _Traits::char_type* __f, 
		     const typename _Traits::char_type* __l) 
    : _M_first(__f), _M_last(__l) {}

  bool operator()(const typename _Traits::char_type& __x) const {
    return find_if((_CharT*)_M_first, (_CharT*)_M_last, 
                   bind1st(_Eq_traits<_Traits>(), __x)) == (_CharT*)_M_last;
  }
};

// ------------------------------------------------------------
// Class _String_base.  

// _String_base is a helper class that makes it it easier to write an
// exception-safe version of basic_string.  The constructor allocates,
// but does not initialize, a block of memory.  The destructor
// deallocates, but does not destroy elements within, a block of
// memory.  The destructor assumes that _M_start either is null, or else
// points to a block of memory that was allocated using _String_base's 
// allocator and whose size is _M_end_of_storage._M_data - _M_start.

// Additionally, _String_base encapsulates the difference between
// old SGI-style allocators and standard-conforming allocators.

template <class _Tp, class _Alloc> class _String_base {
public:
  typedef typename _Alloc_traits<_Tp, _Alloc>::allocator_type allocator_type;
  _Tp*    _M_start;
  _Tp*    _M_finish;
  _STL_alloc_proxy<_Tp*, _Tp, allocator_type> _M_end_of_storage;
                                // Precondition: 0 < __n <= max_size().
  void _M_allocate_block(size_t __n) { 
    if (__n <= max_size()) {
      _M_start  = _M_end_of_storage.allocate(__n);
      _M_finish = _M_start;
      _M_end_of_storage._M_data = _M_start + __n;
    }
    else
      _M_throw_length_error();
  }

  void _M_deallocate_block() 
    { _M_end_of_storage.deallocate(_M_start, _M_end_of_storage._M_data - _M_start); }
  
  size_t max_size() const { return (size_t(-1) / sizeof(_Tp)) - 1; }

  _String_base(const allocator_type& __a)
    : _M_start(0), _M_finish(0), _M_end_of_storage(__a, (_Tp*)0) { }
  
  _String_base(const allocator_type& __a, size_t __n)
    : _M_start(0), _M_finish(0), _M_end_of_storage(__a, (_Tp*)0)
    { _M_allocate_block(__n); }

  ~_String_base() { _M_deallocate_block(); }

  void _M_throw_length_error() const;
  void _M_throw_out_of_range() const;
};

// ------------------------------------------------------------
// Class basic_string.  

// Class invariants:
// (1) [start, finish) is a valid range.
// (2) Each iterator in [start, finish) points to a valid object
//     of type value_type.
// (3) *finish is a valid object of type value_type; in particular,
//     it is value_type().
// (4) [finish + 1, end_of_storage) is a valid range.
// (5) Each iterator in [finish + 1, end_of_storage) points to 
//     unininitialized memory.

// Note one important consequence: a string of length n must manage
// a block of memory whose size is at least n + 1.  

struct _String_reserve_t {};

template <class _CharT, class _Traits, class _Alloc> 
class basic_string : protected _String_base<_CharT,_Alloc> {
private:                        // Protected members inherited from base.
  typedef _String_base<_CharT,_Alloc> _Base;
  typedef basic_string<_CharT, _Traits, _Alloc> _Self;
  // fbp : used to optimize char/wchar_t cases, and to simplify
  // __STL_DEFAULT_CONSTRUCTOR_BUG problem workaround
  typedef typename _Is_integer<_CharT>::_Integral _Char_Is_Integral;
public:
#if defined( __STL_HAS_NAMESPACES )
  __STL_USING_BASE_MEMBER _String_base<_CharT,_Alloc>::_M_allocate_block;
# ifndef __STL_DEBUG
  __STL_USING_BASE_MEMBER _String_base<_CharT,_Alloc>::_M_deallocate_block;
# endif
  __STL_USING_BASE_MEMBER _String_base<_CharT,_Alloc>::_M_throw_length_error;
  __STL_USING_BASE_MEMBER _String_base<_CharT,_Alloc>::_M_throw_out_of_range;

  __STL_USING_BASE_MEMBER _String_base<_CharT,_Alloc>::_M_start;
  __STL_USING_BASE_MEMBER _String_base<_CharT,_Alloc>::_M_finish;
  __STL_USING_BASE_MEMBER _String_base<_CharT,_Alloc>::_M_end_of_storage;
#endif /* __STL_HAS_NAMESPACES */
public:
  typedef _CharT value_type;
  typedef _Traits traits_type;

  typedef value_type* pointer;
  typedef const value_type* const_pointer;
  typedef value_type& reference;
  typedef const value_type& const_reference;
  typedef size_t size_type;
  typedef ptrdiff_t difference_type;

# if defined (__STL_DEBUG)
  typedef _Vec_iter<_CharT, _Nonconst_traits<_CharT> > iterator;
  typedef _Vec_iter<_CharT, _Const_traits<_CharT> > const_iterator;
private:
  mutable __owned_list _M_iter_list;
public:
# else
  typedef const value_type*                const_iterator;
  typedef value_type*                      iterator;
# endif

#if (defined ( __STL_CLASS_PARTIAL_SPECIALIZATION ) && \
  ! defined (__STL_PARTIAL_SPECIALIZATION_BUG)) || \
  defined(CGAL_LIMITED_ITERATOR_TRAITS_SUPPORT)
  typedef __STLPORT_STD::reverse_iterator<const_iterator> const_reverse_iterator;
  typedef __STLPORT_STD::reverse_iterator<iterator> reverse_iterator;
#else /* __STL_CLASS_PARTIAL_SPECIALIZATION */
# if defined (__STL_MSVC50_COMPATIBILITY)
  typedef __STLPORT_STD::reverse_iterator<const_iterator, value_type, const_reference, 
    const_pointer, difference_type>  const_reverse_iterator;
  typedef __STLPORT_STD::reverse_iterator<iterator, value_type, reference, pointer, difference_type>
    reverse_iterator;
# else
  typedef __STLPORT_STD::reverse_iterator<const_iterator, value_type, const_reference, 
    difference_type>  const_reverse_iterator;
  typedef __STLPORT_STD::reverse_iterator<iterator, value_type, reference, difference_type>
    reverse_iterator;
# endif
#endif /* __STL_CLASS_PARTIAL_SPECIALIZATION */

  static const size_type npos;
  typedef _String_reserve_t _Reserve_t;
# ifdef __STL_USE_NATIVE_STRING
  // this typedef is being used for conversions
  typedef __STL_VENDOR_STD::basic_string<_CharT,_Traits, 
    __STL_VENDOR_STD::allocator<_CharT> >  __std_string;
# endif
  
public:                         // Constructor, destructor, assignment.
  typedef typename _Base::allocator_type allocator_type;
  allocator_type get_allocator() const {
    return __STL_CONVERT_ALLOCATOR((const allocator_type&)_M_end_of_storage, _CharT);
  }

  basic_string()
    : _String_base<_CharT,_Alloc>(allocator_type(), 8) { 
            __stl_debug_do(_M_iter_list._Safe_init(&_M_start));
	    _M_terminate_string(); 
  }

  explicit basic_string(const allocator_type& __a)
    : _String_base<_CharT,_Alloc>(__a, 8) { 
    __stl_debug_do(_M_iter_list._Safe_init(&_M_start));
    _M_terminate_string(); 
  }

  basic_string(_Reserve_t, size_t __n,
               const allocator_type& __a = __STL_ALLOC_INSTANCE(allocator_type))
    : _String_base<_CharT,_Alloc>(__a, __n + 1) { 
    __stl_debug_do(_M_iter_list._Safe_init(&_M_start));
    _M_terminate_string(); 
  }

  basic_string(const _Self& __s) : _String_base<_CharT,_Alloc>(__s.get_allocator()) 
    { 
      _M_range_initialize(__s._M_start, __s._M_finish); 
    }

  basic_string(const _Self& __s, size_type __pos, size_type __n = npos,
               const allocator_type& __a = __STL_ALLOC_INSTANCE(allocator_type)) 
    : _String_base<_CharT,_Alloc>(__a) {
    if (__pos > __s.size())
      _M_throw_out_of_range();
    else
      _M_range_initialize(__s._M_start + __pos,
                          __s._M_start + __pos + min(__n, __s.size() - __pos));
  }

  basic_string(const _CharT* __s, size_type __n,
               const allocator_type& __a = __STL_ALLOC_INSTANCE(allocator_type)) 
    : _String_base<_CharT,_Alloc>(__a) 
    { 
      __STL_FIX_LITERAL_BUG(__s)
      _M_range_initialize(__s, __s + __n); 
    }

  basic_string(const _CharT* __s,
               const allocator_type& __a = __STL_ALLOC_INSTANCE(allocator_type))
    : _String_base<_CharT,_Alloc>(__a) 
    {
      __STL_FIX_LITERAL_BUG(__s)
      _M_range_initialize(__s, __s + _Traits::length(__s)); 
    }

  basic_string(size_type __n, _CharT __c,
               const allocator_type& __a = __STL_ALLOC_INSTANCE(allocator_type))
    : _String_base<_CharT,_Alloc>(__a, __n + 1)
  {
    __stl_debug_do(_M_iter_list._Safe_init(&_M_start));
    _M_finish = uninitialized_fill_n(_M_start, __n, __c);
    _M_terminate_string();
  }
  // Check to see if _InputIterator is an integer type.  If so, then
  // it can't be an iterator.
#if defined (__STL_MEMBER_TEMPLATES) && ! defined (__STL_MSVC)
  template <class _InputIterator>
  basic_string(_InputIterator __f, _InputIterator __l,
               const allocator_type & __a = __STL_ALLOC_INSTANCE(allocator_type))
    : _String_base<_CharT,_Alloc>(__a)
  {
    typedef typename _Is_integer<_InputIterator>::_Integral _Integral;
    _M_initialize_dispatch(__f, __l, _Integral());
  }
#else /* __STL_MEMBER_TEMPLATES */

# ifdef __STL_DEBUG
  basic_string(const_iterator __f, const_iterator __l)
    : _String_base<_CharT,_Alloc>(allocator_type())
  {
    _M_range_initialize(_Make_ptr(__f), _Make_ptr(__l));
  }
  basic_string(const_iterator __f, const_iterator __l,
               const allocator_type& __a )
    : _String_base<_CharT,_Alloc>(__a)
  {
    _M_range_initialize(_Make_ptr(__f), _Make_ptr(__l));
  }
# endif

  basic_string(const _CharT* __f, const _CharT* __l,
               const allocator_type& __a = __STL_ALLOC_INSTANCE(allocator_type))
    : _String_base<_CharT,_Alloc>(__a)
  {
    __STL_FIX_LITERAL_BUG(__f)  __STL_FIX_LITERAL_BUG(__l)
    _M_range_initialize(__f, __l);
  }

#endif

# ifdef __STL_USE_NATIVE_STRING
  // these conversion operations still needed for
  // strstream, etc.
  basic_string (const __std_string& __x): _String_base<_CharT,_Alloc>(allocator_type())
    {
      const _CharT* __s = __x.data();
      _M_range_initialize(__s, __s + __x.size()); 
    }
  
  operator __std_string() const { return __std_string(this->data()); }
# endif

  ~basic_string() { destroy(_M_start, _M_finish + 1); }
    
  _Self& operator=(const _Self& __s) {
    if (&__s != this) 
      assign(__s._M_start, __s._M_finish);
    return *this;
  }

  _Self& operator=(const _CharT* __s) { 
    __STL_FIX_LITERAL_BUG(__s)
    return assign(__s, __s + _Traits::length(__s)); 
  }

  _Self& operator=(_CharT __c)
    { return assign(__STATIC_CAST(size_type,1), __c); }

  static _CharT _M_null() {
#   ifndef __STL_DEFAULT_CONSTRUCTOR_BUG
    return _CharT();
#   else
    return (_CharT) 0;
#   endif
  }

private:                        // Helper functions used by constructors
                                // and elsewhere.
  // fbp : simplify integer types (char, wchar)
  void _M_construct_null_aux(_CharT* __p, __false_type) {
    construct(__p);
  }
  void _M_construct_null_aux(_CharT* __p, __true_type) {
    *__p = 0;
  }

  void _M_construct_null(_CharT* __p) {
    _M_construct_null_aux(__p, _Char_Is_Integral());
  }

private:                        
  // Helper functions used by constructors.  It is a severe error for
  // any of them to be called anywhere except from within constructors.

  void _M_terminate_string_aux(__false_type) {
    __STL_TRY {
      _M_construct_null(_M_finish);
    }
    __STL_UNWIND(destroy(_M_start, _M_finish));
  }

  void _M_terminate_string_aux(__true_type) {
    *_M_finish=0;
  }

  void _M_terminate_string() {
    _M_terminate_string_aux(_Char_Is_Integral());
  }

#ifdef __STL_MEMBER_TEMPLATES
    
  template <class _InputIter>
  void _M_range_initialize(_InputIter __f, _InputIter __l,
                           input_iterator_tag) {
    _M_allocate_block(8);
    _M_construct_null(_M_finish);
    __STL_TRY {
      append(__f, __l);
    }
    __STL_UNWIND(destroy(_M_start, _M_finish + 1));
  }

  template <class _ForwardIter>
  void _M_range_initialize(_ForwardIter __f, _ForwardIter __l, 
                           forward_iterator_tag) {
    difference_type __n = 0;
    distance(__f, __l, __n);
    _M_allocate_block(__n + 1);
    _M_finish = uninitialized_copy(__f, __l, _M_start);
    _M_terminate_string();
  }

  template <class _InputIter>
  void _M_range_initialize(_InputIter __f, _InputIter __l) {
    __stl_debug_do(__check_range(__f, __l));
# if defined( __STL_CLASS_PARTIAL_SPECIALIZATION)
    typedef  __STLPORT_STD::iterator_traits<_InputIter> _Iter_traits;
    typedef typename _Iter_traits::iterator_category _Category;
    _M_range_initialize(__f, __l, _Category());
# else
    _M_range_initialize(__f, __l, __ITERATOR_CATEGORY(__f));
# endif
    __stl_debug_do(_M_iter_list._Safe_init(&_M_start));
  }

  template <class _Integer>
  void _M_initialize_dispatch(_Integer __n, _Integer __x, __true_type) {
    _M_allocate_block(__n + 1);
    _M_finish = uninitialized_fill_n(_M_start, __n, __x);
    _M_terminate_string();
    __stl_debug_do(_M_iter_list._Safe_init(&_M_start));
  }

  template <class _InputIter>
  void _M_initialize_dispatch(_InputIter __f, _InputIter __l, __false_type) {
     _M_range_initialize(__f, __l);
  }
    
#else /* __STL_MEMBER_TEMPLATES */

  void _M_range_initialize(const _CharT* __f, const _CharT* __l) {
    __stl_debug_do(__check_range(__f, __l));
    ptrdiff_t __n = __l - __f;
    _M_allocate_block(__n + 1);
    _M_finish = uninitialized_copy(__f, __l, _M_start);
    _M_terminate_string();
    __stl_debug_do(_M_iter_list._Safe_init(&_M_start));
  }

#endif /* __STL_MEMBER_TEMPLATES */

public:                         // Iterators.
# if defined (__STL_DEBUG)
  iterator begin() { return _Make_iterator(_M_start); }
  const_iterator begin() const { return _Make_const_iterator(_M_start); }
  iterator end() { return _Make_iterator(_M_finish); }
  const_iterator end() const { return _Make_const_iterator(_M_finish); }
  void _M_deallocate_block() {
    __stl_debug_do(_M_iter_list._Invalidate_all());
    _Base::_M_deallocate_block();
  }
# else
  iterator begin()             { return _M_start; }
  iterator end()               { return _M_finish; }
  const_iterator begin() const { return _M_start; }
  const_iterator end()   const { return _M_finish; }  
# endif

  reverse_iterator rbegin()             
    { return reverse_iterator(_Make_iterator(_M_finish)); }
  reverse_iterator rend()               
    { return reverse_iterator(_Make_iterator(_M_start)); }
  const_reverse_iterator rbegin() const 
    { return const_reverse_iterator(_Make_const_iterator(_M_finish)); }
  const_reverse_iterator rend()   const 
    { return const_reverse_iterator(_Make_const_iterator(_M_start)); }

public:                         // Size, capacity, etc.
  size_type size() const { return _M_finish - _M_start; }
  size_type length() const { return size(); }

  size_t max_size() const { return _Base::max_size(); }


  void resize(size_type __n, _CharT __c) {
    if (__n <= size())
      erase(begin() + __n, end());
    else
      append(__n - size(), __c);
  }
  void resize(size_type __n) { resize(__n, _M_null()); }

  void reserve(size_type = 0);

  size_type capacity() const { return (_M_end_of_storage._M_data - _M_start) - 1; }

  void clear() {
    if (!empty()) {
      __stl_debug_do(_M_iter_list._Invalidate_all());
      _Traits::assign(*_M_start, _M_null());
      destroy(_M_start+1, _M_finish+1);
      _M_finish = _M_start;
    }
  } 

  bool empty() const { return _M_start == _M_finish; }    

public:                         // Element access.

  const_reference operator[](size_type __n) const
    { return *(_M_start + __n); }
  reference operator[](size_type __n)
    { return *(_M_start + __n); }

  const_reference at(size_type __n) const {
    if (__n >= size())
      _M_throw_out_of_range();
    return *(_M_start + __n);
  }

  reference at(size_type __n) {
    if (__n >= size())
      _M_throw_out_of_range();
    return *(_M_start + __n);
  }

public:                         // Append, operator+=, push_back.

  _Self& operator+=(const _Self& __s) { return append(__s); }
  _Self& operator+=(const _CharT* __s) { __STL_FIX_LITERAL_BUG(__s) return append(__s); }
  _Self& operator+=(_CharT __c) { push_back(__c); return *this; }

  _Self& append(const _Self& __s) 
    { return append(__s._M_start, __s._M_finish); }

  _Self& append(const _Self& __s,
                       size_type __pos, size_type __n)
  {
    if (__pos > __s.size())
      _M_throw_out_of_range();
    return append(__s._M_start + __pos,
                  __s._M_start + __pos + min(__n, __s.size() - __pos));
  }

  _Self& append(const _CharT* __s, size_type __n) 
    { __STL_FIX_LITERAL_BUG(__s) return append(__s, __s+__n); }
  _Self& append(const _CharT* __s) 
    { __STL_FIX_LITERAL_BUG(__s) return append(__s, __s + _Traits::length(__s)); }
  _Self& append(size_type __n, _CharT __c);

#ifdef __STL_MEMBER_TEMPLATES

  // Check to see if _InputIterator is an integer type.  If so, then
  // it can't be an iterator.
  template <class _InputIter>
  _Self& append(_InputIter __first, _InputIter __last) {
    typedef typename _Is_integer<_InputIter>::_Integral _Integral;
    return _M_append_dispatch(__first, __last, _Integral());
  }

#else /* __STL_MEMBER_TEMPLATES */

  _Self& append(const _CharT* __first, const _CharT* __last);

#endif /* __STL_MEMBER_TEMPLATES */

  void push_back(_CharT __c) {
    if (_M_finish + 1 == _M_end_of_storage._M_data)
      reserve(size() + max(size(), __STATIC_CAST(size_type,1)));
    _M_construct_null(_M_finish + 1);
    _Traits::assign(*_M_finish, __c);
    ++_M_finish;
  }

  void pop_back() {
    __stl_debug_do(__invalidate_iterator(&_M_iter_list,end()));
    _Traits::assign(*(_M_finish - 1), _M_null());
    destroy(_M_finish);
    --_M_finish;
  }

private:                        // Helper functions for append.

#ifdef __STL_MEMBER_TEMPLATES

  template <class _InputIter>
  _Self& append(_InputIter __f, _InputIter __last, input_iterator_tag)
  {
	  for ( ; __first != __last ; ++__first)
	    push_back(*__first);
	  return *this;
	}

  template <class _ForwardIter>
  _Self& append(_ForwardIter __first, _ForwardIter __last, 
                       forward_iterator_tag)
# ifndef __STL_INLINE_MEMBER_TEMPLATES
;
# else
  {
    __stl_debug_do(__check_range(__first, __last));
    if (__first != __last) {
	    const size_type __old_size = size();
	    difference_type __n = 0;
	    distance(__first, __last, __n);
	    if (__STATIC_CAST(size_type,__n) > max_size() || __old_size > max_size() - __STATIC_CAST(size_type,__n))
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
# endif /* __STL_INLINE_MEMBER_TEMPLATES */

  template <class _Integer>
  _Self& _M_append_dispatch(_Integer __n, _Integer __x, __true_type) {
    return append((size_type) __n, (_CharT) __x);
  }

  template <class _InputIter>
  _Self& _M_append_dispatch(_InputIter __f, _InputIter __l,
                                   __false_type) {
# if defined ( __STL_CLASS_PARTIAL_SPECIALIZATION )
    typedef typename iterator_traits<_InputIter>::iterator_category _Category;
    return append(__f, __l, _Category());
# else
    return append(__f, __l, __ITERATOR_CATEGORY(__f));
# endif
  }

#endif /* __STL_MEMBER_TEMPLATES */

public:                         // Assign
  
  _Self& assign(const _Self& __s) 
    { return assign(__s._M_start, __s._M_finish); }

  _Self& assign(const _Self& __s, 
                       size_type __pos, size_type __n) {
    if (__pos > __s.size())
      _M_throw_out_of_range();
    return assign(__s._M_start + __pos, 
                  __s._M_start + __pos + min(__n, __s.size() - __pos));
  }

  _Self& assign(const _CharT* __s, size_type __n)
    { return assign(__s, __s + __n); }

  _Self& assign(const _CharT* __s)
    { return assign(__s, __s + _Traits::length(__s)); }

  _Self& assign(size_type __n, _CharT __c);

#ifdef __STL_MEMBER_TEMPLATES

  // Check to see if _InputIterator is an integer type.  If so, then
  // it can't be an iterator.
  template <class _InputIter>
  _Self& assign(_InputIter __first, _InputIter __last) {
    typedef typename _Is_integer<_InputIter>::_Integral _Integral;
    return _M_assign_dispatch(__first, __last, _Integral());
  }

#endif  /* __STL_MEMBER_TEMPLATES */

  _Self& assign(const _CharT* __f, const _CharT* __l);

private:                        // Helper functions for assign.

#ifdef __STL_MEMBER_TEMPLATES

  template <class _Integer>
  _Self& _M_assign_dispatch(_Integer __n, _Integer __x, __true_type) {
    return assign((size_type) __n, (_CharT) __x);
  }

  template <class _InputIter>
  _Self& _M_assign_dispatch(_InputIter __f, _InputIter __l,
                                   __false_type)  {
          __stl_debug_do(__check_range(__f, __l));
          pointer __cur = _M_start;
	  while (__f != __l && __cur != _M_finish) {
	    _Traits::assign(*__cur, *__f);
	    ++__f;
	    ++__cur;
	  }
	  if (__f == __l)
	    erase(_Make_iterator(__cur), end());
	  else
	    append(__f, __l);
	  return *this;
	}

#endif  /* __STL_MEMBER_TEMPLATES */

public:                         // Insert

  _Self& insert(size_type __pos, const _Self& __s) {
    if (__pos > size())
      _M_throw_out_of_range();
    if (size() > max_size() - __s.size())
      _M_throw_length_error();
    insert(begin() + __pos, __s._M_start, __s._M_finish);
    return *this;
  }

  _Self& insert(size_type __pos, const _Self& __s,
                       size_type __beg, size_type __n) {
    if (__pos > size() || __beg > __s.size())
      _M_throw_out_of_range();
    size_type __len = min(__n, __s.size() - __beg);
    if (size() > max_size() - __len)
      _M_throw_length_error();
    insert(begin() + __pos,
           __s._M_start + __beg, __s._M_start + __beg + __len);
    return *this;
  }

  _Self& insert(size_type __pos, const _CharT* __s, size_type __n) {
    __STL_FIX_LITERAL_BUG(__s)
    if (__pos > size())
      _M_throw_out_of_range();
    if (size() > max_size() - __n)
      _M_throw_length_error();
    insert(begin() + __pos, __s, __s + __n);
    return *this;
  }

  _Self& insert(size_type __pos, const _CharT* __s) {
    __STL_FIX_LITERAL_BUG(__s)
    if (__pos > size())
      _M_throw_out_of_range();
    size_type __len = _Traits::length(__s);
    if (size() > max_size() - __len)
      _M_throw_length_error();
    insert(_Make_iterator(_M_start + __pos), __s, __s + __len);
    return *this;
  }
    
  _Self& insert(size_type __pos, size_type __n, _CharT __c) {
    if (__pos > size())
      _M_throw_out_of_range();
    if (size() > max_size() - __n)
      _M_throw_length_error();
    insert(begin() + __pos, __n, __c);
    return *this;
  }

  iterator insert(iterator __p, _CharT __c) {
    __STL_FIX_LITERAL_BUG(__p)
    if (__p == end()) {
      push_back(__c);
      return _Make_iterator(_M_finish - 1);
    }
    else
      return _Make_iterator(_M_insert_aux(_Make_ptr(__p), __c));
  }

  void insert(iterator __p, size_t __n, _CharT __c);

#ifdef __STL_MEMBER_TEMPLATES

  // Check to see if _InputIterator is an integer type.  If so, then
  // it can't be an iterator.
  template <class _InputIter>
  void insert(iterator __p, _InputIter __first, _InputIter __last) {
    typedef typename _Is_integer<_InputIter>::_Integral _Integral;
    _M_insert_dispatch(__p, __first, __last, _Integral());
  }

#else /* __STL_MEMBER_TEMPLATES */

  void insert(iterator __p, const _CharT* __first, const _CharT* __last);
# ifdef __STL_DEBUG
  void insert(iterator __p, const_iterator __first, const_iterator __last) {
    insert(__p, __first._M_iterator, __last._M_iterator); 
  }
# endif
#endif /* __STL_MEMBER_TEMPLATES */

private:                        // Helper functions for insert.

#ifdef __STL_MEMBER_TEMPLATES

  template <class _InputIter>
  void insert(iterator __p, _InputIter __first, _InputIter __last,
	      input_iterator_tag)
  {
	  for ( ; __first != __last; ++__first) {
	    __p = insert(__p, *__first);
	    ++__p;
	  }
	}

  template <class _ForwardIter>
  void insert(iterator __position, _ForwardIter __first, _ForwardIter __last, 
	      forward_iterator_tag)  
# ifndef __STL_INLINE_MEMBER_TEMPLATES
;
# else
{
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
  

  template <class _Integer>
  void _M_insert_dispatch(iterator __p, _Integer __n, _Integer __x,
                          __true_type) {
    insert(__p, (size_type) __n, (_CharT) __x);
  }

  template <class _InputIter>
  void _M_insert_dispatch(iterator __p, _InputIter __first, _InputIter __last,
                          __false_type) {
    __stl_debug_do(__check_range(__first, __last));
# if defined ( __STL_CLASS_PARTIAL_SPECIALIZATION )
    typedef typename iterator_traits<_InputIter>::iterator_category _Category;
    insert(__p, __first, __last, _Category());
# else
    insert(__p, __first, __last, __ITERATOR_CATEGORY(__first));
# endif
  }

  template <class _InputIterator>
  void 
  _M_copy(_InputIterator __first, _InputIterator __last, pointer __result) {
    __stl_debug_do(__check_range(__first, __last));
    for ( ; __first != __last; ++__first, ++__result)
      _Traits::assign(*__result, *__first);
  }

#endif /* __STL_MEMBER_TEMPLATES */

  pointer _M_insert_aux(pointer, _CharT);

  void 
  _M_copy(const _CharT* __first, const _CharT* __last, _CharT* __result) {
    __stl_debug_do(__check_range(__first, __last));
    _Traits::copy(__result, __first, __last - __first);
  }

public:                         // Erase.

  _Self& erase(size_type __pos = 0, size_type __n = npos) {
    if (__pos > size())
      _M_throw_out_of_range();
    erase(begin() + __pos, begin() + __pos + min(__n, size() - __pos));
    return *this;
  }  

  iterator erase(iterator __position) {
    __stl_debug_do(__check_if_owner(&_M_iter_list, __position));
                                // The move includes the terminating _CharT().
    _Traits::move(_Make_ptr(__position), _Make_ptr(__position) + 1, _M_finish - _Make_ptr(__position));
    __stl_debug_do(__invalidate_iterator(&_M_iter_list,end()));
    destroy(_M_finish);
    --_M_finish;
    return __position;
  }

  iterator erase(iterator __first, iterator __last) {
    __stl_debug_do(__check_range(__first, __last)&&__check_if_owner(&_M_iter_list,__first));
    if (__first != __last) {
                                // The move includes the terminating _CharT().
      _Traits::move(_Make_ptr(__first), _Make_ptr(__last), (_M_finish - _Make_ptr(__last)) + 1);
      pointer __new_finish = _M_finish - (__last - __first);
      destroy(__new_finish + 1, _M_finish + 1);
      __stl_debug_do(__invalidate_range(&_M_iter_list, _Make_iterator(__new_finish+1), end()));
      _M_finish = __new_finish;
    }
    return __first;
  }

public:                         // Replace.  (Conceptually equivalent
                                // to erase followed by insert.)
  _Self& replace(size_type __pos, size_type __n, 
                        const _Self& __s) {
    if (__pos > size())
      _M_throw_out_of_range();
    const size_type __len = min(__n, size() - __pos);
    if (size() - __len >= max_size() - __s.size())
      _M_throw_length_error();
    return replace(begin() + __pos, begin() + __pos + __len, 
                   __s._M_start, __s._M_finish);
  }

  _Self& replace(size_type __pos1, size_type __n1,
                        const _Self& __s,
                        size_type __pos2, size_type __n2) {
    if (__pos1 > size() || __pos2 > __s.size())
      _M_throw_out_of_range();
    const size_type __len1 = min(__n1, size() - __pos1);
    const size_type __len2 = min(__n2, __s.size() - __pos2);
    if (size() - __len1 >= max_size() - __len2)
      _M_throw_length_error();
    return replace(begin() + __pos1, begin() + __pos1 + __len1,
                   __s._M_start + __pos2, __s._M_start + __pos2 + __len2);
  }

  _Self& replace(size_type __pos, size_type __n1,
                        const _CharT* __s, size_type __n2) {
    __STL_FIX_LITERAL_BUG(__s)
    if (__pos > size())
      _M_throw_out_of_range();
    const size_type __len = min(__n1, size() - __pos);
    if (__n2 > max_size() || size() - __len >= max_size() - __n2)
      _M_throw_length_error();
    return replace(begin() + __pos, begin() + __pos + __len,
                   __s, __s + __n2);
  }

  _Self& replace(size_type __pos, size_type __n1,
                        const _CharT* __s) {
    __STL_FIX_LITERAL_BUG(__s)
    if (__pos > size())
      _M_throw_out_of_range();
    const size_type __len = min(__n1, size() - __pos);
    const size_type __n2 = _Traits::length(__s);
    if (__n2 > max_size() || size() - __len >= max_size() - __n2)
      _M_throw_length_error();
    return replace(begin() + __pos, begin() + __pos + __len,
                   __s, __s + _Traits::length(__s));
  }

  _Self& replace(size_type __pos, size_type __n1,
                        size_type __n2, _CharT __c) {
    if (__pos > size())
      _M_throw_out_of_range();
    const size_type __len = min(__n1, size() - __pos);
    if (__n2 > max_size() || size() - __len >= max_size() - __n2)
      _M_throw_length_error();
    return replace(begin() + __pos, begin() + __pos + __len, __n2, __c);
  }

  _Self& replace(iterator __first, iterator __last, 
                        const _Self& __s) 
    { return replace(__first, __last, __s._M_start, __s._M_finish); }

  _Self& replace(iterator __first, iterator __last,
                        const _CharT* __s, size_type __n) 
    { __STL_FIX_LITERAL_BUG(__s) return replace(__first, __last, __s, __s + __n); }

  _Self& replace(iterator __first, iterator __last,
                        const _CharT* __s) {
    __STL_FIX_LITERAL_BUG(__s)
    return replace(__first, __last, __s, __s + _Traits::length(__s));
  }

  _Self& replace(iterator __first, iterator __last, 
                        size_type __n, _CharT __c);

  // Check to see if _InputIterator is an integer type.  If so, then
  // it can't be an iterator.
#ifdef __STL_MEMBER_TEMPLATES
  template <class _InputIter>
  _Self& replace(iterator __first, iterator __last,
                        _InputIter __f, _InputIter __l) {
    typedef typename _Is_integer<_InputIter>::_Integral _Integral;
    return _M_replace_dispatch(__first, __last, __f, __l,  _Integral());
  }
#else /* __STL_MEMBER_TEMPLATES */
  _Self& replace(iterator __first, iterator __last,
		 const _CharT* __f, const _CharT* __l);
#endif /* __STL_MEMBER_TEMPLATES */

private:                        // Helper functions for replace.

#ifdef __STL_MEMBER_TEMPLATES

  template <class _Integer>
  _Self& _M_replace_dispatch(iterator __first, iterator __last,
                                    _Integer __n, _Integer __x,
                                    __true_type) {
    return replace(__first, __last, (size_type) __n, (_CharT) __x);
  }

  template <class _InputIter>
  _Self& _M_replace_dispatch(iterator __first, iterator __last,
                                    _InputIter __f, _InputIter __l,
                                    __false_type) {
    __stl_debug_do(__check_if_owner(&_M_iter_list, __first) && __check_range(__first, __last) 
		   && __check_range(__f, __l));
# ifdef __STL_CLASS_PARTIAL_SPECIALIZATION
    typedef typename iterator_traits<_InputIter>::iterator_category _Category;
    return replace(__first, __last, __f, __l, _Category());
# else
    return replace(__first, __last, __f, __l, __ITERATOR_CATEGORY(__f));
# endif
  }

  template <class _InputIter>
  _Self& replace(iterator __first, iterator __last,
                        _InputIter __f, _InputIter __l, input_iterator_tag)  {
	  for ( ; __first != __last && __f != __l; ++__first, ++__f)
	    _Traits::assign(*__first, *__f);

	  if (__f == __l)
	    erase(__first, __last);
	  else
	    insert(__last, __f, __l);
	  return *this;
	}

  template <class _ForwardIter>
  _Self& replace(iterator __first, iterator __last,
                        _ForwardIter __f, _ForwardIter __l, 
                        forward_iterator_tag)  {
	  difference_type __n = 0;
	  distance(__f, __l, __n);
	  const difference_type __len = __last - __first;
	  if (__len >= __n) {
	    _M_copy(__f, __l, _Make_ptr(__first));
	    erase(__first + __n, __last);
	  }
	  else {
	    _ForwardIter __m = __f;
	    advance(__m, __len);
	    _M_copy(__f, __m, _Make_ptr(__first));
	    insert(__last, __m, __l);
	  }
	  return *this;
	}

#endif /* __STL_MEMBER_TEMPLATES */

public:                         // Other modifier member functions.

  size_type copy(_CharT* __s, size_type __n, size_type __pos = 0) const {
    __STL_FIX_LITERAL_BUG(__s)
    if (__pos > size())
      _M_throw_out_of_range();
    const size_type __len = min(__n, size() - __pos);
    _Traits::copy(__s, _M_start + __pos, __len);
    return __len;
  }

  void swap(_Self& __s) {
    __stl_debug_do(_M_iter_list._Swap_owners(__s._M_iter_list));
    __STLPORT_STD::swap(_M_start, __s._M_start);
    __STLPORT_STD::swap(_M_finish, __s._M_finish);
    __STLPORT_STD::swap(_M_end_of_storage, __s._M_end_of_storage);
  }

public:                         // Conversion to C string.

  const _CharT* c_str() const { return _M_start; }
  const _CharT* data()  const { return _M_start; }

public:                         // find.

  size_type find(const _Self& __s, size_type __pos = 0) const 
    { return find(__s._M_start, __pos, __s.size()); }

  size_type find(const _CharT* __s, size_type __pos = 0) const 
    { __STL_FIX_LITERAL_BUG(__s) return find(__s, __pos, _Traits::length(__s)); }

  size_type find(const _CharT* __s, size_type __pos, size_type __n) const;
  size_type find(_CharT __c, size_type __pos = 0) const;

public:                         // rfind.

  size_type rfind(const _Self& __s, size_type __pos = npos) const 
    { return rfind(__s._M_start, __pos, __s.size()); }

  size_type rfind(const _CharT* __s, size_type __pos = npos) const 
    { return rfind(__s, __pos, _Traits::length(__s)); }

  size_type rfind(const _CharT* __s, size_type __pos, size_type __n) const;
  size_type rfind(_CharT __c, size_type __pos = npos) const;

public:                         // find_first_of
  
  size_type find_first_of(const _Self& __s, size_type __pos = 0) const 
    { return find_first_of(__s._M_start, __pos, __s.size()); }

  size_type find_first_of(const _CharT* __s, size_type __pos = 0) const 
    { __STL_FIX_LITERAL_BUG(__s) return find_first_of(__s, __pos, _Traits::length(__s)); }

  size_type find_first_of(const _CharT* __s, size_type __pos, 
                          size_type __n) const;

  size_type find_first_of(_CharT __c, size_type __pos = 0) const 
    { return find(__c, __pos); }

public:                         // find_last_of

  size_type find_last_of(const _Self& __s,
                         size_type __pos = npos) const
    { return find_last_of(__s._M_start, __pos, __s.size()); }

  size_type find_last_of(const _CharT* __s, size_type __pos = npos) const 
    { __STL_FIX_LITERAL_BUG(__s) return find_last_of(__s, __pos, _Traits::length(__s)); }

  size_type find_last_of(const _CharT* __s, size_type __pos, 
                         size_type __n) const;

  size_type find_last_of(_CharT __c, size_type __pos = npos) const {
    return rfind(__c, __pos);
  }

public:                         // find_first_not_of

  size_type find_first_not_of(const _Self& __s, 
                              size_type __pos = 0) const 
    { return find_first_not_of(__s._M_start, __pos, __s.size()); }

  size_type find_first_not_of(const _CharT* __s, size_type __pos = 0) const 
    { __STL_FIX_LITERAL_BUG(__s) return find_first_not_of(__s, __pos, _Traits::length(__s)); }

  size_type find_first_not_of(const _CharT* __s, size_type __pos,
                              size_type __n) const;

  size_type find_first_not_of(_CharT __c, size_type __pos = 0) const;

public:                         // find_last_not_of

  size_type find_last_not_of(const _Self& __s, 
                             size_type __pos = npos) const
    { return find_last_not_of(__s._M_start, __pos, __s.size()); }

  size_type find_last_not_of(const _CharT* __s, size_type __pos = npos) const
    { __STL_FIX_LITERAL_BUG(__s) return find_last_not_of(__s, __pos, _Traits::length(__s)); }

  size_type find_last_not_of(const _CharT* __s, size_type __pos,
                             size_type __n) const;

  size_type find_last_not_of(_CharT __c, size_type __pos = npos) const;

public:                         // Substring.

  _Self substr(size_type __pos = 0, size_type __n = npos) const {
    if (__pos > size())
      _M_throw_out_of_range();
    return _Self(_M_start + __pos, 
                        _M_start + __pos + min(__n, size() - __pos));
  }

public:                         // Compare

  int compare(const _Self& __s) const 
    { return _M_compare(_M_start, _M_finish, __s._M_start, __s._M_finish); }

  int compare(size_type __pos1, size_type __n1,
              const _Self& __s) const {
    if (__pos1 > size())
      _M_throw_out_of_range();
    return _M_compare(_M_start + __pos1, 
                      _M_start + __pos1 + min(__n1, size() - __pos1),
                      __s._M_start, __s._M_finish);
  }
    
  int compare(size_type __pos1, size_type __n1,
              const _Self& __s,
              size_type __pos2, size_type __n2) const {
    if (__pos1 > size() || __pos2 > __s.size())
      _M_throw_out_of_range();
    return _M_compare(_M_start + __pos1, 
                      _M_start + __pos1 + min(__n1, size() - __pos1),
                      __s._M_start + __pos2, 
                      __s._M_start + __pos2 + min(__n2, size() - __pos2));
  }

  int compare(const _CharT* __s) const {
    __STL_FIX_LITERAL_BUG(__s) 
      return _M_compare(_M_start, _M_finish, __s, __s + _Traits::length(__s));
  }

  int compare(size_type __pos1, size_type __n1, const _CharT* __s) const {
    __STL_FIX_LITERAL_BUG(__s)
    if (__pos1 > size())
      _M_throw_out_of_range();
    return _M_compare(_M_start + __pos1, 
                      _M_start + __pos1 + min(__n1, size() - __pos1),
                      __s, __s + _Traits::length(__s));
  }

  int compare(size_type __pos1, size_type __n1, const _CharT* __s,
              size_type __n2) const {
    __STL_FIX_LITERAL_BUG(__s)
    if (__pos1 > size())
      _M_throw_out_of_range();
    return _M_compare(_M_start + __pos1, 
                      _M_start + __pos1 + min(__n1, size() - __pos1),
                      __s, __s + __n2);
  }

public:                        // Helper functions for compare.
  
  static int _M_compare(const _CharT* __f1, const _CharT* __l1,
                        const _CharT* __f2, const _CharT* __l2) {
    const ptrdiff_t __n1 = __l1 - __f1;
    const ptrdiff_t __n2 = __l2 - __f2;
    const int cmp = _Traits::compare(__f1, __f2, min(__n1, __n2));
    return cmp != 0 ? cmp : (__n1 < __n2 ? -1 : (__n1 > __n2 ? 1 : 0));
  }
};


// This is a hook to instantiate STLport exports in a designated DLL
# if defined (__STL_USE_DECLSPEC)
__STL_EXPORT template class __STL_CLASS_DECLSPEC _STL_alloc_proxy<char *,char,allocator<char> >;
__STL_EXPORT template class __STL_CLASS_DECLSPEC _String_base<char, allocator<char> >;
__STL_EXPORT template class __STL_CLASS_DECLSPEC basic_string<char, char_traits<char>, allocator<char> >;
#  if defined (__STL_HAS_WCHAR_T)
__STL_EXPORT template class __STL_CLASS_DECLSPEC _STL_alloc_proxy<wchar_t *,wchar_t,allocator<wchar_t> >;
__STL_EXPORT template class __STL_CLASS_DECLSPEC _String_base<wchar_t, allocator<wchar_t> >;
__STL_EXPORT template class __STL_CLASS_DECLSPEC basic_string<wchar_t, char_traits<wchar_t>, allocator<wchar_t> >;
#  endif
# endif /* __STL_USE_DECLSPEC */

// ------------------------------------------------------------
// Non-member functions.

template <class _CharT, class _Traits, class _Alloc>
inline basic_string<_CharT,_Traits,_Alloc>
operator+(const basic_string<_CharT,_Traits,_Alloc>& __s,
          const basic_string<_CharT,_Traits,_Alloc>& __y)
{
  typedef basic_string<_CharT,_Traits,_Alloc> _Str;
  typedef typename _Str::_Reserve_t _Reserve_t;
# ifdef __GNUC__
  // gcc counts this as a function
  _Str __result  = _Str(_Reserve_t(),__s.size() + __y.size());
# else
  _Str __result(_Reserve_t(), __s.size() + __y.size());
# endif
  __result.append(__s);
  __result.append(__y);
  return __result;
}

template <class _CharT, class _Traits, class _Alloc>
inline basic_string<_CharT,_Traits,_Alloc>
operator+(const _CharT* __s,
          const basic_string<_CharT,_Traits,_Alloc>& __y) {
  __STL_FIX_LITERAL_BUG(__s)
  typedef basic_string<_CharT,_Traits,_Alloc> _Str;
  typedef typename _Str::_Reserve_t _Reserve_t;
  const size_t __n = _Traits::length(__s);
# ifdef __GNUC__
  _Str __result = _Str(_Reserve_t(), __n + __y.size());
# else
  _Str __result(_Reserve_t(), __n + __y.size());
# endif
  __result.append(__s, __s + __n);
  __result.append(__y);
  return __result;
}

template <class _CharT, class _Traits, class _Alloc>
inline basic_string<_CharT,_Traits,_Alloc>
operator+(_CharT __c,
          const basic_string<_CharT,_Traits,_Alloc>& __y) {
  typedef basic_string<_CharT,_Traits,_Alloc> _Str;
  typedef typename _Str::_Reserve_t _Reserve_t;
# ifdef __GNUC__
  _Str __result = _Str(_Reserve_t(), 1 + __y.size());
# else
  _Str __result(_Reserve_t(), 1 + __y.size());
# endif
  __result.push_back(__c);
  __result.append(__y);
  return __result;
}

template <class _CharT, class _Traits, class _Alloc>
inline basic_string<_CharT,_Traits,_Alloc>
operator+(const basic_string<_CharT,_Traits,_Alloc>& __x,
          const _CharT* __s) {
  __STL_FIX_LITERAL_BUG(__s)
  typedef basic_string<_CharT,_Traits,_Alloc> _Str;
  typedef typename _Str::_Reserve_t _Reserve_t;
  const size_t __n = _Traits::length(__s);
# ifdef __GNUC__
  _Str __result = _Str(_Reserve_t(), __x.size() + __n, __x.get_allocator());
# else
  _Str __result(_Reserve_t(), __x.size() + __n, __x.get_allocator());
# endif
  __result.append(__x);
  __result.append(__s, __s + __n);
  return __result;
}

template <class _CharT, class _Traits, class _Alloc>
inline basic_string<_CharT,_Traits,_Alloc>
operator+(const basic_string<_CharT,_Traits,_Alloc>& __x,
          const _CharT __c) {
  typedef basic_string<_CharT,_Traits,_Alloc> _Str;
  typedef typename _Str::_Reserve_t _Reserve_t;
# ifdef __GNUC__
  _Str __result = _Str(_Reserve_t(), __x.size() + 1, __x.get_allocator());
# else
  _Str __result(_Reserve_t(), __x.size() + 1, __x.get_allocator());
# endif
  __result.append(__x);
  __result.push_back(__c);
  return __result;
}

// Operator== and operator!=

template <class _CharT, class _Traits, class _Alloc>
inline bool
operator==(const basic_string<_CharT,_Traits,_Alloc>& __x,
           const basic_string<_CharT,_Traits,_Alloc>& __y) {
  return __x.size() == __y.size() &&
         _Traits::compare(__x.data(), __y.data(), __x.size()) == 0;
}

template <class _CharT, class _Traits, class _Alloc>
inline bool
operator==(const _CharT* __s,
           const basic_string<_CharT,_Traits,_Alloc>& __y) {
  __STL_FIX_LITERAL_BUG(__s)
  size_t __n = _Traits::length(__s);
  return __n == __y.size() && _Traits::compare(__s, __y.data(), __n) == 0;
}

template <class _CharT, class _Traits, class _Alloc>
inline bool
operator==(const basic_string<_CharT,_Traits,_Alloc>& __x,
           const _CharT* __s) {
  __STL_FIX_LITERAL_BUG(__s)
  size_t __n = _Traits::length(__s);
  return __x.size() == __n && _Traits::compare(__x.data(), __s, __n) == 0;
}

// Operator< (and also >, <=, and >=).

template <class _CharT, class _Traits, class _Alloc>
inline bool
operator<(const basic_string<_CharT,_Traits,_Alloc>& __x,
          const basic_string<_CharT,_Traits,_Alloc>& __y) {
  return basic_string<_CharT,_Traits,_Alloc>
    ::_M_compare(_Make_ptr(__x.begin()), _Make_ptr(__x.end()), 
		 _Make_ptr(__y.begin()), _Make_ptr(__y.end())) < 0;
}

template <class _CharT, class _Traits, class _Alloc>
inline bool
operator<(const _CharT* __s,
          const basic_string<_CharT,_Traits,_Alloc>& __y) {
  __STL_FIX_LITERAL_BUG(__s)
  size_t __n = _Traits::length(__s);
  return basic_string<_CharT,_Traits,_Alloc>
    ::_M_compare(__s, __s + __n, _Make_ptr(__y.begin()), _Make_ptr(__y.end())) < 0;
}

template <class _CharT, class _Traits, class _Alloc>
inline bool
operator<(const basic_string<_CharT,_Traits,_Alloc>& __x,
          const _CharT* __s) {
  __STL_FIX_LITERAL_BUG(__s)
  size_t __n = _Traits::length(__s);
  return basic_string<_CharT,_Traits,_Alloc>
    ::_M_compare(_Make_ptr(__x.begin()), _Make_ptr(__x.end()), __s, __s + __n) < 0;
}

#ifdef __STL_USE_SEPARATE_RELOPS_NAMESPACE

template <class _CharT, class _Traits, class _Alloc>
inline bool
operator!=(const basic_string<_CharT,_Traits,_Alloc>& __x,
           const basic_string<_CharT,_Traits,_Alloc>& __y) {
  return !(__x == __y);
}

template <class _CharT, class _Traits, class _Alloc>
inline bool
operator>(const basic_string<_CharT,_Traits,_Alloc>& __x,
          const basic_string<_CharT,_Traits,_Alloc>& __y) {
  return __y < __x;
}

template <class _CharT, class _Traits, class _Alloc>
inline bool
operator<=(const basic_string<_CharT,_Traits,_Alloc>& __x,
           const basic_string<_CharT,_Traits,_Alloc>& __y) {
  return !(__y < __x);
}

template <class _CharT, class _Traits, class _Alloc>
inline bool
operator>=(const basic_string<_CharT,_Traits,_Alloc>& __x,
           const basic_string<_CharT,_Traits,_Alloc>& __y) {
  return !(__x < __y);
}

#endif /* __STL_USE_SEPARATE_RELOPS_NAMESPACE */

template <class _CharT, class _Traits, class _Alloc>
inline bool
operator!=(const _CharT* __s,
           const basic_string<_CharT,_Traits,_Alloc>& __y) {
  __STL_FIX_LITERAL_BUG(__s)
  return !(__s == __y);
}

template <class _CharT, class _Traits, class _Alloc>
inline bool
operator!=(const basic_string<_CharT,_Traits,_Alloc>& __x,
           const _CharT* __s) {
  __STL_FIX_LITERAL_BUG(__s)
  return !(__x == __s);
}

template <class _CharT, class _Traits, class _Alloc>
inline bool
operator>(const _CharT* __s,
          const basic_string<_CharT,_Traits,_Alloc>& __y) {
  __STL_FIX_LITERAL_BUG(__s)
  return __y < __s;
}

template <class _CharT, class _Traits, class _Alloc>
inline bool
operator>(const basic_string<_CharT,_Traits,_Alloc>& __x,
          const _CharT* __s) {
  __STL_FIX_LITERAL_BUG(__s)
  return __s < __x;
}

template <class _CharT, class _Traits, class _Alloc>
inline bool
operator<=(const _CharT* __s,
           const basic_string<_CharT,_Traits,_Alloc>& __y) {
  __STL_FIX_LITERAL_BUG(__s)
  return !(__y < __s);
}

template <class _CharT, class _Traits, class _Alloc>
inline bool
operator<=(const basic_string<_CharT,_Traits,_Alloc>& __x,
           const _CharT* __s) {
  __STL_FIX_LITERAL_BUG(__s)
  return !(__s < __x);
}

template <class _CharT, class _Traits, class _Alloc>
inline bool
operator>=(const _CharT* __s,
           const basic_string<_CharT,_Traits,_Alloc>& __y) {
  __STL_FIX_LITERAL_BUG(__s)
  return !(__s < __y);
}

template <class _CharT, class _Traits, class _Alloc>
inline bool
operator>=(const basic_string<_CharT,_Traits,_Alloc>& __x,
           const _CharT* __s) {
  __STL_FIX_LITERAL_BUG(__s)
  return !(__x < __s);
}


// Swap.

#ifdef __STL_FUNCTION_TMPL_PARTIAL_ORDER

template <class _CharT, class _Traits, class _Alloc>
inline void swap(basic_string<_CharT,_Traits,_Alloc>& __x,
                 basic_string<_CharT,_Traits,_Alloc>& __y) {
  __x.swap(__y);
}

#endif /* __STL_FUNCTION_TMPL_PARTIAL_ORDER */

// I/O.  

#if defined (__STL_USE_NEW_IOSTREAMS)

template <class _CharT, class _Traits, class _Alloc>
basic_ostream<_CharT, _Traits>&
operator<<(basic_ostream<_CharT, _Traits>& __os, 
           const basic_string<_CharT,_Traits,_Alloc>& __s);

template <class _CharT, class _Traits, class _Alloc>
basic_istream<_CharT, _Traits>& 
operator>>(basic_istream<_CharT, _Traits>& __is,
           basic_string<_CharT,_Traits,_Alloc>& __s);

template <class _CharT, class _Traits, class _Alloc>    
basic_istream<_CharT, _Traits>& 
getline(basic_istream<_CharT, _Traits>& __is,
        basic_string<_CharT,_Traits,_Alloc>& __s,
        _CharT __delim);

# if !(defined (__BORLANDC__) && ! defined (__STL_USE_OWN_NAMESPACE))

template <class _CharT, class _Traits, class _Alloc>    
inline basic_istream<_CharT, _Traits>& 
getline(basic_istream<_CharT, _Traits>& __is,
        basic_string<_CharT,_Traits,_Alloc>& __s)
{
  return getline(__is, __s, '\n');
}
# endif

#elif ! defined ( __STL_USE_NO_IOSTREAMS )

template <class _CharT, class _Traits, class _Alloc>
ostream& operator<<(ostream& __os, 
                    const basic_string<_CharT,_Traits,_Alloc>& __s);

template <class _CharT, class _Traits, class _Alloc>
istream& operator>>(istream& __is, basic_string<_CharT,_Traits,_Alloc>& __s);

template <class _CharT, class _Traits, class _Alloc>    
istream& getline(istream& __is,
                 basic_string<_CharT,_Traits,_Alloc>& __s,
                 _CharT __delim);


template <class _CharT, class _Traits, class _Alloc>    
inline istream& 
getline(istream& __is, basic_string<_CharT,_Traits,_Alloc>& __s)
{
  return getline(__is, __s, '\n');
}

#endif /* __STL_USE_NEW_IOSTREAMS */

template <class _CharT, class _Traits, class _Alloc>
void _S_string_copy(const basic_string<_CharT,_Traits,_Alloc>& __s,
                    _CharT* __buf,
                    size_t __n);

__STL_END_NAMESPACE

// cleanup
# if !(defined (__IBMCPP__) || defined (__xlC__))
#  undef _Make_ptr
# endif
#  undef _Make_iterator
#  undef _Make_const_iterator


#if defined(__sgi) && !defined(__GNUC__) && (_MIPS_SIM != _MIPS_SIM_ABI32)
#pragma reset woff 1174
#pragma reset woff 1375
#endif

# if !defined (__STL_LINK_TIME_INSTANTIATION)
#  include <stl_string.c>
# endif

#endif /* __SGI_STL_STRING */


// Local Variables:
// mode:C++
// End:

