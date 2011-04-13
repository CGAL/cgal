// ============================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// ----------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Borland_fixes.h
// chapter       : $CGAL_Chapter: Configuration $
//
// author(s)     : Dimitri Pasechnik <dima@cs.uu.nl>
//
// coordinator   : Utrecht University
// ============================================================================

// Portions of code in this file are
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
 * Permission to use or copy this software for any purpose is hereby
 * granted 
 * without fee, provided the above notices are retained on all copies.
 * Permission to modify the code and to distribute modified code is
 * granted,
 * provided the above notices are retained, and a notice that the code was
 * modified is included with the above copyright notice.
 *
 */

#if defined(__BORLANDC__) && __BORLANDC__ > 0x520
#define CGAL_MASK_FINITE_VALID 1

#include <cstddef>
using std::size_t;
#include <ctime> // for time_t
using std::time_t;
#include <iterator>
#include <functional>
#include <utility>

#include <float>

class Borland_floating_point_initialiser
{
public:
  Borland_floating_point_initialiser();
};

inline Borland_floating_point_initialiser::Borland_floating_point_initialiser()
{
   _control87(MCW_EM, MCW_EM);
}

extern void foogj()
{}

namespace {

Borland_floating_point_initialiser bfpi;

}

namespace std {
// Borland-specific, borrowed from STLPort
template <class _Tp>
struct iterator_traits<_Tp* const> {
  typedef random_access_iterator_tag iterator_category;
  typedef _Tp                         value_type;
  typedef ptrdiff_t                   difference_type;
  typedef const _Tp*                  pointer;
  typedef const _Tp&                  reference;
};
// borrowed from STLPort (HP STL compatibility, stl_iterator_base.h)
template <class _Iter>
inline typename iterator_traits<_Iter>::iterator_category
CGAL__iterator_category(const _Iter&)
{
  typedef typename iterator_traits<_Iter>::iterator_category _Category;
  return _Category();
}

template <class _Iter>
inline typename iterator_traits<_Iter>::difference_type*
CGAL__distance_type(const _Iter&)
{
  typedef typename iterator_traits<_Iter>::difference_type _diff_type;
  return static_cast<_diff_type*>(0);
}

template <class _Iter>
inline typename iterator_traits<_Iter>::value_type*
CGAL__value_type(const _Iter&)
{
  typedef typename iterator_traits<_Iter>::value_type _value_type;
  return static_cast<_value_type*>(0);
}

template <class _Iter>
inline typename iterator_traits<_Iter>::iterator_category
iterator_category(const _Iter& __i) { 
  return CGAL__iterator_category(__i); }


template <class _Iter>
inline typename iterator_traits<_Iter>::difference_type*
distance_type(const _Iter& __i) { return CGAL__distance_type(__i); }

template <class _Iter>
inline typename iterator_traits<_Iter>::value_type*
value_type(const _Iter& __i) { return CGAL__value_type(__i); }

// Strange Borland-specific fix (DVP).
// this fixes a matching problem in algorith.h
template <class _Iter>
inline typename iterator_traits<_Iter>::value_type*
__value_type(_Iter& const)
{
  typedef typename iterator_traits<_Iter>::value_type _value_type;
  return static_cast<_value_type*>(0);
}

// quick (and dirty ?) fix for reverse_bidirectional_iterator
template <class It, class Cat, class T, class Ref = T&, 
  class P = T*, class Dist = ptrdiff_t>
class reverse_bidirectional_iterator : 
  public std::reverse_iterator<It> {};

template <class P_Value,class P_Dist = ptrdiff_t>
class input_iterator
  : public iterator<input_iterator_tag(),P_Value,P_Dist> {};

class output_iterator
  : public iterator<output_iterator_tag(),void *> {};

// borrowed from STLPort
// unary_compose and binary_compose (extensions, not part of the standard).

template <class _Operation1, class _Operation2>
class unary_compose : 
  public unary_function<typename  _Operation2 :: argument_type  ,
                        typename  _Operation1 :: result_type  > {
protected:
  _Operation1 _M_fn1;
  _Operation2 _M_fn2;
public:
  unary_compose(const _Operation1& __x, const _Operation2& __y) 
    : _M_fn1(__x), _M_fn2(__y) {}
  typename _Operation1::result_type
  operator()(const typename _Operation2::argument_type& __x) const {
    return _M_fn1(_M_fn2(__x));
  }
};

template <class _Operation1, class _Operation2>
inline unary_compose<_Operation1,_Operation2> 
compose1(const _Operation1& __fn1, const _Operation2& __fn2)
{
  return unary_compose<_Operation1,_Operation2>(__fn1, __fn2);
}

template <class _Operation1, class _Operation2, class _Operation3>
class binary_compose : 
    public unary_function<typename  _Operation2 :: argument_type  ,
                          typename  _Operation1 :: result_type  > {
protected:
  _Operation1 _M_fn1;
  _Operation2 _M_fn2;
  _Operation3 _M_fn3;
public:
  binary_compose(const _Operation1& __x, const _Operation2& __y, 
                 const _Operation3& __z) 
    : _M_fn1(__x), _M_fn2(__y), _M_fn3(__z) { }
  typename _Operation1::result_type
  operator()(const typename _Operation2::argument_type& __x) const {
    return _M_fn1(_M_fn2(__x), _M_fn3(__x));
  }
};

template <class _Operation1, class _Operation2, class _Operation3>
inline binary_compose<_Operation1, _Operation2, _Operation3> 
compose2(const _Operation1& __fn1, const _Operation2& __fn2, 
         const _Operation3& __fn3)
{
  return binary_compose<_Operation1,_Operation2,_Operation3>
    (__fn1, __fn2, __fn3);
}



// Borrowed from STLPort (and modified)
// copy_n (not part of the C++ standard)

template <class _InputIter, class _Size, class _OutputIter>
inline
std::pair<_InputIter, _OutputIter> __copy_n(_InputIter __first, _Size __count,
                                          _OutputIter __result,
                                          std::input_iterator_tag) {
  for ( ; __count > 0; --__count) {
    *__result = *__first;
    ++__first;
    ++__result;
  }
  return std::pair<_InputIter, _OutputIter>(__first, __result);
}

template <class _RAIter, class _Size, class _OutputIter>
inline std::pair<_RAIter, _OutputIter>
__copy_n(_RAIter __first, _Size __count,
         _OutputIter __result,
         std::random_access_iterator_tag) {
  _RAIter __last = __first + __count;
  return std::pair<_RAIter, _OutputIter>(__last, 
                                       std::copy(__first, __last, __result));
}

template <class _InputIter, class _Size, class _OutputIter>
inline std::pair<_InputIter, _OutputIter>
copy_n(_InputIter __first, _Size __count, _OutputIter __result) {
  return __copy_n(__first, __count, __result,
                std::iterator_category(__first));

  }
  
} // namespace std
#endif
