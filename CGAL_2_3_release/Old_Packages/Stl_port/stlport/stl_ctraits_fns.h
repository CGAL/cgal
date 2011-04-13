/*
 * Copyright (c) 1999
 * Silicon Graphics Computer Systems, Inc.
 *
 * Permission to use, copy, modify, distribute and sell this software
 * and its documentation for any purpose is hereby granted without fee,
 * provided that the above copyright notice appear in all copies and
 * that both that copyright notice and this permission notice appear
 * in supporting documentation.  Silicon Graphics makes no
 * representations about the suitability of this software for any
 * purpose.  It is provided "as is" without express or implied warranty.
 */ 

// WARNING: This is an internal header file, included by other C++
// standard library headers.  You should not attempt to use this header
// file directly.

#ifndef __SGI_STL_INTERNAL_CTRAITS_FUNCTIONS_H
#define __SGI_STL_INTERNAL_CTRAITS_FUNCTIONS_H

// This file contains a few small adapters that allow a character
// traits class to be used as a function object.

__STL_BEGIN_NAMESPACE

template <class _Traits>
struct _Eq_traits
  : public binary_function<typename _Traits::char_type,
                           typename _Traits::char_type,
                           bool>
{
  bool operator()(const typename _Traits::char_type& __x,
                  const typename _Traits::char_type& __y) const
    { return _Traits::eq(__x, __y); }
};

template <class _Traits>
struct _Eq_int_traits
  : public binary_function<typename _Traits::char_type,
                           typename _Traits::int_type,
                           bool>
{
  bool operator()(const typename _Traits::char_type& __x,
                  const typename _Traits::int_type& __y) const
    { return _Traits::eq_int_type(_Traits::to_int_type(__x), __y); }
};

template <class _Traits>
struct _Lt_traits
  : public binary_function<typename _Traits::char_type,
                           typename _Traits::char_type,
                           bool>
{
  bool operator()(const typename _Traits::char_type& __x,
                  const typename _Traits::char_type& __y) const
    { return _Traits::lt(__x, __y); }
};

__STL_END_NAMESPACE

#endif /* __SGI_STL_INTERNAL_CTRAITS_FUNCTIONS_H */

// Local Variables:
// mode:C++
// End:




