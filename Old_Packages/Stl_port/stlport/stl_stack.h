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

#ifndef __SGI_STL_INTERNAL_STACK_H
#define __SGI_STL_INTERNAL_STACK_H

__STL_BEGIN_NAMESPACE

# if !defined ( __STL_LIMITED_DEFAULT_TEMPLATES )
template <class _Tp, class _Sequence = deque<_Tp> >
# elif defined ( __STL_MINIMUM_DEFAULT_TEMPLATE_PARAMS )
# define __STL_STACK_ARGS _Tp
template <class _Tp>
# else
template <class _Tp, class _Sequence>
# endif
class stack {

# ifdef __STL_STACK_ARGS 
  typedef deque<_Tp> _Sequence;
# endif

public:
  typedef typename _Sequence::value_type      value_type;
  typedef typename _Sequence::size_type       size_type;
  typedef          _Sequence                  container_type;

  typedef typename _Sequence::reference       reference;
  typedef typename _Sequence::const_reference const_reference;
protected:
  _Sequence c;
public:
  stack() : c() {}
  explicit stack(const _Sequence& __s) : c(__s) {}

  bool empty() const { return c.empty(); }
  size_type size() const { return c.size(); }
  reference top() { return c.back(); }
  const_reference top() const { return c.back(); }
  void push(const value_type& __x) { c.push_back(__x); }
  void pop() { c.pop_back(); }
  const _Sequence& _Get_c() const { return c; }
};

# ifndef __STL_STACK_ARGS
#  define __STL_STACK_ARGS _Tp, _Sequence
#  define __STL_STACK_HEADER_ARGS class _Tp, class _Sequence
# else
#  define __STL_STACK_HEADER_ARGS class _Tp
# endif

template < __STL_STACK_HEADER_ARGS >
inline bool operator==(const stack< __STL_STACK_ARGS >& __x, const stack< __STL_STACK_ARGS >& __y)
{
  return __x._Get_c() == __y._Get_c();
}

template < __STL_STACK_HEADER_ARGS >
inline bool operator<(const stack< __STL_STACK_ARGS >& __x, const stack< __STL_STACK_ARGS >& __y)
{
  return __x._Get_c() < __y._Get_c();
}

#ifdef __STL_USE_SEPARATE_RELOPS_NAMESPACE

template < __STL_STACK_HEADER_ARGS >
inline bool operator!=(const stack< __STL_STACK_ARGS >& __x, const stack< __STL_STACK_ARGS >& __y)
{
  return !(__x == __y);
}

template < __STL_STACK_HEADER_ARGS >
inline bool operator>(const stack< __STL_STACK_ARGS >& __x, const stack< __STL_STACK_ARGS >& __y)
{
  return __y < __x;
}

template < __STL_STACK_HEADER_ARGS >
inline bool operator<=(const stack< __STL_STACK_ARGS >& __x, const stack< __STL_STACK_ARGS >& __y)
{
  return !(__y < __x);
}

template < __STL_STACK_HEADER_ARGS >
inline bool operator>=(const stack< __STL_STACK_ARGS >& __x, const stack< __STL_STACK_ARGS >& __y)
{
  return !(__x < __y);
}

#endif /* __STL_USE_SEPARATE_RELOPS_NAMESPACE */

__STL_END_NAMESPACE

#  undef __STL_STACK_ARGS
#  undef __STL_STACK_HEADER_ARGS

#endif /* __SGI_STL_INTERNAL_STACK_H */

// Local Variables:
// mode:C++
// End:
