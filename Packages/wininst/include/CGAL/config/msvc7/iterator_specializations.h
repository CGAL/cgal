// ============================================================================
//
// Copyright (c) 1997-2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// file          : include/CGAL/config/msvc7/iterator_specializations.h
// package       : wininst
// author(s)     : 
// release       : 
// release_date  : 
//
// coordinator   : 
//
// ============================================================================

#ifndef CGAL_ITER_VC7
#define CGAL_ITER_VC7
#include <iterator>
#include <CGAL/config/msvc7/stl_iterator_base.h>
namespace std {
template<class C__> inline
typename iterator_traits<C__>::iterator_category
_Iter_cat(const C__&)
  {
    typedef typename iterator_traits<C__>::iterator_category c;
    return c();
  }

template <class _Iter> inline 
typename iterator_traits<_Iter>::difference_type*
  _Dist_type(const _Iter&)
  {
    typedef typename iterator_traits<_Iter>::difference_type _diff_type;
    return static_cast<_diff_type*>(0);
  }

template <class _Iter> inline 
typename iterator_traits<_Iter>::value_type*
  _Val_type(const _Iter&)
{
  typedef typename iterator_traits<_Iter>::value_type _value_type;
  return static_cast<_value_type*>(0);
}


}
#endif // CGAL_ITER_VC7
