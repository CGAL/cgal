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
  
  // af: added the following:
  struct iterator_traits<unsigned int> {
    typedef _Int_iterator_tag iterator_category;
  };
  

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


#if defined (__INTEL_COMPILER)
template <class _Tp>
struct iterator_traits<const _Tp*> {
  typedef random_access_iterator_tag iterator_category;
  typedef _Tp                         value_type;
  typedef ptrdiff_t                   difference_type;
  typedef const _Tp*                  pointer;
  typedef const _Tp&                  reference;
};

template <class _Tp>
struct iterator_traits<_Tp*> {
  typedef random_access_iterator_tag iterator_category;
  typedef _Tp                         value_type;
  typedef ptrdiff_t                   difference_type;
  typedef _Tp*                        pointer;
  typedef _Tp&                        reference;
};
#endif // defined (__INTEL_COMPILER)

}
#endif // CGAL_ITER_VC7
