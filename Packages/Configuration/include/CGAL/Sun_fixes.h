// ======================================================================
//
// Copyright (c) 1997-2002 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-2.4-I-63 $
// release_date  : $CGAL_Date: 2002/03/15 $
//
// file          : include/CGAL/Sun_fixes.h
// package       : Configuration (2.28)
// maintainer    : Geert-Jan Giezeman <geert@cs.uu.nl>
// source        :
// revision      : 1.11
// revision_date : 30 Mar 1998
// author(s)     : Michael Hoffmann (hoffmann@inf.ethz.ch)
//
// coordinator   : Utrecht University
//
// ======================================================================

#ifndef CGAL_SUN_FIXES_H
#define CGAL_SUN_FIXES_H 1

//----------------------------------------------------------------------//
//             enable member templates for stdlib                       //
//----------------------------------------------------------------------//

#include <stdcomp.h>
#undef RWSTD_NO_MEMBER_TEMPLATES
#undef _RWSTD_NO_MEMBER_TEMPLATES

//----------------------------------------------------------------------//
//             fake iterator_traits                                     //
//----------------------------------------------------------------------//

#include <iterator>
namespace std {
  template <class Iterator> struct iterator_traits
  {
    typedef typename Iterator::value_type value_type;
    typedef typename Iterator::difference_type difference_type;
    typedef typename Iterator::pointer pointer;
    typedef typename Iterator::reference reference;
    typedef typename Iterator::iterator_category iterator_category;
  };
  template <class T> struct iterator_traits<T*>
  {
    typedef T value_type;
    typedef ptrdiff_t difference_type;
    typedef T* pointer;
    typedef T& reference;
    typedef random_access_iterator_tag iterator_category;
  };
  template <class T> struct iterator_traits<const T*>
  {
    typedef T value_type;
    typedef ptrdiff_t difference_type;
    typedef const T* pointer;
    typedef const T& reference;
    typedef random_access_iterator_tag iterator_category;
  };
  template <class ForwardIterator>
  inline ptrdiff_t
  distance (ForwardIterator first, ForwardIterator last)
  {
    ptrdiff_t n = 0;
    __distance(first, last, n, 
               iterator_traits<ForwardIterator>::iterator_category());
    return n;
  }

  template < class T >
  inline typename T::value_type*
  __value_type (const T&)
  { return (typename T::value_type*)(0); }

  template < class T >
  inline typename T::difference_type*
  __distance_type(const T&)
  { return (typename T::difference_type*)(0); }

  template < class T >
  inline typename T::iterator_category
  __iterator_category (const T&)
  {
    typename T::iterator_category tmp;
    return tmp;
  }
} // namespace std

#endif // CGAL_SUN_FIXES_H
