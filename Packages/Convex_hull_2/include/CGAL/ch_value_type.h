// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// release       : 
// release_date  : 2000, August 03
//
// file          : include/CGAL/ch_value_type.h
// package       : Convex_hull (3.3)
// maintainer    : Stefan Schirra <stschirr@mpi-sb.mpg.de>
// source        : stl_extensions.lw
// revision      : 3.2.1
// revision_date : 19 Apr 2000
// author(s)     : Stefan Schirra <Stefan.Schirra@@mpi-sb.mpg.de>
//
// coordinator   : MPI, Saarbruecken
// ======================================================================

#ifndef CGAL_CH_VALUE_TYPE_H
#define CGAL_CH_VALUE_TYPE_H

CGAL_BEGIN_NAMESPACE

#ifndef CGAL_CFG_NO_ITERATOR_TRAITS
template <class Iterator>
inline 
typename std::iterator_traits<Iterator>::value_type*
ch_value_type(const Iterator&) 
{
  typedef typename std::iterator_traits<Iterator>::value_type  v_type;
  return static_cast<v_type*>(0);
}
#else
#ifdef _MSC_VER
template<class Iterator> 
inline
Iterator::value_type* 
ch_value_type(const Iterator& it)
{return static_cast<Iterator::value_type*>(0); }

template<class T> 
inline
T*
ch_value_type(const T*) 
{return static_cast<T*>(0); }
#else 
#error To be fixed
#endif // _MSC_VER_
#endif // CGAL_CFG_NO_ITERATOR_TRAITS
CGAL_END_NAMESPACE
#endif // CGAL_CH_VALUE_TYPE_H
