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
// release_date  : 
//
// file          : include/CGAL/ch_value_type.h
// package       : Convex_hull_2 
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
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
