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
// file          : include/CGAL/gnu_istream_iterator_value_type_fix.h
// package       : Convex_hull_2
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
// coordinator   : MPI, Saarbruecken
// ======================================================================


#ifndef GNU_ISTREAM_ITERATOR_VALUE_TYPE_FIX_H
#define GNU_ISTREAM_ITERATOR_VALUE_TYPE_FIX_H

#ifdef STL_GCC
#include <iterator>
template <class T, class Distance> 
inline 
T* 
value_type(const istream_iterator<T, Distance>&) 
{ return (T*)(0); }
#endif // STL_GCC

#endif // GNU_ISTREAM_ITERATOR_VALUE_TYPE_FIX_H
