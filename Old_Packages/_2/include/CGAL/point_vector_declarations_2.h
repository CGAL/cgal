// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// 
// release       : 
// release_date  : 
// 
// file          : include/CGAL/point_vector_declarations_2.h
// package       : _2
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_POINT_VECTOR_DECLARATIONS_2_H
#define CGAL_POINT_VECTOR_DECLARATIONS_2_H

#include <CGAL/user_classes.h>

CGAL_BEGIN_NAMESPACE

template < class R >
inline
Vector_2<R>
point_to_vector_conversion(const Point_2<R>& p);

template < class R >
inline
Point_2<R>
vector_to_point_conversion(const Vector_2<R>& v);

template < class R >
inline
Point_2<R>
operator+(const Point_2<R>& p, const Vector_2<R>& v);

template < class R >
inline
Point_2<R>
operator-(const Point_2<R>& p, const Vector_2<R>& v);

class Origin;

template < class R >
inline
Point_2<R>
operator+(const Origin& , const Vector_2<R>& v);

template < class R >
inline
Point_2<R>
operator-(const Origin& , const Vector_2<R>& v);

template < class R >
inline
Vector_2<R>
operator-(const Point_2<R>& p, const Point_2<R>& q);

template < class R >
inline
Vector_2<R>
operator-(const Point_2<R>& p, const Origin& );

template < class R >
inline
Vector_2<R>
operator-(const Origin& , const Point_2<R>& p);
CGAL_END_NAMESPACE

#endif // CGAL_POINT_VECTOR_DECLARATIONS_2_H
