
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
// release_date  : 2000, October 15
// 
// source        : PV_decl.fw
// file          : include/CGAL/point_vector_declarations_3.h
// package       : _3 (3.9)
// maintainer    : Stefan Schirra <stschirr@mpi-sb.mpg.de>
// revision      : 3.9
// revision_date : 15 Oct 2000 
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_POINT_VECTOR_DECLARATIONS_3_H
#define CGAL_POINT_VECTOR_DECLARATIONS_3_H

#include <CGAL/user_classes.h>

CGAL_BEGIN_NAMESPACE


template < class R >
inline
Vector_3<R>
point_to_vector_conversion(const Point_3<R>& p);

template < class R >
inline
Point_3<R>
vector_to_point_conversion(const Vector_3<R>& v);

template < class R >
inline
Point_3<R>
operator+(const Point_3<R>& p, const Vector_3<R>& v);

template < class R >
inline
Point_3<R>
operator-(const Point_3<R>& p, const Vector_3<R>& v);

template < class R >
inline
Point_3<R>
operator+(const Origin& , const Vector_3<R>& v);

template < class R >
inline
Point_3<R>
operator-(const Origin& , const Vector_3<R>& v);

template < class R >
inline
Vector_3<R>
operator-(const Point_3<R>& p, const Point_3<R>& q);

template < class R >
inline
Vector_3<R>
operator-(const Point_3<R>& p, const Origin& );

template < class R >
inline
Vector_3<R>
operator-(const Origin& , const Point_3<R>& p);
CGAL_END_NAMESPACE


#endif // CGAL_POINT_VECTOR_DECLARATIONS_3_H
