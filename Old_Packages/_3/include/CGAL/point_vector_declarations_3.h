
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
// file          : include/CGAL/point_vector_declarations_3.h
// package       : _3
// revision      : $Revision$
// revision_date : $Date$
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
point_to_vector_conversion(const Point_3<R>& p)
{ return Vector_3<R>(p); }

template < class R >
inline
Point_3<R>
vector_to_point_conversion(const Vector_3<R>& v)
{ return Point_3<R>(v); }

template < class R >
inline
Point_3<R>
operator+(const Point_3<R>& p, const Vector_3<R>& v)
{
  typedef typename  R::Point_3_base  RPoint_3;
  typedef typename  R::Vector_3_base  RVector_3;
  return Point_3<R>(static_cast<const RPoint_3&>(p) +
                    static_cast<const RVector_3&>(v)) ;
}

template < class R >
inline
Point_3<R>
operator-(const Point_3<R>& p, const Vector_3<R>& v)
{
  typedef typename  R::Point_3_base  RPoint_3;
  typedef typename  R::Vector_3_base  RVector_3;
  return Point_3<R>(static_cast<const RPoint_3&>(p) -
                    static_cast<const RVector_3&>(v));
}

template < class R >
inline
Point_3<R>
operator+(const Origin& , const Vector_3<R>& v)
{ return vector_to_point_conversion(v) ; }

template < class R >
inline
Point_3<R>
operator-(const Origin& , const Vector_3<R>& v)
{ return vector_to_point_conversion(-v) ; }

template < class R >
inline
Vector_3<R>
operator-(const Point_3<R>& p, const Point_3<R>& q)
{
  typedef typename  R::Point_3_base  RPoint_3;
  return Vector_3<R>(static_cast<const RPoint_3&>(p) -
                     static_cast<const RPoint_3&>(q));
}

template < class R >
inline
Vector_3<R>
operator-(const Point_3<R>& p, const Origin& )
{ return point_to_vector_conversion(p) ; }

template < class R >
inline
Vector_3<R>
operator-(const Origin& , const Point_3<R>& p)
{
  typedef typename  R::Point_3_base  RPoint_3;
  return Vector_3<R>(ORIGIN - static_cast<const RPoint_3&>(p));
}

CGAL_END_NAMESPACE

#endif // CGAL_POINT_VECTOR_DECLARATIONS_3_H
