// ============================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
// 
// release       :
// release_date  :
// 
// source        : basic_constructions_3.fw
// file          : basic_constructions_3.h
// revision      : 2.4
// revision_date : 24 Aug 1999 
// author(s)     : Stefan Schirra <Stefan.Schirra@mpi-sb.mpg.de>
//
// coordinator   : MPI, Saarbruecken
// ============================================================================


#ifndef CGAL_BASIC_CONSTRUCTIONS_3_H
#define CGAL_BASIC_CONSTRUCTIONS_3_H

#ifdef CGAL_HOMOGENEOUS_H
#include <CGAL/basic_constructionsH3.h>
#endif // CGAL_HOMOGENEOUS_H

#ifdef CGAL_CARTESIAN_H
#include <CGAL/basic_constructionsC3.h>
#endif // CGAL_CARTESIAN_H

CGAL_BEGIN_NAMESPACE

template < class R >
inline
Point_3<R>
midpoint( const Point_3<R>& p,
               const Point_3<R>& q )
{
    typedef typename R::Point_3_base  RPoint_3;
    return( midpoint( CGAL_static_cast( const RPoint_3& , p),
                           CGAL_static_cast( const RPoint_3& , q)));
}
template < class R >
inline
Point_3<R>
circumcenter( const Point_3<R>& p,
                   const Point_3<R>& q,
                   const Point_3<R>& r,
                   const Point_3<R>& s)
{
    typedef typename R::Point_3_base  RPoint_3;
    return( circumcenter( CGAL_static_cast( const RPoint_3& , p),
                               CGAL_static_cast( const RPoint_3& , q),
                               CGAL_static_cast( const RPoint_3& , r),
                               CGAL_static_cast( const RPoint_3& , s)));
}
CGAL_END_NAMESPACE


#endif // CGAL_BASIC_CONSTRUCTIONS_3_H
