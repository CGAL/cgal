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
// file          : basic_constructions_3.h
// package       : _3
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 
#ifndef CGAL_BASIC_CONSTRUCTIONS_3_H
#define CGAL_BASIC_CONSTRUCTIONS_3_H

CGAL_BEGIN_NAMESPACE

template < class R >
inline
Point_3<R>
midpoint( const Point_3<R>& p,
          const Point_3<R>& q )
{
    typedef typename R::Point_3_base  RPoint_3;
    return( midpoint( static_cast<const RPoint_3&>(p),
                      static_cast<const RPoint_3&>(q)));
}

template < class R >
inline
Point_3<R>
centroid( const Point_3<R>& p,
          const Point_3<R>& q,
          const Point_3<R>& r,
          const Point_3<R>& s)
{
    typedef typename R::Point_3_base  RPoint_3;
    return( centroid( static_cast<const RPoint_3&>(p),
                      static_cast<const RPoint_3&>(q),
                      static_cast<const RPoint_3&>(r),
                      static_cast<const RPoint_3&>(s)));
}

template < class R >
inline
Point_3<R>
centroid( const Point_3<R>& p,
          const Point_3<R>& q,
          const Point_3<R>& r)
{
    typedef typename R::Point_3_base  RPoint_3;
    return( centroid( static_cast<const RPoint_3&>(p),
                      static_cast<const RPoint_3&>(q),
                      static_cast<const RPoint_3&>(r)));
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
    return( circumcenter( static_cast<const RPoint_3&>(p),
                          static_cast<const RPoint_3&>(q),
                          static_cast<const RPoint_3&>(r),
                          static_cast<const RPoint_3&>(s)));
}

CGAL_END_NAMESPACE

#endif // CGAL_BASIC_CONSTRUCTIONS_3_H
