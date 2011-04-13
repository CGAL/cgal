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
// file          : basic_constructions_2.h
// package       : _2
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Sven Schoenherr
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================

#ifndef CGAL_BASIC_CONSTRUCTIONS_2_H
#define CGAL_BASIC_CONSTRUCTIONS_2_H

CGAL_BEGIN_NAMESPACE

template < class R >
inline
Point_2<R>
midpoint( const Point_2<R>& p,
          const Point_2<R>& q )
{
    typedef typename R::Point_2_base  RPoint_2;
    return( midpoint( static_cast<const RPoint_2&>(p),
                      static_cast<const RPoint_2&>(q)));
}

template < class R >
inline
Point_2<R>
circumcenter( const Point_2<R>& p,
              const Point_2<R>& q,
              const Point_2<R>& r)
{
    typedef typename R::Point_2_base  RPoint_2;
    return( circumcenter( static_cast<const RPoint_2&>(p),
                          static_cast<const RPoint_2&>(q),
                          static_cast<const RPoint_2&>(r)));
}

template < class R >
inline
Point_2<R>
centroid( const Point_2<R>& p,
          const Point_2<R>& q,
          const Point_2<R>& r)
{
    typedef typename R::Point_2_base  RPoint_2;
    return( centroid( static_cast<const RPoint_2&>(p),
                      static_cast<const RPoint_2&>(q),
                      static_cast<const RPoint_2&>(r)));
}

template < class R >
inline
Point_2<R>
centroid( const Point_2<R>& p,
          const Point_2<R>& q,
          const Point_2<R>& r,
          const Point_2<R>& s)
{
    typedef typename R::Point_2_base  RPoint_2;
    return( centroid( static_cast<const RPoint_2&>(p),
                      static_cast<const RPoint_2&>(q),
                      static_cast<const RPoint_2&>(r),
                      static_cast<const RPoint_2&>(s)));
}

CGAL_END_NAMESPACE

#endif // CGAL_BASIC_CONSTRUCTIONS_2_H
