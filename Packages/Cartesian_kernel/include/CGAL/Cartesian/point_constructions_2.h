// ======================================================================
//
// Copyright (c) 2000 The CGAL Consortium
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
// file          : include/CGAL/Cartesian/point_constructions_2.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve Bronnimann
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_POINT_CONSTRUCTIONS_2_H
#define CGAL_CARTESIAN_POINT_CONSTRUCTIONS_2_H

#include <CGAL/Cartesian/Point_2.h>

CGAL_BEGIN_NAMESPACE

template < class R >
inline
PointC2<R>
midpoint(const PointC2<R> &p,
         const PointC2<R> &q )
{
  typename R::FT x, y;
  midpointC2(p.x(), p.y(), q.x(), q.y(), x, y);
  return PointC2<R>(x, y);
}

template < class R >
inline
PointC2<R>
circumcenter(const PointC2<R> &p,
             const PointC2<R> &q,
             const PointC2<R> &r)
{
  typename R::FT x, y;
  circumcenterC2(p.x(), p.y(), q.x(), q.y(), r.x(), r.y(), x, y);
  return PointC2<R>(x, y);
}

template < class R >
inline
PointC2<R>
centroid(const PointC2<R> &p,
         const PointC2<R> &q,
         const PointC2<R> &r)
{
  typename R::FT x, y;
  centroidC2(p.x(), p.y(), q.x(), q.y(), r.x(), r.y(), x, y);
  return PointC2<R>(x, y);
}

template < class R >
inline
PointC2<R>
centroid(const PointC2<R> &p,
         const PointC2<R> &q,
         const PointC2<R> &r,
         const PointC2<R> &s)
{
  typename R::FT x, y;
  centroidC2(p.x(), p.y(), q.x(), q.y(), r.x(), r.y(), s.x(), s.y(), x, y);
  return PointC2<R>(x, y);
}

template < class R >
inline
PointC2<R>
line_get_point(const LineC2<R> &l, int i)
{
  typename R::FT x, y;
  line_get_pointC2(l.a(), l.b(), l.c(), i, x, y);
  return PointC2<R>(x, y);
}

template < class R >
inline
PointC2<R>
line_project_point(const LineC2<R> &l,
                   const PointC2<R> &p)
{
  typename R::FT x, y;
  line_project_pointC2(l.a(), l.b(), l.c(), p.x(), p.y(), x, y);
  return PointC2<R>(x, y);
}

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_POINT_CONSTRUCTIONS_2_H
