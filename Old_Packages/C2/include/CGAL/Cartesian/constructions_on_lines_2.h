// ============================================================================
//
// Copyright (c) 1998, 1999 The CGAL Consortium
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
// file          : include/CGAL/Cartesian/constructions_on_lines_2.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve Bronnimann
//
// coordinator   : INRIA Sophia-Antipolis (Herve.Bronnimann@sophia.inria.fr)
//
// ============================================================================

#ifndef CGAL_CARTESIAN_CONSTRUCTIONS_ON_LINES_2_H
#define CGAL_CARTESIAN_CONSTRUCTIONS_ON_LINES_2_H

#include <CGAL/Cartesian/Point_2.h>
#include <CGAL/Cartesian/Line_2.h>

CGAL_BEGIN_NAMESPACE

template < class R >
inline
LineC2<R CGAL_CTAG>
line_from_points( PointC2<R CGAL_CTAG> const& p,
                  PointC2<R CGAL_CTAG> const& q)
{
  typename R::FT a,b,c;
  line_from_pointsC2(p.x(),p.y(),q.x(),q.y(),a,b,c);
  return LineC2<R CGAL_CTAG>(a,b,c);
}

template < class R >
inline
LineC2<R CGAL_CTAG>
line_from_point_direction( PointC2<R CGAL_CTAG> const& p,
                           DirectionC2<R CGAL_CTAG> const& d)
{
  typename R::FT a,b,c;
  line_from_point_directionC2(p.x(),p.y(),d.dx(),d.dy(),a,b,c);
  return LineC2<R CGAL_CTAG>(a,b,c);
}

template < class R >
inline
typename R::FT
line_y_at_x( LineC2<R CGAL_CTAG> const& l,
             typename R::FT const& x)
{
  return line_y_at_xC2(l.a(), l.b(), l.c(), x);
}

template < class R >
inline
typename R::FT
line_x_at_y( LineC2<R CGAL_CTAG> const& l,
             typename R::FT const& y)
{
  return line_y_at_xC2(l.b(), l.a(), l.c(), y);
}

template < class R >
inline
PointC2<R>
line_get_point( LineC2<R CGAL_CTAG> const& l, int i)
{
  typename R::FT x, y;
  line_get_pointC2(l.a(), l.b(), l.c(), i, x, y);
  return PointC2<R>(x,y);
}

template < class R >
inline
PointC2<R>
line_project_point( LineC2<R CGAL_CTAG> const& l,
                    PointC2<R CGAL_CTAG> const& p)
{
  typename R::FT x, y;
  line_project_pointC2(l.a(), l.b(), l.c(), p.x(), p.y());
  return PointC2<R>(x,y);
}

template < class R >
inline
Oriented_side
side_of_oriented_line( LineC2<R CGAL_CTAG> const& l,
                       PointC2<R CGAL_CTAG> const& p)
{
  return side_of_oriented_lineC2(l.a(), l.b(), l.c(), p.x(), p.y());
}

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_CONSTRUCTIONS_ON_LINES_2_H
