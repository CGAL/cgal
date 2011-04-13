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
// file          : include/CGAL/Cartesian/ft_constructions_2.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve Bronnimann
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_FT_CONSTRUCTIONS_2_H
#define CGAL_CARTESIAN_FT_CONSTRUCTIONS_2_H

#include <CGAL/Cartesian/redefine_names_2.h>

CGAL_BEGIN_NAMESPACE

template < class R >
inline
typename R::FT
squared_distance(const PointC2<R CGAL_CTAG> &p,
                 const PointC2<R CGAL_CTAG> &q)
{
  return squared_distanceC2(p.x(), p.y(), q.x(), q.y());
}

template < class R >
inline
typename R::FT
scaled_distance_to_line(const LineC2<R CGAL_CTAG> &l,
                        const PointC2<R CGAL_CTAG> &p)
{
  return scaled_distance_to_lineC2(l.a(), l.b(), l.c(), p.x(), p.y());
}

template < class R >
inline
typename R::FT
scaled_distance_to_line(const PointC2<R CGAL_CTAG> &p,
                        const PointC2<R CGAL_CTAG> &q,
                        const PointC2<R CGAL_CTAG> &r)
{
  return scaled_distance_to_lineC2(p.x(), p.y(), q.x(), q.y(), r.x(), r.y());
}

template < class R >
inline
typename R::FT
line_y_at_x(const LineC2<R CGAL_CTAG> &l, const typename R::FT &x)
{
  return line_y_at_xC2(l.a(), l.b(), l.c(), x);
}

template < class R >
inline
typename R::FT
line_x_at_y(const LineC2<R CGAL_CTAG> &l, const typename R::FT &y)
{
  return line_y_at_xC2(l.b(), l.a(), l.c(), y);
}

template < class R >
inline
typename R::FT
squared_radius(const PointC2<R CGAL_CTAG> &p,
                     const PointC2<R CGAL_CTAG> &q,
                     const PointC2<R CGAL_CTAG> &r)
{
  return squared_radiusC2(p.x(), p.y(), q.x(), q.y(), r.x(), r.y());
}

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_FT_CONSTRUCTIONS_2_H
