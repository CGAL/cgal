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

CGAL_BEGIN_NAMESPACE

template < class K >
inline
typename K::FT
squared_distance(const PointC2<K> &p,
                 const PointC2<K> &q)
{
  return squared_distanceC2(p.x(), p.y(), q.x(), q.y());
}

template < class K >
inline
typename K::FT
scaled_distance_to_line(const LineC2<K> &l,
                        const PointC2<K> &p)
{
  return scaled_distance_to_lineC2(l.a(), l.b(), l.c(), p.x(), p.y());
}

template < class K >
inline
typename K::FT
scaled_distance_to_line(const PointC2<K> &p,
                        const PointC2<K> &q,
                        const PointC2<K> &r)
{
  return scaled_distance_to_lineC2(p.x(), p.y(), q.x(), q.y(), r.x(), r.y());
}

template < class K >
inline
typename K::FT
line_y_at_x(const LineC2<K> &l, const typename K::FT &x)
{
  return line_y_at_xC2(l.a(), l.b(), l.c(), x);
}

template < class K >
inline
typename K::FT
line_x_at_y(const LineC2<K> &l, const typename K::FT &y)
{
  return line_y_at_xC2(l.b(), l.a(), l.c(), y);
}

template < class K >
inline
typename K::FT
squared_radius(const PointC2<K> &p, const PointC2<K> &q, const PointC2<K> &r)
{
  return K().compute_squared_radius_2_object()(p, q, r);
}

template < class K >
inline
typename K::FT
squared_radius(const PointC2<K> &p, const PointC2<K> &q)
{
  return K().compute_squared_radius_2_object()(p, q);
}

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_FT_CONSTRUCTIONS_2_H
