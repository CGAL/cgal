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
// file          : include/CGAL/Cartesian/line_constructions_2.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve Bronnimann
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_LINE_CONSTRUCTIONS_2_H
#define CGAL_CARTESIAN_LINE_CONSTRUCTIONS_2_H

#include <CGAL/Cartesian/Point_2.h>
#include <CGAL/Cartesian/Line_2.h>

CGAL_BEGIN_NAMESPACE

template < class K >
inline
LineC2<K>
line_from_points(const PointC2<K> &p,
                 const PointC2<K> &q)
{
  return K().construct_line_2_object()(p, q);
}

template < class K >
inline
LineC2<K>
line_from_point_direction(const PointC2<K> &p,
                          const DirectionC2<K> &d)
{
  return K().construct_line_2_object()(p, d);
}

template < class K >
inline
LineC2<K>
bisector(const PointC2<K> &p,
         const PointC2<K> &q)
{
  return K().construct_bisector_2_object()(p, q);
}

template < class K >
inline
LineC2<K>
perpendicular_through_point(const LineC2<K> &l,
                            const PointC2<K> &p)
{
  typename K::FT a, b, c;
  perpendicular_through_pointC2(l.a(), l.b(), p.x(), p.y(), a, b, c);
  return LineC2<K>(a, b, c);
}

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_LINE_CONSTRUCTIONS_2_H
