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

template < class K >
inline
PointC2<K>
midpoint(const PointC2<K> &p,
         const PointC2<K> &q )
{
  return K().construct_midpoint_2_object()(p, q);
}

template < class K >
inline
PointC2<K>
circumcenter(const PointC2<K> &p,
             const PointC2<K> &q,
             const PointC2<K> &r)
{
  return K().construct_circumcenter_2_object()(p, q, r);
}

template < class K >
inline
PointC2<K>
centroid(const PointC2<K> &p,
         const PointC2<K> &q,
         const PointC2<K> &r)
{
  return K().construct_centroid_2_object()(p, q, r);
}

template < class K >
inline
PointC2<K>
centroid(const PointC2<K> &p,
         const PointC2<K> &q,
         const PointC2<K> &r,
         const PointC2<K> &s)
{
  return K().construct_centroid_2_object()(p, q, r, s);
}

template < class K >
inline
PointC2<K>
line_get_point(const LineC2<K> &l, int i)
{
  typename K::FT x, y;
  line_get_pointC2(l.a(), l.b(), l.c(), i, x, y);
  return PointC2<K>(x, y);
}

template < class K >
inline
PointC2<K>
line_project_point(const LineC2<K> &l,
                   const PointC2<K> &p)
{
  typename K::FT x, y;
  line_project_pointC2(l.a(), l.b(), l.c(), p.x(), p.y(), x, y);
  return PointC2<K>(x, y);
}

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_POINT_CONSTRUCTIONS_2_H
