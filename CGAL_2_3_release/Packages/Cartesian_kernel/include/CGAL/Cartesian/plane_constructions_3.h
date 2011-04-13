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
// file          : include/CGAL/Cartesian/plane_constructions_3.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve Bronnimann
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_PLANE_CONSTRUCTIONS_3_H
#define CGAL_CARTESIAN_PLANE_CONSTRUCTIONS_3_H

#include <CGAL/Cartesian/redefine_names_3.h>
#include <CGAL/Cartesian/Point_3.h>
#include <CGAL/constructions/kernel_ftC3.h>

CGAL_BEGIN_NAMESPACE

template <class R>
CGAL_KERNEL_LARGE_INLINE
PlaneC3<R CGAL_CTAG>
plane_from_points(const PointC3<R CGAL_CTAG> &p,
                  const PointC3<R CGAL_CTAG> &q,
		  const PointC3<R CGAL_CTAG> &r)
{
  typename R::FT a, b, c, d;
  plane_from_pointsC3(p.x(), p.y(), p.z(),
                      q.x(), q.y(), q.z(),
                      r.x(), r.y(), r.z(),
                      a, b, c, d);
  return PlaneC3<R CGAL_CTAG>(a, b, c, d);
}

template <class R>
CGAL_KERNEL_LARGE_INLINE
PlaneC3<R CGAL_CTAG>
plane_from_point_direction(const PointC3<R CGAL_CTAG> &p,
                           const DirectionC3<R CGAL_CTAG> &d)
{
  typename R::FT A, B, C, D;
  plane_from_point_directionC3(p.x(), p.y(), p.z(), d.dx(), d.dy(), d.dz(),
                               A, B, C, D);
  return PlaneC3<R CGAL_CTAG>(A, B, C, D);
}

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_PLANE_CONSTRUCTIONS_3_H
