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
// release       : 4.3
// release_date  :  6 Apr 2000
//
// file          : include/CGAL/Cartesian/constructions_on_planes_3.h
// package       : C3 (4.3)
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve Bronnimann
// coordinator   : INRIA Sophia-Antipolis (Herve.Bronnimann@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_CONSTRUCTIONS_ON_PLANES_3_H
#define CGAL_CARTESIAN_CONSTRUCTIONS_ON_PLANES_3_H

#include <CGAL/Cartesian/redefine_names_3.h>
#include <CGAL/Cartesian/Point_3.h>
#include <CGAL/constructions/kernel_ftC3.h>

CGAL_BEGIN_NAMESPACE

template <class R>
CGAL_KERNEL_LARGE_INLINE
PlaneC3<R CGAL_CTAG>
plane_from_points(const PointC3<R CGAL_CTAG>& p,
                  const PointC3<R CGAL_CTAG>& q,
		  const PointC3<R CGAL_CTAG>& r)
{
  typename R::FT a,b,c,d;
  plane_from_pointsC3(p.x(),p.y(),p.z(),
                      q.x(),q.y(),q.z(),
                      r.x(),r.y(),r.z(),
                      a,b,c,d);
  return PlaneC3<R CGAL_CTAG>(a,b,c,d);
}

template <class R>
CGAL_KERNEL_LARGE_INLINE
PlaneC3<R CGAL_CTAG>
plane_from_point_direction(const PointC3<R CGAL_CTAG>& p,
                           const DirectionC3<R CGAL_CTAG>& d)
{
  typename R::FT A, B, C, D;
  plane_from_point_directionC3(p.x(),p.y(),p.z(),d.dx(),d.dy(),d.dz(),
                               A,B,C,D);
  return PlaneC3<R CGAL_CTAG>(A,B,C,D);
}

template <class R>
CGAL_KERNEL_LARGE_INLINE
PointC3<R CGAL_CTAG>
point_on_plane(const PlaneC3<R CGAL_CTAG>& p)
{
  typename R::FT x, y, z;
  point_on_planeC3(p.a(),p.b(),p.c(),p.d(),x,y,z);
  return PointC3<R CGAL_CTAG>(x,y,z);
}

template <class R>
CGAL_KERNEL_LARGE_INLINE
PointC3<R CGAL_CTAG>
projection_plane(const PointC3<R CGAL_CTAG>& p,
                 const PlaneC3<R CGAL_CTAG>& h)
{
  typename R::FT x,y,z;
  projection_planeC3(h.a(),h.b(),h.c(),h.d(),
                     p.x(),p.y(),p.z(),
                     x,y,z);
  return PointC3<R CGAL_CTAG>(x,y,z);
}

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_CONSTRUCTIONS_ON_PLANES_3_H
