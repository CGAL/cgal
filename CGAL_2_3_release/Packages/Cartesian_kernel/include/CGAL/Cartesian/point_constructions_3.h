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
// file          : include/CGAL/Cartesian/point_constructions_3.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve Bronnimann
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_POINT_CONSTRUCTIONS_3_H
#define CGAL_CARTESIAN_POINT_CONSTRUCTIONS_3_H

#include <CGAL/Cartesian/redefine_names_3.h>
#include <CGAL/Cartesian/Point_3.h>
#include <CGAL/constructions/kernel_ftC3.h>

CGAL_BEGIN_NAMESPACE

template < class R >
CGAL_KERNEL_MEDIUM_INLINE
PointC3<R CGAL_CTAG>
midpoint(const PointC3<R CGAL_CTAG> &p,
         const PointC3<R CGAL_CTAG> &q)
{
  typename R::FT x, y, z;
  midpointC3(p.x(), p.y(), p.z(), q.x(), q.y(), q.z(), x, y, z);
  return PointC3<R CGAL_CTAG>(x, y, z);
}

template < class R >
PointC3<R CGAL_CTAG>
centroid(const PointC3<R CGAL_CTAG> &p,
         const PointC3<R CGAL_CTAG> &q,
         const PointC3<R CGAL_CTAG> &r,
         const PointC3<R CGAL_CTAG> &s)
{
  typename R::FT x, y, z;
  centroidC3(p.x(), p.y(), p.z(),
             q.x(), q.y(), q.z(),
             r.x(), r.y(), r.z(),
             s.x(), s.y(), s.z(),
             x, y, z);
  return PointC3<R CGAL_CTAG>(x, y, z);
}

template < class R >
PointC3<R CGAL_CTAG>
centroid(const PointC3<R CGAL_CTAG> &p,
         const PointC3<R CGAL_CTAG> &q,
         const PointC3<R CGAL_CTAG> &r)
{
  typename R::FT x, y, z;
  centroidC3(p.x(), p.y(), p.z(),
             q.x(), q.y(), q.z(),
             r.x(), r.y(), r.z(),
             x, y, z);
  return PointC3<R CGAL_CTAG>(x, y, z);
}

template < class R >
PointC3<R CGAL_CTAG>
circumcenter(const PointC3<R CGAL_CTAG> &p,
             const PointC3<R CGAL_CTAG> &q,
             const PointC3<R CGAL_CTAG> &r,
             const PointC3<R CGAL_CTAG> &s)
{
  typename R::FT x, y, z;
  circumcenterC3(p.x(), p.y(), p.z(),
                 q.x(), q.y(), q.z(),
                 r.x(), r.y(), r.z(),
                 s.x(), s.y(), s.z(),
                 x, y, z);
  return PointC3<R CGAL_CTAG>(x, y, z);
}

template < class R >
PointC3<R CGAL_CTAG>
circumcenter(const PointC3<R CGAL_CTAG> &p,
             const PointC3<R CGAL_CTAG> &q,
             const PointC3<R CGAL_CTAG> &r)
{
  typename R::FT x, y, z;
  circumcenterC3(p.x(), p.y(), p.z(),
                 q.x(), q.y(), q.z(),
                 r.x(), r.y(), r.z(),
                 x, y, z);
  return PointC3<R CGAL_CTAG>(x, y, z);
}

template <class R>
CGAL_KERNEL_LARGE_INLINE
PointC3<R CGAL_CTAG>
point_on_line(int i, const LineC3<R CGAL_CTAG> &l)
{
  typename R::FT x, y, z;
  point_on_lineC3(l.point().x(), l.point().y(), l.point().z(),
                  l.direction().dx(), l.direction().dy(), l.direction().dz(),
                  i, x, y, z);
  return PointC3<R CGAL_CTAG>(x, y, z);
}

template <class R>
CGAL_KERNEL_LARGE_INLINE
PointC3<R CGAL_CTAG>
projection_line(const PointC3<R CGAL_CTAG> &p, const LineC3<R CGAL_CTAG> &l)
{
  typename R::FT x, y, z;
  projection_lineC3(p.x(), p.y(), p.z(),
		    l.point().x(), l.point().y(), l.point().z(),
                    l.direction().dx(), l.direction().dy(), l.direction().dz(),
                    x, y, z);
  return PointC3<R CGAL_CTAG>(x, y, z);
}

template <class R>
CGAL_KERNEL_LARGE_INLINE
PointC3<R CGAL_CTAG>
point_on_plane(const PlaneC3<R CGAL_CTAG> &p)
{
  typename R::FT x, y, z;
  point_on_planeC3(p.a(), p.b(), p.c(), p.d(), x, y, z);
  return PointC3<R CGAL_CTAG>(x, y, z);
}

template <class R>
CGAL_KERNEL_LARGE_INLINE
PointC3<R CGAL_CTAG>
projection_plane(const PointC3<R CGAL_CTAG> &p,
                 const PlaneC3<R CGAL_CTAG> &h)
{
  typename R::FT x, y, z;
  projection_planeC3(h.a(), h.b(), h.c(), h.d(),
                     p.x(), p.y(), p.z(),
                     x, y, z);
  return PointC3<R CGAL_CTAG>(x, y, z);
}

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_POINT_CONSTRUCTIONS_3_H
