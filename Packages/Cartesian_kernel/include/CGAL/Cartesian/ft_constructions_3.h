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
// file          : include/CGAL/Cartesian/ft_constructions_3.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve Bronnimann
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_FT_CONSTRUCTIONS_3_H
#define CGAL_CARTESIAN_FT_CONSTRUCTIONS_3_H

#include <CGAL/Cartesian/redefine_names_3.h>
#include <CGAL/Cartesian/Point_3.h>
#include <CGAL/Cartesian/Vector_3.h>
#include <CGAL/Cartesian/Plane_3.h>
#include <CGAL/constructions/kernel_ftC3.h>

CGAL_BEGIN_NAMESPACE

template < class R >
inline
typename R::FT
squared_distance(const PointC3<R CGAL_CTAG> &p,
                 const PointC3<R CGAL_CTAG> &q)
{
  return squared_distanceC3(p.x(), p.y(), p.z(), q.x(), q.y(), q.z());
}

template < class R >
inline
typename R::FT
scaled_distance_to_plane(const PlaneC3<R CGAL_CTAG> &h,
                         const PointC3<R CGAL_CTAG> &p)
{
  return scaled_distance_to_planeC3(h.a(), h.b(), h.c(), h.d(),
                                    p.x(), p.y(), p.z());
}

template < class R >
inline
typename R::FT
scaled_distance_to_plane(const PointC3<R CGAL_CTAG> &hp,
                         const PointC3<R CGAL_CTAG> &hq,
                         const PointC3<R CGAL_CTAG> &hr,
                         const PointC3<R CGAL_CTAG> &p)
{
  return scaled_distance_to_planeC3(hp.x(), hp.y(), hp.z(),
                                    hq.x(), hq.y(), hq.z(),
                                    hr.x(), hr.y(), hr.z(),
                                    p.x(), p.y(), p.z());
}

template < class R >
inline
typename R::FT
squared_radius(const PointC3<R CGAL_CTAG> &p,
                     const PointC3<R CGAL_CTAG> &q,
                     const PointC3<R CGAL_CTAG> &r)
{
  return squared_radiusC3(p.x(), p.y(), p.z(),
	                  q.x(), q.y(), q.z(),
			  r.x(), r.y(), r.z());
}

template < class R >
inline
typename R::FT
squared_radius(const PointC3<R CGAL_CTAG> &p,
                     const PointC3<R CGAL_CTAG> &q,
                     const PointC3<R CGAL_CTAG> &r,
                     const PointC3<R CGAL_CTAG> &s)
{
  return squared_radiusC3(p.x(), p.y(), p.z(),
	                  q.x(), q.y(), q.z(),
			  r.x(), r.y(), r.z(),
			  s.x(), s.y(), s.z());
}

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_FT_CONSTRUCTIONS_3_H
