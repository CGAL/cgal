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

#include <CGAL/Cartesian/Point_3.h>
#include <CGAL/Cartesian/Vector_3.h>
#include <CGAL/Cartesian/Plane_3.h>
#include <CGAL/constructions/kernel_ftC3.h>

CGAL_BEGIN_NAMESPACE

template < class K >
inline
typename K::FT
squared_distance(const PointC3<K> &p,
                 const PointC3<K> &q)
{
  return squared_distanceC3(p.x(), p.y(), p.z(), q.x(), q.y(), q.z());
}

template < class K >
inline
typename K::FT
scaled_distance_to_plane(const PlaneC3<K> &h,
                         const PointC3<K> &p)
{
  return scaled_distance_to_planeC3(h.a(), h.b(), h.c(), h.d(),
                                    p.x(), p.y(), p.z());
}

template < class K >
inline
typename K::FT
scaled_distance_to_plane(const PointC3<K> &hp,
                         const PointC3<K> &hq,
                         const PointC3<K> &hr,
                         const PointC3<K> &p)
{
  return scaled_distance_to_planeC3(hp.x(), hp.y(), hp.z(),
                                    hq.x(), hq.y(), hq.z(),
                                    hr.x(), hr.y(), hr.z(),
                                    p.x(), p.y(), p.z());
}

template < class K >
inline
typename K::FT
squared_radius(const PointC3<K> &p, const PointC3<K> &q,
	       const PointC3<K> &r, const PointC3<K> &s)
{
  return K().compute_squared_radius_3_object()(p, q, r, s);
}

template < class K >
inline
typename K::FT
squared_radius(const PointC3<K> &p, const PointC3<K> &q, const PointC3<K> &r)
{
  return K().compute_squared_radius_3_object()(p, q, r);
}

template < class K >
inline
typename K::FT
squared_radius(const PointC3<K> &p, const PointC3<K> &q)
{
  return K().compute_squared_radius_3_object()(p, q);
}

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_FT_CONSTRUCTIONS_3_H
