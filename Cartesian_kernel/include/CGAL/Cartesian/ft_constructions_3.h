// Copyright (c) 2000
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Herve Bronnimann

#ifndef CGAL_CARTESIAN_FT_CONSTRUCTIONS_3_H
#define CGAL_CARTESIAN_FT_CONSTRUCTIONS_3_H

#include <CGAL/Cartesian/Point_3.h>
#include <CGAL/Cartesian/Vector_3.h>
#include <CGAL/Cartesian/Plane_3.h>
#include <CGAL/constructions/kernel_ftC3.h>

namespace CGAL {

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

} //namespace CGAL

#endif // CGAL_CARTESIAN_FT_CONSTRUCTIONS_3_H
