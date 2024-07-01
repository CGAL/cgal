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

#ifndef CGAL_CARTESIAN_POINT_CONSTRUCTIONS_3_H
#define CGAL_CARTESIAN_POINT_CONSTRUCTIONS_3_H

#include <CGAL/Cartesian/Point_3.h>
#include <CGAL/constructions/kernel_ftC3.h>

namespace CGAL {

template <class K>
CGAL_KERNEL_LARGE_INLINE
typename K::Point_3
point_on_plane(const PlaneC3<K> &p)
{
  typename K::FT x, y, z;
  point_on_planeC3(p.a(), p.b(), p.c(), p.d(), x, y, z);
  return K().construct_point_3_object()(x, y, z);
}

template <class K>
CGAL_KERNEL_LARGE_INLINE
typename K::Point_3
projection_plane(const typename K::Point_3 &p,
                 const PlaneC3<K> &h)
{
  typename K::FT x, y, z;
  projection_planeC3(h.a(), h.b(), h.c(), h.d(),
                     p.x(), p.y(), p.z(),
                     x, y, z);
  return K().construct_point_3_object()(x, y, z);
}

} //namespace CGAL

#endif // CGAL_CARTESIAN_POINT_CONSTRUCTIONS_3_H
