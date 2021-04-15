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
// Author(s)     : Andreas Fabri, Herve Bronnimann

#ifndef CGAL_CARTESIAN_PREDICATES_ON_PLANES_3_H
#define CGAL_CARTESIAN_PREDICATES_ON_PLANES_3_H

#include <CGAL/predicates/kernel_ftC3.h>

namespace CGAL {

template < class K >
inline
typename K::Oriented_side
side_of_oriented_plane(const PlaneC3<K> &h,
                       const PointC3<K> &p)
{
  return side_of_oriented_planeC3(h.a(), h.b(), h.c(), h.d(),
                                  p.x(), p.y(), p.z());
}

template < class K >
inline
typename K::Boolean
equal_plane(const PlaneC3<K> &h, const PlaneC3<K> &p)
{
  return equal_planeC3(h.a(), h.b(), h.c(), h.d(),
                       p.a(), p.b(), p.c(), p.d());
}

} //namespace CGAL

#endif // CGAL_CARTESIAN_PREDICATES_ON_PLANES_3_H
