// Copyright (c) 2000  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
// 
//
// Author(s)     : Herve Bronnimann, Sylvain Pion

#ifndef CGAL_CARTESIAN_PLANE_CONSTRUCTIONS_3_H
#define CGAL_CARTESIAN_PLANE_CONSTRUCTIONS_3_H

#include <CGAL/Cartesian/Point_3.h>
#include <CGAL/Cartesian/Direction_3.h>
#include <CGAL/constructions/kernel_ftC3.h>

namespace CGAL {

template <class R_>
class PlaneC3;

template <class R>
CGAL_KERNEL_LARGE_INLINE
PlaneC3<R>
plane_from_points(const typename R::Point_3 &p,
                  const typename R::Point_3 &q,
		  const typename R::Point_3 &r)
{
  typename R::FT a, b, c, d;
  plane_from_pointsC3(p.x(), p.y(), p.z(),
                      q.x(), q.y(), q.z(),
                      r.x(), r.y(), r.z(),
                      a, b, c, d);
  return PlaneC3<R>(a, b, c, d);
}

template <class R>
CGAL_KERNEL_LARGE_INLINE
PlaneC3<R>
plane_from_point_direction(const typename R::Point_3 &p,
                           const typename R::Direction_3 &d)
{
  typename R::FT A, B, C, D;
  plane_from_point_directionC3(p.x(), p.y(), p.z(), d.dx(), d.dy(), d.dz(),
                               A, B, C, D);
  return PlaneC3<R>(A, B, C, D);
}

} //namespace CGAL

#endif // CGAL_CARTESIAN_PLANE_CONSTRUCTIONS_3_H
