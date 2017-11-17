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
