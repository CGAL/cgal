// Copyright (c) 2018  GeometryFactory Sarl (France).
// All rights reserved.
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
// Author(s)     : Andreas Fabri

#ifndef CGAL_INTERNAL_INTERSECTIONS_3_TETRAHEDRON_3_DO_INTERSECT_H
#define CGAL_INTERNAL_INTERSECTIONS_3_TETRAHEDRON_3_DO_INTERSECT_H

#include <CGAL/Intersections_3/Ray_3_Triangle_3.h>
#include <CGAL/Intersections_3/Sphere_3_Triangle_3.h>
#include <CGAL/Intersections_3/Plane_3_Triangle_3.h>

namespace CGAL {

namespace Intersections {

namespace internal {

template<typename K, class Unbounded>
typename K::Boolean
do_intersect_tetrahedron_unbounded(const typename K::Tetrahedron_3& tet,
                                   const Unbounded& unb,
                                   const K& k) {
  typedef typename K::Triangle_3 Triangle;
  typedef typename K::Boolean Boolean;
  Boolean result = false;
  for (int i = 0; i < 4; ++i)
  {
    const Boolean b = do_intersect(unb,
                                   Triangle(tet[i],
                                            tet[(i+1)%4],
                                            tet[(i+2)%4]),
                                   k);
    if(certainly(b)) return b;
    if(is_indeterminate(b)) result = b;
  }
  return result;
}



template<typename K>
typename K::Boolean
do_intersect(const typename K::Plane_3& unb,
             const typename K::Tetrahedron_3& tet,
             const K& k) {
  return do_intersect_tetrahedron_unbounded(tet, unb, k);
}

template<typename K>
typename K::Boolean
do_intersect(const typename K::Tetrahedron_3& tet,
             const typename K::Plane_3& unb,
             const K& k) {
  return do_intersect_tetrahedron_unbounded(tet, unb, k);
}


template<typename K>
typename K::Boolean
do_intersect(const typename K::Line_3& unb,
             const typename K::Tetrahedron_3& tet,
             const K& k) {
  return do_intersect_tetrahedron_unbounded(tet, unb, k);
}

template<typename K>
typename K::Boolean
do_intersect(const typename K::Tetrahedron_3& tet,
             const typename K::Line_3& unb,
             const K& k) {
  return do_intersect_tetrahedron_unbounded(tet, unb, k);
}


template<typename K>
typename K::Boolean
do_intersect(const typename K::Ray_3& unb,
             const typename K::Tetrahedron_3& tet,
             const K& k) {
  return do_intersect_tetrahedron_unbounded(tet, unb, k);
}

template<typename K>
typename K::Boolean
do_intersect(const typename K::Tetrahedron_3& tet,
             const typename K::Ray_3& unb,
             const K& k) {
  return do_intersect_tetrahedron_unbounded(tet, unb, k);
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_3_TETRAHEDRON_3_DO_INTERSECT_H
