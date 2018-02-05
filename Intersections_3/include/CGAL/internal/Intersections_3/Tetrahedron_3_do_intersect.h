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

#include <CGAL/Tetrahedron_3.h>
#include <CGAL/Plane_3.h>
#include <CGAL/Line_3.h>
#include <CGAL/Ray_3.h>

namespace CGAL {

namespace internal {

template<typename K, class Unbounded>
bool do_intersect_tetrahedron_unbounded(const typename K::Tetrahedron_3& tet,
                                        const Unbounded& unb,
                                        const K& k) {
    if (do_intersect(unb, Triangle(tet[0], tet[1], tet[2]), k)) return true;
    if (do_intersect(unb, Triangle(tet[0], tet[1], tet[3]), k)) return true;
    if (do_intersect(unb, Triangle(tet[0], tet[2], tet[3]), k)) return true;
    if (do_intersect(unb, Triangle(tet[1], tet[2], tet[3]), k)) return true;
    return false;
}



template<typename K>
bool do_intersect(const typename K::Tetrahedron_3& tet,
                  const typename K::Plane_3& unb,
                  const K& k) {
  return do_intersect_tetrahedron_unbounded(tet, unb, k);
}


template<typename K>
bool do_intersect(const typename K::Tetrahedron_3& tet,
                  const typename K::Line_3& unb,
                  const K& k) {
  return do_intersect_tetrahedron_unbounded(tet, unb, k);
}


template<typename K>
bool do_intersect(const typename K::Tetrahedron_3& tet,
                  const typename K::Ray_3& unb,
                  const K& k) {
  return do_intersect_tetrahedron_unbounded(tet, unb, k);
}
  
} // namespace internal

    
template<typename K>
bool do_intersect(const CGAL::Tetrahedron_3<K>& a,
                  const CGAL::Plane_3<K>& b) {
  return K().do_intersect_3_object()(a, b);
}

  
template<typename K>
bool do_intersect(const CGAL::Plane_3<K>& b,
                  const CGAL::Tetrahedron_3<K>& a) {
  return K().do_intersect_3_object()(a, b);
}

  
template<typename K>
bool do_intersect(const CGAL::Tetrahedron_3<K>& a,
                  const CGAL::Line_3<K>& b) {
  return K().do_intersect_3_object()(a, b);
}

  
template<typename K>
bool do_intersect(const CGAL::Line_3<K>& b,
                  const CGAL::Tetrahedron_3<K>& a) {
  return K().do_intersect_3_object()(a, b);
}

  
template<typename K>
bool do_intersect(const CGAL::Tetrahedron_3<K>& a,
                  const CGAL::Ray_3<K>& b) {
  return K().do_intersect_3_object()(a, b);
}

  
template<typename K>
bool do_intersect(const CGAL::Ray_3<K>& b,
                  const CGAL::Tetrahedron_3<K>& a) {
  return K().do_intersect_3_object()(a, b);
}

} // namespace CGAL

#endif CGAL_INTERNAL_INTERSECTIONS_3_TETRAHEDRON_3_DO_INTERSECT_H
