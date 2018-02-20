// Copyright (c) 2008  INRIA Sophia-Antipolis (France), ETH Zurich (Switzerland).
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
// Author(s)     : Camille Wormser, Jane Tournois, Pierre Alliez


#ifndef CGAL_INTERNAL_INTERSECTIONS_3_ISO_CUBOID_3_SPHERE_3_DO_INTERSECT_H
#define CGAL_INTERNAL_INTERSECTIONS_3_ISO_CUBOID_3_SPHERE_3_DO_INTERSECT_H

#include <CGAL/Sphere_3.h>
#include <CGAL/Iso_cuboid_3.h>
#include <CGALL/Intersections_3/internal/Bbox_3_Sphere_3_do_intersect.h>

namespace CGAL {
  
namespace Intersections {

namespace internal {

    template <class K>
    bool do_intersect(const typename K::Iso_cuboid_3& ic,
                      const typename K::Sphere_3& sphere,
                      const K&)
    {
      return do_intersect_sphere_box_3(sphere, ic, K());
    }

    template <class K>
    bool do_intersect(const typename K::Sphere_3& sphere,
                      const typename K::Iso_cuboid_3& ic,
                      const K&)
    {
      return do_intersect_sphere_box_3(sphere, ic, K());
    }
  
} // namespace internal
} // namespace Intersections
  
template<typename K>
bool do_intersect(const Iso_cuboid_3<K>& a,
                  const Sphere_3<K>& b) {
  return K().do_intersect_3_object()(a, b);
}

template<typename K>
bool do_intersect(const Sphere_3<K>& a,
                  const Iso_cuboid_3<K>& b) {
  return K().do_intersect_3_object()(a, b);
}


} //namespace CGAL

#endif  // CGAL_INTERNAL_INTERSECTIONS_3_ISO_CUBOID_3_SPHERE_3_DO_INTERSECT_H
