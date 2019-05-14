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


#ifndef CGAL_INTERNAL_INTERSECTIONS_3_BBOX_3_PLANE_3_DO_INTERSECT_H
#define CGAL_INTERNAL_INTERSECTIONS_3_BBOX_3_PLANE_3_DO_INTERSECT_H

#include <CGAL/Plane_3.h>
#include <CGAL/Bbox_3.h>

// Opcode like

namespace CGAL {

namespace internal {

  template <class K>
  Uncertain<bool> get_min_max(const typename K::Vector_3& p,
    const CGAL::Bbox_3& bbox,
    typename K::Point_3& p_min,
    typename K::Point_3& p_max)
  {
    if(certainly(p.x() > 0)) {
      if(certainly(p.y() > 0)) {
	if(certainly(p.z() > 0)) { 
            p_min = typename K::Point_3(bbox.xmin(), bbox.ymin(),bbox.zmin());
            p_max = typename K::Point_3(bbox.xmax(), bbox.ymax(),bbox.zmax());
            return true;
        } else if(certainly(p.z() <= 0)) {							    
          p_min = typename K::Point_3(bbox.xmin(), bbox.ymin(),bbox.zmax());
          p_max = typename K::Point_3(bbox.xmax(), bbox.ymax(),bbox.zmin());
          return true;
        }
      } else if(certainly(p.y() <= 0)) {
        if(certainly(p.z() > 0)) {
          p_min = typename K::Point_3(bbox.xmin(), bbox.ymax(),bbox.zmin());
          p_max = typename K::Point_3(bbox.xmax(), bbox.ymin(),bbox.zmax());
          return true;
        } else if(certainly(p.z() <= 0)) {
          p_min = typename K::Point_3(bbox.xmin(), bbox.ymax(),bbox.zmax());
          p_max = typename K::Point_3(bbox.xmax(), bbox.ymin(),bbox.zmin());
          return true;
        }
      }
    }
    else if(certainly(p.x() <= 0)) {
      if(certainly(p.y() > 0)) {
	if(certainly(p.z() > 0)) { 
          p_min = typename K::Point_3(bbox.xmax(), bbox.ymin(),bbox.zmin());
          p_max = typename K::Point_3(bbox.xmin(), bbox.ymax(),bbox.zmax());
          return true;}
	else if(certainly(p.z() <= 0)) {
          p_min = typename K::Point_3(bbox.xmax(), bbox.ymin(),bbox.zmax());
          p_max = typename K::Point_3(bbox.xmin(), bbox.ymax(),bbox.zmin());
          return true;
        }
      }
      else if(certainly(p.y() <= 0)) {
	if(certainly(p.z() > 0)) { 
          p_min = typename K::Point_3(bbox.xmax(), bbox.ymax(),bbox.zmin());
          p_max = typename K::Point_3(bbox.xmin(), bbox.ymin(),bbox.zmax());
          return true;
        }
	else if(certainly(p.z() <= 0)) {
          p_min = typename K::Point_3(bbox.xmax(), bbox.ymax(),bbox.zmax());
          p_max = typename K::Point_3(bbox.xmin(), bbox.ymin(),bbox.zmin());
          return true;
        }
      }
    }
    return  Uncertain<bool>::indeterminate();
  }

  template <class K>
  bool do_intersect(const typename K::Plane_3& plane,
    const CGAL::Bbox_3& bbox,
    const K&)
  {
    typedef typename K::Point_3 Point_3;
    Point_3 p_max, p_min;
    Uncertain<bool> b = get_min_max<K>(plane.orthogonal_vector(), bbox, p_min, p_max);
    if(is_certain(b)){
      return ! (plane.oriented_side(p_max) == ON_NEGATIVE_SIDE ||
                plane.oriented_side(p_min) == ON_POSITIVE_SIDE);
    }
    CGAL::Oriented_side side = plane.oriented_side(Point_3(bbox.xmin(), bbox.ymin(),bbox.zmin()));
    if(side == ON_ORIENTED_BOUNDARY) return true;
    CGAL::Oriented_side s = plane.oriented_side(Point_3(bbox.xmax(), bbox.ymax(),bbox.zmax()));
    if(s != side) return true;
    s = plane.oriented_side(Point_3(bbox.xmin(), bbox.ymin(),bbox.zmax()));
    if(s != side) return true;
    s = plane.oriented_side(Point_3(bbox.xmax(), bbox.ymax(),bbox.zmin()));
    if(s != side) return true;
    s = plane.oriented_side(Point_3(bbox.xmin(), bbox.ymax(),bbox.zmin()));
    if(s != side) return true;
     s = plane.oriented_side(Point_3(bbox.xmax(), bbox.ymin(),bbox.zmax()));
    if(s != side) return true;
     s = plane.oriented_side(Point_3(bbox.xmin(), bbox.ymax(),bbox.zmax()));
    if(s != side) return true;
     s = plane.oriented_side(Point_3(bbox.xmax(), bbox.ymin(),bbox.zmin()));
    if(s != side) return true;
    return false;
  }
  
  template <class K>
  bool do_intersect(const CGAL::Bbox_3& bbox, const typename K::Plane_3& plane,
                    const K&) {
    return do_intersect(plane, bbox, K());
  }


} // namespace internal

template<typename K>
bool do_intersect(const CGAL::Bbox_3 a,
                  const Plane_3<K>& b) {
  return K().do_intersect_3_object()(a, b);
}

template<typename K>
bool do_intersect(const Plane_3<K>& a,
                  const CGAL::Bbox_3& b) {
  return K().do_intersect_3_object()(a, b);
}


} //namespace CGAL

#endif  // CGAL_INTERNAL_INTERSECTIONS_3_BBOX_3_PLANE_3_DO_INTERSECT_H
