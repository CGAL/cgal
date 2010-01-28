// Copyright (c) 1997  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Geert-Jan Giezeman <geert@cs.uu.nl>


#ifndef CGAL_BBOX_INTERSECTION_3_H
#define CGAL_BBOX_INTERSECTION_3_H

#include <CGAL/Bbox_3.h>
#include <CGAL/Object.h>

CGAL_BEGIN_NAMESPACE

// This function intersects a bbox with a ray, line or segment 
// Its essentially a copy of the function that was in Bbox_3_intersections.cpp
// But it must be a template function since the original kernel must be 
// taken into account. (Michael.Hemmer@sophia.inria.fr)
template <class R_cd> 
Object
intersection_bl(const Bbox_3 &box,
        double lpx, double lpy, double lpz,
        double ldx, double ldy, double ldz,
        bool min_infinite, bool max_infinite)
{
  Point_3<R_cd> ref_point(lpx,lpy, lpz);
  Vector_3<R_cd> dir(ldx, ldy, ldz);
  double seg_min = 0.0, seg_max = 1.0;
// first on x value
  if (dir.x() == 0.0) { 
    if (ref_point.x() < box.xmin())
      return Object();
    if (ref_point.x() > box.xmax())
      return Object();
  } else {
    double newmin, newmax;
    if (dir.x() > 0.0) {
      newmin = (box.xmin()-ref_point.x())/dir.x();
      newmax = (box.xmax()-ref_point.x())/dir.x();
    } else {
      newmin = (box.xmax()-ref_point.x())/dir.x();
      newmax = (box.xmin()-ref_point.x())/dir.x();
    }
    if (min_infinite) {
      min_infinite = false;
      seg_min = newmin;
    } else {
      if (newmin > seg_min)
        seg_min = newmin;
    }
    if (max_infinite) {
      max_infinite = false;
      seg_max = newmax;
    } else {
      if (newmax < seg_max)
        seg_max = newmax;
    }
    if (seg_max < seg_min)
      return Object();
  }
// now on y value
  if (dir.y() == 0.0) {
    if (ref_point.y() < box.ymin())
      return Object();
    if (ref_point.y() > box.ymax())
      return Object();
  } else {
    double newmin, newmax;
    if (dir.y() > 0.0) {
      newmin = (box.ymin()-ref_point.y())/dir.y();
      newmax = (box.ymax()-ref_point.y())/dir.y();
    } else {
      newmin = (box.ymax()-ref_point.y())/dir.y();
      newmax = (box.ymin()-ref_point.y())/dir.y();
    }
    if (min_infinite) {
      min_infinite = false;
      seg_min = newmin;
    } else {
      if (newmin > seg_min)
        seg_min = newmin;
    }
    if (max_infinite) {
      max_infinite = false;
      seg_max = newmax;
    } else {
      if (newmax < seg_max)
        seg_max = newmax;
    }
    if (seg_max < seg_min)
      return Object();
  }
// now on z value
  if (dir.z() == 0.0) {
    if (ref_point.z() < box.zmin())
      return Object();
    if (ref_point.z() > box.zmax())
      return Object();
  } else {
    double newmin, newmax;
    if (dir.z() > 0.0) {
      newmin = (box.zmin()-ref_point.z())/dir.z();
      newmax = (box.zmax()-ref_point.z())/dir.z();
    } else {
      newmin = (box.zmax()-ref_point.z())/dir.z();
      newmax = (box.zmin()-ref_point.z())/dir.z();
    }
    if (min_infinite) {
      min_infinite = false;
      seg_min = newmin;
    } else {
      if (newmin > seg_min)
        seg_min = newmin;
    }
    if (max_infinite) {
      max_infinite = false;
      seg_max = newmax;
    } else {
      if (newmax < seg_max)
        seg_max = newmax;
    }
    if (seg_max < seg_min)
      return Object();
  }
  if (min_infinite || max_infinite) {
    seg_max = 0.0;
    CGAL_kernel_assertion_msg(true,
        "Zero direction vector of line detected.");
  }
  if (seg_max == seg_min)
    return make_object(ref_point + dir*seg_max);
  return make_object(Segment_3<R_cd>(
                         ref_point + dir*seg_min, ref_point + dir*seg_max));
}


CGAL_END_NAMESPACE

#endif // CGAL_BBOX_INTERSECTION_3_H
