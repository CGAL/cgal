// Copyright (c) 1997  
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
// 
//
// Author(s)     : Geert-Jan Giezeman <geert@cs.uu.nl>


#ifndef CGAL_BBOX_INTERSECTION_3_H
#define CGAL_BBOX_INTERSECTION_3_H

#include <CGAL/Bbox_3.h>
#include <CGAL/Object.h>

namespace CGAL {

// This function intersects a bbox with a ray, line or segment 
// Its essentially a copy of the function that was in Bbox_3_intersections.cpp
// But it must be a template function since the original kernel must be 
// taken into account. (Michael.Hemmer@sophia.inria.fr)
template <class K> 
#if CGAL_INTERSECTION_VERSION < 2
CGAL::Object
#else
typename boost::optional< boost::variant< 
                            typename K::Segment_3, 
                            typename K::Point_3 > >
#endif
intersection_bl(const Bbox_3 &box,
        double lpx, double lpy, double lpz,
        double ldx, double ldy, double ldz,
        bool min_infinite, bool max_infinite)
{
  typedef 
#if CGAL_INTERSECTION_VERSION < 2
    CGAL::Object
#else
    typename 
    boost::optional< 
      boost::variant< typename K::Segment_3, 
                      typename K::Point_3 > >
#endif
    result_type;

  double seg_min = 0.0, seg_max = 1.0;
  // first on x value
  if (ldx == 0.0) { 
    if (lpx < box.xmin())
      return result_type();
    if (lpx > box.xmax())
      return result_type();
  } else {
    double newmin, newmax;
    if (ldx > 0.0) {
      newmin = (box.xmin()-lpx)/ldx;
      newmax = (box.xmax()-lpx)/ldx;
    } else {
      newmin = (box.xmax()-lpx)/ldx;
      newmax = (box.xmin()-lpx)/ldx;
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
      return result_type();
  }
// now on y value
  if (ldy == 0.0) {
    if (lpy < box.ymin())
      return result_type();
    if (lpy > box.ymax())
      return result_type();
  } else {
    double newmin, newmax;
    if (ldy > 0.0) {
      newmin = (box.ymin()-lpy)/ldy;
      newmax = (box.ymax()-lpy)/ldy;
    } else {
      newmin = (box.ymax()-lpy)/ldy;
      newmax = (box.ymin()-lpy)/ldy;
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
      return result_type();
  }
// now on z value
  if (ldz == 0.0) {
    if (lpz < box.zmin())
      return result_type();
    if (lpz > box.zmax())
      return result_type();
  } else {
    double newmin, newmax;
    if (ldz > 0.0) {
      newmin = (box.zmin()-lpz)/ldz;
      newmax = (box.zmax()-lpz)/ldz;
    } else {
      newmin = (box.zmax()-lpz)/ldz;
      newmax = (box.zmin()-lpz)/ldz;
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
      return result_type();
  }
  if (min_infinite || max_infinite) {
    seg_max = 0.0;
    CGAL_kernel_assertion_msg(true,
        "Zero direction vector of line detected.");
  }

  typedef typename K::FT FT;
  typedef typename K::Point_3 Point_3;   
  typedef typename K::Vector_3  Vector_3; 
  typedef typename K::Segment_3 Segment_3; 
  
  Point_3 ref_point = Point_3( FT(lpx), FT(lpy), FT(lpz));
  Vector_3 dir = Vector_3( FT(ldx), FT(ldy), FT(ldz));
  
  if (seg_max == seg_min) {
#if CGAL_INTERSECTION_VERSION < 2
    return make_object(ref_point + dir * FT(seg_max));
#else
    return result_type(ref_point + dir * FT(seg_max));
#endif
  } 
#if CGAL_INTERSECTION_VERSION < 2
  return make_object(
      Segment_3(ref_point + dir*FT(seg_min), ref_point + dir*FT(seg_max)));
#else
  return result_type(
      Segment_3(ref_point + dir*FT(seg_min), ref_point + dir*FT(seg_max)));
#endif
}


} //namespace CGAL

#endif // CGAL_BBOX_INTERSECTION_3_H
