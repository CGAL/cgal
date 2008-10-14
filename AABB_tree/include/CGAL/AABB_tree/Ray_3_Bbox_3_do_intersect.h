// Copyright (c) 2008  INRIA Sophia-Antipolis (France), ETHZ (Suisse).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s)     : Camille Wormser, Pierre Alliez

#ifndef CGAL_RAY_3_BBOX_3_DO_INTERSECT_H
#define CGAL_RAY_3_BBOX_3_DO_INTERSECT_H

#include <CGAL/Ray_3.h>

// inspired from http://cag.csail.mit.edu/~amy/papers/box-jgt.pdf

CGAL_BEGIN_NAMESPACE

namespace CGALi {

  template <class K>
  bool do_intersect(const typename K::Ray_3& ray, 
    const CGAL::Bbox_3& bbox,
    const K& kernel)
  {
    typedef typename K::FT FT;
    typedef typename K::Point_3 Point;
    typedef typename K::Vector_3 Vector;

    Point parameters[2];
    parameters[0] = Point(bbox.xmin(), bbox.ymin(), bbox.zmin());
    parameters[1] = Point(bbox.xmax(), bbox.ymax(), bbox.zmax());


    const Point source = ray.source();

    if((parameters[0].x() <= source.x()) && (source.x() <= parameters[1].x()) && 
      (parameters[0].y() <= source.y()) && (source.y() <= parameters[1].y()) &&
      (parameters[0].z() <= source.z()) && (source.z() <= parameters[1].z()))
      return true;

    const Vector direction = ray.to_vector();
    // CAREFUL, when 1.0/0.0, this works only with doubles, not with filtered kernels
    const Vector inv_direction(1.0/direction.x(),
      1.0/direction.y(),
      1.0/direction.z()); 
    const int sign_x = inv_direction.x() < (FT)0.0;
    const int sign_y = inv_direction.y() < (FT)0.0;
    const int sign_z = inv_direction.z() < (FT)0.0;

    FT tmin = (parameters[sign_x].x() - source.x()) * inv_direction.x();
    FT tmax = (parameters[1-sign_x].x() - source.x()) * inv_direction.x();

    if(tmax < (FT)0.0)
      return false;
    if(tmin < (FT)0.0)
      tmin = (FT)0.0;

    const FT tymin = (parameters[sign_y].y() - source.y()) * inv_direction.y();
    const FT tymax = (parameters[1-sign_y].y() - source.y()) * inv_direction.y();

    if(tmin > tymax || tymin > tmax) 
      return false;

    if(tymin > tmin)
      tmin = tymin;

    if(tymax < tmax)
      tmax = tymax;

    FT tzmin = (parameters[sign_z].z() - source.z()) * inv_direction.z();
    FT tzmax = (parameters[1-sign_z].z() - source.z()) * inv_direction.z();

    if(tmin > tzmax || tzmin > tmax) 
      return false;

    return true;
  }

} // namespace CGALi

template <class K>
bool do_intersect(const CGAL::Ray_3<K>& ray, 
		  const CGAL::Bbox_3& bbox)
{
  return CGALi::do_intersect(ray, bbox, K());
}

template <class K>
bool do_intersect(const CGAL::Bbox_3& bbox, 
		  const CGAL::Ray_3<K>& ray)
{
  return CGALi::do_intersect(ray, bbox, K());
}

CGAL_END_NAMESPACE

#endif  // CGAL_RAY_3_BBOX_3_DO_INTERSECT_H


