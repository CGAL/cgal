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
// Author(s)     : Camille Wormser, Jane Tournois, Pierre Alliez

#ifndef CGAL_TRIANGLE_3_BBOX_3_DO_INTERSECT_H
#define CGAL_TRIANGLE_3_BBOX_3_DO_INTERSECT_H

#include <CGAL/Triangle_3.h>
#include <CGAL/Bbox_3.h>

// Fast Triangle-Cuboid intersection test, following Tomas Akenine-Moeller description.
// The code looks slightly different from his code because we avoid the translation at 
// a minimal cost (and we use C++ ;).

#include <CGAL/number_utils.h>
#include "Plane_3_Bbox_3_do_intersect.h"

CGAL_BEGIN_NAMESPACE

namespace CGALi {

  template <class K>
  inline
    bool do_bbox_intersect(const typename K::Triangle_3& triangle, 
    const CGAL::Bbox_3& bbox)
  {
    const typename K::Point_3& p = triangle.vertex(0), 
      q = triangle.vertex(1),
      r = triangle.vertex(2);
    for(int i = 0; i < 3; ++i) {
      if(p[i] <= q[i]) {
	if(q[i] <= r[i]) { // pqr
	  if(bbox.max(i) < p[i] || bbox.min(i) > r[i])
	    return false;
	}
	else {
	  if(p[i] <= r[i]) { // prq
	    if(bbox.max(i) < p[i] || bbox.min(i) > q[i])
	      return false;
	  }
	  else { // rpq
	    if(bbox.max(i) < r[i] || bbox.min(i) > q[i])
	      return false;
	  }
	}
      }
      else {
	if(p[i] <= r[i]) { // qpr
	  if(bbox.max(i) < q[i] || bbox.min(i) > r[i])
	    return false;
	}
	else {
	  if(q[i] <= r[i]) { // qrp
	    if(bbox.max(i) < q[i] || bbox.min(i) > p[i])
	      return false;
	  }
	  else { // rqp
	    if(bbox.max(i) < r[i] || bbox.min(i) > p[i])
	      return false;
	  }
	}
      }
    }
    return true;
  }

  template <class K, int AXE, int SIDE>
  inline
    bool do_axis_intersect(const typename K::Triangle_3& triangle,
    const typename K::Vector_3* sides, 
    const CGAL::Bbox_3& bbox)
  {
    const typename K::Point_3& j = triangle.vertex(SIDE);
    const typename K::Point_3& k = triangle.vertex((SIDE+2)%3);

    typename K::FT t_min = AXE==0? -sides[SIDE].z()*j.y() + sides[SIDE].y()*j.z():
      AXE==1?  sides[SIDE].z()*j.x() - sides[SIDE].x()*j.z(): 
      -sides[SIDE].y()*j.x() + sides[SIDE].x()*j.y();

    typename K::FT t_max = AXE==0? -sides[SIDE].z()*k.y() + sides[SIDE].y()*k.z():
      AXE==1?  sides[SIDE].z()*k.x() - sides[SIDE].x()*k.z(): 
      -sides[SIDE].y()*k.x() + sides[SIDE].x()*k.y();

    typename K::FT tmp;
    if(t_max < t_min)
    {
      tmp = t_min;
      t_min = t_max;
      t_max = tmp;
    }
    typename K::Point_3 p_min, p_max;
    get_min_max<K, AXE>(AXE==0? 0: AXE==1? sides[SIDE].z(): -sides[SIDE].y(),
      AXE==0? -sides[SIDE].z(): AXE==1? 0: sides[SIDE].x(),
      AXE==0? sides[SIDE].y(): AXE==1? -sides[SIDE].x(): 0,
      bbox, p_min, p_max);

    if((AXE==0? -sides[SIDE].z()*p_min.y() + sides[SIDE].y()*p_min.z() :
      AXE==1?  sides[SIDE].z()*p_min.x() - sides[SIDE].x()*p_min.z() :
      -sides[SIDE].y()*p_min.x() + sides[SIDE].x()*p_min.y() ) > t_max ||
      (AXE==0? -sides[SIDE].z()*p_max.y() + sides[SIDE].y()*p_max.z() :
      AXE==1?  sides[SIDE].z()*p_max.x() - sides[SIDE].x()*p_max.z() :
      -sides[SIDE].y()*p_max.x() + sides[SIDE].x()*p_max.y() ) < t_min )
      return false;
    return true;
  }

  // AXE is the axe such that p is orthogonal to it.
  // if you do not know it, or if it does not exist, 
  // use get_min_max without the AXE template parameter
  // available in _plane_is_cuboid_do_intersect.h
  template <class K, int AXE>
  inline
    void get_min_max(const typename K::FT& px,
    const typename K::FT& py,
    const typename K::FT& pz,
    const CGAL::Bbox_3& c,
    typename K::Point_3& p_min, 
    typename K::Point_3& p_max)
  {
    if(AXE == 0 || px > 0) {
      if(AXE == 1 || py > 0) {
	if(AXE == 2 || pz > 0) { p_min = typename K::Point_3(c.xmin(), c.ymin(),c.zmin()); 
	p_max = typename K::Point_3(c.xmax(), c.ymax(),c.zmax());}
	else {							     p_min = typename K::Point_3(c.xmin(), c.ymin(),c.zmax()); 
	p_max = typename K::Point_3(c.xmax(), c.ymax(),c.zmin());}
      }
      else {
	if(AXE == 2 || pz > 0) { p_min = typename K::Point_3(c.xmin(), c.ymax(),c.zmin()); 
	p_max = typename K::Point_3(c.xmax(), c.ymin(),c.zmax());}
	else {					         p_min = typename K::Point_3(c.xmin(), c.ymax(),c.zmax()); 
	p_max = typename K::Point_3(c.xmax(), c.ymin(),c.zmin());}
      }
    }
    else {
      if(AXE == 1 || py > 0) {
	if(AXE == 2 || pz > 0) { p_min = typename K::Point_3(c.xmax(), c.ymin(),c.zmin()); 
	p_max = typename K::Point_3(c.xmin(), c.ymax(),c.zmax());}
	else {					         p_min = typename K::Point_3(c.xmax(), c.ymin(),c.zmax()); 
	p_max = typename K::Point_3(c.xmin(), c.ymax(),c.zmin());}
      }
      else {
	if(AXE == 2 || pz > 0) { p_min = typename K::Point_3(c.xmax(), c.ymax(),c.zmin()); 
	p_max = typename K::Point_3(c.xmin(), c.ymin(),c.zmax());}
	else {					         p_min = typename K::Point_3(c.xmax(), c.ymax(),c.zmax()); 
	p_max = typename K::Point_3(c.xmin(), c.ymin(),c.zmin());}
      }
    }
  }

  // assumes that the intersection with the supporting plane has
  // already been checked.
  template <class K>
  bool do_intersect(const typename K::Triangle_3& triangle, 
    const CGAL::Bbox_3& bbox,
    const K& kernel)
  {
    if(! do_bbox_intersect<K>(triangle, bbox))
      return false;

    typename K::Vector_3 sides[3];
    sides[0] = triangle[1] - triangle[0];
    sides[1] = triangle[2] - triangle[1];
    sides[2] = triangle[0] - triangle[2];

    if(! do_axis_intersect<K,0,0>(triangle, sides, bbox))
      return false;
    if(! do_axis_intersect<K,0,1>(triangle, sides, bbox))
      return false;
    if(! do_axis_intersect<K,0,2>(triangle, sides, bbox))
      return false;
    if(! do_axis_intersect<K,1,0>(triangle, sides, bbox))
      return false;
    if(! do_axis_intersect<K,1,1>(triangle, sides, bbox))
      return false;
    if(! do_axis_intersect<K,1,2>(triangle, sides, bbox))
      return false;
    if(! do_axis_intersect<K,2,0>(triangle, sides, bbox))
      return false;
    if(! do_axis_intersect<K,2,1>(triangle, sides, bbox))
      return false;
    if(! do_axis_intersect<K,2,2>(triangle, sides, bbox))
      return false;
    return true; 
  }

} // namespace CGALi

template <class K>
bool do_intersect(const CGAL::Triangle_3<K>& triangle, 
		  const CGAL::Bbox_3& bbox)
{
  return typename K::Do_intersect_3()(triangle, bbox);
}

template <class K>
bool do_intersect(const CGAL::Bbox_3& bbox, 
		  const CGAL::Triangle_3<K>& triangle)
{
  return typename K::Do_intersect_3()(triangle, bbox);
}

CGAL_END_NAMESPACE

#endif  // CGAL_TRIANGLE_3_BBOX_3_DO_INTERSECT_H
