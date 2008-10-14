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

#ifndef CGAL_SPHERE_3_BBOX_DO_INTERSECT_H
#define CGAL_SPHERE_3_BBOX_DO_INTERSECT_H

#include <CGAL/Sphere_3.h>
#include <CGAL/Bbox_3.h>

#include <CGAL/number_utils.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

  // assumes that the intersection with the supporting plane has
  // already been checked.
  template <class K>
  bool do_intersect(const typename K::Sphere_3& sphere, 
    const CGAL::Bbox_3& bbox,
    const K& kernel)
  {
    typename K::FT d, distance = 0;
    for(int i = 0; i < 3; ++i)
    {
      if(sphere.center()[i] < bbox.min(i))
      {
	d = bbox.min(i) - sphere.center()[i];
	distance += d*d;
      }
      else if(sphere.center()[i] > bbox.max(i))
      {
	d = sphere.center()[i] - bbox.max(i);
	distance += d*d;
      }
    }
    return distance <= sphere.squared_radius();
  }

} // namespace CGALi

template <class K>
bool do_intersect(const CGAL::Sphere_3<K>& sphere, 
		  const CGAL::Bbox_3& bbox)
{
  return typename K::Do_intersect_3()(sphere, bbox);
}

template <class K>
bool do_intersect(const CGAL::Bbox_3& bbox, 
		  const CGAL::Sphere_3<K>& sphere)
{
  return typename K::Do_intersect_3()(sphere, bbox);
}

CGAL_END_NAMESPACE

#endif  // CGAL_SPHERE_3_BBOX_DO_INTERSECT_H
