// Copyright (c) 2019 GeometryFactory(France).
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
// Author(s)     : Maxime Gimeno
//

#ifndef CGAL_INTERSECTIONS_3_INTERNAL_ISO_CUBOID_3_TRIANGLE_3_INTERSECTION_H
#define CGAL_INTERSECTIONS_3_INTERNAL_ISO_CUBOID_3_TRIANGLE_3_INTERSECTION_H

#include <CGAL/kernel_basic.h>
#include <CGAL/intersections.h>

#include <vector>
#include <list>

namespace CGAL {

namespace Intersections {

namespace internal {

//only work for convex polygons, but in here that's always the case
template<class K>
void clip_poly_halfspace(
    const std::vector<typename K::Point_3>& input,
    const typename K::Plane_3& pl,
    std::vector<typename K::Point_3>& output)
{
  if(input.empty())
    return;

  typedef typename K::Point_3 Point;
  typedef typename K::Plane_3 Plane;
  typedef typename K::Segment_3 S;

  typedef typename Intersection_traits<K,
      CGAL::Plane_3<K>,
      CGAL::Segment_3<K> >::result_type SP_type;

  std::list<Point> p_list(input.begin(), input.end());
  auto it = p_list.begin();
  //corefine with plane.
  while(it != p_list.end())
  {
    Point p = *it;
    ++it;
    if(it == p_list.end())
      break;
    if(do_intersect(S(p, *it), pl))
    {
      SP_type inter = typename K::Intersect_3()(S(p, *it), pl);
      if(inter)
      {
        Point* p_inter = boost::get<Point>(&*inter);
        if(p_inter
           && *p_inter != p
           && *p_inter != *it)
          p_list.insert(it, *p_inter);
      }
    }
  }

  if(input.size() >2)
  {
    Point p2(p_list.front()),
        p1(p_list.back());
    S seg(p1, p2);
    if(do_intersect(seg, pl))
    {
      SP_type inter = typename K::Intersect_3()(seg, pl);
      if(inter)
      {
        Point* p_inter = boost::get<Point>(&*inter);
        if(p_inter
           && *p_inter != p1
           && *p_inter != p2)
          p_list.push_back(*p_inter);
      }
    }
  }
  //remove all points on positive side

  for(auto p_it = p_list.begin();
      p_it != p_list.end();)
  {
    if(pl.has_on_positive_side(*p_it))
      p_it = p_list.erase(p_it);
    else
      ++p_it;
  }
  for(auto p : p_list)
    output.push_back(p);

}

template <class K>
typename Intersection_traits<K, typename K::Iso_cuboid_3, typename K::Triangle_3>::result_type
intersection(
    const typename K::Iso_cuboid_3 &cub,
    const typename K::Triangle_3 &tr,
    const K&)
{
  typedef typename K::Point_3 Point;
  typedef typename K::Segment_3 Segment;
  typedef typename K::Triangle_3 Triangle;
  typedef typename K::Plane_3 Plane_3;
  typedef std::vector<Point> Poly;

  typedef typename Intersection_traits<K,
      CGAL::Iso_cuboid_3<K>,
      CGAL::Triangle_3<K> >::result_type Res_type;

  //Lazy implem: clip 6 times the input triangle.
  Plane_3 planes[6];
  planes[0] = Plane_3(cub.vertex(0),
                      cub.vertex(1),
                      cub.vertex(5));

  planes[1] = Plane_3(cub.vertex(0),
                      cub.vertex(4),
                      cub.vertex(3));

  planes[2]=Plane_3(cub.vertex(0),
                    cub.vertex(3),
                    cub.vertex(1));

  planes[3] = Plane_3(cub.vertex(7),
                      cub.vertex(6),
                      cub.vertex(1));

  planes[4] = Plane_3(cub.vertex(7),
                      cub.vertex(3),
                      cub.vertex(4));

  planes[5] = Plane_3(cub.vertex(7),
                      cub.vertex(4),
                      cub.vertex(6));

  std::vector<Point> poly;
  poly.push_back(tr.vertex(0));
  poly.push_back(tr.vertex(1));
  poly.push_back(tr.vertex(2));

  for (int i = 0; i < 6; ++i)
  {
    Poly clipped;
    clip_poly_halfspace<K>(poly, planes[i], clipped);
    poly = clipped;
    for(auto p : poly)
      std::cout<<p<<std::endl;
  }

  switch(poly.size())
  {
  case 0:
    return Res_type();
    break;
  case 1:
  {
    Point res = poly.front();
    return Res_type(std::forward<Point>(res));
  }
    break;
  case 2:
  {
    Segment res = Segment(poly.front(), poly.back());
    return Res_type(std::forward<Segment>(res));
  }
    break;
  case 3:
  {

    Triangle res = Triangle (poly[0], poly[1], poly[2]);
    return Res_type(std::forward<Triangle>(res));
  }
    break;
  default:
  {
    return Res_type(std::forward<Poly>(poly));
  }
    break;
  }
}

template <class K>
typename Intersection_traits<K, typename K::Iso_cuboid_3, typename K::Triangle_3>::result_type
intersection(
    const typename K::Triangle_3 &tr,
    const typename K::Iso_cuboid_3 &cub,
    const K& k)
{
  return intersection(cub, tr, k);
}
}}}

#endif // CGAL_INTERSECTIONS_3_INTERNAL_ISO_CUBOID_3_TRIANGLE_3_INTERSECTION_H
