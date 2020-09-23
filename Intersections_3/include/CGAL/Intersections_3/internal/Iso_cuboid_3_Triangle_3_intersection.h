// Copyright (c) 2019 GeometryFactory(France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Maxime Gimeno
//

#ifndef CGAL_INTERSECTIONS_3_INTERNAL_ISO_CUBOID_3_TRIANGLE_3_INTERSECTION_H
#define CGAL_INTERSECTIONS_3_INTERNAL_ISO_CUBOID_3_TRIANGLE_3_INTERSECTION_H

#include <CGAL/kernel_basic.h>
#include <CGAL/intersections.h>

#include <iterator>
#include <list>
#include <vector>

namespace CGAL {

namespace Intersections {

namespace internal {

//only work for convex polygons, but in here that's always the case
template<class K>
void clip_poly_halfspace(
    std::vector<typename K::Point_3>& polygon,
    const typename K::Plane_3& pl,
    const K& k)
{
  if(polygon.empty())
    return;

  typedef typename K::Point_3 Point;
  typedef typename K::Plane_3 Plane;
  typedef typename K::Segment_3 Segment;

  typedef typename Intersection_traits<K,
      Plane,
      CGAL::Segment_3<K> >::result_type SP_type;

  // Keep in memory which points we are going to delete later (newer intersection points
  // by construction will not be deleted)
  std::list<std::pair<Point, bool> > p_list;
  for(const Point& p : polygon)
    p_list.emplace_back(p, pl.has_on_positive_side(p));

  //corefine with plane.
  auto it = p_list.begin();
  while(it != p_list.end())
  {
    const Point& p1 = (it++)->first;
    if(it == p_list.end())
      break;

    const Point& p2 = it->first;
    const Segment seg = k.construct_segment_3_object()(p1, p2);

    if(do_intersect(seg, pl))
    {
      SP_type inter = k.intersect_3_object()(seg, pl);
      if(inter)
      {
        Point* p_inter = boost::get<Point>(&*inter);
        if(p_inter
           && !(k.equal_3_object()(*p_inter, p1))
           && !(k.equal_3_object()(*p_inter, p2)))
        {
          // 'false' because we know the intersection is by construction not on the positive side of the plane
          p_list.insert(it, std::make_pair(*p_inter, false));
        }
      }
    }
  }

  if(polygon.size() > 2)
  {
    const Point& p2 = p_list.front().first;
    const Point& p1 = p_list.back().first;
    const Segment seg(p1, p2);

    if(do_intersect(seg, pl))
    {
      SP_type inter = typename K::Intersect_3()(seg, pl);
      if(inter)
      {
        Point* p_inter = boost::get<Point>(&*inter);
        if(p_inter
           && !(k.equal_3_object()(*p_inter, p1))
           && !(k.equal_3_object()(*p_inter, p2)))
        {
          // 'false' because we know the intersection is by construction not on the positive side of the plane
          p_list.emplace_back(*p_inter, false);
        }
      }
    }
  }

  //remove all points on positive side
  for(auto it = p_list.begin(); it != p_list.end();)
  {
    if(it->second)
      it = p_list.erase(it);
    else
      ++it;
  }

  // Update the polygon
  polygon.clear();
  for(const auto& pr : p_list)
    polygon.push_back(pr.first);
}

template <class K>
typename Intersection_traits<K, typename K::Iso_cuboid_3, typename K::Triangle_3>::result_type
intersection(
    const typename K::Iso_cuboid_3 &cub,
    const typename K::Triangle_3 &tr,
    const K& k)
{
  typedef typename K::Point_3 Point;
  typedef typename K::Segment_3 Segment;
  typedef typename K::Triangle_3 Triangle;
  typedef typename K::Plane_3 Plane;
  typedef std::vector<Point> Poly;

  typedef typename Intersection_traits<K,
      CGAL::Iso_cuboid_3<K>,
      CGAL::Triangle_3<K> >::result_type Res_type;

  //Lazy implem: clip 6 times the input triangle.
  Plane planes[6];
  planes[0] = Plane(cub.vertex(0),
                    cub.vertex(1),
                    cub.vertex(5));

  planes[1] = Plane(cub.vertex(0),
                    cub.vertex(4),
                    cub.vertex(3));

  planes[2] = Plane(cub.vertex(0),
                    cub.vertex(3),
                    cub.vertex(1));

  planes[3] = Plane(cub.vertex(7),
                    cub.vertex(6),
                    cub.vertex(1));

  planes[4] = Plane(cub.vertex(7),
                    cub.vertex(3),
                    cub.vertex(4));

  planes[5] = Plane(cub.vertex(7),
                    cub.vertex(4),
                    cub.vertex(6));

  std::vector<Point> poly;
  poly.push_back(tr.vertex(0));
  poly.push_back(tr.vertex(1));
  poly.push_back(tr.vertex(2));

  for (int i = 0; i < 6; ++i)
    clip_poly_halfspace<K>(poly, planes[i], k);

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

    Triangle res = Triangle(poly[0], poly[1], poly[2]);
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

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERSECTIONS_3_INTERNAL_ISO_CUBOID_3_TRIANGLE_3_INTERSECTION_H
