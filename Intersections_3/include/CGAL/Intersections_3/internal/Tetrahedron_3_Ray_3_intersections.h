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

#ifndef CGAL_INTERNAL_INTERSECTIONS_3_TETRAHEDRON_3_RAY_3_INTERSECTION_H
#define CGAL_INTERNAL_INTERSECTIONS_3_TETRAHEDRON_3_RAY_3_INTERSECTION_H

#include <CGAL/kernel_basic.h>
#include <CGAL/intersections.h>

namespace CGAL {

namespace Intersections {

namespace internal {
//Tetrahedron_3 Ray_3
template <class K>
typename Intersection_traits<K, typename K::Tetrahedron_3, typename K::Ray_3>::result_type
intersection(
    const typename K::Tetrahedron_3 &tet,
    const typename K::Ray_3 &ray,
    const K&)
{
  typedef typename Intersection_traits<K,
      typename K::Tetrahedron_3,
      typename K::Ray_3>::result_type Result_type;

  typedef typename Intersection_traits<K,
      typename K::Triangle_3,
      typename K::Ray_3>::result_type Inter_type;

  int res_id = -1;

  Inter_type tr_seg[4];
  for(std::size_t i = 0; i < 4; ++i)
  {
   const typename K::Triangle_3 triangle(tet.vertex((i+1)%4),
                           tet.vertex((i+2)%4),
                           tet.vertex((i+3)%4));
    tr_seg[i] = CGAL::intersection(ray, triangle);
    if(tr_seg[i])
      if( boost::get<typename K::Ray_3>(&*tr_seg[i]) != nullptr)
        res_id = i;
  }
  //if there is a segment in the intersections, then we return it
  if(res_id !=-1)
    return tr_seg[res_id];

  //else if there is only 1 intersection
  std::vector<typename K::Point_3> res_points;
  res_points.reserve(4);
  res_id = -1;
  for(std::size_t i = 0; i< 4; ++i)
  {
    if(tr_seg[i])
    {
      if (const typename K::Point_3* p = boost::get<typename K::Point_3>(&*tr_seg[i]))
      {
        if(res_points.empty())
        {
          res_id = i;
        }
        else {
          if(*p != res_points.front())
          {
            res_id = -1;
          }
        }
        res_points.push_back(*p);
      }
    }
  }
  if(res_id != -1)
  {
    //If source is inside tet : return a segment
    if(tet.has_on_bounded_side(ray.source()))
    {
      typename K::Ray_3 result(ray.source(), res_points.front());
      return Result_type(std::forward<typename K::Ray_3>(result));
    }
    //else segment of intersections:
    return tr_seg[res_id];
  }
  //else, we return a segment of the 2 intersection points (the most far away, in case of inexact)
  typename K::FT max_dist = 0;
  std::size_t res_id_2 = -1;
  std::vector<std::vector<typename K::FT> > sq_distances(res_points.size());
  for(std::size_t i = 0; i< res_points.size(); ++i)
  {
    auto p1 = res_points[i];
    for(auto p2 : res_points)
    {
      sq_distances[i].push_back(CGAL::squared_distance(p1, p2));
      if(sq_distances[i].back() > max_dist)
      {
        res_id = i;
        res_id_2 = sq_distances[i].size()-1;
        max_dist = sq_distances[i].back();
      }
    }
  }
  CGAL_assertion(res_id != -1);
  CGAL_assertion(res_id_2 != -1);
  CGAL_assertion(max_dist >0 );

  typename K::Ray_3 res_seg(res_points[res_id], res_points[res_id_2]);

  return Result_type(std::forward<typename K::Ray_3>(res_seg));
}

template <class K>
typename Intersection_traits<K, typename K::Tetrahedron_3, typename K::Ray_3>::result_type
intersection(
    const typename K::Ray_3 &seg,
    const typename K::Tetrahedron_3 &tet,
    const K& k)
{
  return intersection(tet, seg, k);
}

}}}
#endif // CGAL_INTERNAL_INTERSECTIONS_3_TETRAHEDRON_3_RAY_3_INTERSECTION_H
