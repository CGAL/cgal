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


#ifndef CGAL_INTERNAL_INTERSECTIONS_3_TETRAHEDRON_3_PLANE_3_INTERSECTIONS_H
#define CGAL_INTERNAL_INTERSECTIONS_3_TETRAHEDRON_3_PLANE_3_INTERSECTIONS_H

#include <CGAL/kernel_basic.h>
#include <CGAL/intersections.h>
#include <CGAL/Intersections_3/internal/Triangle_3_Plane_3_do_intersect.h>
#include <CGAL/Intersections_3/internal/tetrahedron_intersection_helpers.h>
#include <set>
namespace CGAL {

namespace Intersections {

namespace internal {


//Tetrahedron_3 Segment_3
template <class K>
typename Intersection_traits<K, typename K::Tetrahedron_3, typename K::Plane_3>::result_type
intersection(
    const typename K::Tetrahedron_3 &tet,
    const typename K::Plane_3 &pl,
    const K&)
{
  typedef typename Intersection_traits<K,
      typename K::Tetrahedron_3,
      typename K::Plane_3>::result_type Result_type;

  typedef typename Intersection_traits<K,
      typename K::Triangle_3,
      typename K::Plane_3>::result_type Inter_type;

  typedef typename K::Segment_3 Segment_3;
  Inter_type intersections[4];
  std::size_t seg_id = -1,
      p_id = -1;
  std::vector<typename K::Point_3> points;
  std::vector<Segment_3> segments;
  for(std::size_t i = 0; i < 4; ++i)
  {
   const typename K::Triangle_3 triangle(tet.vertex((i+1)%4),
                           tet.vertex((i+2)%4),
                           tet.vertex((i+3)%4));
    intersections[i] = typename K::Intersect_3()(pl, triangle);
    if(intersections[i]){
      if(const typename K::Triangle_3* tr = boost::get<typename K::Triangle_3>(&*intersections[i]))
      {
        typename K::Triangle_3 res = *tr;
        return Result_type(std::forward<typename K::Triangle_3>(res));
      }
      else if( const Segment_3* s
               = boost::get<Segment_3>(&*intersections[i]))
      {
        segments.push_back(*s);
        seg_id = i;
      }
      else if( const typename K::Point_3* p
               = boost::get<typename K::Point_3>(&*intersections[i]))
      {
        points.push_back(*p);
        p_id = i;
      }
    }
  }
  CGAL_assertion(segments.size() != 1);

  switch(segments.size())
  {
  case 0:
  {
    if(p_id == -1)
      return Result_type();
    else
    {
      typename K::Point_3 p
                    = *boost::get<typename K::Point_3>(&*intersections[p_id]);

      return Result_type(std::forward<typename  K::Point_3>(p));
    }
  }
    break;
  case 2:
  {
    return Result_type(std::forward<typename  K::Segment_3>(segments.back()));
  }
    break;
  default:
  {
    std::set<typename K::Point_3> all_points;
    for (auto s : segments)
    {
      all_points.insert(s.source());
      all_points.insert(s.target());
    }
    if(all_points.size() == 3)
    {
      auto p_it = all_points.begin();
      ++p_it;
      typename K::Point_3 mid_point = *p_it;
      ++p_it;
      typename K::Triangle_3 result(*all_points.begin(), mid_point, *p_it );
      return Result_type(std::forward<typename K::Triangle_3>(result));
    }
    else //size = 4
    {
      std::list<Segment_3> segs;
      for(auto s : segments)
        segs.push_back(s);
      std::list<typename K::Point_3> tmp;
      fill_points_list(segs, tmp);
      std::vector<typename K::Point_3> res;
      for( auto p : tmp)
        res.push_back(p);
      return Result_type(std::forward<std::vector<typename K::Point_3> >(res));
    }
  }
    break;
  }
  return Result_type();
}

template <class K>
typename Intersection_traits<K, typename K::Tetrahedron_3, typename K::Plane_3>::result_type
intersection(
    const typename K::Plane_3 &pl,
    const typename K::Tetrahedron_3 &tet,
    const K& k)
{
  return intersection(tet, pl, k);
}

}}}
#endif // CGAL_INTERNAL_INTERSECTIONS_3_TETRAHEDRON_3_PLANE_3_INTERSECTIONS_H
