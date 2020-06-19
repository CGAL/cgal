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


//Tetrahedron_3 Plane_3
template <class K>
typename Intersection_traits<K, typename K::Tetrahedron_3, typename K::Plane_3>::result_type
intersection(
    const typename K::Tetrahedron_3 &tet,
    const typename K::Plane_3 &pl,
    const K&)
{
  typedef typename Intersection_traits<K,
      typename K::Tetrahedron_3,
      typename K::Plane_3>::result_type result_type;

  typedef typename Intersection_traits<K,
      typename K::Triangle_3,
      typename K::Plane_3>::result_type Inter_type;

  typedef typename K::Segment_3 Segment_3;
  Inter_type intersections[4];
  int p_id = -1;
  std::vector<typename K::Point_3> points;
  std::vector<Segment_3> segments;
  for(int i = 0; i < 4; ++i)
  {
    const typename K::Triangle_3 triangle(tet.vertex((i+1)%4),
                                          tet.vertex((i+2)%4),
                                          tet.vertex((i+3)%4));
    intersections[i] = typename K::Intersect_3()(pl, triangle);
    if(intersections[i]){
      if(const typename K::Triangle_3* tr = boost::get<typename K::Triangle_3>(&*intersections[i]))
      {
        typename K::Triangle_3 res = *tr;
        return result_type(std::forward<typename K::Triangle_3>(res));
      }
      else if( const Segment_3* s
               = boost::get<Segment_3>(&*intersections[i]))
      {
        segments.push_back(*s);
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
      return result_type();
    else
    {
      typename K::Point_3 p
          = *boost::get<typename K::Point_3>(&*intersections[p_id]);

      return result_type(std::forward<typename  K::Point_3>(p));
    }
  }
    break;
  case 2:
  {
    return result_type(std::forward<typename  K::Segment_3>(segments.back()));
  }
    break;
  default:
  {
    std::set<typename K::Point_3> all_points;
    for (const auto& s : segments)
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
      return result_type(std::forward<typename K::Triangle_3>(result));
    }
    else //size = 4
    {
      std::list<Segment_3> segs(segments.begin(), segments.end());
      std::list<typename K::Point_3> tmp;
      fill_points_list(segs, tmp);
      std::vector<typename K::Point_3> res;
      for( const auto& p : tmp)
        res.push_back(p);
      return result_type(std::forward<std::vector<typename K::Point_3> >(res));
    }
  }
    break;
  }
  CGAL_assertion(false);
  return result_type();
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
