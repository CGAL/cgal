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

#ifndef CGAL_INTERNAL_INTERSECTIONS_3_TETRAHEDRON_LINES_INTERSECTIONS_3_H
#define CGAL_INTERNAL_INTERSECTIONS_3_TETRAHEDRON_LINES_INTERSECTIONS_3_H

#include <vector>
#include <CGAL/kernel_basic.h>
#include <CGAL/intersections.h>

namespace CGAL {
namespace Intersections {
namespace internal {

template<class K, class O, class Derived>
struct Tetrahedron_lines_intersection_3_base
{
  typedef typename Intersection_traits<K,
  CGAL::Tetrahedron_3<K>,
  O>::result_type Result_type;

  typedef typename Intersection_traits<K,
  CGAL::Triangle_3<K>,
  O>::result_type Inter_type;

  const typename K::Tetrahedron_3& tet;
  const O& o;
  std::vector<typename K::Point_3> res_points;
  Result_type output;

  Tetrahedron_lines_intersection_3_base(const typename K::Tetrahedron_3& tet,
                                        const O& o):tet(tet), o(o)
  {}

  bool all_inside_test()
  {
    return false;
  }

  bool are_extremities_inside_test()
  {
    return false;
  }

  void do_procede()
  {
    if(static_cast<Derived*>(this)->all_inside_test())
      return;

    int res_id = -1;

    Inter_type tr_seg[4];
    for(int i = 0; i < 4; ++i)
    {
      const typename K::Triangle_3 triangle(tet.vertex((i+1)%4),
                                            tet.vertex((i+2)%4),
                                            tet.vertex((i+3)%4));
      if(do_intersect(o, triangle))
      {
        tr_seg[i] = typename K::Intersect_3()(o, triangle);
        if( boost::get<typename K::Segment_3>(&*tr_seg[i]) != nullptr)
        {
          res_id = i;
          break;
        }
      }
    }

    //if there is a segment in the intersections, then we return it
    if(res_id !=-1)
    {
      output = tr_seg[res_id];
      return;
    }

    //else if there is only 1 intersection
    res_points.reserve(4);
    res_id = -1;

    for(int i = 0; i< 4; ++i)
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

    if(res_points.empty())
    {
      return;
    }

    if(res_id != -1)
    {
      if(static_cast<Derived*>(this)->are_extremities_inside_test())
        return;
      //else point and segment entirely not inside or one intersection on boundary:
      output = tr_seg[res_id];
      return;
    }

    //else, we return a segment of the 2 intersection points (the most far away, in case of inexact)
    typename K::FT max_dist = 0;
    std::size_t res_id_2 = -1;
    std::vector<std::vector<typename K::FT> > sq_distances(res_points.size());

    for(std::size_t i = 0; i< res_points.size(); ++i)
    {
      const auto& p1 = res_points[i];
      for(const auto& p2 : res_points)
      {
        sq_distances[i].push_back(CGAL::squared_distance(p1, p2));
        if(sq_distances[i].back() > max_dist)
        {
          res_id = static_cast<int>(i);
          res_id_2 = sq_distances[i].size()-1;
          max_dist = sq_distances[i].back();
        }
      }
    }

    CGAL_assertion(res_id != -1);
    CGAL_assertion(res_id_2 != static_cast<std::size_t>(-1));
    CGAL_assertion(max_dist >0 );

    typename K::Segment_3 res_seg(res_points[res_id], res_points[res_id_2]);
    output = Result_type(std::forward<typename K::Segment_3>(res_seg));
    return;
  }
};

}
}
} //end namespaces
#endif // CGAL_INTERNAL_INTERSECTIONS_3_TETRAHEDRON_LINES_INTERSECTIONS_3_H
