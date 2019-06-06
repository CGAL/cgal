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

#ifndef CGAL_INTERNAL_INTERSECTIONS_3_ISO_CUBOID_3_PLANE_3_INTERSECTION_H
#define CGAL_INTERNAL_INTERSECTIONS_3_ISO_CUBOID_3_PLANE_3_INTERSECTION_H
#include <CGAL/kernel_basic.h>
#include <CGAL/intersections.h>
#include <CGAL/Intersections_3/internal/tetrahedron_intersection_helpers.h>
#include <CGAL/Intersections_3/Iso_cuboid_3_Segment_3.h>
#include <CGAL/Intersections_3/Plane_3_Plane_3.h>

#include <set>
namespace CGAL {

namespace Intersections {

namespace internal {

template<typename Point>
void filter_points(const std::vector<Point>& input,
                     std::vector<Point>& output)
{
  std::set<Point> tmp;
  for( auto p : input)
    tmp.insert(p);
  for(auto p: tmp)
    output.push_back(p);
}

//Tetrahedron_3 Line_3
template <class K>
typename Intersection_traits<K, typename K::Iso_cuboid_3, typename K::Plane_3>::result_type
intersection(
    const typename K::Iso_cuboid_3 &cub,
    const typename K::Plane_3 &pl,
    const K&)
{
  typedef typename K::Point_3 Point_3;
  typedef typename K::Segment_3 Segment_3;
  typedef typename K::Line_3 Line_3;
  typedef typename K::Plane_3 Plane_3;
  typedef std::vector<Point_3> Poly;

  typedef typename Intersection_traits<K,
      CGAL::Iso_cuboid_3<K>,
      CGAL::Plane_3<K> >::result_type Result_type;

  typedef typename Intersection_traits<K,
      CGAL::Segment_3<K>,
      CGAL::Plane_3<K> >::result_type Inter_type;

  std::vector<Segment_3> edges;
  edges.reserve(12);

  //get all edges of cub
  for(int i=0; i< 4; ++i)
  {
    edges.push_back(Segment_3(cub.vertex(i), cub.vertex((i+1)%4)));
  }

  for(int i=0; i < 4; ++i)
  {
    edges.push_back(Segment_3(cub.vertex(i+4), cub.vertex((i+1)%4+4)));
  }

  for(int i=0; i < 4; ++i)
  {
    edges.push_back(Segment_3(cub.vertex(i), cub.vertex((i+1)%4+4)));
  }
  //get all intersections between pl and cub edges
  std::vector<Segment_3> segments;
  std::vector<Point_3> points;

  for(int i=0; i < 12; ++i)
  {
    Inter_type inter = typename K::Intersect_3()(pl, edges[i]);
    if(inter){
      if(const Segment_3* seg = boost::get<Segment_3>(&*inter))
      {
        segments.push_back(*seg);
      }
      else if(const Point_3* p = boost::get<Point_3>(&*inter))
      {
        points.push_back(*p);
      }
    }
  }

  switch(segments.size())
  {
  case 1: //adj to an edge
  {
    return Result_type(std::forward<Segment_3>(segments.front()));
  }
    break;

  case 2: //intersects diagonally
  {
    Poly res(4);
    Segment_3 front(segments.front()),
        back(segments.back());
    res[0] = front.target();
    if((front.target() - front.source())
       * (back.target() - back.source()) > 0)
    {
      res[1] = back.target();
      res[2] = back.source();
    }
    else
    {
      res[1] = back.source();
      res[2] = back.target();
    }
    res[3] = front.source();

    return Result_type(std::forward<Poly>(res));
  }
    break;
  case 4: // intersects a face
  {
    Poly res;
    res.reserve(4);
    std::list<Point_3> tmp;
    std::list<Segment_3> seg_list;
    for(auto s : segments)
      seg_list.push_back(s);
    fill_points_list(seg_list, tmp);
    for(auto p : tmp)
      res.push_back(p);
    return Result_type(std::forward<Poly>(res));
  }
    break;
  default:
    break;
  }

  Poly filtered_points;
  filter_points(points, filtered_points);
  if(filtered_points.empty())
    return Result_type();
  //adjacent to a vertex
  if(filtered_points.size() == 1)
  {
    return Result_type(std::forward<Point_3>(filtered_points.front()));
  }
  else
  {
    //get intersections between pl and each face -> line. Foreach line, creates segment with points. Then use helper_function to recover polygon.
    typedef typename Intersection_traits<K,
        CGAL::Plane_3<K>,
        CGAL::Plane_3<K> >::result_type Pl_pl_type;

    std::vector<Line_3> plane_intersections;
    Pl_pl_type pl_inter = CGAL::intersection(pl, Plane_3(cub.vertex(0),
                                                         cub.vertex(1),
                                                         cub.vertex(5)));
    if(const Line_3* line = boost::get<Line_3>(&*pl_inter)){
      plane_intersections.push_back(*line);
    }
    pl_inter = CGAL::intersection(pl, Plane_3(cub.vertex(0),
                                              cub.vertex(3),
                                              cub.vertex(4)));
    if(const Line_3* line = boost::get<Line_3>(&*pl_inter)){
      plane_intersections.push_back(*line);
    }
    pl_inter = CGAL::intersection(pl, Plane_3(cub.vertex(0),
                                              cub.vertex(1),
                                              cub.vertex(3)));
    if(const Line_3* line = boost::get<Line_3>(&*pl_inter)){
      plane_intersections.push_back(*line);
    }
    pl_inter = CGAL::intersection(pl, Plane_3(cub.vertex(7),
                                              cub.vertex(6),
                                              cub.vertex(1)));
    if(const Line_3* line = boost::get<Line_3>(&*pl_inter)){
      plane_intersections.push_back(*line);
    }
    pl_inter = CGAL::intersection(pl, Plane_3(cub.vertex(7),
                                              cub.vertex(4),
                                              cub.vertex(3)));
    if(const Line_3* line = boost::get<Line_3>(&*pl_inter)){
      plane_intersections.push_back(*line);
    }
    pl_inter = CGAL::intersection(pl, Plane_3(cub.vertex(7),
                                              cub.vertex(6),
                                              cub.vertex(4)));
    if(const Line_3* line = boost::get<Line_3>(&*pl_inter)){
      plane_intersections.push_back(*line);
    }

    std::list<Segment_3> tmp_segs;
    for(auto line : plane_intersections)
    {
      bool first_found = false;
      Point_3 first_p;
      for(auto p : filtered_points)
      {
        if(line.has_on(p))
        {
          if(!first_found)
          {
            first_found = true;
            first_p = p;
          }
          else
          {
            tmp_segs.push_back(Segment_3(first_p, p));
            break;
          }
        }
      }
    }
    if(tmp_segs.size() < 3)
      return Result_type();
    std::list<Point_3> tmp_pts;
    fill_points_list(tmp_segs,tmp_pts);
    Poly res;
    for(auto p : tmp_pts)
      res.push_back(p);
    if(res.size() == 3){
      typename K::Triangle_3 tr(res[0], res[1], res[2]);
      return Result_type(std::forward<typename K::Triangle_3>(tr));
    }
    else
    {
      return Result_type(std::forward<Poly>(res));
    }
  }
}

template <class K>
typename Intersection_traits<K, typename K::Iso_cuboid_3, typename K::Plane_3>::result_type
intersection(
    const typename K::Plane_3 &pl,
    const typename K::Iso_cuboid_3 &cub,
    const K& k)
{
  return intersection(cub, pl, k);
}

}}}

#endif // CGAL_INTERNAL_INTERSECTIONS_3_ISO_CUBOID_3_PLANE_3_INTERSECTION_H
