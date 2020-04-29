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
void filter_points(std::vector<Point>& input,
                   std::vector<Point>& output)
{
  std::sort(input.begin(), input.end());
  auto last = std::unique(input.begin(), input.end());
  input.erase(last, input.end());
  output = std::move(input);
}

//Tetrahedron_3 Line_3
template <class K>
typename Intersection_traits<K, typename K::Iso_cuboid_3, typename K::Plane_3>::result_type
intersection(
    const typename K::Iso_cuboid_3 &cub,
    const typename K::Plane_3 &pl,
    const K& k)
{
  typedef typename K::Point_3 Point_3;
  typedef typename K::Segment_3 Segment_3;
  typedef typename K::Line_3 Line_3;
  typedef typename K::Plane_3 Plane_3;
  typedef std::vector<Point_3> Poly;

  typedef typename Intersection_traits<K, CGAL::Iso_cuboid_3<K>,
      CGAL::Plane_3<K> >::result_type result_type;

  typedef typename Intersection_traits<K, CGAL::Segment_3<K>,
      CGAL::Plane_3<K> >::result_type Inter_type;

  std::vector<Segment_3> edges;
  edges.reserve(12);

  //get all edges of cub
  for(int i=0; i< 4; ++i)
    edges.emplace_back(cub.vertex(i), cub.vertex((i+1)%4));
  for(int i=0; i < 4; ++i)
    edges.emplace_back(cub.vertex(i+4), cub.vertex((i+1)%4+4));
  for(int i=0; i < 4; ++i)
    edges.emplace_back(cub.vertex(i), cub.vertex((i+1)%4+4));

  //get all intersections between pl and cub edges
  std::vector<Segment_3> segments;
  std::vector<Point_3> points;

  for(int i=0; i < 12; ++i)
  {
    // Intersect_3 checks the orientation of the segment's extremities to avoid actually computing
    // the intersection if possible
    Inter_type inter = typename K::Intersect_3()(pl, edges[i]);
    if(inter)
    {
      if(const Segment_3* seg = boost::get<Segment_3>(&*inter))
        segments.push_back(*seg);
      else if(const Point_3* p = boost::get<Point_3>(&*inter))
        points.push_back(*p);
    }
  }

  if(segments.empty() && points.empty())
    return result_type();

  switch(segments.size())
  {
    case 0:
    //points dealt with later
    break;
    case 1:
    {
      //adj to an edge
      if(points.size() == 4)
      {
        return result_type(std::forward<Segment_3>(segments.front()));
      }
      //plane intersecting through an edge (not 2)
      else
      {
        Poly res(4);
        const Segment_3& entry_seg(segments.front());
        Point_3 p1, p2;
        bool p1_found(false),
             p2_found(false);

        for(const Point_3& p : points)
        {
          if(!k.equal_3_object()(entry_seg.source(), p)
             && ! k.equal_3_object()(entry_seg.target(), p))
          {
            if(!p1_found)
            {
              p1 = p;
              p1_found = true;
            }
            else {
              p2 = p;
              p2_found = true;
              break;
            }
          }
        }
        CGAL_USE(p2_found);
        CGAL_assertion(p1_found && p2_found);
        res[0] = entry_seg.target();
        res[1] = p2;
        res[2] = p1;
        res[3] = entry_seg.source();

        typename K::Coplanar_orientation_3 coplanar_orientation =
          k.coplanar_orientation_3_object();

        if( coplanar_orientation(res[0], res[1], res[2])
            != coplanar_orientation(res[0], res[1], res[3]))
        {
          std::swap(res[1], res[2]);
        }

        return result_type(std::forward<Poly>(res));
      }
    }
      break;

    case 2: //intersects diagonally
    {
      Poly res(4);
      Segment_3 &front(segments.front()),
          &back(segments.back());
      res[0] = front.target();
      res[1] = back.target();
      res[2] = back.source();
      res[3] = front.source();
      typename K::Coplanar_orientation_3 coplanar_orientation =
        k.coplanar_orientation_3_object();

      if( coplanar_orientation(res[0], res[1], res[2])
          != coplanar_orientation(res[0], res[1], res[3]))
      {
        std::swap(res[1], res[2]);
      }

      return result_type(std::forward<Poly>(res));
    }
    break;
  case 4: // intersects a face
  {
    Poly res;
    res.reserve(4);
    typename K::Has_on_3 has_on;
    if(has_on(pl, cub.vertex(0))
       && has_on(pl, cub.vertex(5))
       && has_on(pl, cub.vertex(4)))
    {
      res.push_back(cub.vertex(0));
      res.push_back(cub.vertex(5));
      res.push_back(cub.vertex(4));
      res.push_back(cub.vertex(3));
    }
    else if(has_on(pl, cub.vertex(0))
            && has_on(pl, cub.vertex(1))
            && has_on(pl, cub.vertex(6)))
    {
      res.push_back(cub.vertex(0));
      res.push_back(cub.vertex(1));
      res.push_back(cub.vertex(6));
      res.push_back(cub.vertex(5));

    }
    else if(has_on(pl, cub.vertex(1))
            && has_on(pl, cub.vertex(2))
            && has_on(pl, cub.vertex(7)))
    {
      res.push_back(cub.vertex(1));
      res.push_back(cub.vertex(2));
      res.push_back(cub.vertex(7));
      res.push_back(cub.vertex(6));
    }
    else if(has_on(pl, cub.vertex(2))
            && has_on(pl, cub.vertex(3))
            && has_on(pl, cub.vertex(4)))
    {
      res.push_back(cub.vertex(2));
      res.push_back(cub.vertex(3));
      res.push_back(cub.vertex(4));
      res.push_back(cub.vertex(7));
    }
    else if(has_on(pl, cub.vertex(6))
            && has_on(pl, cub.vertex(7))
            && has_on(pl, cub.vertex(4)))
    {
      res.push_back(cub.vertex(6));
      res.push_back(cub.vertex(7));
      res.push_back(cub.vertex(4));
      res.push_back(cub.vertex(5));
    }
    else if(has_on(pl, cub.vertex(2))
            && has_on(pl, cub.vertex(1))
            && has_on(pl, cub.vertex(0)))
    {
      res.push_back(cub.vertex(2));
      res.push_back(cub.vertex(1));
      res.push_back(cub.vertex(0));
      res.push_back(cub.vertex(3));
    }
    return result_type(std::forward<Poly>(res));
  }
    break;
  default:
    break;
  }

  Poly filtered_points;
  filter_points(points, filtered_points);

  //adjacent to a vertex
  if(filtered_points.size() == 1)
  {
    return result_type(std::forward<Point_3>(filtered_points.front()));
  }

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
  for(const auto& line : plane_intersections)
  {
    bool first_found = false;
    Point_3 first_p;
    typename K::Has_on_3 has_on;
    for(const auto  &p : filtered_points)
    {
      if(has_on(line, p))
      {
        if(!first_found)
        {
          first_found = true;
          first_p = p;
        }
        else
        {
          tmp_segs.emplace_back(first_p, p);
          break;
        }
      }
    }
  }

  if(tmp_segs.size() < 3)
    return result_type();

  std::list<Point_3> tmp_pts;
  fill_points_list(tmp_segs,tmp_pts);

  Poly res;
  for(const auto& p : tmp_pts)
    res.push_back(p);

  if(res.size() == 3){
    typename K::Triangle_3 tr(res[0], res[1], res[2]);
    return result_type(std::forward<typename K::Triangle_3>(tr));
  }
  else
  {
    return result_type(std::forward<Poly>(res));
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
