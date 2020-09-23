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


#ifndef CGAL_INTERNAL_INTERSECTIONS_3_TETRAHEDRON_3_TRIANGLE_3_INTERSECTIONS_H
#define CGAL_INTERNAL_INTERSECTIONS_3_TETRAHEDRON_3_TRIANGLE_3_INTERSECTIONS_H

#include <CGAL/kernel_basic.h>
#include <CGAL/intersections.h>
#include <CGAL/Intersections_3/internal/Triangle_3_Triangle_3_intersection.h>
#include <CGAL/Intersections_3/internal/Tetrahedron_3_Plane_3_intersection.h>
#include <CGAL/Intersections_3/internal/tetrahedron_intersection_helpers.h>
namespace CGAL {

namespace Intersections {

namespace internal {

//Tetrahedron_3 Segment_3
template <class K>
typename Intersection_traits<K, typename K::Tetrahedron_3, typename K::Triangle_3>::result_type
intersection(
    const typename K::Tetrahedron_3 &tet,
    const typename K::Triangle_3 &tr,
    const K& k)
{
  typedef typename Intersection_traits<K,
      typename K::Tetrahedron_3,
      typename K::Triangle_3>::result_type result_type;

  typedef typename Intersection_traits<K,
      typename K::Triangle_3,
      typename K::Triangle_3>::result_type Inter_type;

  typedef typename K::Segment_3 Segment_3;
  typedef typename K::Point_3 Point_3;
  typedef typename K::Triangle_3 Triangle_3;
  typedef std::vector<Point_3> Poly;

  std::vector<Point_3> inside_points;
  for(int i = 0; i< 3; ++i)
  {
    if(tet.has_on_bounded_side(tr.vertex(i))
       || tet.has_on_boundary(tr.vertex(i))){
      inside_points.push_back(tr.vertex(i));
    }
  }
  switch(inside_points.size())
  {
  case 0:
  {
    Inter_type intersections[4];
    std::vector<Segment_3> segments;
    std::vector<std::size_t> seg_ids;
    std::vector<Point_3> points;
    for(std::size_t i = 0; i < 4; ++i)
    {
     const typename K::Triangle_3 triangle(tet.vertex((i+1)%4),
                             tet.vertex((i+2)%4),
                             tet.vertex((i+3)%4));
      intersections[i] = typename K::Intersect_3()(tr, triangle);
      if(intersections[i]){
        //a face is inside the input tr
        if(const Triangle_3* t = boost::get<typename K::Triangle_3>(&*intersections[i]))
        {
          Triangle_3 res = *t;
          return result_type(std::forward<Triangle_3>(res));
        }
        //get segs and pts to construct poly
        else if( const Segment_3* s
                 = boost::get<Segment_3>(&*intersections[i]))
        {
          segments.push_back(*s);
          seg_ids.push_back(i);
        }
        else if( const typename K::Point_3* p
                 = boost::get<typename K::Point_3>(&*intersections[i]))
        {
          points.push_back(*p);
        }
        //if poly : then the input is in a supporting plane of a face, return the poly.
        else if( const Poly* p
                 = boost::get<Poly>(&*intersections[i]))
        {
          Poly res = *p;
          return result_type(std::forward<Poly>(res));
        }
      }
    }
    if(segments.size() > 1)
    {
      std::vector<Segment_3> filtered;
      filter_segments(segments, filtered);
      segments = filtered;
    }
    //if there are several segments, then we need to compute the polygone.
    if(segments.size() > 1)
    {
      std::list<Point_3> tmp;
      fill_segments_infos(segments,tmp, tr);
      Poly res;
      res.reserve(4);
      for(const auto& p : tmp)
        res.push_back(p);
      return result_type(std::forward<Poly>(res));
    }
    //else it must be adjacent to an vertex, so we return the point
    else if(segments.size() == 1)
    {
      //adjacency to an edge, return resulting segment.
      return result_type(std::forward<Segment_3>(segments.front()));
    }
    else
    {
      //no segment = adjacency to an vertex or an edge : return result point
      return result_type(std::forward<Point_3>(points.front()));

    }
  }
    break;
  case 1:
  case 2:
  {
    //tricky cases
    Inter_type intersections[4];
    std::vector<typename K::Point_3> points;
    std::vector<Segment_3> segments;
    for(std::size_t i = 0; i < 4; ++i)
    {
     const typename K::Triangle_3 triangle(tet.vertex((i+1)%4),
                             tet.vertex((i+2)%4),
                             tet.vertex((i+3)%4));
      intersections[i] = typename K::Intersect_3()(tr, triangle);
      if(intersections[i]){
        if(const Triangle_3* t = boost::get<typename K::Triangle_3>(&*intersections[i]))
        {
          Triangle_3 res = *t;
          return result_type(std::forward<Triangle_3>(res));
        }
        //get segs and pts to construct poly
        else if( const Segment_3* s
                 = boost::get<Segment_3>(&*intersections[i]))
        {
          segments.push_back(*s);
        }
        else if( const typename K::Point_3* p
                 = boost::get<typename K::Point_3>(&*intersections[i]))
        {
          points.push_back(*p);
        }
        //if poly : then the input is in a supporting plane of a face, return the poly.
        else if( const Poly* p
                 = boost::get<Poly>(&*intersections[i]))
        {
          Poly res = *p;
          return result_type(std::forward<Poly>(res));
        }
      }
    }
    if(segments.empty())
    {
      //then there is only one point of contact. Return it:
      return result_type(std::forward<Point_3>(points.front()));
    }

    if(segments.size() > 1)
    {
      std::vector<Segment_3> filtered;
      filter_segments(segments, filtered);
      segments = filtered;
    }

    switch(segments.size())
    {
    case 1:
    {
      bool return_solo_seg = true;
      //only one intersection, a triangle edge is one of the tet edges, and
      //the 3rd point is outside. This is the only intersection.
      for(const auto& p : inside_points)
      {
        if(!tet.has_on_boundary(p))
        {
          return_solo_seg = false;
          break;
        }
      }
      if(return_solo_seg)
      {
        return result_type(std::forward<Segment_3>(segments.front()));
      }

      if(inside_points.size() == 1)
      {
        Triangle_3 res(inside_points.front(), segments.front().source(),
                       segments.front().target());
        return result_type(std::forward<Triangle_3>(res));
      }
      else //size 2
      {
        Poly res(4);
        res[0] = inside_points.front();
        res[1] = inside_points.back();
        if((inside_points.front() - inside_points.back()) *
           (segments.front().source() - segments.front().target()) > 0)
        {
          res[2] = segments.front().target();
          res[3] = segments.front().source();
        }
        else
        {
          res[3] = segments.front().target();
          res[2] = segments.front().source();
        }
        return result_type(std::forward<Poly>(res));
      }
      }
      break;
    case 2:
    case 3:
    {
      std::list<Segment_3> segs(segments.begin(), segments.end());
      std::list<Point_3> tmp;
      fill_points_list(segs, tmp);
      if(inside_points.size() == 1)
      {
        Poly res;
        res.reserve(4);
        res.push_back(inside_points.front());
        for(const auto& p : tmp)
          res.push_back(p);
        return result_type(std::forward<Poly>(res));
      }
      else //size 2
      {
        Poly res;
        res.reserve(5);
        typename K::Compute_scalar_product_3 scalar =
          k.compute_scalar_product_3_object();
        typename K::Construct_vector_3 construct_vector =
          k.construct_vector_3_object();

        if(scalar(construct_vector(inside_points.front(), inside_points.back()),
           construct_vector(tmp.front(), tmp.back())) > 0)
        {
          res.push_back(inside_points.front());
        }
        res.push_back(inside_points.back());
        for(const auto& p : tmp)
          res.push_back(p);
        return result_type(std::forward<Poly>(res));
      }
    }
      break;
    default:
      //3 faces max if a point or more are inside tetrahedron
      break;
    }
  }
    break;

  case 3:
  {
    //triangle entirely inside tetra : return input triangle
    return result_type(tr);
  }
    break;
  default:
    //never happens (only 3 pts in a tr)
    break;
  }
  return result_type();
}

template <class K>
typename Intersection_traits<K, typename K::Tetrahedron_3, typename K::Triangle_3>::result_type
intersection(
    const typename K::Triangle_3 &pl,
    const typename K::Tetrahedron_3 &tet,
    const K& k)
{
  return intersection(tet, pl, k);
}

}}}
#endif // CGAL_INTERNAL_INTERSECTIONS_3_TETRAHEDRON_3_TRIANGLE_3_INTERSECTIONS_H
