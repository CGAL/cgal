// Copyright (c) 2003  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Philippe Guigue

#ifndef CGAL_INTERNAL_INTERSECTIONS_PLANE_3_TRIANGLE_3_INTERSECTION_H
#define CGAL_INTERNAL_INTERSECTIONS_PLANE_3_TRIANGLE_3_INTERSECTION_H

#include <CGAL/Intersection_traits_3.h>
#include <CGAL/Intersections_3/internal/Line_3_Plane_3_intersection.h>

#include <CGAL/enum.h>
#include <CGAL/kernel_assertions.h>
#include <iterator>

namespace CGAL {
namespace Intersections {
namespace internal {

template <class K>
inline
typename Intersection_traits<K, typename K::Plane_3, typename K::Triangle_3>::result_type
intersection(const typename K::Plane_3& plane,
             const typename K::Triangle_3& tri,
             const K& k)
{
  typedef typename Intersection_traits<K, typename K::Plane_3, typename K::Line_3>::result_type pl_res;

  typename K::Construct_vertex_3 vertex_on = k.construct_vertex_3_object();

  CGAL::Oriented_side or0 = plane.oriented_side(vertex_on(tri, 0));
  CGAL::Oriented_side or1 = plane.oriented_side(vertex_on(tri, 1));
  CGAL::Oriented_side or2 = plane.oriented_side(vertex_on(tri, 2));

  if(or0 == ON_ORIENTED_BOUNDARY)
  {
    if(or1 == ON_ORIENTED_BOUNDARY)
    {
      if(or2 == ON_ORIENTED_BOUNDARY)
        return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Triangle_3>(tri);
      else
        return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Triangle_3>(
                k.construct_segment_3_object()(tri.vertex(0), tri.vertex(1)));
    }
    else
    {
      if(or2 == ON_ORIENTED_BOUNDARY)
      {
        return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Triangle_3>(
                k.construct_segment_3_object()(tri.vertex(0), tri.vertex(2)));
      }
      else
      {
        if(or1 == or2)
        {
          return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Triangle_3>(tri.vertex(0));
        }
        else
        {
          pl_res v = internal::intersection(plane, k.construct_line_3_object()(tri.vertex(1), tri.vertex(2)), k);
          const typename K::Point_3* p = intersect_get<typename K::Point_3>(v);
          CGAL_kernel_assertion(p!=nullptr);
          return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Triangle_3>(
                   k.construct_segment_3_object()(*p, tri.vertex(0)));
        }
      }
    }
  }

  if(or1 == ON_ORIENTED_BOUNDARY)
  {
    if(or2 == ON_ORIENTED_BOUNDARY)
      return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Triangle_3>(
               k.construct_segment_3_object()(tri.vertex(1), tri.vertex(2)));
    if(or2 == or0)
    {
      return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Triangle_3>(tri.vertex(1));
    }
    else
    {
      pl_res v = intersection(plane, k.construct_line_3_object()(tri.vertex(0), tri.vertex(2)), k);
      const typename K::Point_3* p = intersect_get<typename K::Point_3>(v);
      CGAL_kernel_assertion(p!=nullptr);
      return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Triangle_3>(
               k.construct_segment_3_object()(*p, tri.vertex(1)));
    }
  }

  if(or2 == ON_ORIENTED_BOUNDARY)
  {
    if(or1 == or0)
    {
      return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Triangle_3>(tri.vertex(2));
    }
    else
    {
      pl_res v = intersection(plane, k.construct_line_3_object()(tri.vertex(0),tri.vertex(1)), k);
      const typename K::Point_3* p = intersect_get<typename K::Point_3>(v);
      CGAL_kernel_assertion(p!=nullptr);
      return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Triangle_3>(
               k.construct_segment_3_object()(*p,tri.vertex(2)));
    }
  }

  //triangle vertices are not in the plane
  std::vector<typename K::Point_3> pts;
  pts.reserve(2);
  if(or0 != or1)
  {
    pl_res v = intersection(plane, k.construct_line_3_object()(tri.vertex(0),tri.vertex(1)), k);
    const typename K::Point_3* pt_ptr = intersect_get<typename K::Point_3>(v);
    CGAL_kernel_assertion(pt_ptr != nullptr);
    pts.push_back(*pt_ptr);
  }

  if(or0 != or2)
  {
    pl_res v = intersection(plane, k.construct_line_3_object()(tri.vertex(0),tri.vertex(2)), k);
    const typename K::Point_3* pt_ptr = intersect_get<typename K::Point_3>(v);
    CGAL_kernel_assertion(pt_ptr != nullptr);
    pts.push_back(*pt_ptr);
  }

  if(or1 != or2)
  {
    pl_res v = intersection(plane, k.construct_line_3_object()(tri.vertex(1),tri.vertex(2)), k);
    const typename K::Point_3* pt_ptr = intersect_get<typename K::Point_3>(v);
    CGAL_kernel_assertion(pt_ptr != nullptr);
    pts.push_back(*pt_ptr);
  }

  if(pts.empty())
    return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Triangle_3>();

  CGAL_kernel_assertion(pts.size() == 2);

  return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Triangle_3>(
           k.construct_segment_3_object()(*pts.begin(), *std::prev(pts.end())));
}

template <class K>
inline
typename Intersection_traits<K, typename K::Triangle_3, typename K::Plane_3>::result_type
intersection(const typename K::Triangle_3& triangle,
             const typename K::Plane_3& plane,
             const K& k)
{
  return intersection(plane, triangle, k);
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_PLANE_3_TRIANGLE_3_INTERSECTION_H
