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
//                 Mael Rouxel-Labbé
//

#ifndef CGAL_INTERNAL_INTERSECTIONS_3_TETRAHEDRON_3_TRIANGLE_3_INTERSECTIONS_H
#define CGAL_INTERNAL_INTERSECTIONS_3_TETRAHEDRON_3_TRIANGLE_3_INTERSECTIONS_H

#include <CGAL/Intersections_3/internal/Line_3_Plane_3_intersection.h>

#include <CGAL/kernel_basic.h>

#include <algorithm>
#include <iterator>
#include <list>
#include <utility>
#include <vector>
#include <bitset>

namespace CGAL {
namespace Intersections {
namespace internal {

template <class K>
typename Intersection_traits<K, typename K::Tetrahedron_3, typename K::Triangle_3>::result_type
intersection(const typename K::Tetrahedron_3& tet,
             const typename K::Triangle_3& tr,
             const K& k)
{
  typedef typename Intersection_traits<K, typename K::Tetrahedron_3, typename K::Triangle_3>::result_type result_type;

  CGAL_precondition(!tet.is_degenerate());
  CGAL_precondition(!tr.is_degenerate());

  typedef typename K::Point_3 Point_3;
  typedef typename K::Plane_3 Plane_3;

  typename K::Construct_plane_3 plane = k.construct_plane_3_object();
  typename K::Construct_vertex_3 vertex = k.construct_vertex_3_object();
  typename K::Construct_triangle_3 triangle = k.construct_triangle_3_object();
  typename K::Construct_segment_3 segment = k.construct_segment_3_object();
  typename K::Construct_line_3 line = k.construct_line_3_object();
  typename K::Oriented_side_3 oriented_side = k.oriented_side_3_object();
  typename K::Orientation_3 orientation = k.orientation_3_object();


  std::vector<Point_3> res = { vertex(tr,0), vertex(tr,1), vertex(tr,2) };
  std::vector<std::bitset<4>> supporting_planes(3); // bitset used to indicate when a point is on a plane

  // iteratively clip `tr` with the halfspaces whose intersection form `tet`
  static constexpr std::array<int8_t, 12> vids = { 1,2,3, 0,3,2, 0,1,3, 1,0,2 };
  const bool tet_ori_positive = (orientation(tet)==POSITIVE);
  for (int pid=0; pid<4; ++pid)
  {
    Plane_3 pl = tet_ori_positive
               ? plane(vertex(tet, vids[pid*3]), vertex(tet, vids[pid*3+2]),vertex(tet, vids[pid*3+1]))
               : plane(vertex(tet, vids[pid*3]), vertex(tet, vids[pid*3+1]),vertex(tet, vids[pid*3+2]));
    CGAL_assertion(oriented_side(pl, vertex(tet,pid))==ON_POSITIVE_SIDE);

    std::vector<Point_3> current;
    std::vector<std::bitset<4>> current_sp;
    std::vector<Oriented_side> orientations(res.size());
    for (std::size_t i=0; i<res.size(); ++i)
    {
      orientations[i]=oriented_side(pl, res[i]);
      if (orientations[i]==ON_ORIENTED_BOUNDARY)
      {
        supporting_planes[i].set(pid);
        //workaround for kernels with inexact constructions
        //--
        if (supporting_planes[i].count()==3)
        {
          for (int b=0; i<4; ++b)
          {
            if (!supporting_planes[i].test(b))
            {
              res[i] = vertex(tet, b);
              break;
            }
          }
        }
        //--
      }
    }

    for (std::size_t i=0; i<res.size(); ++i)
    {
      const bool test_segment = i!=1 || res.size()!=2;
      std::size_t j = (i+1)%res.size();
      switch(orientations[j])
      {
        case ON_POSITIVE_SIDE:
          if (test_segment && orientations[i]==ON_NEGATIVE_SIDE)
          {
            current_sp.push_back(supporting_planes[i] & supporting_planes[j]);
            current_sp.back().set(pid);
            if (current_sp.back().count()==3)
            {
              for (int b=0; i<4; ++b)
                if (!current_sp.back().test(b))
                {
                  current.push_back(vertex(tet, b));
                  break;
                }
            }
            else
              current.push_back(*CGAL::Intersections::internal::intersection_point(pl, line(res[i], res[j]), k));
          }
          current.push_back(res[j]);
          current_sp.push_back(supporting_planes[j]);
        break;
        case ON_NEGATIVE_SIDE:
          if (test_segment && orientations[i]==ON_POSITIVE_SIDE)
          {
            current_sp.push_back(supporting_planes[i] & supporting_planes[j]);
            current_sp.back().set(pid);
            if (current_sp.back().count()==3)
            {
              for (int b=0; i<4; ++b)
                if (!current_sp.back().test(b))
                {
                  current.push_back(vertex(tet, b));
                  break;
                }
            }
            else
              current.push_back(*CGAL::Intersections::internal::intersection_point(pl, line(res[i], res[j]), k));
          }
        break;
        default:
        {
          CGAL_assertion(supporting_planes[j].test(pid));
          current.push_back(res[j]);
          current_sp.push_back(supporting_planes[j]);
        }
      }
    }
    res.swap(current);
    supporting_planes.swap(current_sp);

    if (res.empty())
      return std::nullopt;
  }

  switch(res.size())
  {
    case 1:
      return result_type(res[0]);
    case 2:
      return result_type(segment(res[0], res[1]));
    case 3:
      return result_type(triangle(res[0], res[1], res[2]));
    default:
      return result_type(res);
  }
}

template <class K>
typename Intersection_traits<K, typename K::Tetrahedron_3, typename K::Triangle_3>::result_type
intersection(const typename K::Triangle_3& pl,
             const typename K::Tetrahedron_3& tet,
             const K& k)
{
  return intersection(tet, pl, k);
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_3_TETRAHEDRON_3_TRIANGLE_3_INTERSECTIONS_H
