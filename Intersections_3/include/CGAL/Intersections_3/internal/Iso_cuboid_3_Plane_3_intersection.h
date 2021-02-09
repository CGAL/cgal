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
#include <CGAL/utility.h>

#include <CGAL/Intersections_3/internal/tetrahedron_intersection_helpers.h>
#include <CGAL/Intersections_3/Iso_cuboid_3_Segment_3.h>
#include <CGAL/Intersections_3/Plane_3_Plane_3.h>

#include <set>

namespace CGAL {
namespace Intersections {
namespace internal {

template <class Geom_traits, class Plane_3, class Point_3>
int
inter_pt_index(int i, int j,
               const Plane_3& plane,
               std::vector<Point_3>& points,
               std::map<std::pair<int, int>, int>& id_map)
{
  std::pair<std::map<std::pair<int, int>, int>::iterator, bool> res =
    id_map.insert(std::make_pair(make_sorted_pair(i, j),
                                 static_cast<int> (points.size())));
  if (res.second)
    points.push_back(typename Geom_traits::Construct_plane_line_intersection_point_3()
                    (plane, points[i], points[j]));

  return res.first->second;
}

//Iso_cuboid_3 Plane_3
template <class K>
typename Intersection_traits<K, typename K::Iso_cuboid_3, typename K::Plane_3>::result_type
intersection(
             const typename K::Iso_cuboid_3& cub,
             const typename K::Plane_3& plane,
             const K& k)
{
  typedef typename K::Point_3 Point_3;
  typedef std::vector<Point_3> Poly;
  typedef typename Intersection_traits<K, typename K::Iso_cuboid_3, typename K::Plane_3>::result_type result_type;

  std::vector<Point_3> corners(8);
  corners.reserve(14); // 8 corners + up to 6 polygon points
  corners[0] = cub[0];
  corners[1] = cub[3];
  corners[2] = cub[2];
  corners[3] = cub[1];
  corners[4] = cub[5];
  corners[5] = cub[4];
  corners[6] = cub[7];
  corners[7] = cub[6];

  std::array<CGAL::Oriented_side, 8> orientations = { {
      plane.oriented_side(corners[0]),
      plane.oriented_side(corners[1]),
      plane.oriented_side(corners[2]),
      plane.oriented_side(corners[3]),
      plane.oriented_side(corners[4]),
      plane.oriented_side(corners[5]),
      plane.oriented_side(corners[6]),
      plane.oriented_side(corners[7])
    } };

  // description of faces of the bbox
  std::array<int, 24> face_indices =
    { { 0, 1, 2, 3,
        2, 1, 5, 6,
        3, 2, 6, 7,
        1, 0, 4, 5,
        4, 0, 3, 7,
        6, 5, 4, 7 } };

  std::map<std::pair<int, int>, int> id_map;
  bool all_in = true;
  bool all_out = true;

  std::vector<int> next(14, -1);
  std::vector<int> prev(14, -1);

  int start = -1;
  int solo = -1;
  // for each face of the bbox, we look for intersection of the plane with its edges
  std::vector<int> ids;
  for (int i = 0; i < 6; ++i)
    {
      ids.clear();
      for (int k = 0; k < 4; ++k)
        {

          int current_id = face_indices[4 * i + k];
          int next_id = face_indices[4 * i + (k + 1) % 4];

          switch (orientations[current_id])
            {
            case ON_NEGATIVE_SIDE:
              {
                all_out = false;
                // check for intersection of the edge
                if (orientations[next_id] == ON_POSITIVE_SIDE)
                  {
                    ids.push_back(
                                  inter_pt_index<K>(current_id, next_id, plane, corners, id_map));
                  }
                break;
              }
            case ON_POSITIVE_SIDE:
              {
                all_in = false;
                // check for intersection of the edge
                if (orientations[next_id] == ON_NEGATIVE_SIDE)
                  {
                    ids.push_back(
                                  inter_pt_index<K>(current_id, next_id, plane, corners, id_map));
                  }
                break;
              }
            case ON_ORIENTED_BOUNDARY:
              {
                all_in = all_out = false;
                ids.push_back(current_id);
              }
            }
        }
      if (ids.size() == 4){
        std::vector<Point_3> res(4);
        res[0] = corners[ids[0]];
        res[1] = corners[ids[1]];
        res[2] = corners[ids[2]];
        res[3] = corners[ids[3]];
        return result_type(std::forward<Poly>(res));
      } else
        {
          if (ids.size() == 2)
            {
              if (start == -1) start = ids[0];
              if (next[ids[0]] == -1) {
                next[ids[0]] = ids[1];
              }
              else {
                prev[ids[0]] = ids[1];
              }
              if (next[ids[1]] == -1) {
                next[ids[1]] = ids[0];
              }
              else {
                prev[ids[1]] = ids[0];
              }

            }
          else
            if (ids.size() == 1)
              solo = ids[0];
        }
    }

  if (all_in || all_out) return boost::none;
  if (start == -1) return { result_type(corners[solo]) };

  int pre = -1;
  int cur = start;
  std::vector<Point_3> res;
  res.reserve(6);
  do {
    res.push_back(corners[cur]);
    int n = next[cur] == pre ? prev[cur] : next[cur];
    if (n == -1 || n == start){
      if(res.size() == 2){
        typename K::Segment_3 seg(res[0], res[1]);
        return result_type(std::forward<typename K::Segment_3>(seg));
      }
      if(res.size() == 3){
        typename K::Triangle_3 tr(res[0], res[1], res[2]);
        return result_type(std::forward<typename K::Triangle_3>(tr));
      }
      return result_type(std::forward<Poly>(res));;
    }
    pre = cur;
    cur = n;
  } while (true);
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
