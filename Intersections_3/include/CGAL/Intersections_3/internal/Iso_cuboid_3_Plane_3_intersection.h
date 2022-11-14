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
// Author(s)     : Sebastien Loriot and Andreas Fabri
//

#ifndef CGAL_INTERNAL_INTERSECTIONS_3_ISO_CUBOID_3_PLANE_3_INTERSECTION_H
#define CGAL_INTERNAL_INTERSECTIONS_3_ISO_CUBOID_3_PLANE_3_INTERSECTION_H

#include <CGAL/enum.h>
#include <CGAL/kernel_basic.h>
#include <CGAL/intersections.h>
#include <CGAL/utility.h>

#include <boost/none.hpp>

#include <vector>
#include <array>

namespace CGAL {
namespace Intersections {
namespace internal {

//Iso_cuboid_3 Plane_3
template <class K>
typename Intersection_traits<K, typename K::Iso_cuboid_3, typename K::Plane_3>::result_type
intersection(const typename K::Iso_cuboid_3& cub,
             const typename K::Plane_3& plane,
             const K& k)
{
  typedef typename K::Point_3 Point_3;
  typedef typename Intersection_traits<K, typename K::Iso_cuboid_3, typename K::Plane_3>::result_type result_type;
  typename K::Oriented_side_3 oriented_side = k.oriented_side_3_object();

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

  const std::array<CGAL::Oriented_side, 8> orientations  { {
      oriented_side(plane, corners[0]),
      oriented_side(plane, corners[1]),
      oriented_side(plane, corners[2]),
      oriented_side(plane, corners[3]),
      oriented_side(plane, corners[4]),
      oriented_side(plane, corners[5]),
      oriented_side(plane, corners[6]),
      oriented_side(plane, corners[7])
    } };

  // description of faces of the bbox
  static constexpr std::array<int, 24> face_indices
    { { 0, 1, 2, 3,
        2, 1, 5, 6,
        3, 2, 6, 7,
        1, 0, 4, 5,
        4, 0, 3, 7,
        6, 5, 4, 7 } };

  static constexpr std::array<int, 24> edge_indices
    { { 0,  1,  2, 3,
        1,  4,  5, 6,
        2,  6,  7, 8,
        0,  9, 10, 4,
        9,  3,  8, 11,
        5, 10, 11, 7 } };

  std::array<int, 12> edge_ipt_id;
  edge_ipt_id.fill(-1);

  auto inter_pt_index =
    [&k, &plane, &corners, &edge_ipt_id](int i, int j, int edge_id)
  {
    if (edge_ipt_id[edge_id]==-1)
    {
      edge_ipt_id[edge_id] = static_cast<int> (corners.size());
      corners.push_back(k.construct_plane_line_intersection_point_3_object()(plane, corners[i], corners[j]));
    }

    return edge_ipt_id[edge_id];
  };

  bool all_in = true;
  bool all_out = true;

  std::vector<std::array<int,2> > neighbor_ids(14, {-1,-1});

  auto add_neighbor = [&neighbor_ids](int i, int j)
  {
    if (neighbor_ids[i][0] == -1 ) {
      neighbor_ids[i][0] = j;
    }
    else {
      if (neighbor_ids[i][0]!=j && neighbor_ids[i][1]==-1)
      {
        neighbor_ids[i][1] = j;
      }
    }
  };

  int start_id = -1;
  int solo_id = -1;
  // for each face of the bbox, we look for intersection of the plane with its edges
  std::vector<int> ids;
  for (int i = 0; i < 6; ++i)
  {
    ids.clear();
    for (int k = 0; k < 4; ++k)
    {

      int current_id = face_indices[4 * i + k];
      int next_id = face_indices[4 * i + (k + 1) % 4];
      int edge_id = edge_indices[4 * i + k];

      switch (orientations[current_id])
      {
        case ON_NEGATIVE_SIDE:
        {
          all_out = false;
          // check for intersection of the edge
          if (orientations[next_id] == ON_POSITIVE_SIDE)
          {
            ids.push_back(inter_pt_index(current_id, next_id, edge_id));
          }
          break;
        }
        case ON_POSITIVE_SIDE:
        {
          all_in = false;
          // check for intersection of the edge
          if (orientations[next_id] == ON_NEGATIVE_SIDE)
          {
            ids.push_back(inter_pt_index(current_id, next_id, edge_id));
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

    switch (ids.size())
    {
      case 4:
      {
        std::vector<Point_3> res({ corners[ids[0]],
                                   corners[ids[1]],
                                   corners[ids[2]],
                                   corners[ids[3]] });
        return result_type(res);
      }
      case 2:
      {
        if (start_id == -1) start_id = ids[0];
        add_neighbor(ids[0], ids[1]);
        add_neighbor(ids[1], ids[0]);
        break;
      }
      case 1:
        solo_id = ids[0];
      default:
        break;
    }
  }

  if (all_in || all_out)
    return boost::none;
  if (start_id == -1)
    return { result_type(corners[solo_id]) };

  int prv_id = -1;
  int cur_id = start_id;
  std::vector<Point_3> res;
  res.reserve(6);

  for (;;) {
    res.push_back(corners[cur_id]);
    int nxt_id = neighbor_ids[cur_id][0] == prv_id
               ? neighbor_ids[cur_id][1]
               : neighbor_ids[cur_id][0];
    if (nxt_id == -1 || nxt_id == start_id) {
      if(res.size() == 2) {
        typename K::Segment_3 seg(res[0], res[1]);
        return result_type(seg);
      }
      if(res.size() == 3) {
        typename K::Triangle_3 tr(res[0], res[1], res[2]);
        return result_type(tr);
      }
      return result_type(res);
    }
    prv_id = cur_id;
    cur_id = nxt_id;
  }
}

template <class K>
typename Intersection_traits<K, typename K::Iso_cuboid_3, typename K::Plane_3>::result_type
intersection(const typename K::Plane_3& pl,
             const typename K::Iso_cuboid_3& cub,
             const K& k)
{
  return intersection(cub, pl, k);
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_3_ISO_CUBOID_3_PLANE_3_INTERSECTION_H
