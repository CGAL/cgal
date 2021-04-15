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
    const typename K::Plane_3 &plane,
    const K& k)
{
 typedef typename K::Point_3 Point_3;
  typedef typename K::Triangle_3 Triangle_3;
 typedef typename Intersection_traits<K, typename K::Tetrahedron_3, typename K::Plane_3>::result_type result_type;
  typename K::Oriented_side_3 oriented_side = k.oriented_side_3_object();

  std::vector<Point_3> corners(4);
  corners.reserve(8); // 4 corners + up to 4 polygon points
  corners[0] = tet[0];
  corners[1] = tet[1];
  corners[2] = tet[2];
  corners[3] = tet[3];

  const std::array<CGAL::Oriented_side, 4> orientations  { {
      oriented_side(plane, corners[0]),
      oriented_side(plane, corners[1]),
      oriented_side(plane, corners[2]),
      oriented_side(plane, corners[3])
    } };

  // description of faces of the bbox
  constexpr std::array<int, 12> face_indices
    { { 0, 1, 2,
        0, 1, 3,
        1, 2, 3,
        2, 0, 3 } };

  constexpr std::array<int, 12> edge_indices
    { { 0,  1,  2,
        0,  3,  5,
        1,  4,  3,
        2,  5,  4 } };

  std::array<int, 12> edge_ipt_id;
  edge_ipt_id.fill(-1);

  auto inter_pt_index =
    [&plane, &corners, &edge_ipt_id](int i, int j, int edge_id)
  {
    if (edge_ipt_id[edge_id]==-1)
    {
      edge_ipt_id[edge_id] = static_cast<int> (corners.size());
      corners.push_back(typename K::Construct_plane_line_intersection_point_3()
                      (plane, corners[i], corners[j]));
    }

    return edge_ipt_id[edge_id];
  };

  bool all_in = true;
  bool all_out = true;

  std::vector<std::array<int,2>> neighbor_ids(8, {-1,-1});

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
  for (int i = 0; i < 4; ++i)
  {
    ids.clear();
    for (int k = 0; k < 3; ++k)
    {

      int current_id = face_indices[3 * i + k];
      int next_id = face_indices[3 * i + (k + 1) % 3];
      int edge_id = edge_indices[3 * i + k];

      switch (orientations[current_id])
      {
        case ON_NEGATIVE_SIDE:
        {
          all_out = false;
          // check for intersection of the edge
          if (orientations[next_id] == ON_POSITIVE_SIDE)
          {
            ids.push_back(
                          inter_pt_index(current_id, next_id, edge_id));
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
      case 3:
      {
        Triangle_3 res(corners[ids[0]],
                       corners[ids[1]],
                       corners[ids[2]]);
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

  if (all_in || all_out) return boost::none;
  if (start_id == -1) return { result_type(corners[solo_id]) };

  int prv_id = -1;
  int cur_id = start_id;
  std::vector<Point_3> res;
  res.reserve(4);
  do {
    res.push_back(corners[cur_id]);
    int nxt_id = neighbor_ids[cur_id][0] == prv_id
               ? neighbor_ids[cur_id][1]
               : neighbor_ids[cur_id][0];
    if (nxt_id == -1 || nxt_id == start_id){
      if(res.size() == 2){
        typename K::Segment_3 seg(res[0], res[1]);
        return result_type(seg);
      }
      if(res.size() == 3){
        typename K::Triangle_3 tr(res[0], res[1], res[2]);
        return result_type(tr);
      }
      return result_type(res);
    }
    prv_id = cur_id;
    cur_id = nxt_id;
  } while (true);

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
