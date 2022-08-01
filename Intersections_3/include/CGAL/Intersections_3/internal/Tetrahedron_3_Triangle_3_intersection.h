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
  typedef typename K::Segment_3 Segment_3;
  typedef typename K::Plane_3 Plane_3;
  typedef typename K::Triangle_3 Triangle_3;
  typedef typename K::Line_3 Line_3;
  typedef typename K::Point_3 Point_3;

  typename K::Construct_plane_3 plane = k.construct_plane_3_object();
  typename K::Construct_vertex_3 vertex = k.construct_vertex_3_object();
  typename K::Construct_triangle_3 triangle = k.construct_triangle_3_object();
  typename K::Construct_segment_3 segment = k.construct_segment_3_object();
  typename K::Construct_line_3 line = k.construct_line_3_object();
  typename K::Oriented_side_3 oriented_side = k.oriented_side_3_object();


  std::vector<Point_3> res = { vertex(tr,0), vertex(tr,1), vertex(tr,2) };

  for (int i=0; i<4; ++i)
  {
    Plane_3 pl = plane(vertex(tet, (1+i)%4), vertex(tet, (2+i)%4),vertex(tet, (3+i)%4));
    if (oriented_side(pl, vertex(tet,i))!=ON_POSITIVE_SIDE) // TODO: this should be precomputed based on the tetra orientation
      pl = pl.opposite();

    std::vector<Point_3> current;
    std::vector<Oriented_side> orientations(res.size());
    for (std::size_t i=0; i<res.size(); ++i)
      orientations[i]=oriented_side(pl, res[i]);

    for (std::size_t i=0; i<res.size(); ++i)
    {
      const bool test_segment = i!=1 || res.size()!=2;
      std::size_t j = (i+1)%res.size();
      switch(orientations[j])
      {
        case ON_POSITIVE_SIDE:
          if (test_segment && orientations[i]==ON_NEGATIVE_SIDE)
            current.push_back(*CGAL::Intersections::internal::intersection_point(pl, line(res[i], res[j]), k));
          current.push_back(res[j]);
        break;
        case ON_NEGATIVE_SIDE:
          if (test_segment && orientations[i]==ON_POSITIVE_SIDE)
            current.push_back(*CGAL::Intersections::internal::intersection_point(pl, line(res[i], res[j]), k));
        break;
        default:
          current.push_back(res[j]);
      }
    }
    res.swap(current);
    if (res.empty())
      return boost::none;
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
