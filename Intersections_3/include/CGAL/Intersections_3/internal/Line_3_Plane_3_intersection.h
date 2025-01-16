// Copyright (c) 2008  INRIA Sophia-Antipolis (France), ETH Zurich (Switzerland).
// Copyright (c) 2010, 2014  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Geert-Jan Giezeman,
//                 Sylvain Pion

#ifndef CGAL_INTERNAL_INTERSECTIONS_3_LINE_3_PLANE_3_INTERSECTION_H
#define CGAL_INTERNAL_INTERSECTIONS_3_LINE_3_PLANE_3_INTERSECTION_H

#include <CGAL/Intersection_traits_3.h>

#include <CGAL/wmult.h>

namespace CGAL {
namespace Intersections {
namespace internal {

template <class K>
std::optional<typename K::Point_3>
intersection_point(const typename K::Plane_3& plane,
                   const typename K::Line_3& line,
                   const K& /*k*/)
{
  typedef typename K::Point_3 Point_3;
  typedef typename K::Direction_3 Direction_3;
  typedef typename K::RT RT;

  const Point_3& line_pt = line.point();
  const Direction_3& line_dir = line.direction();

  RT num = plane.a()*line_pt.hx() + plane.b()*line_pt.hy()
           + plane.c()*line_pt.hz() + wmult_hw((K*)0, plane.d(), line_pt);
  RT den = plane.a()*line_dir.dx() + plane.b()*line_dir.dy()
           + plane.c()*line_dir.dz();
  if(den == 0)
    return std::nullopt;

  return std::make_optional(Point_3(den*line_pt.hx()-num*line_dir.dx(),
                                      den*line_pt.hy()-num*line_dir.dy(),
                                      den*line_pt.hz()-num*line_dir.dz(),
                                      wmult_hw((K*)0, den, line_pt)));
}

template <class K>
typename Intersection_traits<K, typename K::Plane_3, typename K::Line_3>::result_type
intersection(const typename K::Plane_3& plane,
             const typename K::Line_3& line,
             const K& /*k*/)
{
  typedef typename K::Point_3 Point_3;
  typedef typename K::Direction_3 Direction_3;
  typedef typename K::RT RT;

  const Point_3& line_pt = line.point();
  const Direction_3& line_dir = line.direction();

  RT num = plane.a()*line_pt.hx() + plane.b()*line_pt.hy()
           + plane.c()*line_pt.hz() + wmult_hw((K*)0, plane.d(), line_pt);
  RT den = plane.a()*line_dir.dx() + plane.b()*line_dir.dy()
           + plane.c()*line_dir.dz();

  if (den == 0) {
    if (num == 0) {
      // all line
      return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Line_3>(line);
    } else {
      // no intersection
      return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Line_3>();
    }
  }
  return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Line_3>(
           Point_3{den*line_pt.hx() - num*line_dir.dx(),
                   den*line_pt.hy() - num*line_dir.dy(),
                   den*line_pt.hz() - num*line_dir.dz(),
                   wmult_hw((K*)0, den, line_pt)});
}

template <class K>
inline
typename Intersection_traits<K, typename K::Plane_3, typename K::Line_3>::result_type
intersection(const typename K::Line_3& line,
             const typename K::Plane_3& plane,
             const K& k)
{
  return intersection(plane, line, k);
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_3_LINE_3_PLANE_3_INTERSECTION_H
