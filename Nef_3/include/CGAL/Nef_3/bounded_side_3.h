// Copyright (c) 1997-2002  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael Seel    <seel@mpi-sb.mpg.de>
//                 Miguel Granados <granados@mpi-sb.mpg.de>
//                 Susan Hert      <hert@mpi-sb.mpg.de>
//                 Lutz Kettner    <kettner@mpi-sb.mpg.de>
#ifndef CGAL_BOUNDED_SIDE_3_H
#define CGAL_BOUNDED_SIDE_3_H

#include <CGAL/license/Nef_3.h>


#include <CGAL/basic.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/Iterator_project.h>
#include <CGAL/representation_tags.h>
#include <vector>

#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 17
#include <CGAL/Nef_2/debug.h>

namespace CGAL {

template <class Point_2, class Point_3>
Point_2 point_3_get_x_y_point_2(const Point_3& p, const Homogeneous_tag&) {
  return( Point_2(p.hx(), p.hy(), p.hw()) );
}

template <class Point_2, class Point_3>
Point_2 point_3_get_y_z_point_2(const Point_3& p, const Homogeneous_tag&) {
  return( Point_2(p.hy(), p.hz(), p.hw()) );
}

template <class Point_2, class Point_3>
Point_2 point_3_get_z_x_point_2(const Point_3& p, const Homogeneous_tag&) {
  return( Point_2(p.hz(), p.hx(), p.hw()) );
}

template <class Point_2, class Point_3>
Point_2 point_3_get_x_y_point_2(const Point_3& p, const Cartesian_tag&) {
  return( Point_2(p.x(), p.y()) );
}

template <class Point_2, class Point_3>
Point_2 point_3_get_y_z_point_2(const Point_3& p, const Cartesian_tag&) {
  return( Point_2(p.y(), p.z()) );
}

template <class Point_2, class Point_3>
Point_2 point_3_get_z_x_point_2(const Point_3& p, const Cartesian_tag&) {
  return( Point_2(p.z(), p.x()) );
}

template <class Point_2, class Point_3, class R>
Point_2 point_3_get_x_y_point_2(const Point_3& p) {
    return point_3_get_x_y_point_2<Point_2,Point_3>(p,typename R::Kernel_tag());
}

template <class Point_2, class Point_3, class R>
Point_2 point_3_get_y_z_point_2(const Point_3& p) {
    return point_3_get_y_z_point_2<Point_2,Point_3>(p,typename R::Kernel_tag());
}

template <class Point_2, class Point_3, class R>
Point_2 point_3_get_z_x_point_2(const Point_3& p) {
    return point_3_get_z_x_point_2<Point_2,Point_3>(p,typename R::Kernel_tag());
}

template <class IteratorForward, class R>
Bounded_side bounded_side_3(IteratorForward first,
                            IteratorForward last,
                            const Point_3<R>& point,
                            typename R::Plane_3 plane = typename R::Plane_3(0,0,0,0)) {
  typedef typename R::Point_2 Point_2;
  typedef typename R::Point_3 Point_3;
  typedef typename R::Plane_3 Plane_3;

  if(plane == Plane_3(0,0,0,0)) {
    // TO TEST: code never tested
    IteratorForward p(first);
    Point_3 p0(*(p++));
    CGAL_assertion(p != last);
    Point_3 p1(*(p++));
    CGAL_assertion(p != last);
    Point_3 p2(*(p++));
    plane = Plane_3(p0, p1, p2);

    /* since we just need to project the points to a non-perpendicular plane
       we don't need to care about the plane orientation */
  }

  typename R::Non_zero_dimension_3 non_zero_dimension_3;
  int dir = non_zero_dimension_3(plane.orthogonal_vector());


  CGAL_assertion(!plane.is_degenerate());
  Point_2 (*t)(const Point_3&);

  if(dir == 0){
    t = &point_3_get_y_z_point_2< Point_2, Point_3, R>;
  }else if(dir == 1){
    t = &point_3_get_z_x_point_2< Point_2, Point_3, R>;
  }else{
    t = &point_3_get_x_y_point_2< Point_2, Point_3, R>;
  }

  std::vector< Point_2> points;
  CGAL_NEF_TRACEN("facet:");
  for( ; first != last; ++first ) {
    CGAL_NEF_TRACEN(t(*first)<<" "<<*first);
    points.push_back( t(*first));
  }
  Bounded_side side = bounded_side_2( points.begin(), points.end(), t(point));
  return side;
}

} //namespace CGAL


#endif // CGAL_BOUNDED_SIDE_3_H
