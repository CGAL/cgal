// ============================================================================
//
// Copyright (c) 1997-2002 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: $
// release_date  : $CGAL_Date: $
//
// file          : include/CGAL/Nef_3/bounded_side_3.h
// package       : Nef_3
// chapter       : 3D-Nef Polyhedra
//
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Michael Seel    <seel@mpi-sb.mpg.de>
//                 Miguel Granados <granados@mpi-sb.mpg.de>
//                 Susan Hert      <hert@mpi-sb.mpg.de>
//                 Lutz Kettner    <kettner@mpi-sb.mpg.de>
// maintainer    : Susan Hert      <hert@mpi-sb.mpg.de>
//                 Lutz Kettner    <kettner@mpi-sb.mpg.de>
// coordinator   : MPI Saarbruecken
//
// 
// ============================================================================
#ifndef CGAL_BOUNDED_SIDE_3_H
#define CGAL_BOUNDED_SIDE_3_H

#include <CGAL/basic.h>
#include <CGAL/Polygon_2_algorithms.h>

CGAL_BEGIN_NAMESPACE

template <class Point_2, class Point_3> 
Point_2 point_3_get_x_y_point_2(Point_3 p) {
  return( Point_2(p.hx(), p.hy(), p.hw()) );
}

template <class Point_2, class Point_3> 
Point_2 point_3_get_x_z_point_2(Point_3 p) {
  return( Point_2(p.hx(), p.hz(), p.hw()) );
}

template <class Point_2, class Point_3> 
Point_2 point_3_get_y_z_point_2(Point_3 p) {
  return( Point_2(p.hy(), p.hz(), p.hw()) );
}

template <class ForwardIterator, class R>
Bounded_side bounded_side_3(ForwardIterator first,
			    ForwardIterator last,
			    const Point_3<R>& point,
			    Plane_3<R> plane = Plane_3<R>()) {
  if(plane == Plane_3<R>()) {
    CGAL_assertion(last-first > 3);
    /* we need at least 3 points do discover the original plane */
    ForwardIterator p(first);
    Point_3<R> p0(*(p++)), p1(*(p++)), p2(*(p++));
    plane = Plane_3<R>(p0, p1, p2);
    /* since we just need to project points to a non-perpendicular plane
       we don't need to care about the plane orientation */
  }
  Point_2<R> (*t)(Point_3<R>);
  Direction_3<R> pd(plane.orthogonal_vector()), pyz(1,0,0), pxz(0,1,0);
  if(pd == pyz || pd == -pyz)
    /* the plane is parallel to the YZ plane */
    t = &point_3_get_y_z_point_2<Point_2<R>, Point_3<R> >;
  else if(pd == pxz || pd ==- pxz)
    /* the plane is parallel to the XZ plane */
    t = &point_3_get_x_z_point_2<Point_2<R>, Point_3<R> >;
  else {
    CGAL_assertion(cross_product(pd.vector(),Vector_3<R>(0,0,1))==NULL_VECTOR);
    /* the plane is not perpendicular to the XY plane */
    t = &point_3_get_x_y_point_2<Point_2<R>, Point_3<R> >;
  }

  // TODO: implement an interator projector instead of make a copy
  std::vector<Point_2<R> > points;
  points.reserve(last-first);
  CGAL_For_all(first, last)
    points.push_back(t(*first));
  return bounded_side_2(points.begin(), points.end(), t(point));

  /*
  CGAL_For_all(first, last) {
    if(pd == pyz || pd == -pyz)
      points.push_back
	(point_3_get_y_z_point_2<Point_2<R>, Point_3<R> >(*first));
    else if(pd == pxz || pd ==- pxz)
      points.push_back
	(point_3_get_x_z_point_2<Point_2<R>, Point_3<R> >(*first));
     else {
       CGAL_assertion
	(cross_product(Vector_3<R>(0,0,1), pd.vector()) == NULL_VECTOR);
       points.push_back
	 (point_3_get_x_y_point_2<Point_2<R>, Point_3<R> >(*first));
     }
  }
  Point_2<R> p;
  if(pd == pyz || pd == -pyz)
    p = point_3_get_y_z_point_2<Point_2<R>, Point_3<R> >(point);
  else if(pd == pxz || pd ==- pxz)
    p = point_3_get_x_z_point_2<Point_2<R>, Point_3<R> >(point);
  else {
    CGAL_assertion
      (cross_product(Vector_3<R>(0,0,1), pd.vector()) == NULL_VECTOR);
    p = point_3_get_x_y_point_2<Point_2<R>, Point_3<R> >(point);
  }
  return bounded_side_2(points.begin(), points.end(), p);
  */
}

CGAL_END_NAMESPACE

#endif // CGAL_BOUNDED_SIDE_3_H
