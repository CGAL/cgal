// Copyright (c) 1997-2002  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
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
#include <vector>

#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 17
#include <CGAL/Nef_2/debug.h>

namespace CGAL {

template <class Point_2, class Point_3> 
Point_2 point_3_get_x_y_point_2(Point_3 p) {
  return( Point_2(p.hx(), p.hy(), p.hw()) );
}

template <class Point_2, class Point_3> 
Point_2 point_3_get_y_z_point_2(Point_3 p) {
  return( Point_2(p.hy(), p.hz(), p.hw()) );
}

template <class Point_2, class Point_3> 
Point_2 point_3_get_z_x_point_2(Point_3 p) {
  return( Point_2(p.hz(), p.hx(), p.hw()) );
}

template <class IteratorForward, class R>
Bounded_side bounded_side_3(IteratorForward first,
			    IteratorForward last,
                            const Point_3<R>& point,
                            typename R::Plane_3 plane = typename R::Plane_3(0,0,0,0)) {
  typedef typename R::Point_2 Point_2;
  typedef typename R::Point_3 Point_3;
  typedef typename R::Vector_3 Vector_3;
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



  CGAL_assertion(!plane.is_degenerate());
  Point_2 (*t)(Point_3);
  Vector_3 pv(plane.orthogonal_vector()), pxy(0,0,1), pyz(1,0,0), pzx(0,1,0);
  CGAL_NEF_TRACEN("pv*pxz: "<<pv*pzx);
  CGAL_NEF_TRACEN("pv*pyz: "<<pv*pyz);
  CGAL_NEF_TRACEN("pv*pxy: "<<pv*pxy);
  if( !CGAL_NTS is_zero(pv*pzx) )
    /* the plane is not perpendicular to the ZX plane */
    t = &point_3_get_z_x_point_2< Point_2, Point_3>;
  else if( !CGAL_NTS is_zero(pv*pyz) )
    /* the plane is not perpendicular to the YZ plane */
    t = &point_3_get_y_z_point_2< Point_2, Point_3>;
  else {
    CGAL_assertion( !CGAL_NTS is_zero(pv*pxy) );
    /* the plane is not perpendicular to the XY plane */
    t = &point_3_get_x_y_point_2< Point_2, Point_3>;
  }

  std::vector< Point_2> points;
  CGAL_NEF_TRACEN("facet:");
  for( ; first != last; ++first ) {
    CGAL_NEF_TRACEN(t(*first)<<" "<<*first);
    points.push_back( t(*first));
  }
  Bounded_side side = bounded_side_2( points.begin(), points.end(), t(point));
  points.clear();
  return side;  
}

} //namespace CGAL

#ifdef WRONG_IMPLEMENTATION
/* The following code is wrong since Proyector_.. structures must not return
   references to temporal objects */
template < class Point_2, class Point_3> 
struct Project_XY {
  typedef Point_3                  argument_type;
  typedef Point_2                  result_type;
  Point_2 operator()( Point_3& p) const { 
    return Point_2(p.hx(), p.hy(), p.hw());
  }
  const Point_2 operator()( const Point_3& p) const { 
    return Point_2(p.hx(), p.hy(), p.hw());
  }
};

template < class Point_2, class Point_3> 
struct Project_YZ {
  typedef Point_3                  argument_type;
  typedef Point_2                  result_type;
  Point_2 operator()( Point_3& p) const { 
    return Point_2(p.hy(), p.hz(), p.hw());
  }
  const Point_2 operator()( const Point_3& p) const { 
    return Point_2(p.hy(), p.hz(), p.hw());
  }
};

template < class Point_2, class Point_3> 
struct Project_XZ {
  typedef Point_3                  argument_type;
  typedef Point_2                  result_type;
  Point_2 operator()( Point_3& p) const { 
    return Point_2(p.hx(), p.hz(), p.hw());
  }
  const Point_2 operator()( const Point_3& p) const { 
    return Point_2(p.hx(), p.hz(), p.hw());
  }
};

template <class IC, class R>
Bounded_side bounded_side_3(IC first,
			    IC last,
			    const Point_3<R>& point,
			    Plane_3<R> plane = Plane_3<R>(0,0,0,0)) {

  typedef typename R::Point_2 Point_2;
  typedef typename R::Point_3 Point_3;
  typedef typename R::Vector_3 Vector_3;
  typedef typename R::Plane_3 Plane_3;

  CGAL_assertion( !CGAL::is_empty_range( first, last));

  if(plane == Plane_3(0,0,0,0)) {
    Vector_3 hv;
    normal_vector_newell_3( first, last, hv);
    plane = Plane_3( *first, Vector_3(hv));
  }
  CGAL_assertion(!plane.is_degenerate());
  Vector_3 pd(plane.orthogonal_vector()), pyz(1,0,0), pxz(0,1,0);
  if(pd == pyz || pd == -pyz) {
    /* the plane is parallel to the YZ plane */
    typedef Project_YZ< Point_2, Point_3>                  Project_YZ;
    typedef Iterator_project< IC, Project_YZ> Iterator_YZ;
    Project_YZ project;
    Point_2 p = project(point);
    Iterator_YZ pfirst(first), plast(last);
    return bounded_side_2(pfirst, plast, p);
  }
  else if(pd == pxz || pd ==- pxz) {
    /* the plane is parallel to the XZ plane */
    typedef Project_XZ< Point_2, Point_3>                  Project_XZ;
    typedef Iterator_project< IC, Project_XZ> Iterator_XZ;
    Project_XZ project;
    Point_2 p = project(point);
    Iterator_XZ pfirst(first), plast(last);
    return bounded_side_2(pfirst, plast, p);
  }
  else {
    CGAL_assertion(cross_product(pd.vector(), Vector_3(0,0,1)) == NULL_VECTOR);
    /* the plane is not perpendicular to the XY plane */
    typedef Project_XY< Point_2, Point_3>                  Project_XY;
    typedef Iterator_project< IC, Project_XY> Iterator_XY;
    Project_XY project;
    Point_2 p = project(point);
    Iterator_XY pfirst(first), plast(last);
    return bounded_side_2(pfirst, plast, p);
  }
}
#endif // WRONG_IMPLEMENTATION

#endif // CGAL_BOUNDED_SIDE_3_H
