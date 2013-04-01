// Copyright (c) 1999  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Andreas Fabri
 

#ifndef CGAL__TEST_FURTHER_FCT_POINT_3_PLANE_3_H
#define CGAL__TEST_FURTHER_FCT_POINT_3_PLANE_3_H

template <class R>
bool
_test_further_fct_point_plane_3(const R& )
{
 std::cout << "Testing further functions Point_3 Plane_3" ;

 typedef typename  R::RT    RT;

 CGAL::Point_3<R> p0( RT(0), RT(0), RT(0) );
 CGAL::Point_3<R> px( RT(1), RT(0), RT(0) );
 CGAL::Point_3<R> py( RT(0), RT(1), RT(0) );
 CGAL::Point_3<R> pz( RT(0), RT(0), RT(1) );

 CGAL::Point_3<R> q1( RT(-1), RT(-1), RT(-1) );
 CGAL::Point_3<R> q2( RT(-4), RT(-3), RT(-4) );
 CGAL::Point_3<R> q3( RT(-12), RT(13), RT(-4) );

 CGAL::Plane_3<R> hxzy(px, pz, py);
 CGAL::Plane_3<R> hxyz(px,py, pz);

 CGAL::Plane_3<R> hxy(p0, px, py);
 CGAL::Plane_3<R> hyx(p0, py, px);

 std::cout << '.';

 assert( CGAL::compare_signed_distance_to_plane(hxyz, px, py) == CGAL::EQUAL );
 assert( CGAL::compare_signed_distance_to_plane(px, py, pz, px, py) == CGAL::EQUAL );

 assert( CGAL::compare_signed_distance_to_plane(hxzy, px,py) == CGAL::EQUAL );
 assert( CGAL::compare_signed_distance_to_plane(px, pz, py, px, py) == CGAL::EQUAL );


 assert( CGAL::compare_signed_distance_to_plane(hxzy, q1,q2) == CGAL::SMALLER );
 assert( CGAL::compare_signed_distance_to_plane(px, pz, py, q1, q2) == CGAL::SMALLER );

 assert( CGAL::compare_signed_distance_to_plane(hxzy, q1,q2) == CGAL::SMALLER );
 assert( CGAL::compare_signed_distance_to_plane(px, pz, py, q1, q2) == CGAL::SMALLER );

 assert( CGAL::compare_signed_distance_to_plane(hxyz, q1,q2) == CGAL::LARGER );
 assert( CGAL::compare_signed_distance_to_plane(px, py, pz, q1, q2) == CGAL::LARGER );
 
 assert( CGAL::compare_signed_distance_to_plane(hxy, q2, q3) == CGAL::EQUAL );
 assert( CGAL::compare_signed_distance_to_plane(p0, px, py, q2, q3) == CGAL::EQUAL );

 assert( CGAL::compare_signed_distance_to_plane(hyx, q2, q3) == CGAL::EQUAL );
 assert( CGAL::compare_signed_distance_to_plane(p0, py, px, q2, q3) == CGAL::EQUAL );

 std::cout << "done" << std::endl;
 return true;
}
#endif // CGAL__TEST_FURTHER_FCT_POINT_3_PLANE_3_H
