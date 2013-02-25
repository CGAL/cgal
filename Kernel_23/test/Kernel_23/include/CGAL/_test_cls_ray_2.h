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
// Author(s)     : Stefan Schirra
 

#ifndef CGAL__TEST_CLS_RAY_2_H
#define CGAL__TEST_CLS_RAY_2_H

template <class R>
bool
_test_cls_ray_2(const R& )
{
 std::cout << "Testing class Ray_2";

 typedef typename  R::RT    RT;

 typename R::Ray_2  ir;
 CGAL::Ray_2<R> r0;

 RT  n2 = 2;
 RT  n3 = 3;
 RT  n4 = 4;
 RT  n5 = 5;
 RT  n8 = 8;
 RT  n9 = 9;
 RT n10 =10;

 CGAL::Point_2<R> p1( n2, n8, n2);
 CGAL::Point_2<R> p2( n10, n4, n2);
 CGAL::Point_2<R> p3( n9, n9, n3);
 CGAL::Point_2<R> p4( n10, n8, n2);
 CGAL::Direction_2<R> d12( p2 - p1);
 CGAL::Direction_2<R> d24( p4 - p2);
 CGAL::Line_2<R> l12( p1, p2);
 CGAL::Line_2<R> l24( p2, p4);
 CGAL::Vector_2<R> v12( p2 - p1);
 CGAL::Vector_2<R> v24( p4 - p2);

 CGAL::Ray_2<R> r1( p1, p2);
 CGAL::Ray_2<R> r2( p1, d12);
 CGAL::Ray_2<R> r3( p3, p2);
 CGAL::Ray_2<R> r4( p2, d24);
 CGAL::Ray_2<R> r5( p2, p4);
 CGAL::Ray_2<R> r6( p1, p4);
 CGAL::Ray_2<R> r7( p1, l12);
 CGAL::Ray_2<R> r8( p2, l24);
 CGAL::Ray_2<R> r7v( p1, v12);
 CGAL::Ray_2<R> r8v( p2, v24);
 r0 = r3;

 assert(   CGAL::parallel(r1, r2) );
 assert(   CGAL::parallel(r1, r3) );
 assert(   CGAL::parallel(r4, r5) );
 assert( ! CGAL::parallel(r1, r6) );

 std::cout << '.';

 assert( r1 == r1 );
 assert( r2 == r1 );
 assert( r4 == r5 );
 assert( r0 == r3 );
 assert( r6 != r1 );
 assert( r1 != r3 );
 assert( r1 != r5 );

 assert( r2.source() == p1 );
 assert( r0.source() == p3 );
 assert( r4.source() == r4.point(0) );

 assert( r1.direction() == d12 );
 assert( r2.direction() == d12 );
 assert( r3.direction() == r1.direction() );
 assert( r3.to_vector().direction() == r1.to_vector().direction() );

 std::cout << '.';

 assert( r3.supporting_line() == r1.supporting_line() );
 assert( r5.supporting_line() == CGAL::Line_2<R>( p2, p4 ) );

 assert( r4.opposite() == CGAL::Ray_2<R>( p2, -d24 ) );
 assert( r1.opposite() == CGAL::Ray_2<R>( p1, -d12 ) );
 assert( r2.opposite().opposite() == r2 );
 assert( r2 == r7 );
 assert( r4 == r8 );
 assert( r2 == r7v );
 assert( r4 == r8v );

 assert( r6.is_horizontal() );
 assert( ! r0.is_horizontal() );
 assert( r5.is_vertical() );
 assert( ! r5.is_horizontal() );

 assert( r1.has_on( p1 ) );
 assert( r1.has_on( p2 ) );
 assert( r1.has_on( p3 ) );
 assert( r3.opposite().has_on( p1 ) );
 assert( ! r1.has_on( p4 ) );
 assert( ! r0.has_on( CGAL::Point_2<R>( n8, n5, n8 )) );
 assert( r4.has_on( r4.point(7)) );
 assert( r3.collinear_has_on( r3.point(7)) );
 assert( r1.collinear_has_on( p3) );
 assert( ! r3.collinear_has_on( p1 ) );

 std::cout << '.';

 assert( CGAL::Ray_2<R>( p1, p1).is_degenerate() );
 assert( ! r0.is_degenerate() );
 assert( CGAL::parallel(r0, r0) );

 std::cout << "done" << std::endl;
 return true;
}
#endif // CGAL__TEST_CLS_RAY_2_H
