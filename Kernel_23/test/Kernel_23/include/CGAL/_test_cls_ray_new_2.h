// Copyright (c) 1999
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later
//
//
// Author(s)     : Stefan Schirra


#ifndef CGAL__TEST_CLS_RAY_NEW_2_H
#define CGAL__TEST_CLS_RAY_NEW_2_H

template <class R>
bool
_test_cls_ray_new_2(const R& )
{
 std::cout << "Testing class Ray_2";

 typedef typename  R::RT    RT;
 typedef typename  R::FT    FT;

 typedef typename  R::Point_2 Point_2;
 typedef typename  R::Vector_2 Vector_2;
 typedef typename  R::Direction_2 Direction_2;

 typedef typename  R::Line_2 Line_2;
 typedef typename  R::Ray_2 Ray_2;

 typename R::Construct_point_2 construct_point;
 typename R::Construct_vector_2 construct_vector;

 Ray_2  ir;
 Ray_2 r0; // af was: CGAL::Ray_2<R> r0;

 RT  n2 = 2;
 RT  n3 = 3;
 RT  n4 = 4;
 RT  n5 = 5;
 RT  n8 = 8;
 RT  n9 = 9;
 RT n10 =10;

 Point_2 p1 = construct_point( n2, n8, n2);
 Point_2 p2 = construct_point( n10, n4, n2);
 Point_2 p3 = construct_point( n9, n9, n3);
 Point_2 p4 = construct_point( n10, n8, n2);
 Direction_2 d12(construct_vector( p1, p2));
 Direction_2 d24(construct_vector( p2, p4));
 Line_2 l12( p1, p2);
 Line_2 l24( p2, p4);
 Vector_2 v12(construct_vector( p1, p2));
 Vector_2 v24(construct_vector( p2, p4));

 Ray_2 r1( p1, p2);
 Ray_2 r2( p1, d12);
 Ray_2 r3( p3, p2);
 Ray_2 r4( p2, d24);
 Ray_2 r5( p2, p4);
 Ray_2 r6( p1, p4);
 Ray_2 r7( p1, l12);
 Ray_2 r8( p2, l24);
 Ray_2 r7v( p1, v12);
 Ray_2 r8v( p2, v24);
 r0 = r3;

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
 assert( r5.supporting_line() == Line_2( p2, p4 ) );

 assert( r4.opposite() == Ray_2( p2, -d24 ) );
 assert( r1.opposite() == Ray_2( p1, -d12 ) );
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
 assert( ! r0.has_on( construct_point( n8, n5, n8 )) );
 assert( r4.has_on( r4.point(7)) );
 assert( r3.collinear_has_on( r3.point(7)) );
 assert( r3.collinear_has_on( r3.point(FT(42)/FT(13))) );
 assert( r1.collinear_has_on( p3) );
 assert( ! r3.collinear_has_on( p1 ) );

 std::cout << '.';

 assert( Ray_2( p1, p1).is_degenerate() );
 assert( ! r0.is_degenerate() );

 std::cout << "done" << std::endl;
 return true;
}
#endif // CGAL__TEST_CLS_RAY_2_H
