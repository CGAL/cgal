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


#ifndef CGAL__TEST_CLS_TRIANGLE_NEW_2_H
#define CGAL__TEST_CLS_TRIANGLE_NEW_2_H

template <class R>
bool
_test_cls_triangle_new_2(const R& )
{
 std::cout << "Testing class Triangle_2";

 typedef typename  R::RT    RT;
 typedef typename  R::FT    FT;

 typedef typename R::Point_2 Point_2;
 typedef typename R::Triangle_2 Triangle_2;

 typename R::Construct_point_2 construct_point;
 typename R::Triangle_2 it;
 Triangle_2 t0(it); // af:  CGAL::Triangle_2<R> t0(it);

 RT n0 =  0;
 RT n1 =  1;
 RT n2 =  2;
 RT n3 =  3;
 RT n4 =  4;
 RT n5 =  5;
 RT n6 =  6;
 RT n8 =  8;
 RT n9 =  9;
 RT n10= 10;
 RT n12= 12;
 RT n21= 21;

 Point_2 p1 = construct_point( n6, n6, n6);    // ( 1, 1)
 Point_2 p2 = construct_point( n6, n9, n3);    // ( 2, 3)
 Point_2 p3 = construct_point( n6, n10, n2);   // ( 3, 5)
 Point_2 p4 = construct_point( n5, n4, n1);    // ( 5, 4)
 Point_2 p5 = construct_point( n21, n9, n3);   // ( 7, 3)
 Point_2 p6 = construct_point( n8, n4, n2);    // ( 4, 2)
 Point_2 p7 = construct_point( n4, n0);        // ( 4, 0)
 Point_2 p8 = construct_point(-n12,-n8,-n2);   // ( 6, 4)
 Point_2 p9 = construct_point( n9, n9, n3);    // ( 3, 3)

 Triangle_2 t1( p1, p3, p5);
 Triangle_2 t2( p3, p1, p5);
 Triangle_2 t3( p7, p8, p9);
 Triangle_2 t4( p3, p5, p1);
 Triangle_2 t5( p5, p1, p3);
 t0 = t3;

 std::cout << '.';

 assert( t1 == t1 );
 assert( t4 == t1 );
 assert( t1 == t5 );
 assert( t0 == t3 );
 assert( t1 != t2 );
 assert( t3 != t2 );

 assert( t3.vertex(0) == p7 );
 assert( t3.vertex(1) == p8 );
 assert( t3.vertex(2) == p9 );
 assert( t3.vertex(3) == p7 );
 assert( t3.vertex(4) == p8 );
 assert( t3.vertex(5) == p9 );
 assert( t2[5] == t2.vertex(5) );
 assert( t2[6] == t2.vertex(6) );

 assert( t1.orientation() == CGAL::NEGATIVE );
 assert( t2.orientation() == CGAL::POSITIVE );
 assert( t0.orientation() == CGAL::POSITIVE );

 std::cout << '.';

 assert( t1.oriented_side( p9 ) == CGAL::ON_NEGATIVE_SIDE );
 assert( t1.oriented_side( p7 ) == CGAL::ON_POSITIVE_SIDE );
 assert( t1.oriented_side( p8 ) == CGAL::ON_POSITIVE_SIDE );
 assert( t1.oriented_side( p6 ) == CGAL::ON_ORIENTED_BOUNDARY );
 assert( t2.oriented_side( p8 ) == CGAL::ON_NEGATIVE_SIDE );
 assert( t2.oriented_side( p9 ) == CGAL::ON_POSITIVE_SIDE );
 assert( t2.oriented_side( p6 ) == CGAL::ON_ORIENTED_BOUNDARY );
 assert( t2.oriented_side( p3 ) == CGAL::ON_ORIENTED_BOUNDARY );

 assert( t1.bounded_side( p9 ) == CGAL::ON_BOUNDED_SIDE );
 assert( t1.bounded_side( p7 ) == CGAL::ON_UNBOUNDED_SIDE );
 assert( t1.bounded_side( p2 ) == CGAL::ON_BOUNDARY );
 assert( t2.bounded_side( p9 ) == CGAL::ON_BOUNDED_SIDE );
 assert( t2.bounded_side( p7 ) == CGAL::ON_UNBOUNDED_SIDE );
 assert( t2.bounded_side( p2 ) == CGAL::ON_BOUNDARY );
 assert( t2.bounded_side( p1 ) == CGAL::ON_BOUNDARY );
 assert( t2.bounded_side( p5 ) == CGAL::ON_BOUNDARY );

 assert( t1.opposite().has_on_positive_side( p9 ) );
 assert( t1.has_on_positive_side( p8 ) );
 assert( t3.has_on_negative_side( p2 ) );
 assert( t2.has_on_boundary( p1 ) );
 assert( t2.has_on_boundary( p2 ) );
 assert( t2.has_on_boundary( p3 ) );
 assert( t2.has_on_boundary( p4 ) );
 assert( t2.has_on_boundary( p5 ) );
 assert( t2.has_on_boundary( p6 ) );
 assert( t1.has_on_bounded_side( construct_point( n6, n8, n2)) );
 assert( t1.has_on_unbounded_side( construct_point( -n4, n8, n6)) );

 std::cout << '.';

 assert( t1.opposite() == t2 );
 assert( t3 == t3.opposite().opposite() );

 Triangle_2 tdeg1( p1, p7, p7);
 Triangle_2 tdeg2( p6, p6, p6);
 assert( tdeg1.orientation() == CGAL::ZERO );
 assert( tdeg2.orientation() == CGAL::ZERO );
 assert( tdeg1.is_degenerate() );
 assert( tdeg2.is_degenerate() );

 std::cout << '.';

 assert( tdeg1.area() == FT(0) );
 assert( tdeg2.area() == FT(0) );
 assert( t1.area() == FT(-10) );
 assert( t1.area() == -t2.area() );
 assert( t1.area() == t4.area() );
 assert( t2.area() == -t5.area() );
 assert( t3.area() == FT(5) );

 std::cout << "done" << std::endl;
 return true;
}
#endif // CGAL__TEST_CLS_TRIANGLE_NEW_2_H
