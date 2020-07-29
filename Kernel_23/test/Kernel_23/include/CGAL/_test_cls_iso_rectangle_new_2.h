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


#ifndef CGAL__TEST_CLS_ISO_RECTANGLE_NEW_2_H
#define CGAL__TEST_CLS_ISO_RECTANGLE_NEW_2_H

template <class R>
bool
_test_cls_iso_rectangle_new_2(const R& )
{
 std::cout << "Testing class Iso_rectangle_2";

 typedef typename  R::RT    RT;
 typedef typename  R::FT    FT;

 typedef typename  R::Point_2 Point_2;
 typedef typename  R::Iso_rectangle_2 Iso_rectangle_2;

 typename R::Construct_point_2 construct_point;
 typename R::Iso_rectangle_2 ir;
 Iso_rectangle_2 r0(ir); // CGAL::Iso_rectangle_2<R>  r0(ir);

 RT n1 =  1;
 RT n2 =  2;
 RT n3 =  3;
 RT n4 =  4;
 RT n5 =  5;
 RT n6 =  6;
 RT n7 =  7;
 RT n8 =  8;
 RT n12= 12;

 Point_2 p1 = construct_point( n5, n5, n5);    // ( 1, 1)
 Point_2 p2 = construct_point( n2, n8, n2);    // ( 1, 4)
 Point_2 p3 = construct_point( n12, n12, n3);  // ( 4, 4)
 Point_2 p4 = construct_point( n5, n4, n1);    // ( 5, 4)
 Point_2 p5 = construct_point( n4, n1);        // ( 4, 1)
 Point_2 p6 = construct_point( n8, n4, n2);    // ( 4, 2)
 Point_2 p7 = construct_point( n6, n3, n2);    // ( 3, 1.5)
 Point_2 p8 = construct_point( n4, n6, n2);    // ( 2, 3)
 Point_2 p9 = construct_point(-n3, n7);        // (-3, 7)

 Iso_rectangle_2 r1( p1, p3);
 Iso_rectangle_2 r2( p3, p1);
 Iso_rectangle_2 r3( p2, p5);
 Iso_rectangle_2 r4( p5, p2);
 Iso_rectangle_2 r5( p9, p2);
 Iso_rectangle_2 r6( r2 );
 Iso_rectangle_2 r7( n3, n3, n12, n12, n3);
 Iso_rectangle_2 r8( n4, n2, n6, n8, n2);
 Iso_rectangle_2 r9( n2, n1, n3, n4, n1);
 Iso_rectangle_2 r10( n2, n1, n3, n4);

 Iso_rectangle_2 r11(p1, p1, p1, p1);
 Iso_rectangle_2 r12(p1, p1, p2, p2);
 Iso_rectangle_2 r13(p1, p2, p1, p2);
 Iso_rectangle_2 r14(p3, p4, p1, p2);

 r0 = r1;

 assert( r1 == r1 );
 assert( r0 == r1 );
 assert( r1 == r2 );
 assert( r1 == r3 );
 assert( r1 == r4 );
 assert( r2 == r6 );
 assert( r2 != r5 );
 assert( r7 == r1 );
 assert( r8 == r9 );
 assert( r9 == r10 );

 std::cout << '.';

 assert( r1.vertex(0) == p1 );
 assert( r1.vertex(1) == p5 );
 assert( r1.vertex(2) == p3 );
 assert( r1.vertex(3) == p2 );
 assert( r1.vertex(4) == p1 );
 assert( r3.vertex(0) == p1 );
 assert( r3.vertex(1) == p5 );
 assert( r3.vertex(2) == p3 );
 assert( r3.vertex(3) == p2 );
 assert( r3.vertex(4) == p1 );
 assert( r2[0] == r2.vertex(0) );
 assert( r2[1] == r2.vertex(1) );
 assert( r2[3] == r2.vertex(3) );
 assert( r2[4] == r2.vertex(0) );

 std::cout << '.';

 assert( (r4.min)() == p1 );
 assert( (r1.min)() == p1 );
 assert( (r5.min)() != p9 );
 assert( (r2.max)() == p3 );

 std::cout << '.';

 assert( r1.min_coord(0) == r1.xmin() );
 assert( r1.min_coord(1) == r1.ymin() );
 assert( r2.max_coord(0) == r2.xmax() );
 assert( r2.max_coord(1) == r2.ymax() );

 std::cout << '.';

 assert( r1.bounded_side( p8 ) == CGAL::ON_BOUNDED_SIDE );
 assert( r2.bounded_side( p7 ) == CGAL::ON_BOUNDED_SIDE );
 assert( r3.bounded_side( p9 ) == CGAL::ON_UNBOUNDED_SIDE );
 assert( r1.bounded_side( p4 ) == CGAL::ON_UNBOUNDED_SIDE );
 assert( r4.bounded_side( p6 ) == CGAL::ON_BOUNDARY );
 assert( r4.bounded_side( p1 ) == CGAL::ON_BOUNDARY );

 assert( r5.has_on_boundary( p2 ) );
 assert( r4.has_on_boundary( p2 ) );
 assert( r2.has_on_bounded_side( p7 ) );
 assert( r4.has_on_unbounded_side( p9 ) );

 std::cout << '.';

 assert( Iso_rectangle_2( p1, p1 ).is_degenerate() );
 assert( Iso_rectangle_2( p1, p2 ).is_degenerate() );
 assert( Iso_rectangle_2( p3, p4 ).is_degenerate() );

 std::cout << '.';

 assert( Iso_rectangle_2( p1, p1 ).area() == FT(0) );
 assert( Iso_rectangle_2( p1, p2 ).area() == FT(0) );
 assert( Iso_rectangle_2( p3, p4 ).area() == FT(0) );
 assert( Iso_rectangle_2( p1, p3 ).area() == FT(9) );
 assert( Iso_rectangle_2( p3, p1 ).area() == FT(9) );
 assert( Iso_rectangle_2( p1, p7 ).area() == FT(1) );
 assert( Iso_rectangle_2( p9, p3 ).area() == FT(21) );

 std::cout << "done" << std::endl;
 return true;
}
#endif // CGAL__TEST_CLS_ISO_RECTANGLE_NEW_2_H
