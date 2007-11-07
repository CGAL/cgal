// Copyright (c) 1999  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Stefan Schirra
 

#ifndef CGAL__TEST_CLS_TRIANGLE_2_H
#define CGAL__TEST_CLS_TRIANGLE_2_H

#include <CGAL/Bbox_2.h>
#include <CGAL/Testsuite/assert.h>

template <class R>
bool
_test_cls_triangle_2(const R& )
{
 std::cout << "Testing class Triangle_2";

 typedef typename  R::RT    RT;
 typedef typename  R::FT    FT;

 typename R::Triangle_2 it;
 CGAL::Triangle_2<R> t0(it);

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

 CGAL::Point_2<R> p1( n6, n6, n6);    // ( 1, 1)
 CGAL::Point_2<R> p2( n6, n9, n3);    // ( 2, 3)
 CGAL::Point_2<R> p3( n6, n10, n2);   // ( 3, 5)
 CGAL::Point_2<R> p4( n5, n4, n1);    // ( 5, 4)
 CGAL::Point_2<R> p5( n21, n9, n3);   // ( 7, 3)
 CGAL::Point_2<R> p6( n8, n4, n2);    // ( 4, 2)
 CGAL::Point_2<R> p7( n4, n0);        // ( 4, 0)
 CGAL::Point_2<R> p8(-n12,-n8,-n2);   // ( 6, 4)
 CGAL::Point_2<R> p9( n9, n9, n3);    // ( 3, 3)

 CGAL::Triangle_2<R> t1( p1, p3, p5);
 CGAL::Triangle_2<R> t2( p3, p1, p5);
 CGAL::Triangle_2<R> t3( p7, p8, p9);
 CGAL::Triangle_2<R> t4( p3, p5, p1);
 CGAL::Triangle_2<R> t5( p5, p1, p3);
 t0 = t3;

 std::cout << '.';

 CGAL_test_assert( t1 == t1 );
 CGAL_test_assert( t4 == t1 );
 CGAL_test_assert( t1 == t5 );
 CGAL_test_assert( t0 == t3 );
 CGAL_test_assert( t1 != t2 );
 CGAL_test_assert( t3 != t2 );

 CGAL_test_assert( t3.vertex(0) == p7 );
 CGAL_test_assert( t3.vertex(1) == p8 );
 CGAL_test_assert( t3.vertex(2) == p9 );
 CGAL_test_assert( t3.vertex(3) == p7 );
 CGAL_test_assert( t3.vertex(4) == p8 );
 CGAL_test_assert( t3.vertex(5) == p9 );
 CGAL_test_assert( t2[5] == t2.vertex(5) );
 CGAL_test_assert( t2[6] == t2.vertex(6) );

 CGAL_test_assert( t1.orientation() == CGAL::NEGATIVE );
 CGAL_test_assert( t2.orientation() == CGAL::POSITIVE );
 CGAL_test_assert( t0.orientation() == CGAL::POSITIVE );

 std::cout << '.';

 CGAL_test_assert( t1.oriented_side( p9 ) == CGAL::ON_NEGATIVE_SIDE );
 CGAL_test_assert( t1.oriented_side( p7 ) == CGAL::ON_POSITIVE_SIDE );
 CGAL_test_assert( t1.oriented_side( p8 ) == CGAL::ON_POSITIVE_SIDE );
 CGAL_test_assert( t1.oriented_side( p6 ) == CGAL::ON_ORIENTED_BOUNDARY );
 CGAL_test_assert( t2.oriented_side( p8 ) == CGAL::ON_NEGATIVE_SIDE );
 CGAL_test_assert( t2.oriented_side( p9 ) == CGAL::ON_POSITIVE_SIDE );
 CGAL_test_assert( t2.oriented_side( p6 ) == CGAL::ON_ORIENTED_BOUNDARY );
 CGAL_test_assert( t2.oriented_side( p3 ) == CGAL::ON_ORIENTED_BOUNDARY );

 CGAL_test_assert( t1.bounded_side( p9 ) == CGAL::ON_BOUNDED_SIDE );
 CGAL_test_assert( t1.bounded_side( p7 ) == CGAL::ON_UNBOUNDED_SIDE );
 CGAL_test_assert( t1.bounded_side( p2 ) == CGAL::ON_BOUNDARY );
 CGAL_test_assert( t2.bounded_side( p9 ) == CGAL::ON_BOUNDED_SIDE );
 CGAL_test_assert( t2.bounded_side( p7 ) == CGAL::ON_UNBOUNDED_SIDE );
 CGAL_test_assert( t2.bounded_side( p2 ) == CGAL::ON_BOUNDARY );
 CGAL_test_assert( t2.bounded_side( p1 ) == CGAL::ON_BOUNDARY );
 CGAL_test_assert( t2.bounded_side( p5 ) == CGAL::ON_BOUNDARY );

 CGAL_test_assert( t1.opposite().has_on_positive_side( p9 ) );
 CGAL_test_assert( t1.has_on_positive_side( p8 ) );
 CGAL_test_assert( t3.has_on_negative_side( p2 ) );
 CGAL_test_assert( t2.has_on_boundary( p1 ) );
 CGAL_test_assert( t2.has_on_boundary( p2 ) );
 CGAL_test_assert( t2.has_on_boundary( p3 ) );
 CGAL_test_assert( t2.has_on_boundary( p4 ) );
 CGAL_test_assert( t2.has_on_boundary( p5 ) );
 CGAL_test_assert( t2.has_on_boundary( p6 ) );
 CGAL_test_assert( t1.has_on_bounded_side( CGAL::Point_2<R>( n6, n8, n2)) );
 CGAL_test_assert( t1.has_on_unbounded_side( CGAL::Point_2<R>( -n4, n8, n6)) );

 std::cout << '.';

 CGAL_test_assert( t1.opposite() == t2 );
 CGAL_test_assert( t3 == t3.opposite().opposite() );

 CGAL::Triangle_2<R> tdeg1( p1, p7, p7);
 CGAL::Triangle_2<R> tdeg2( p6, p6, p6);
 CGAL_test_assert( tdeg1.orientation() == CGAL::ZERO );
 CGAL_test_assert( tdeg2.orientation() == CGAL::ZERO );
 CGAL_test_assert( tdeg1.is_degenerate() );
 CGAL_test_assert( tdeg2.is_degenerate() );

 std::cout << '.';

 CGAL_test_assert( tdeg1.area() == FT(0) );
 CGAL_test_assert( tdeg2.area() == FT(0) );
 CGAL_test_assert( t1.area() == FT(-10) );
 CGAL_test_assert( t1.area() == -t2.area() );
 CGAL_test_assert( t1.area() == t4.area() );
 CGAL_test_assert( t2.area() == -t5.area() );
 CGAL_test_assert( t3.area() == FT(5) );

 std::cout << '.';

 CGAL::Bbox_2 bb = t1.bbox();
 CGAL_test_assert(bb.xmin() <= 1.0);
 CGAL_test_assert(bb.xmax() >= 7.0);
 CGAL_test_assert(bb.ymin() <= 1.0);
 CGAL_test_assert(bb.ymax() >= 5.0);

 std::cout << "done" << std::endl;
 return true;
}
#endif // CGAL__TEST_CLS_TRIANGLE_2_H
