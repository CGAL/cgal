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


#ifndef CGAL__TEST_CLS_CIRCLE_2_H
#define CGAL__TEST_CLS_CIRCLE_2_H

#include <CGAL/Bbox_2.h>
#include <CGAL/Testsuite/assert.h>

template <class R>
bool
_test_cls_circle_2(const R& )
{
 std::cout << "Testing class Circle_2";

 typedef typename  R::RT    RT;
 typedef typename  R::FT    FT;

 typename R::Circle_2  ic;
 CGAL::Circle_2<R> c0;

 RT n0 =  0;
 RT n1 = 16;
 RT n2 = -4;
 RT n3 =  2;
 RT n4 =  5;
 RT n5 = 10;
 RT n6 =  6;
 RT n8 =  9;
 RT n9 = 36;
 RT n10=  3;
 RT n11=-11;
 RT n12=200;
 RT n13= 13;

 CGAL::Point_2<R> p0( n1, n2, -n2);  // ( 4, -1)
 CGAL::Point_2<R> p1( n6, n8, n10);  // ( 2,  3)
 CGAL::Point_2<R> p2( n2, n0,  n2);  // ( 1,  0)
 CGAL::Point_2<R> p3( n5, n5,  n4);  // ( 2,  2)
 CGAL::Point_2<R> p4( n0, n2,  n2);  // ( 0,  1)

 CGAL::Vector_2<R> vx = p2 - CGAL::ORIGIN;
 CGAL::Vector_2<R> vy = p4 - CGAL::ORIGIN;
 CGAL::Vector_2<R> v1 = p1 - CGAL::ORIGIN;

 CGAL::Circle_2<R> c1( p0, p1, p2);
 CGAL::Circle_2<R> c2( p0, p1, p3);
 CGAL::Circle_2<R> c3( p1, p0, p2);
 const CGAL::Circle_2<R> c4( p3, FT( n9 ));      // n9 = (n6)^2
 CGAL::Vector_2<R> vx6 = vx * n6;
 CGAL::Vector_2<R> vy6 = vy * n6;
 CGAL::Circle_2<R> c5( p3 - vx6, p3 + vx6, p3 + vy6);
 CGAL::Circle_2<R> c6( c3 );
 CGAL::Circle_2<R> c7( p3, n9, CGAL::POSITIVE);
 CGAL::Circle_2<R> c8( p3, n9, CGAL::NEGATIVE);
 CGAL::Circle_2<R> cc( p3 - vx6, p3 + vx6);
 CGAL::Circle_2<R> cp( p3 - vx6, p3 + vx6, CGAL::POSITIVE);
 CGAL::Circle_2<R> cn( p3 - vx6, p3 + vx6, CGAL::NEGATIVE);
 c0 = c3;

 CGAL_test_assert( c1 == c1 );
 CGAL_test_assert( c1 != c2 );
 CGAL_test_assert( c3 == c0 );
 CGAL_test_assert( c0 == c3 );
 CGAL_test_assert( c3 == c6 );
 CGAL_test_assert( c7 != c8 );
 CGAL_test_assert( c4 == c7 );
 CGAL_test_assert( cc == cp );
 CGAL_test_assert( cn != cp );
 CGAL_test_assert( cc != c8 );
 CGAL_test_assert( cc == c7 );

 CGAL_test_assert( c5.center() == p3 );
 CGAL_test_assert( cc.center() == p3 );
 CGAL_test_assert( c5.squared_radius() == FT( n9 ) );
 CGAL_test_assert( c4.squared_radius() == cc.squared_radius() );
 CGAL_test_assert( c4 == c5 );
 CGAL_test_assert( c4 == c7 );
 CGAL_test_assert( c4 != c8 );
 CGAL_test_assert( cn == cp.opposite() );
 CGAL_test_assert( c7.opposite() == c8 );
 CGAL_test_assert( c8.opposite() == c7 );
 CGAL_test_assert( c1.opposite() == c3 );
 CGAL_test_assert( c3.opposite() == c1 );
 CGAL_test_assert( c7.orientation() == CGAL::POSITIVE );
 CGAL_test_assert( c8.orientation() == CGAL::NEGATIVE );
 CGAL_test_assert( c5.orientation() == CGAL::POSITIVE );
 CGAL_test_assert( cc.orientation() == CGAL::POSITIVE );
 CGAL_test_assert( cp.orientation() == CGAL::POSITIVE );
 CGAL_test_assert( cn.orientation() == CGAL::NEGATIVE );

 std::cout << '.';

 CGAL_test_assert( c4.center() == p3 );
 CGAL_test_assert( c5.center() == p3 );
 CGAL_test_assert( c4.squared_radius() == FT( n9 ) );
 CGAL_test_assert( c5.squared_radius() == FT( n9 ) );
 CGAL_test_assert( c8.squared_radius() == FT( n9 ) );

 CGAL_test_assert( c7.bounded_side( p3 + vx*n2 ) == CGAL::ON_BOUNDED_SIDE );
 CGAL_test_assert( c7.bounded_side( p3 + vy*n11 ) == CGAL::ON_UNBOUNDED_SIDE );
 CGAL_test_assert( c7.bounded_side( p3 - vy6 ) == CGAL::ON_BOUNDARY );
 CGAL_test_assert( c8.bounded_side( p3 + vx*n2 ) == CGAL::ON_BOUNDED_SIDE );
 CGAL_test_assert( c8.bounded_side( p3 + vy*n11 ) == CGAL::ON_UNBOUNDED_SIDE );
 CGAL_test_assert( c8.bounded_side( p3 - vy6 ) == CGAL::ON_BOUNDARY );
 CGAL_test_assert( cc.has_on_boundary( p3 + vy6) );
 CGAL_test_assert( cc.has_on_boundary( p3 - vx6) );

 std::cout << '.';

 CGAL::Aff_transformation_2<R>
          rotate1(CGAL::ROTATION,CGAL::Direction_2<R>(n11,n13),-n2,n12),
          rotate2(CGAL::ROTATION,CGAL::Direction_2<R>(-n8, n9),-n2,n12),
          rotate3(CGAL::ROTATION,CGAL::Direction_2<R>( n5,-n1),-n2,n12),
          rotate4(CGAL::ROTATION,CGAL::Direction_2<R>(-n5,-n11),-n2,n12);
 CGAL::Point_2<R> p6 = p2.transform( rotate1 );
 CGAL::Point_2<R> p7 = p2.transform( rotate2 );
 CGAL::Point_2<R> p8 = p2.transform( rotate3 );
 CGAL::Point_2<R> p9 = p2.transform( rotate4 );
 p6 = p6 + v1;  // v1 = ( 2, 3 )
 p7 = p7 + v1;
 p8 = p8 + v1;
 p9 = p9 + v1;
 CGAL::Circle_2<R> c10 (p6, p8, p7);
 CGAL_test_assert( c10.orientation() == CGAL::POSITIVE );
 CGAL_test_assert( c10.opposite().orientation() == CGAL::NEGATIVE );

 CGAL_test_assert( c10.oriented_side(c10.center() ) == CGAL::ON_POSITIVE_SIDE );
 CGAL_test_assert( c10.oriented_side(CGAL::ORIGIN + v1 + vx/n2 ) \
         == CGAL::ON_POSITIVE_SIDE );
 CGAL_test_assert( c10.oriented_side(CGAL::ORIGIN + v1 + vx*n2 ) \
         == CGAL::ON_NEGATIVE_SIDE );
 CGAL_test_assert( c10.oriented_side(p9 ) == CGAL::ON_ORIENTED_BOUNDARY );
 CGAL_test_assert( c10.has_on_boundary(p9) );
 CGAL_test_assert( c10.has_on_boundary(p4 + v1) );
 CGAL::Point_2<R> p11( n4, n4, n3) ; // (2.5, 2.5)
 CGAL::Point_2<R> p12( n5, n5, n3) ; // ( 5 ,  5 )
 CGAL_test_assert( c10.has_on_bounded_side( p11 ) );
 CGAL_test_assert( ! c10.has_on_bounded_side( p12 ) );
 CGAL_test_assert( c10.has_on_unbounded_side( p12 ) );
 CGAL_test_assert( c10.has_on_positive_side( p11 ) );
 CGAL_test_assert( c10.has_on_negative_side( p12 ) );
 CGAL_test_assert( c10.opposite().has_on_negative_side( p11 ) );
 CGAL_test_assert( c10.opposite().has_on_positive_side( p12 ) );
 CGAL_test_assert( c10.has_on_boundary( p6 ) );
 CGAL_test_assert( c10.has_on_boundary( p8 ) );

 std::cout << '.';

 CGAL::Circle_2<R> c11( p0 );
 CGAL::Circle_2<R> c12( p0, CGAL::POSITIVE );
 CGAL::Circle_2<R> c13( p0, CGAL::NEGATIVE );
 CGAL_test_assert( c11.orientation() == CGAL::POSITIVE );
 CGAL_test_assert( c12.orientation() == CGAL::POSITIVE );
 CGAL_test_assert( c13.orientation() == CGAL::NEGATIVE );
 CGAL_test_assert( c11.is_degenerate() );
 CGAL_test_assert( c12.is_degenerate() );
 CGAL_test_assert( c13.is_degenerate() );

 std::cout << '.';

 CGAL::Bbox_2 bb = c4.bbox();
 CGAL_test_assert(bb.xmin() <= -4.0);
 CGAL_test_assert(bb.xmax() >= 8.0);
 CGAL_test_assert(bb.ymin() <= -4.0);
 CGAL_test_assert(bb.ymax() >= 8.0);

 std::cout << "done" << std::endl;
 return true;
}

#endif // CGAL__TEST_CLS_CIRCLE_2_H
