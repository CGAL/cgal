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


#ifndef CGAL__TEST_CLS_CIRCLE_NEW_2_H
#define CGAL__TEST_CLS_CIRCLE_NEW_2_H

template <class R>
bool
_test_cls_circle_new_2(const R& )
{
 std::cout << "Testing class Circle_2";

 typedef typename  R::RT    RT;
 typedef typename  R::FT    FT;

 typedef typename  R::Point_2 Point_2;
 typedef typename  R::Vector_2 Vector_2;
 typedef typename  R::Direction_2 Direction_2;
 typedef typename  R::Circle_2 Circle_2;
 typedef typename  R::Aff_transformation_2 Aff_transformation_2;

 typename R::Construct_vector_2 construct_vector;
 typename R::Construct_point_2 construct_point;
 typename R::Construct_translated_point_2 construct_translated_point;

 typename R::Circle_2  ic;
 Circle_2 c0; // af: CGAL::Circle_2<R> c0;

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

 Point_2 p0 = construct_point( n1, n2, -n2);  // ( 4, -1)
 Point_2 p1 = construct_point( n6, n8, n10);  // ( 2,  3)
 Point_2 p2 = construct_point( n2, n0,  n2);  // ( 1,  0)
 Point_2 p3 = construct_point( n5, n5,  n4);  // ( 2,  2)
 Point_2 p4 = construct_point( n0, n2,  n2);  // ( 0,  1)

 Vector_2 vx = construct_vector(construct_point(CGAL::ORIGIN), p2);
 Vector_2 vy = construct_vector(construct_point(CGAL::ORIGIN), p4);
 Vector_2 v1 = construct_vector(construct_point(CGAL::ORIGIN), p1);

 Circle_2 c1( p0, p1, p2);
 Circle_2 c2( p0, p1, p3);
 Circle_2 c3( p1, p0, p2);
 Circle_2 c4( p3, FT( n9 ));      // n9 = (n6)^2

 Vector_2 vx6 = vx * n6;
 Vector_2 vy6 = vy * n6;
 Circle_2 c5( construct_translated_point(p3, - vx6), construct_translated_point(p3, vx6), construct_translated_point(p3 , vy6));
 Circle_2 c6( c3 );
 Circle_2 c7( p3, n9, CGAL::POSITIVE);
 Circle_2 c8( p3, n9, CGAL::NEGATIVE);
 Circle_2 cc( construct_translated_point(p3, - vx6), construct_translated_point(p3, vx6));
 Circle_2 cp( construct_translated_point(p3, - vx6), construct_translated_point(p3, vx6), CGAL::POSITIVE);
 Circle_2 cn( construct_translated_point(p3, - vx6), construct_translated_point(p3,  vx6), CGAL::NEGATIVE);
 c0 = c3;

 assert( c1 == c1 );
 assert( c1 != c2 );
 assert( c3 == c0 );
 assert( c0 == c3 );
 assert( c3 == c6 );
 assert( c7 != c8 );
 assert( c4 == c7 );
 assert( cc == cp );
 assert( cn != cp );
 assert( cc != c8 );
 assert( cc == c7 );

 assert( c5.center() == p3 );
 assert( cc.center() == p3 );
 assert( c5.squared_radius() == FT( n9 ) );
 assert( c4.squared_radius() == cc.squared_radius() );
 assert( c4 == c5 );
 assert( c4 == c7 );
 assert( c4 != c8 );
 assert( cn == cp.opposite() );
 assert( c7.opposite() == c8 );
 assert( c8.opposite() == c7 );
 assert( c1.opposite() == c3 );
 assert( c3.opposite() == c1 );
 assert( c7.orientation() == CGAL::POSITIVE );
 assert( c8.orientation() == CGAL::NEGATIVE );
 assert( c5.orientation() == CGAL::POSITIVE );
 assert( cc.orientation() == CGAL::POSITIVE );
 assert( cp.orientation() == CGAL::POSITIVE );
 assert( cn.orientation() == CGAL::NEGATIVE );

 std::cout << '.';

 assert( c4.center() == p3 );
 assert( c5.center() == p3 );
 assert( c4.squared_radius() == FT( n9 ) );
 assert( c5.squared_radius() == FT( n9 ) );
 assert( c8.squared_radius() == FT( n9 ) );

 assert( c7.bounded_side( construct_translated_point(p3, vx*n2) ) == CGAL::ON_BOUNDED_SIDE );
 assert( c7.bounded_side( construct_translated_point(p3, vy*n11) ) == CGAL::ON_UNBOUNDED_SIDE );
 assert( c7.bounded_side( construct_translated_point(p3, - vy6) ) == CGAL::ON_BOUNDARY );
 assert( c8.bounded_side( construct_translated_point(p3, vx*n2) ) == CGAL::ON_BOUNDED_SIDE );
 assert( c8.bounded_side( construct_translated_point(p3, vy*n11) ) == CGAL::ON_UNBOUNDED_SIDE );
 assert( c8.bounded_side( construct_translated_point(p3, - vy6) ) == CGAL::ON_BOUNDARY );
 assert( cc.has_on_boundary( construct_translated_point(p3, vy6)) );
 assert( cc.has_on_boundary( construct_translated_point(p3, - vx6)) );

 std::cout << '.';

 Aff_transformation_2
          rotate1(CGAL::ROTATION,Direction_2(n11,n13),-n2,n12),
          rotate2(CGAL::ROTATION,Direction_2(-n8, n9),-n2,n12),
          rotate3(CGAL::ROTATION,Direction_2( n5,-n1),-n2,n12),
          rotate4(CGAL::ROTATION,Direction_2(-n5,-n11),-n2,n12);
 // af: Point_2 p6 = p2.transform( rotate1 );
 Point_2 p6 = rotate1(p2); // af: try this instead
 Point_2 p7 = rotate2(p2 );
 Point_2 p8 = rotate3(p2 );
 Point_2 p9 = rotate4(p2 );
 p6 =construct_translated_point( p6, v1);  // v1 = ( 2, 3 )
 p7 = construct_translated_point(p7, v1);
 p8 = construct_translated_point(p8, v1);
 p9 = construct_translated_point(p9, v1);
 Circle_2 c10 (p6, p8, p7);
 assert( c10.orientation() == CGAL::POSITIVE );
 assert( c10.opposite().orientation() == CGAL::NEGATIVE );

 assert( c10.oriented_side(c10.center() ) == CGAL::ON_POSITIVE_SIDE );
 assert( c10.oriented_side(construct_translated_point( construct_translated_point( construct_point(CGAL::ORIGIN), v1), vx/n2) ) \
         == CGAL::ON_POSITIVE_SIDE );
 assert( c10.oriented_side(construct_translated_point( construct_translated_point( construct_point(CGAL::ORIGIN), v1), vx*n2) ) \
         == CGAL::ON_NEGATIVE_SIDE );
 assert( c10.oriented_side(p9 ) == CGAL::ON_ORIENTED_BOUNDARY );
 assert( c10.has_on_boundary(p9) );
 assert( c10.has_on_boundary(construct_translated_point(p4, v1)) );
 Point_2 p11 = construct_point( n4, n4, n3) ; // (2.5, 2.5)
 Point_2 p12 = construct_point( n5, n5, n3) ; // ( 5 ,  5 )
 assert( c10.has_on_bounded_side( p11 ) );
 assert( ! c10.has_on_bounded_side( p12 ) );
 assert( c10.has_on_unbounded_side( p12 ) );
 assert( c10.has_on_positive_side( p11 ) );
 assert( c10.has_on_negative_side( p12 ) );
 assert( c10.opposite().has_on_negative_side( p11 ) );
 assert( c10.opposite().has_on_positive_side( p12 ) );
 assert( c10.has_on_boundary( p6 ) );
 assert( c10.has_on_boundary( p8 ) );

 std::cout << '.';

 Circle_2 c11( p0 );
 Circle_2 c12( p0, CGAL::POSITIVE );
 Circle_2 c13( p0, CGAL::NEGATIVE );
 assert( c11.orientation() == CGAL::POSITIVE );
 assert( c12.orientation() == CGAL::POSITIVE );
 assert( c13.orientation() == CGAL::NEGATIVE );
 assert( c11.is_degenerate() );
 assert( c12.is_degenerate() );
 assert( c13.is_degenerate() );

 std::cout << "done" << std::endl;
 return true;
}
#endif // CGAL__TEST_CLS_CIRCLE_NEW_2_H
