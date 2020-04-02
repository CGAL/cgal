// Copyright (c) 2001
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
// Author(s)     : Sylvain Pion


#ifndef CGAL__TEST_FCT_POINT_2_SEGMENT_2_H
#define CGAL__TEST_FCT_POINT_2_SEGMENT_2_H


template <class R>
bool
_test_fct_point_segment_2(const R& )
{
 std::cout << "Testing functions Point_2 Segment_2" ;

 typedef typename  R::RT    RT;

 RT n0 =  0;
 RT n1 =  1;
 RT n2 =  2;
 RT n3 =  3;
 RT n4 =  4;
 RT n5 =  5;
 RT n6 =  6;
 RT n7 =  7;
 RT n8 =  8;
 RT n9 =  9;
 RT n10= 10;
 RT n12= 12;
 RT n20= 20;

 CGAL::Point_2<R>  p1( n2, n8, n2);   // ( 1, 4)
 CGAL::Point_2<R>  p2( n4, n4, n2);   // ( 2, 2)
 CGAL::Point_2<R>  p3(n12,n10, n2);   // ( 6, 5)
 CGAL::Point_2<R>  p4( n7, n3);       // ( 7, 3)
 CGAL::Point_2<R>  p5( n9, n0, n3);   // ( 3, 0)
 CGAL::Point_2<R>  p6(n12,n20, n4);   // ( 3, 5)
 CGAL::Point_2<R>  p7(n12, n3, n3);   // ( 4, 1)
 CGAL::Point_2<R>  p8( n8, n4);       // ( 8, 4)
 CGAL::Point_2<R>  p9(-n4, n8, n4);   // (-1, 2)
 CGAL::Point_2<R>  p10(n10, n4, n2);  // ( 5, 2)

 CGAL::Point_2<R>  p11( n1, -n5,-n1); // (-1, 5)
 CGAL::Point_2<R>  p0( CGAL::ORIGIN ); // ( 0, 0)
 CGAL::Point_2<R>  p13( n8, n8, n4);  // ( 2, 2)
 CGAL::Point_2<R>  p14( n5,-n1);      // ( 5,-1)
 CGAL::Point_2<R>  p16(n12, n9, n3);  // ( 4, 3)
 CGAL::Point_2<R>  p17( n0,n1);       // ( 0, 1)
 CGAL::Point_2<R>  p18( n3,n2);       // ( 3, 2)
 CGAL::Point_2<R>  p19( n4, n6, n2);  // ( 2, 3)

 CGAL::Point_2<R>  p20( n8, n8, n2);  // ( 4, 4)
 CGAL::Point_2<R>  p21( n5, n0);      // ( 5, 0)
 CGAL::Point_2<R>  p22( n5, n1);      // ( 5, 1)
 CGAL::Point_2<R>  p23( n8, n20, n2); // ( 4, 10)

 CGAL::Segment_2<R>   s15( p1, p5 );
 CGAL::Segment_2<R>   s16( p1, p6 );
 CGAL::Segment_2<R>   s23( p2, p3 );
 CGAL::Segment_2<R>   s58( p5, p8 );
 CGAL::Segment_2<R>   s716( p7, p16 ); // vertical
 CGAL::Segment_2<R>   s1617( p16, p17);
 CGAL::Segment_2<R>   s013( p0, p13 );
 CGAL::Segment_2<R>   s910( p9, p10 );
 CGAL::Segment_2<R>   s114( p1, p14 );
 CGAL::Segment_2<R>   s121( p1, p21 );
 CGAL::Segment_2<R>   s122( p1, p22 );
 CGAL::Segment_2<R>   s2023( p20, p23 );
 CGAL::Segment_2<R>   s1620( p16, p20 );

 std::cout << '.';

 assert( CGAL::compare_y_at_x( p6, s23 ) == CGAL::LARGER );
 assert( CGAL::compare_y_at_x( p6, s23.opposite() ) == CGAL::LARGER );
 assert( CGAL::compare_y_at_x( p10, s23 ) == CGAL::SMALLER );
 assert( CGAL::compare_y_at_x( p17, s910 ) == CGAL::SMALLER );
 assert( CGAL::compare_y_at_x( p8, s58 ) == CGAL::EQUAL );
 assert( CGAL::compare_y_at_x( p2, s1617 ) == CGAL::EQUAL );

 assert( CGAL::compare_y_at_x( p7, s716 ) == CGAL::EQUAL );
 assert( CGAL::compare_y_at_x( p16, s716 ) == CGAL::EQUAL );
 assert( CGAL::compare_y_at_x( p20, s716 ) == CGAL::LARGER );

 assert( CGAL::compare_y_at_x( p5, s23, s16 ) == CGAL::SMALLER );
 assert( CGAL::compare_y_at_x( p2, s013, s15 ) == CGAL::EQUAL );

 assert( CGAL::compare_y_at_x( p7, s716, s114 ) == CGAL::LARGER );
 assert( CGAL::compare_y_at_x( p7, s114, s716 ) == CGAL::SMALLER );
 assert( CGAL::compare_y_at_x( p16, s716, s114 ) == CGAL::LARGER );
 assert( CGAL::compare_y_at_x( p16, s114, s716 ) == CGAL::SMALLER );

 assert( CGAL::compare_y_at_x( p7, s716, s121 ) == CGAL::EQUAL );
 assert( CGAL::compare_y_at_x( p7, s121, s716 ) == CGAL::EQUAL );
 assert( CGAL::compare_y_at_x( p16, s716, s121 ) == CGAL::EQUAL );
 assert( CGAL::compare_y_at_x( p16, s121, s716 ) == CGAL::EQUAL );

 assert( CGAL::compare_y_at_x( p7, s716, s122 ) == CGAL::EQUAL );
 assert( CGAL::compare_y_at_x( p7, s122, s716 ) == CGAL::EQUAL );
 assert( CGAL::compare_y_at_x( p16, s716, s122 ) == CGAL::EQUAL );
 assert( CGAL::compare_y_at_x( p16, s122, s716 ) == CGAL::EQUAL );

 assert( CGAL::compare_y_at_x( p16, s2023, s716 ) == CGAL::LARGER );
 assert( CGAL::compare_y_at_x( p16, s716, s2023 ) == CGAL::SMALLER );

 assert( CGAL::compare_y_at_x( p16, s1620, s716 ) == CGAL::EQUAL );
 assert( CGAL::compare_y_at_x( p16, s716, s1620 ) == CGAL::EQUAL );

 std::cout << "done" << std::endl;
 return true;
}
#endif // CGAL__TEST_FCT_POINT_2_SEGMENT_2_H
