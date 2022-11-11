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


#ifndef CGAL__TEST_CLS_SEGMENT_NEW_2_H
#define CGAL__TEST_CLS_SEGMENT_NEW_2_H

template <class R>
bool
_test_cls_segment_new_2(const R& )
{
 std::cout << "Testing class Segment_2";

 typedef typename  R::RT    RT;
 typedef typename  R::FT    FT;

 typedef typename  R::Point_2 Point_2;
 typedef typename  R::Vector_2 Vector_2;
 typedef typename  R::Direction_2 Direction_2;

 typedef typename  R::Segment_2 Segment_2;
 typedef typename  R::Line_2 Line_2;

 typename R::Construct_vector_2 construct_vector;
 typename R::Construct_point_2 construct_point;

 typename R::Segment_2  is;
 Segment_2 s0; // af: CGAL::Segment_2<R> s0;

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

 Segment_2 s1 ( p1, p2 );
 Segment_2 s2 ( p2, p1 );
 Segment_2 s3 ( p2, p3 );
 Segment_2 s4 ( p2, p4 );
 Segment_2 s5 ( p4, p1 );
 Segment_2 s6 ( s3 );
 s0 = s2;

 Vector_2 v1 = s1.to_vector();
 assert( v1 == construct_vector(p1, p2) );

 assert( s5 == s5 );
 assert( s6 == s3 );
 assert( s2 == s0 );
 assert( s0 == s2 );
 assert( s1 != s5 );
 assert( s1 != s2 );
 assert( s3 != s2 );

 std::cout << '.';

 assert( s1.source() == p1 );
 assert( s2.source() == p2 );
 assert( s6.source() == p2 );
 assert( s4.target() == p4 );
 assert( s5.target() == p1 );

 assert( (s1.min)() == p1 );
 assert( (s3.min)() == p3 );
 assert( (s4.min)() == p2 );
 assert( (s4.max)() == p4 );
 assert( (s5.max)() == p4 );

 assert( s3.vertex(0) == p2 );
 assert( s3.vertex(1) == p3 );
 assert( s3.vertex(2) == p2 );
 assert( s4.point(8) == s4.vertex(8) );
 assert( s1.point(3) == s1.vertex(3) );
 assert( s5[0] == s5.vertex(0) );
 assert( s6[1] == s6.vertex(1) );

 std::cout << '.';

 assert( s1.squared_length() == FT( 20 ) );
 assert( s5.squared_length() == FT( 16 ) );
 assert( s4.direction() == Direction_2(construct_vector( p2, p4 )) );
 assert( s2.opposite() == s1 );
 assert( s1.opposite() == s2 );
 assert( s1.supporting_line() == Line_2( p1, p2 ) );
 assert( s3.supporting_line() == Line_2( p2, p3 ) );
 assert( ! s1.is_horizontal() );
 assert( ! s1.is_vertical() );
 assert( s5.is_horizontal() );
 assert( s4.is_vertical() );

 std::cout << '.';

 assert( s1.has_on( p1 ) );
 assert( s1.has_on( p2 ) );
 assert( s1.has_on( p3 ) );
 assert( s2.has_on( p3 ) );
 assert( ! s4.has_on( p3 ) );
 assert( s1.collinear_has_on( p3 ) );
 assert( s2.collinear_has_on( p1 ) );
 assert( ! s3.collinear_has_on( p1 ) );
 assert( s3.collinear_has_on( construct_point( n8, n5, n2)) );

 assert( Segment_2( p3, p3).is_degenerate() );


 std::cout << "done" << std::endl;
 return true;
}
#endif // CGAL__TEST_CLS_SEGMENT_NEW_2_H
