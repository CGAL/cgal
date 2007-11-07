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
 

#ifndef CGAL__TEST_CLS_SEGMENT_3_H
#define CGAL__TEST_CLS_SEGMENT_3_H

template <class R>
bool
_test_cls_segment_3(const R& )
{
 std::cout << "Testing class Segment_3" ;

 typedef typename  R::RT    RT;
 typedef typename  R::FT    FT;

 typename R::Segment_3 is;
 CGAL::Segment_3<R>  s1(is);

 RT  n1 =  7;
 RT  n2 = 21;
 RT  n3 = 14;
 RT  n4 =-10;
 RT  n5 =  5;
 RT  n6 = 20;
 RT  n8 =  3;

 CGAL::Point_3<R> p1( n1, n2, n3, n1);
 CGAL::Point_3<R> p2( n4, n5, n6, n5);
 CGAL::Point_3<R> p3( n2, n8, n2, n8);

 CGAL::Segment_3<R> s2( p1, p2 );
 CGAL::Segment_3<R> s3( p2, p1 );
 CGAL::Segment_3<R> s4( s2 );
 s1 = s4;

 CGAL_test_assert( CGAL::parallel(s2, s3) );

 CGAL::Vector_3<R> v0(p1, p2);
 CGAL_test_assert( v0 == s2.to_vector() );

 CGAL_test_assert( s1 == s1 );
 CGAL_test_assert( s4 == s2 );
 CGAL_test_assert( s1 == s4 );
 CGAL_test_assert( s1 == s2 );
 CGAL_test_assert( s2 != s3 );

 CGAL::Segment_3<R> s5 (p3, p3 + (p1 - p3) + (p1 - p3) );
 CGAL_test_assert( s5.has_on( p1 ) );
 CGAL_test_assert( s5.has_on( p3 ) );
 CGAL_test_assert( s2.has_on( p2 ) );

 CGAL_test_assert( ! CGAL::parallel(s1, s5) );

 std::cout <<'.';

 CGAL_test_assert( s5.source() == p3 );
 CGAL_test_assert( s5.target() == p1 + (p1 - p3) );
 CGAL_test_assert( s2.source() != s3.source() );
 CGAL_test_assert( s2.target() != s3.target() );

 std::cout <<'.';

 CGAL_test_assert( (s2.min)() == p2 );
 CGAL_test_assert( (s2.max)() == p1 );
 CGAL_test_assert( (s2.min)() == (s3.min)() );
 CGAL_test_assert( (s2.max)() == (s3.max)() );
 CGAL_test_assert( (s5.max)() != (s5.min)() );
 CGAL_test_assert( (s5.max)() == (s5.opposite().max)() );
 CGAL_test_assert( s5.vertex(0) == s5.source() );
 CGAL_test_assert( s2.vertex(1) == s2.target() );
 CGAL_test_assert( s2.vertex(1) == (s2.min)() );
 CGAL_test_assert( s2[1] == s1[1] );
 CGAL_test_assert( s2[1] == s3[0] );

 std::cout << '.';

 CGAL_test_assert( s2.squared_length() == FT( RT(17) ) );
 CGAL_test_assert( s2.direction() == CGAL::Direction_3<R>(s2.target() - s2.source() ));
 CGAL_test_assert( s2.direction() == s3.opposite().direction() );

 CGAL_test_assert( s1.supporting_line() == s2.supporting_line() );
 CGAL::Line_3<R> lin(p1,p2);
 CGAL_test_assert( s2.supporting_line() == lin );

 CGAL::Segment_3<R> sdeg(p3,p3);
 CGAL_test_assert( sdeg.is_degenerate() );
 CGAL_test_assert( CGAL::parallel(sdeg, sdeg) );

 std::cout << '.';

 CGAL::Bbox_3 bb = s2.bbox();
 CGAL_test_assert(bb.xmin() <= -2.0);
 CGAL_test_assert(bb.xmax() >= 1.0);
 CGAL_test_assert(bb.ymin() <= 1.0);
 CGAL_test_assert(bb.ymax() >= 3.0);
 CGAL_test_assert(bb.zmin() <= 2.0);
 CGAL_test_assert(bb.zmax() >= 4.0);

 std::cout << "done" << std::endl;
 return true;
}
#endif // CGAL__TEST_CLS_SEGMENT_3_H
