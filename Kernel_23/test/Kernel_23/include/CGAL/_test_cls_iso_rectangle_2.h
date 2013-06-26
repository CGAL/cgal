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


#ifndef CGAL__TEST_CLS_ISO_RECTANGLE_2_H
#define CGAL__TEST_CLS_ISO_RECTANGLE_2_H

#include <CGAL/Bbox_2.h>
#include <cassert>

template <class R>
bool
_test_cls_iso_rectangle_2(const R& )
{
 std::cout << "Testing class Iso_rectangle_2";

 typedef typename  R::RT    RT;
 typedef typename  R::FT    FT;

 typename R::Iso_rectangle_2 ir;
 CGAL::Iso_rectangle_2<R>  r0(ir);

 RT n1 =  1;
 RT n2 =  2;
 RT n3 =  3;
 RT n4 =  4;
 RT n5 =  5;
 RT n6 =  6;
 RT n7 =  7;
 RT n8 =  8;
 RT n12= 12;

 CGAL::Point_2<R> p1( n5, n5, n5);    // ( 1, 1)
 CGAL::Point_2<R> p2( n2, n8, n2);    // ( 1, 4)
 CGAL::Point_2<R> p3( n12, n12, n3);  // ( 4, 4)
 CGAL::Point_2<R> p4( n5, n4, n1);    // ( 5, 4)
 CGAL::Point_2<R> p5( n4, n1);        // ( 4, 1)
 CGAL::Point_2<R> p6( n8, n4, n2);    // ( 4, 2)
 CGAL::Point_2<R> p7( n6, n3, n2);    // ( 3, 1.5)
 CGAL::Point_2<R> p8( n4, n6, n2);    // ( 2, 3)
 CGAL::Point_2<R> p9(-n3, n7);        // (-3, 7)

 const CGAL::Iso_rectangle_2<R> r1( p1, p3);
 CGAL::Iso_rectangle_2<R> r1_( p1, p3, 0);
 CGAL::Iso_rectangle_2<R> r2( p3, p1);
 CGAL::Iso_rectangle_2<R> r3( p2, p5);
 CGAL::Iso_rectangle_2<R> r4( p5, p2);
 CGAL::Iso_rectangle_2<R> r5( p9, p2);
 CGAL::Iso_rectangle_2<R> r6( r2 );
 CGAL::Iso_rectangle_2<R> r7( n3, n3, n12, n12, n3);
 CGAL::Iso_rectangle_2<R> r8( n4, n2, n6, n8, n2);
 CGAL::Iso_rectangle_2<R> r9( n2, n1, n3, n4, n1);
 CGAL::Iso_rectangle_2<R> r10( n2, n1, n3, n4);

 CGAL::Iso_rectangle_2<R> r11(p1, p1, p1, p1);
 CGAL::Iso_rectangle_2<R> r12(p1, p1, p2, p2);
 CGAL::Iso_rectangle_2<R> r13(p1, p2, p1, p2);
 CGAL::Iso_rectangle_2<R> r14(p3, p4, p1, p2);

 CGAL::Iso_rectangle_2<R> r15( r14.bbox() );
 typename R::Iso_rectangle_2 r16( r15.bbox() );

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

 assert( CGAL::Iso_rectangle_2<R>( p1, p1 ).is_degenerate() );
 assert( CGAL::Iso_rectangle_2<R>( p1, p2 ).is_degenerate() );
 assert( CGAL::Iso_rectangle_2<R>( p3, p4 ).is_degenerate() );

 std::cout << '.';

 assert( CGAL::Iso_rectangle_2<R>( p1, p1 ).area() == FT(0) );
 assert( CGAL::Iso_rectangle_2<R>( p1, p2 ).area() == FT(0) );
 assert( CGAL::Iso_rectangle_2<R>( p3, p4 ).area() == FT(0) );
 assert( CGAL::Iso_rectangle_2<R>( p1, p3 ).area() == FT(9) );
 assert( CGAL::Iso_rectangle_2<R>( p3, p1 ).area() == FT(9) );
 assert( CGAL::Iso_rectangle_2<R>( p1, p7 ).area() == FT(1) );
 assert( CGAL::Iso_rectangle_2<R>( p9, p3 ).area() == FT(21) );

 std::cout << '.';

 CGAL::Bbox_2 bb = r1.bbox();
 assert(bb.xmin() <= 1.0);
 assert(bb.xmax() >= 4.0);
 assert(bb.ymin() <= 1.0);
 assert(bb.ymax() >= 4.0);

 std::cout << "done" << std::endl;
 return true;
}
#endif // CGAL__TEST_CLS_ISO_RECTANGLE_2_H
