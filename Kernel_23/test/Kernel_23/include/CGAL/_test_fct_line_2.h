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
// Author(s)     : Susan Hert, Sylvain Pion


#ifndef CGAL__TEST_FCT_LINE_2_H
#define CGAL__TEST_FCT_LINE_2_H

template <class R>
bool
_test_fct_line_sqrt_2(const R&)
{
 typedef typename  R::Point_2  Point_2;
 typedef typename  R::Line_2   Line_2;

 // bisector of 2 lines (uses sqrt()...)
 Point_2 q0(0, 0, 1);
 Point_2 q1(1, 0, 1);
 Point_2 q2(0, 1, 1);
 Point_2 q3(1, 1, 1);
 Point_2 q4(2, 0, 1);

 Line_2 ql1 (q0, q1);
 Line_2 ql2 (q0, q2);
 Line_2 ql3 (q0, q3);
 Line_2 ql4 (q0, q4);
 Line_2 ql5 (q1, q0);
 Line_2 bl3 = CGAL::bisector(ql1, ql2);

 assert( bl3 == ql3 );
 assert( CGAL::bisector(ql4, ql2) == ql3 );
 assert( CGAL::bisector(ql1, ql5) == ql1 );

 return true;
}

template <class R>
bool
_test_fct_line_2(const R& )
{
 std::cout << "Testing functions Line_2" ;

 typedef typename  R::RT       RT;
 typedef typename  R::Point_2  Point_2;
 typedef typename  R::Line_2   Line_2;

 Point_2 p1 ( RT(18), RT(12), RT(3) );  // ( 6, 4)
 Point_2 p2 ( RT(18), RT(15), RT(3) );  // ( 6, 5)
 Point_2 p3 ( RT(18), RT( 9), RT(3) );  // ( 6, 3)
 Point_2 p4 ( RT(28), RT(40), RT(4) );  // ( 7,10)
 Point_2 p5 ( RT(12), RT(40), RT(4) );  // ( 3,10)
 Point_2 p6 ( RT(28), RT(12), RT(4) );  // ( 7, 3)
 Point_2 p7 ( RT(18), RT( 6), RT(3) );  // ( 6, 2)
 Point_2 p8 ( RT(24), RT( 9), RT(3) );  // ( 8, 3)
 Point_2 p9 ( RT( 6), RT(10), RT(1) );  // ( 6,10)
 Point_2 p10( RT( 8), RT( 5), RT(1) );  // ( 8, 5)
 Point_2 p11( RT( 7), RT( 5), RT(1) );  // ( 7, 5)
 Point_2 p12( RT( 0), RT( 4), RT(1) );  // ( 0, 4)

 // vertical lines
 Line_2 l1(p1, p2);
 Line_2 l2(p3, p2);
 Line_2 l3(p4, p6);

 assert( CGAL::compare_slope(l1,l2) == CGAL::EQUAL );
 assert( CGAL::compare_slope(l1,l3) == CGAL::EQUAL );
 assert( CGAL::compare_slope(l3,l1) == CGAL::EQUAL );

 std::cout <<'.';

 // horizontal lines
 Line_2 l4(p3, p8);
 Line_2 l5(p4, p9);
 assert( CGAL::compare_slope(l4, l5) == CGAL::EQUAL );
 assert( CGAL::compare_slope(l3, l4) == CGAL::LARGER );
 assert( CGAL::compare_slope(l4, l3) == CGAL::SMALLER );

 std::cout <<'.';

 // parallel lines
 Line_2 l5a(p6, p7);
 Line_2 l5b(p11, p1);
 assert( CGAL::compare_slope(l5a, l5b) == CGAL::EQUAL );
 assert( CGAL::parallel(l5a, l5b) );

 // two positive slopes
 Line_2 l6(p2, p4);
 Line_2 l7(p2, p6);
 Line_2 l8(p7, p10);
 assert( CGAL::compare_slope(l6, l6) == CGAL::EQUAL );
 assert( CGAL::compare_slope(l6, l7) == CGAL::LARGER );
 assert( CGAL::compare_slope(l7, l6) == CGAL::SMALLER );
 assert( CGAL::compare_slope(l6, l8) == CGAL::LARGER );
 assert( CGAL::compare_slope(l8, l6) == CGAL::SMALLER );
 assert(   CGAL::parallel(l6, l6) );
 assert( ! CGAL::parallel(l6, l7) );
 assert( ! CGAL::parallel(l7, l6) );
 assert( ! CGAL::parallel(l6, l8) );
 assert( ! CGAL::parallel(l8, l6) );

 // vertical and positive slope
 assert( CGAL::compare_slope(l1, l6) == CGAL::LARGER );
 assert( CGAL::compare_slope(l6, l1) == CGAL::SMALLER );
 assert( ! CGAL::parallel(l1, l6) );
 assert( ! CGAL::parallel(l6, l1) );

 // horizontal and positive slope
 assert( CGAL::compare_slope(l5, l6) == CGAL::SMALLER );
 assert( CGAL::compare_slope(l6, l5) == CGAL::LARGER );
 assert( ! CGAL::parallel(l5, l6) );
 assert( ! CGAL::parallel(l6, l5) );

 std::cout <<'.';

 // two negative slopes
 Line_2 l9 (p4, p8);
 Line_2 l10(p9, p8);
 Line_2 l11(p5, p3);

 assert( CGAL::compare_slope(l9, l10) == CGAL::SMALLER );
 assert( CGAL::compare_slope(l10, l9) == CGAL::LARGER );
 assert( CGAL::compare_slope(l11, l10) == CGAL::LARGER );
 assert( ! CGAL::parallel(l9, l10) );
 assert( ! CGAL::parallel(l10, l9) );
 assert( ! CGAL::parallel(l11, l10) );

 // vertical and negative slope
 assert( CGAL::compare_slope(l2, l9) == CGAL::LARGER );
 assert( CGAL::compare_slope(l9, l2) == CGAL::SMALLER );
 assert( ! CGAL::parallel(l2, l9) );
 assert( ! CGAL::parallel(l9, l2) );

 // horizontal and negative slope
 assert( CGAL::compare_slope(l5, l9) == CGAL::LARGER );
 assert( CGAL::compare_slope(l9, l5) == CGAL::SMALLER );

 std::cout <<'.';

 // positive and negative slope
 assert( CGAL::compare_slope(l6, l9) == CGAL::LARGER );
 assert( CGAL::compare_slope(l9, l7) == CGAL::SMALLER );
 assert( ! CGAL::parallel(l6, l9) );
 assert( ! CGAL::parallel(l9, l7) );

 std::cout <<'.';

 // bisector construction
 Line_2 bl1 = CGAL::bisector(p2, p3);
 Line_2 bl2 = CGAL::bisector(p3, p2);
 assert(bl1 == Line_2(p12, p1));
 assert(bl2 == Line_2(p1, p12));
 assert(bl1.oriented_side(p2) == CGAL::ON_POSITIVE_SIDE);
 assert( CGAL::parallel(bl1, bl2) );

 // More tests, that require sqrt() or use approx.
  _test_fct_line_sqrt_2(R());

 std::cout << "done" << std::endl;
 return true;
}

#endif // CGAL__TEST_FCT_LINE_2_H
