// Copyright (c) 2003
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

#ifndef CGAL__TEST_FCT_PLANE_3_H
#define CGAL__TEST_FCT_PLANE_3_H

template <class R>
bool
_test_fct_plane_sqrt_3(const R&)
{
 typedef typename  R::Point_3  Point_3;
 typedef typename  R::Plane_3  Plane_3;

 // bisector of 2 planes
 Point_3 q0(0, 0, 0, 1);
 Point_3 q1(1, 0, 0, 1);
 Point_3 q2(0, 1, 0, 1);
 Point_3 q3(1, 1, 0, 1);
 Point_3 q4(2, 0, 0, 1);
 Point_3 q5(0, 0, 1, 1);

 Plane_3 ql1 (q0, q1, q5);
 Plane_3 ql2 (q0, q2, q5);
 Plane_3 ql3 (q0, q3, q5);
 Plane_3 ql4 (q0, q4, q5);
 Plane_3 ql5 (q1, q0, q5);
 Plane_3 bl3 = CGAL::bisector(ql1, ql2);

 assert( bl3 == ql3 );
 assert( CGAL::bisector(ql4, ql2) == ql3 );
 assert( CGAL::bisector(ql1, ql5) == ql1 );

 return true;
}

template <class R>
bool
_test_fct_plane_3(const R& )
{
 std::cout << "Testing functions Plane_3" ;

 typedef typename  R::RT       RT;
 typedef typename  R::Point_3  Point_3;
 typedef typename  R::Plane_3  Plane_3;

 Point_3 p1 ( RT(18), RT(12), RT(3), RT(3) );  // ( 6, 4, 1)
 Point_3 p2 ( RT(18), RT(15), RT(3), RT(3) );  // ( 6, 5, 1)
 Point_3 p3 ( RT(18), RT( 9), RT(3), RT(3) );  // ( 6, 3, 1)
 Point_3 p4 ( RT(28), RT(40), RT(4), RT(4) );  // ( 7,10, 1)
 Point_3 p5 ( RT(12), RT(40), RT(4), RT(4) );  // ( 3,10, 1)
 Point_3 p6 ( RT(28), RT(12), RT(4), RT(4) );  // ( 7, 3, 1)
 Point_3 p7 ( RT(18), RT( 6), RT(3), RT(3) );  // ( 6, 2, 1)
 Point_3 p8 ( RT(24), RT( 9), RT(3), RT(3) );  // ( 8, 3, 1)
 Point_3 p9 ( RT( 6), RT(10), RT(1), RT(1) );  // ( 6,10, 1)
 Point_3 p10( RT( 8), RT( 5), RT(1), RT(1) );  // ( 8, 5, 1)
 Point_3 p11( RT( 7), RT( 5), RT(1), RT(1) );  // ( 7, 5, 1)
 Point_3 p12( RT( 0), RT( 4), RT(1), RT(1) );  // ( 0, 4, 1)
 Point_3 p13( RT(18), RT(12), RT(0), RT(3) );  // ( 6, 4, 0)

 // bisector construction
 Plane_3 h1 = CGAL::bisector(p2, p3);
 Plane_3 h2 = CGAL::bisector(p3, p2);
 Plane_3 h3 (p12, p1, p13);
 Plane_3 h4 (p1, p12, p13);
 assert(h1 == h3);
 assert(h2 == h4);
 assert(h1.oriented_side(p2) == CGAL::ON_POSITIVE_SIDE);

 Plane_3 h5 (p7, p8, p9);
 assert(   CGAL::parallel(h1, h2) );
 assert( ! CGAL::parallel(h1, h5) );

 // More tests, that require sqrt() or use approx.
 _test_fct_plane_sqrt_3(R());

 std::cout << "done" << std::endl;
 return true;
}

#endif // CGAL__TEST_FCT_PLANE_3_H
