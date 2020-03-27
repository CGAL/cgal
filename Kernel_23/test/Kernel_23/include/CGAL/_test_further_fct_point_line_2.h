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


#ifndef CGAL__TEST_FURTHER_FCT_POINT_2_LINE_2_H
#define CGAL__TEST_FURTHER_FCT_POINT_2_LINE_2_H

template <class R>
bool
_test_further_fct_point_line_2(const R& )
{
 std::cout << "Testing further functions Point_2 Line_2" ;

 typedef typename  R::RT    RT;

 CGAL::Point_2<R> p1( RT(18), RT(12), RT(3) );  // ( 6, 4)
 CGAL::Point_2<R> p2( RT(18), RT(15), RT(3) );  // ( 6, 5)
 CGAL::Point_2<R> p3( RT(18), RT( 9), RT(3) );  // ( 6, 3)
 CGAL::Point_2<R> p4( RT(28), RT(40), RT(4) );  // ( 7,10)
 CGAL::Point_2<R> p5( RT(12), RT(-4), RT(4) );  // ( 3,-1)
 CGAL::Point_2<R> p6( RT(28), RT(12), RT(4) );  // ( 7, 3)
 CGAL::Point_2<R> p7( RT(18), RT( 6), RT(3) );  // ( 6, 2)
 CGAL::Point_2<R> p8( RT(24), RT( 9), RT(3) );  // ( 8, 3)
 CGAL::Point_2<R> p9( RT( 6), RT(10), RT(1) );  // ( 6,10)

 CGAL::Line_2<R> l12(p1,p2);
 CGAL::Line_2<R> l21(p2,p1);
 CGAL::Line_2<R> l18(p1,p8);
 CGAL::Line_2<R> l64(p6,p4);
 CGAL::Line_2<R> l68(p6,p8);
 CGAL::Line_2<R> l61(p6,p1);
 CGAL::Line_2<R> l46(p4,p6);

 std::cout << '.';

 assert( CGAL::compare_signed_distance_to_line(l12, p3,p9) == CGAL::EQUAL );
 assert( CGAL::compare_signed_distance_to_line(l12, p3,p8) == CGAL::LARGER );
 assert( CGAL::compare_signed_distance_to_line(l21, p3,p8) == CGAL::SMALLER );
 assert( CGAL::compare_signed_distance_to_line(l18, p9,p8) == CGAL::LARGER );
 assert( CGAL::compare_signed_distance_to_line(l18, p9,p9) == CGAL::EQUAL );
 assert( CGAL::compare_signed_distance_to_line(l18, p5,p4) == CGAL::SMALLER );
 assert( CGAL::compare_signed_distance_to_line(l18, p4,p5) == CGAL::LARGER );
 assert( CGAL::compare_signed_distance_to_line(l18.opposite(),p4,p5) ==
                                                              CGAL::SMALLER );
 assert( CGAL::compare_signed_distance_to_line(l64, p1,p2) == CGAL::EQUAL );
 assert( CGAL::compare_signed_distance_to_line(l64, p1,p9) == CGAL::EQUAL );
 assert( CGAL::compare_signed_distance_to_line(l68, p7,p9) == CGAL::SMALLER );
 assert( CGAL::compare_signed_distance_to_line(l68, p2,p9) == CGAL::SMALLER );
 assert( CGAL::compare_signed_distance_to_line(l68, p2,p3) == CGAL::LARGER );
 assert( CGAL::compare_signed_distance_to_line(l61, p1,p6) == CGAL::EQUAL );
 assert( CGAL::compare_signed_distance_to_line(l61, p6,p1) == CGAL::EQUAL );
 assert( CGAL::compare_signed_distance_to_line(l61, p4,p1) == CGAL::SMALLER );
 assert( CGAL::compare_signed_distance_to_line(l61, p5,p6) == CGAL::LARGER );
 assert( CGAL::compare_signed_distance_to_line(l46, p8,p6) == CGAL::LARGER );
 assert( CGAL::compare_signed_distance_to_line(l64, p8,p6) == CGAL::SMALLER );

 std::cout << '.';

 assert( CGAL::has_larger_signed_distance_to_line(l12, p3,p8) );
 assert( CGAL::has_larger_signed_distance_to_line(l18, p9,p8) );
 assert( CGAL::has_larger_signed_distance_to_line(l61, p1,p4) );
 assert( CGAL::has_larger_signed_distance_to_line(l64, p6,p8) );

 std::cout << '.';

 assert( CGAL::has_smaller_signed_distance_to_line(l18, p5,p4) );
 assert( CGAL::has_smaller_signed_distance_to_line(l21, p3,p8) );
 assert( CGAL::has_smaller_signed_distance_to_line(l61, p6,p5) );
 assert( CGAL::has_smaller_signed_distance_to_line(l68, p2,p9) );
 assert( CGAL::has_smaller_signed_distance_to_line(l64, p8,p6) );

 std::cout << "done" << std::endl;
 return true;
}
#endif // CGAL__TEST_FURTHER_FCT_POINT_2_LINE_2_H
