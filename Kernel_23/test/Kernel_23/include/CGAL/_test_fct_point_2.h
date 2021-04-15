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


#ifndef CGAL__TEST_FCT_POINT_2_H
#define CGAL__TEST_FCT_POINT_2_H

#include <cassert>
#include <iostream>

template <class R>
bool
_test_fct_point_2(const R& )
{
 std::cout << "Testing functions Point_2" ;

 typedef typename  R::RT    RT;
 typedef typename  R::FT    FT;

 CGAL::Point_2<R> p1( RT(18), RT(12), RT(3) );  // ( 6, 4)
 CGAL::Point_2<R> p2( RT(18), RT(15), RT(3) );  // ( 6, 5)
 CGAL::Point_2<R> p3( RT(18), RT( 9), RT(3) );  // ( 6, 3)
 CGAL::Point_2<R> p4( RT(28), RT(40), RT(4) );  // ( 7,10)
 CGAL::Point_2<R> p5( RT(12), RT(-4), RT(4) );  // ( 3,-1)
 CGAL::Point_2<R> p6( RT(28), RT(12), RT(4) );  // ( 7, 3)
 CGAL::Point_2<R> p7( RT(18), RT( 6), RT(3) );  // ( 6, 2)
 CGAL::Point_2<R> p8( RT(24), RT( 9), RT(3) );  // ( 8, 3)
 CGAL::Point_2<R> p9( RT( 6), RT(10), RT(1) );  // ( 6,10)

 assert( !CGAL::less_x(p1,p1) );
 assert( !CGAL::less_x(p1,p2) );
 assert(  CGAL::less_x(p1,p4) );
 assert( !CGAL::less_x(p4,p1) );

 assert( !CGAL::less_y(p1,p1) );
 assert( !CGAL::less_y(p3,p6) );
 assert( !CGAL::less_y(p1,p3) );
 assert(  CGAL::less_y(p3,p1) );

 assert( CGAL::compare_xy(p1,p2) == CGAL::SMALLER );
 assert( CGAL::compare_xy(p3,p2) == CGAL::SMALLER );
 assert( CGAL::compare_xy(p3,p1) == CGAL::SMALLER );
 assert( CGAL::compare_xy(p3,p2) == CGAL::SMALLER );
 assert( CGAL::compare_xy(p2,p1) == CGAL::LARGER );
 assert( CGAL::compare_xy(p2,p3) == CGAL::LARGER );
 assert( CGAL::compare_xy(p4,p3) == CGAL::LARGER );
 assert( CGAL::compare_xy(p4,p4) == CGAL::EQUAL );

 assert( CGAL::compare_lexicographically(p1,p2) == CGAL::SMALLER );
 assert( CGAL::compare_lexicographically(p3,p2) == CGAL::SMALLER );
 assert( CGAL::compare_lexicographically(p3,p1) == CGAL::SMALLER );
 assert( CGAL::compare_lexicographically(p2,p1) == CGAL::LARGER );
 assert( CGAL::compare_lexicographically(p2,p3) == CGAL::LARGER );
 assert( CGAL::compare_lexicographically(p4,p3) == CGAL::LARGER );
 assert( CGAL::compare_lexicographically(p4,p4) == CGAL::EQUAL );

 assert( CGAL::lexicographically_xy_smaller_or_equal(p1,p1) );
 assert( CGAL::lexicographically_xy_smaller_or_equal(p3,p1) );
 assert( CGAL::lexicographically_xy_smaller_or_equal(p3,p2) );
 assert( CGAL::lexicographically_xy_smaller_or_equal(p3,p4) );

 assert( !CGAL::lexicographically_xy_smaller(p3,p3) );
 assert(  CGAL::lexicographically_xy_smaller(p3,p2) );
 assert( !CGAL::lexicographically_xy_smaller(p4,p3) );

 assert( CGAL::lexicographically_xy_larger(p2,p1) );
 assert( CGAL::lexicographically_xy_larger(p1,p3) );
 assert( CGAL::lexicographically_xy_larger(p2,p3) );
 assert( CGAL::lexicographically_xy_larger(p4,p3) );

 assert( CGAL::lexicographically_xy_larger_or_equal(p3,p3) );
 assert( CGAL::lexicographically_xy_larger_or_equal(p2,p3) );
 assert( CGAL::lexicographically_xy_larger_or_equal(p4,p3) );

 std::cout <<'.';

 CGAL::Point_2<R> pe0( RT(1), RT(0) );
 CGAL::Point_2<R> pe1( RT(0), RT(1) );

 assert( CGAL::orientation( CGAL::Point_2<R>(CGAL::ORIGIN), pe0, pe1 ) \
                           == CGAL::POSITIVE);

 assert( CGAL::orientation( p1, p2, p3) == CGAL::COLLINEAR );
 assert( CGAL::orientation( p1, p2, p4) == CGAL::RIGHT_TURN );
 assert( CGAL::orientation( p2, p1, p4) == CGAL::LEFT_TURN );
 assert( CGAL::orientation( p5, p4, p3) == CGAL::RIGHT_TURN );
 assert( CGAL::orientation( p2, p4, p6) == CGAL::RIGHT_TURN );
 assert( CGAL::orientation( p6, p4, p2) == CGAL::LEFT_TURN );
 assert( CGAL::orientation( p4, p6, p2) == CGAL::RIGHT_TURN );
 assert( CGAL::orientation( p5, p6, p7) == CGAL::COLLINEAR );
 assert( CGAL::orientation( p6, p5, p7) == CGAL::COLLINEAR );

 assert( CGAL::collinear( p1, p2, p3 ) );
 assert( CGAL::collinear( p1, p2, p7 ) );
 assert( CGAL::collinear( p6, p5, p7 ) );
 assert( CGAL::collinear( p1, p2, p3 ) );
 assert( !CGAL::collinear( p1, p2, p4 ) );
 assert( CGAL::collinear( p6, p6, p3 ) );

 assert( CGAL::left_turn( p1, p4, p2 ) );
 assert( CGAL::left_turn( p6, p4, p2 ) );

 assert( CGAL::right_turn( p4, p6, p2 ) );
 assert( CGAL::right_turn( p1, p2, p4 ) );

 std::cout << '.';

 assert( CGAL::are_ordered_along_line( p5, p7, p6 ) );   // p7 between p6 and p5
 assert( CGAL::are_ordered_along_line( p6, p7, p5 ) );
 assert( CGAL::are_ordered_along_line( p2, p1, p3 ) );
 assert( !CGAL::are_ordered_along_line( p7, p6, p5 ) );
 assert( !CGAL::are_ordered_along_line( p7, p5, p6 ) );
 assert( !CGAL::are_ordered_along_line( p7, p4, p6 ) );
 assert( !CGAL::are_ordered_along_line( p2, p4, p6 ) );
 assert( CGAL::collinear_are_ordered_along_line( p5, p7, p6 ) );
 assert( CGAL::collinear_are_ordered_along_line( p6, p7, p5 ) );
 assert( CGAL::collinear_are_ordered_along_line( p2, p1, p3 ) );
 assert( !CGAL::collinear_are_ordered_along_line( p7, p6, p5 ) );

 assert( CGAL::collinear_are_ordered_along_line( p7, p7, p7 ) );
 assert( !CGAL::collinear_are_ordered_along_line( p5, p6, p5 ) );
 assert( !CGAL::collinear_are_ordered_along_line( p1, p3, p1 ) );
 assert( !CGAL::collinear_are_ordered_along_line( p3, p6, p3 ) );

 assert( CGAL::are_strictly_ordered_along_line( p5, p7, p6 ) );
 assert( CGAL::are_strictly_ordered_along_line( p6, p7, p5 ) );
 assert( CGAL::are_strictly_ordered_along_line( p2, p1, p3 ) );
 assert( !CGAL::are_strictly_ordered_along_line( p7, p6, p5 ) );
 assert( !CGAL::are_strictly_ordered_along_line( p7, p5, p6 ) );
 assert( !CGAL::are_strictly_ordered_along_line( p7, p4, p6 ) );
 assert( !CGAL::are_strictly_ordered_along_line( p2, p4, p6 ) );
 assert( CGAL::collinear_are_strictly_ordered_along_line( p5, p7, p6 ) );
 assert( CGAL::collinear_are_strictly_ordered_along_line( p6, p7, p5 ) );
 assert( CGAL::collinear_are_strictly_ordered_along_line( p2, p1, p3 ) );
 assert( !CGAL::collinear_are_strictly_ordered_along_line( p7, p6, p5 ) );

 assert( !CGAL::are_strictly_ordered_along_line( p7, p7, p6 ) );
 assert( !CGAL::collinear_are_strictly_ordered_along_line( p5, p6, p6 ) );
 assert( CGAL::are_ordered_along_line( p7, p7, p6 ) );
 assert( CGAL::collinear_are_ordered_along_line( p5, p6, p6 ) );
 assert( CGAL::are_ordered_along_line( p6, p6, p6 ) );
 assert( CGAL::collinear_are_ordered_along_line( p5, p5, p5 ) );
 assert( !CGAL::are_strictly_ordered_along_line( p6, p6, p6 ) );
 assert( !CGAL::collinear_are_strictly_ordered_along_line( p5, p5, p5 ) );

 std::cout << '.';

 assert( CGAL::compare_yx(p1,p2) == CGAL::SMALLER );
 assert( CGAL::compare_yx(p2,p1) == CGAL::LARGER );
 assert( CGAL::compare_yx(p2,p2) == CGAL::EQUAL );
 assert( CGAL::compare_yx(p2,p4) == CGAL::SMALLER );
 assert( CGAL::compare_yx(p3,p6) == CGAL::SMALLER );
 assert( CGAL::compare_yx(p6,p3) == CGAL::LARGER );
 assert( CGAL::compare_yx(p4,p9) == CGAL::LARGER );

 assert( CGAL::lexicographically_yx_smaller(p1,p2) );
 assert( CGAL::lexicographically_yx_smaller(p2,p4) );
 assert( CGAL::lexicographically_yx_smaller(p3,p6) );
 assert( CGAL::lexicographically_yx_smaller(p2,p4) );
 assert( CGAL::lexicographically_yx_smaller(p8,p1) );
 assert(!CGAL::lexicographically_yx_smaller(p1,p8) );

 assert( CGAL::lexicographically_yx_larger(p2,p1) );
 assert( CGAL::lexicographically_yx_larger(p6,p3) );
 assert( CGAL::lexicographically_yx_larger(p4,p9) );
 assert( CGAL::lexicographically_yx_larger(p6,p7) );
 assert(!CGAL::lexicographically_yx_larger(p6,p8) );

 assert( CGAL::lexicographically_yx_smaller_or_equal(p1,p2) );
 assert( CGAL::lexicographically_yx_smaller_or_equal(p3,p6) );
 assert( CGAL::lexicographically_yx_smaller_or_equal(p4,p4) );
 assert( CGAL::lexicographically_yx_smaller_or_equal(p5,p3) );

 assert(!CGAL::lexicographically_yx_larger_or_equal(p5,p3) );
 assert( CGAL::lexicographically_yx_larger_or_equal(p3,p5) );
 assert( CGAL::lexicographically_yx_larger_or_equal(p2,p7) );
 assert( CGAL::lexicographically_yx_larger_or_equal(p7,p7) );
 assert( CGAL::lexicographically_yx_larger_or_equal(p8,p3) );

 std::cout << '.';

 assert( CGAL::compare_distance_to_point(p3,p2,p1) == CGAL::LARGER );
 assert( CGAL::compare_distance_to_point(p1,p5,p1) == CGAL::LARGER );
 assert( CGAL::compare_distance_to_point(p4,p6,p5) == CGAL::SMALLER );
 assert( CGAL::compare_distance_to_point(p4,p9,p1) == CGAL::SMALLER );
 assert( CGAL::compare_distance_to_point(p8,p3,p3) == CGAL::EQUAL );
 assert( CGAL::compare_distance_to_point(p2,p3,p3) == CGAL::EQUAL );
 assert( CGAL::compare_distance_to_point(p2,p4,p9) == CGAL::LARGER );

 assert( CGAL::compare_distance(p3,p2,p1) == CGAL::LARGER );
 assert( CGAL::compare_distance(p1,p5,p1) == CGAL::LARGER );
 assert( CGAL::compare_distance(p4,p6,p5) == CGAL::SMALLER );
 assert( CGAL::compare_distance(p4,p9,p1) == CGAL::SMALLER );
 assert( CGAL::compare_distance(p8,p3,p3) == CGAL::EQUAL );
 assert( CGAL::compare_distance(p2,p3,p3) == CGAL::EQUAL );
 assert( CGAL::compare_distance(p2,p4,p9) == CGAL::LARGER );

 assert( CGAL::compare_distance(p3,p2,p1,p3) == CGAL::LARGER );
 assert( CGAL::compare_distance(p1,p5,p1,p1) == CGAL::LARGER );
 assert( CGAL::compare_distance(p4,p6,p5,p4) == CGAL::SMALLER );
 assert( CGAL::compare_distance(p4,p9,p1,p4) == CGAL::SMALLER );
 assert( CGAL::compare_distance(p8,p3,p3,p8) == CGAL::EQUAL );
 assert( CGAL::compare_distance(p2,p3,p3,p2) == CGAL::EQUAL );
 assert( CGAL::compare_distance(p2,p4,p9,p2) == CGAL::LARGER );

 assert( CGAL::compare_squared_distance(p3,p2,CGAL::squared_distance(p1,p3)) == CGAL::LARGER );
 assert( CGAL::compare_squared_distance(p1,p5,CGAL::squared_distance(p1,p1)) == CGAL::LARGER );
 assert( CGAL::compare_squared_distance(p4,p6,CGAL::squared_distance(p5,p4)) == CGAL::SMALLER );
 assert( CGAL::compare_squared_distance(p4,p9,CGAL::squared_distance(p1,p4)) == CGAL::SMALLER );
 assert( CGAL::compare_squared_distance(p8,p3,CGAL::squared_distance(p3,p8)) == CGAL::EQUAL );
 assert( CGAL::compare_squared_distance(p2,p3,CGAL::squared_distance(p3,p2)) == CGAL::EQUAL );
 assert( CGAL::compare_squared_distance(p2,p4,CGAL::squared_distance(p9,p2)) == CGAL::LARGER );

 assert( CGAL::has_larger_distance_to_point(p3,p2,p1) );
 assert( CGAL::has_larger_distance_to_point(p1,p5,p2) );
 assert( CGAL::has_larger_distance_to_point(p1,p5,p1) );
 assert(!CGAL::has_larger_distance_to_point(p1,p5,p5) );

 assert( CGAL::has_smaller_distance_to_point(p4,p9,p1) );
 assert( CGAL::has_smaller_distance_to_point(p4,p6,p5) );
 assert( CGAL::has_smaller_distance_to_point(p9,p9,p6) );
 assert(!CGAL::has_smaller_distance_to_point(p8,p3,p3) );

 std::cout << '.';

 assert( CGAL::compare_signed_distance_to_line(p1,p2, p3,p9) == CGAL::EQUAL );
 assert( CGAL::compare_signed_distance_to_line(p1,p2, p3,p8) == CGAL::LARGER );
 assert( CGAL::compare_signed_distance_to_line(p2,p1, p3,p8) == CGAL::SMALLER );
 assert( CGAL::compare_signed_distance_to_line(p1,p8, p9,p8) == CGAL::LARGER );
 assert( CGAL::compare_signed_distance_to_line(p1,p8, p9,p9) == CGAL::EQUAL );
 assert( CGAL::compare_signed_distance_to_line(p1,p8, p5,p4) == CGAL::SMALLER );
 assert( CGAL::compare_signed_distance_to_line(p1,p8, p4,p5) == CGAL::LARGER );
 assert( CGAL::compare_signed_distance_to_line(p8,p1, p4,p5) == CGAL::SMALLER );
 assert( CGAL::compare_signed_distance_to_line(p6,p4, p1,p2) == CGAL::EQUAL );
 assert( CGAL::compare_signed_distance_to_line(p6,p4, p1,p9) == CGAL::EQUAL );
 assert( CGAL::compare_signed_distance_to_line(p6,p8, p7,p9) == CGAL::SMALLER );
 assert( CGAL::compare_signed_distance_to_line(p6,p8, p2,p9) == CGAL::SMALLER );
 assert( CGAL::compare_signed_distance_to_line(p6,p8, p2,p3) == CGAL::LARGER );
 assert( CGAL::compare_signed_distance_to_line(p6,p1, p1,p6) == CGAL::EQUAL );
 assert( CGAL::compare_signed_distance_to_line(p6,p1, p6,p1) == CGAL::EQUAL );
 assert( CGAL::compare_signed_distance_to_line(p6,p1, p4,p1) == CGAL::SMALLER );
 assert( CGAL::compare_signed_distance_to_line(p6,p1, p5,p6) == CGAL::LARGER );
 assert( CGAL::compare_signed_distance_to_line(p4,p6, p8,p6) == CGAL::LARGER );
 assert( CGAL::compare_signed_distance_to_line(p6,p4, p8,p6) == CGAL::SMALLER );

 assert( CGAL::has_larger_signed_distance_to_line(p1,p2, p3,p8) );
 assert( CGAL::has_larger_signed_distance_to_line(p1,p8, p9,p8) );
 assert( CGAL::has_larger_signed_distance_to_line(p6,p1, p1,p4) );
 assert( CGAL::has_larger_signed_distance_to_line(p6,p4, p6,p8) );

 assert( CGAL::has_smaller_signed_distance_to_line(p1,p8, p5,p4) );
 assert( CGAL::has_smaller_signed_distance_to_line(p2,p1, p3,p8) );
 assert( CGAL::has_smaller_signed_distance_to_line(p6,p1, p6,p5) );
 assert( CGAL::has_smaller_signed_distance_to_line(p6,p8, p2,p9) );
 assert( CGAL::has_smaller_signed_distance_to_line(p6,p4, p8,p6) );

 std::cout << '.';

 assert( CGAL::area(p1, p2, p3) == FT(0) );
 assert( CGAL::area(p2, p3, p6) == FT(1) );
 assert( CGAL::area(p2, p6, p3) == FT(-1) );

 std::cout << '.';

 assert(CGAL::l_infinity_distance(p1,p2) == FT(1));
 assert(CGAL::l_infinity_distance(p1,p8) == FT(2));
 std::cout << "done" << std::endl;
 return true;
}

#endif // CGAL__TEST_FCT_POINT_2_H
