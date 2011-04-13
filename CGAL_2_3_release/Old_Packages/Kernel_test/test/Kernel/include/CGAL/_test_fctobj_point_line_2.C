// ============================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
// 
// release       :
// release_date  :
// 
// file          : _test_fctobj_point_line_2.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#ifndef CGAL__TEST_FCTOBJ_POINT_LINE_2_C
#define CGAL__TEST_FCTOBJ_POINT_LINE_2_C

#include <CGAL/predicate_classes_2.h>


template <class Point, class Line>
bool
_test_fctobj_point_line_2(const Point&, const Line& );


template <class Point, class Line>
bool
_test_fctobj_point_line_2(const Point&, const Line& )
{
 std::cout << "Testing function objects Point_2 and Line_2";

 typedef typename  Point::R::RT    RT;

  Point p1( RT(18), RT(12), RT(3) );  // ( 6, 4)
  Point p2( RT(18), RT(15), RT(3) );  // ( 6, 5)
  Point p3( RT(18), RT( 9), RT(3) );  // ( 6, 3)
  Point p4( RT(28), RT(40), RT(4) );  // ( 7,10)
  Point p5( RT(12), RT(-4), RT(4) );  // ( 3,-1)
  Point p6( RT(28), RT(12), RT(4) );  // ( 7, 3)
  Point p7( RT(18), RT( 6), RT(3) );  // ( 6, 2)
  Point p8( RT(24), RT( 9), RT(3) );  // ( 8, 3)
  Point p9( RT( 6), RT(10), RT(1) );  // ( 6,10)
 
  CGAL::Compare_xy_2< Point> compare_xy;
  CGAL::Less_xy_2< Point>    lexicographically_xy_smaller;
 
  assert( compare_xy(p1,p2) == CGAL::SMALLER );
  assert( compare_xy(p3,p2) == CGAL::SMALLER );
  assert( compare_xy(p3,p1) == CGAL::SMALLER );
  assert( compare_xy(p3,p2) == CGAL::SMALLER );
  assert( compare_xy(p2,p1) == CGAL::LARGER );
  assert( compare_xy(p2,p3) == CGAL::LARGER );
  assert( compare_xy(p4,p3) == CGAL::LARGER );
  assert( compare_xy(p4,p4) == CGAL::EQUAL );
 
  assert( !lexicographically_xy_smaller(p3,p3) );
  assert(  lexicographically_xy_smaller(p3,p2) );
  assert( !lexicographically_xy_smaller(p4,p3) );
 
  std::cout <<'.';
  CGAL::Compare_yx_2< Point> compare_yx;
  CGAL::Less_yx_2< Point>    lexicographically_yx_smaller;
 
  assert( compare_yx(p1,p2) == CGAL::SMALLER );
  assert( compare_yx(p2,p1) == CGAL::LARGER );
  assert( compare_yx(p2,p2) == CGAL::EQUAL );
  assert( compare_yx(p2,p4) == CGAL::SMALLER );
  assert( compare_yx(p3,p6) == CGAL::SMALLER );
  assert( compare_yx(p6,p3) == CGAL::LARGER );
  assert( compare_yx(p4,p9) == CGAL::LARGER );
 
  assert( lexicographically_yx_smaller(p1,p2) );
  assert( lexicographically_yx_smaller(p2,p4) );
  assert( lexicographically_yx_smaller(p3,p6) );
  assert( lexicographically_yx_smaller(p2,p4) );
  assert( lexicographically_yx_smaller(p8,p1) );
  assert(!lexicographically_yx_smaller(p1,p8) );
 
  std::cout << '.';
 
  CGAL::Orientation_2< Point>   orientation;
  CGAL::Collinear_2< Point>     collinear;
  CGAL::Left_turn_2< Point>     left_turn;
  CGAL::Right_turn_2< Point>    right_turn;
 
  Point pe0( RT(1), RT(0) );
  Point pe1( RT(0), RT(1) );
 
  assert( orientation( Point(CGAL::ORIGIN), pe0, pe1)== CGAL::POSITIVE);
 
  assert( orientation( p1, p2, p3) == CGAL::COLLINEAR );
  assert( orientation( p1, p2, p4) == CGAL::RIGHT_TURN );
  assert( orientation( p2, p1, p4) == CGAL::LEFT_TURN );
  assert( orientation( p5, p4, p3) == CGAL::RIGHT_TURN );
  assert( orientation( p2, p4, p6) == CGAL::RIGHT_TURN );
  assert( orientation( p6, p4, p2) == CGAL::LEFT_TURN );
  assert( orientation( p4, p6, p2) == CGAL::RIGHT_TURN );
  assert( orientation( p5, p6, p7) == CGAL::COLLINEAR );
  assert( orientation( p6, p5, p7) == CGAL::COLLINEAR );
 
  assert( collinear( p1, p2, p3 ) );
  assert( collinear( p1, p2, p7 ) );
  assert( collinear( p6, p5, p7 ) );
  assert( collinear( p1, p2, p3 ) );
  assert( !collinear( p1, p2, p4 ) );
  assert( collinear( p6, p6, p3 ) );
 
  assert( left_turn( p1, p4, p2 ) );
  assert( left_turn( p6, p4, p2 ) );
 
  assert( right_turn( p4, p6, p2 ) );
  assert( right_turn( p1, p2, p4 ) );
 
  std::cout << '.';
 
  CGAL::Are_ordered_along_line_2< Point>
       are_ordered_along_line;
 
  CGAL::Collinear_are_ordered_along_line_2< Point>
       collinear_are_ordered_along_line;
 
  CGAL::Are_strictly_ordered_along_line_2< Point>
       are_strictly_ordered_along_line;
 
  CGAL::Collinear_are_strictly_ordered_along_line_2< Point>
       collinear_are_strictly_ordered_along_line;
 
  assert( are_ordered_along_line( p5, p7, p6 ) );   // p7 between p6 and p5
  assert( are_ordered_along_line( p6, p7, p5 ) );
  assert( are_ordered_along_line( p2, p1, p3 ) );
  assert( !are_ordered_along_line( p7, p6, p5 ) );
  assert( !are_ordered_along_line( p7, p5, p6 ) );
  assert( !are_ordered_along_line( p7, p4, p6 ) );
  assert( !are_ordered_along_line( p2, p4, p6 ) );
  assert( collinear_are_ordered_along_line( p5, p7, p6 ) );
  assert( collinear_are_ordered_along_line( p6, p7, p5 ) );
  assert( collinear_are_ordered_along_line( p2, p1, p3 ) );
  assert( !collinear_are_ordered_along_line( p7, p6, p5 ) );
 
  assert( collinear_are_ordered_along_line( p7, p7, p7 ) );
  assert( !collinear_are_ordered_along_line( p5, p6, p5 ) );
  assert( !collinear_are_ordered_along_line( p1, p3, p1 ) );
  assert( !collinear_are_ordered_along_line( p3, p6, p3 ) );
 
  assert( are_strictly_ordered_along_line( p5, p7, p6 ) );
  assert( are_strictly_ordered_along_line( p6, p7, p5 ) );
  assert( are_strictly_ordered_along_line( p2, p1, p3 ) );
  assert( !are_strictly_ordered_along_line( p7, p6, p5 ) );
  assert( !are_strictly_ordered_along_line( p7, p5, p6 ) );
  assert( !are_strictly_ordered_along_line( p7, p4, p6 ) );
  assert( !are_strictly_ordered_along_line( p2, p4, p6 ) );
  assert( collinear_are_strictly_ordered_along_line( p5, p7, p6 ) );
  assert( collinear_are_strictly_ordered_along_line( p6, p7, p5 ) );
  assert( collinear_are_strictly_ordered_along_line( p2, p1, p3 ) );
  assert( !collinear_are_strictly_ordered_along_line( p7, p6, p5 ) );
 
  assert( !are_strictly_ordered_along_line( p7, p7, p6 ) );
  assert( !collinear_are_strictly_ordered_along_line( p5, p6, p6 ) );
  assert( are_ordered_along_line( p7, p7, p6 ) );
  assert( collinear_are_ordered_along_line( p5, p6, p6 ) );
  assert( are_ordered_along_line( p6, p6, p6 ) );
  assert( collinear_are_ordered_along_line( p5, p5, p5 ) );
  assert( !are_strictly_ordered_along_line( p6, p6, p6 ) );
  assert( !collinear_are_strictly_ordered_along_line( p5, p5, p5 ) );
 
  std::cout << '.';
 
  CGAL::Less_distance_to_point_2< Point>     has_smaller_distance_to_point;
  assert( has_smaller_distance_to_point(p4, p9,p1) );
  assert( has_smaller_distance_to_point(p4, p6,p5) );
  assert( has_smaller_distance_to_point(p9, p9,p6) );
  assert(!has_smaller_distance_to_point(p8, p3,p3) );
 
  std::cout << '.';
 
 
  CGAL::Compare_signed_distance_to_implicit_line_2< Point>
       compare_signed_distance_to_line_p1_p2(p1,p2);
  assert( compare_signed_distance_to_line_p1_p2( p3,p9) == CGAL::EQUAL );
  assert( compare_signed_distance_to_line_p1_p2( p3,p8) == CGAL::LARGER );
 
  CGAL::Compare_signed_distance_to_implicit_line_2< Point>
       compare_signed_distance_to_line_p2_p1(p2,p1);
  assert( compare_signed_distance_to_line_p2_p1( p3,p8) == CGAL::SMALLER );
 
  CGAL::Compare_signed_distance_to_implicit_line_2< Point>
       compare_signed_distance_to_line_p1_p8(p1,p8);
  assert( compare_signed_distance_to_line_p1_p8( p9,p8) == CGAL::LARGER );
  assert( compare_signed_distance_to_line_p1_p8( p9,p9) == CGAL::EQUAL );
  assert( compare_signed_distance_to_line_p1_p8( p5,p4) == CGAL::SMALLER );
  assert( compare_signed_distance_to_line_p1_p8( p4,p5) == CGAL::LARGER );
 
  CGAL::Compare_signed_distance_to_implicit_line_2< Point>
       compare_signed_distance_to_line_p8_p1(p8,p1);
  assert( compare_signed_distance_to_line_p8_p1( p4,p5) == CGAL::SMALLER );
 
  CGAL::Compare_signed_distance_to_implicit_line_2< Point>
       compare_signed_distance_to_line_p6_p4(p6,p4);
  assert( compare_signed_distance_to_line_p6_p4( p1,p2) == CGAL::EQUAL );
  assert( compare_signed_distance_to_line_p6_p4( p1,p9) == CGAL::EQUAL );
  assert( compare_signed_distance_to_line_p6_p4( p8,p6) == CGAL::SMALLER );
 
  CGAL::Compare_signed_distance_to_implicit_line_2< Point>
       compare_signed_distance_to_line_p6_p8(p6,p8);
  assert( compare_signed_distance_to_line_p6_p8( p7,p9) == CGAL::SMALLER );
  assert( compare_signed_distance_to_line_p6_p8( p2,p9) == CGAL::SMALLER );
  assert( compare_signed_distance_to_line_p6_p8( p2,p3) == CGAL::LARGER );
 
  CGAL::Compare_signed_distance_to_implicit_line_2< Point>
       compare_signed_distance_to_line_p6_p1(p6,p1);
  assert( compare_signed_distance_to_line_p6_p1( p1,p6) == CGAL::EQUAL );
  assert( compare_signed_distance_to_line_p6_p1( p6,p1) == CGAL::EQUAL );
  assert( compare_signed_distance_to_line_p6_p1( p4,p1) == CGAL::SMALLER );
  assert( compare_signed_distance_to_line_p6_p1( p5,p6) == CGAL::LARGER );
 
  CGAL::Compare_signed_distance_to_implicit_line_2< Point>
       compare_signed_distance_to_line_p4_p6(p4,p6);
  assert( compare_signed_distance_to_line_p4_p6( p8,p6) == CGAL::LARGER );
 
 
  CGAL::Less_signed_distance_to_implicit_line_2< Point>
       has_smaller_signed_distance_to_line_p1_p8(p1,p8);
  assert( has_smaller_signed_distance_to_line_p1_p8( p5,p4) );
 
  CGAL::Less_signed_distance_to_implicit_line_2< Point>
       has_smaller_signed_distance_to_line_p2_p1(p2,p1);
  assert( has_smaller_signed_distance_to_line_p2_p1( p3,p8) );
 
  CGAL::Less_signed_distance_to_implicit_line_2< Point>
       has_smaller_signed_distance_to_line_p6_p1(p6,p1);
  assert( has_smaller_signed_distance_to_line_p6_p1( p6,p5) );
 
  CGAL::Less_signed_distance_to_implicit_line_2< Point>
       has_smaller_signed_distance_to_line_p6_p8(p6,p8);
  assert( has_smaller_signed_distance_to_line_p6_p8( p2,p9) );
 
  CGAL::Less_signed_distance_to_implicit_line_2< Point>
       has_smaller_signed_distance_to_line_p6_p4(p6,p4);
  assert( has_smaller_signed_distance_to_line_p6_p4( p8,p6) );
 
  std::cout << '.';
 
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
 
  Point  q1( n2, n8, n2);   // ( 1, 4)
  Point  q2( n4, n4, n2);   // ( 2, 2)
  Point  q3(n12,n10, n2);   // ( 6, 5)
  Point  q4( n7, n3);       // ( 7, 3)
  Point  q5( n9, n0, n3);   // ( 3, 0)
  Point  q6(n12,n20, n4);   // ( 3, 5)
  Point  q7(n12, n3, n3);   // ( 4, 1)
  Point  q8( n8, n4);       // ( 8, 4)
  Point  q9(-n4, n8, n4);   // (-1, 2)
  Point  q10(n10, n4, n2);  // ( 5, 2)
 
  Point  q11( n1, -n5,-n1); // (-1, 5)
  Point  q0( CGAL::ORIGIN ); // ( 0, 0)
  Point  q13( n8, n8, n4);  // ( 4, 4)
  Point  q14( n5,-n1);      // ( 5,-1)
  Point  q16(n12, n9, n3);  // ( 4, 3)
  Point  q17( n0,n1);       // ( 0, 1)
 
  Line   l14( q1, q4 );
  Line   l23( q2, q3 );
  Line   l67( q6, q7 );
  Line   l58( q5, q8 );
  Line   l1617( q16, q17);
  Line   l1114( q11, q14);
  Line   l1716( q17, q16);
  Line   l1411( q14, q11);
  Line   l013( q0, q13 );
  Line   l910( q9, q10 );
 
  std::cout << '.';
 
  CGAL::Compare_x_2< Point>                       compare_x;
  assert( compare_x( q16, q14 ) == CGAL::SMALLER );
 
  CGAL::Compare_x_implicit_point_2< Point, Line>  compare_xip;
  assert( compare_xip( q9, l14, l23) == CGAL::SMALLER );
  assert( compare_xip( q8, l14, l23) == CGAL::LARGER );
  assert( compare_xip( q2, l1617, l910) == CGAL::EQUAL );
  assert( compare_xip( q2, l1716, l910) == CGAL::EQUAL );
  assert( compare_xip( q2, l1114, l013) == CGAL::EQUAL );
  assert( compare_xip( q2, l1411, l013) == CGAL::EQUAL );
  assert( compare_xip( q2, l1716, l013) == CGAL::EQUAL );
 
  CGAL::Compare_y_implicit_point_2< Point, Line>  compare_yip;
  assert( compare_yip( q6, l14, l23 ) == CGAL::LARGER );
  assert( compare_yip( q9, l14, l23 ) == CGAL::SMALLER );
 
  std::cout << '.';
 
  CGAL::Compare_x_implicit_points_same_line_2< Line>
                                                 compare_xip_sl;
  assert( compare_xip_sl( l14, l23, l58 ) == CGAL::SMALLER);
  assert( compare_xip_sl( l14, l58, l23 ) == CGAL::LARGER);
  assert( compare_xip_sl( l14, l58, l58 ) == CGAL::EQUAL);
  assert( compare_xip_sl( l1114, l013, l910 ) == CGAL::EQUAL);
  assert( compare_xip_sl( l1617, l910, l013 ) == CGAL::EQUAL);
 
  CGAL::Compare_y_implicit_points_same_line_2< Line>
                                                 compare_yip_sl;
  assert( compare_yip_sl( l14, l58, l23 ) == CGAL::SMALLER);
  assert( compare_yip_sl( l14, l23, l58 ) == CGAL::LARGER);
  assert( compare_yip_sl( l14, l58, l58 ) == CGAL::EQUAL);
  assert( compare_yip_sl( l1114, l013, l910 ) == CGAL::EQUAL);
  assert( compare_yip_sl( l1617, l910, l013 ) == CGAL::EQUAL);
 
  CGAL::Compare_x_implicit_points_2< Line>       compare_xips;
  assert( compare_xips( l14, l23, l67, l58 ) == CGAL::SMALLER);
  assert( compare_xips( l67, l58, l23, l14 ) == CGAL::LARGER);
  assert( compare_xips( l1114, l1617, l910, l013 ) == CGAL::EQUAL);
 
  CGAL::Compare_y_implicit_points_2< Line>       compare_yips;
  assert( compare_yips( l14, l23, l67, l58 ) == CGAL::LARGER);
  assert( compare_yips( l67, l58, l23, l14 ) == CGAL::SMALLER);
  assert( compare_yips( l1114, l1617, l910, l013 ) == CGAL::EQUAL);
 
  std::cout << '.';
 
  CGAL::Compare_y_at_point_2< Point, Line>       compare_y_at_x;
  assert( compare_y_at_x( q6, l23 ) == CGAL::LARGER );
  assert( compare_y_at_x( q6, l23.opposite() ) == CGAL::LARGER );
  assert( compare_y_at_x( q10, l23 ) == CGAL::SMALLER );
  assert( compare_y_at_x( q9, l23 ) == CGAL::LARGER );
  assert( compare_y_at_x( q17, l910 ) == CGAL::SMALLER );
  assert( compare_y_at_x( q0, l23 ) == CGAL::SMALLER );
  assert( compare_y_at_x( q8, l58 ) == CGAL::EQUAL );
  assert( compare_y_at_x( q2, l1617 ) == CGAL::EQUAL );
 
  CGAL::Compare_y_at_implicit_point_2< Line>     compare_y_at_ip;
  assert( compare_y_at_ip( l14, l23, l58 ) == CGAL::LARGER );
  assert( compare_y_at_ip( l67, l58, l23 ) == CGAL::SMALLER );
  assert( compare_y_at_ip( l910, l1716, l1114) == CGAL::EQUAL);
 
  CGAL::Compare_y_of_lines_at_implicit_point_2< Line>
                                                compare_y_of_l_ip;
  assert( compare_y_of_l_ip( l14, l23, l58, l67 ) == CGAL::SMALLER );
  assert( compare_y_of_l_ip( l14, l23, l67, l58 ) == CGAL::LARGER );
  assert( compare_y_of_l_ip( l14, l23, l1411, l1114 ) == CGAL::EQUAL );
  assert( compare_y_of_l_ip( l1617, l013, l910, l67 ) == CGAL::SMALLER);
  assert( compare_y_of_l_ip( l1617, l013, l67, l910 ) == CGAL::LARGER);
  assert( compare_y_of_l_ip( l1617, l013, l1114, l910 ) == CGAL::EQUAL);
  assert( compare_y_of_l_ip( l1617, l013, l910, l1114 ) == CGAL::EQUAL);
 
 

 std::cout << "done" << std::endl;
 return true;
}
#endif // _TEST_FCTOBJ_POINT_LINE_2_C
