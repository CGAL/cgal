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


#ifndef CGAL__TEST_FCT_POINT_3_H
#define CGAL__TEST_FCT_POINT_3_H

// Accessory function testing functions that require sqrt().
// Doesn't instantiate anything if RT doesn't support sqrt().
template <class R>
bool
_test_fct_point_sqrt_3(const R&, CGAL::Tag_false)
{
//  bool UNTESTED_STUFF_BECAUSE_SQRT_IS_NOT_SUPPORTED;
  std::cout << std::endl
            << "NOTE : FT doesn't support sqrt(),"
               " hence some functions are not tested." << std::endl;
  return true;
}

template <class R>
bool
_test_fct_point_sqrt_3(const R&, CGAL::Tag_true)
{
 typedef typename  R::RT       RT;

 // unit_normal of three points
 CGAL::Point_3<R> pe0( RT(1), RT(0), RT(0) );
 CGAL::Point_3<R> pe1( RT(0), RT(1), RT(0) );
 CGAL::Vector_3<R> res( RT(0), RT(0), RT(1) );
 assert( CGAL::unit_normal( CGAL::Point_3<R>(CGAL::ORIGIN), pe0, pe1) \
                           == res);
 return true;
}

template <class R>
bool
_test_fct_point_3(const R& )
{
 std::cout << "Testing functions Point_3" ;

 typedef typename  R::RT    RT;
 typedef typename  R::FT    FT;

 CGAL::Point_3<R> p1(RT(18), RT(15), RT(-21), RT(3) ); // 6,5,-7
 CGAL::Point_3<R> p2(RT(18), RT(15), RT( 12), RT(3) ); // 6,5,4
 CGAL::Point_3<R> p3(RT(18), RT(12), RT(-21), RT(3) ); // 6,4,-7
 CGAL::Point_3<R> p4(RT(28), RT(40), RT( 20), RT(4) ); // 7,10,5
 CGAL::Point_3<R> p5(RT(12), RT(-4), RT(-20), RT(4) ); // 3,-1,-5

 assert( CGAL::compare_xyz(p1,p2) == CGAL::SMALLER );
 assert( CGAL::compare_xyz(p3,p2) == CGAL::SMALLER );
 assert( CGAL::compare_xyz(p3,p1) == CGAL::SMALLER );
 assert( CGAL::compare_xyz(p3,p2) == CGAL::SMALLER );
 assert( CGAL::compare_xyz(p2,p1) == CGAL::LARGER );
 assert( CGAL::compare_xyz(p2,p3) == CGAL::LARGER );
 assert( CGAL::compare_xyz(p4,p3) == CGAL::LARGER );
 assert( CGAL::compare_xyz(p4,p4) == CGAL::EQUAL );

 assert( CGAL::compare_lexicographically(p1,p2) == CGAL::SMALLER );
 assert( CGAL::compare_lexicographically(p3,p2) == CGAL::SMALLER );
 assert( CGAL::compare_lexicographically(p3,p1) == CGAL::SMALLER );
 assert( CGAL::compare_lexicographically(p2,p1) == CGAL::LARGER );
 assert( CGAL::compare_lexicographically(p2,p3) == CGAL::LARGER );
 assert( CGAL::compare_lexicographically(p4,p3) == CGAL::LARGER );
 assert( CGAL::compare_lexicographically(p4,p4) == CGAL::EQUAL );

 assert( CGAL::lexicographically_xyz_smaller_or_equal(p1,p1) );
 assert( CGAL::lexicographically_xyz_smaller_or_equal(p3,p1) );
 assert( CGAL::lexicographically_xyz_smaller_or_equal(p3,p2) );
 assert( CGAL::lexicographically_xyz_smaller_or_equal(p3,p4) );

 assert( !CGAL::lexicographically_xyz_smaller(p3,p3) );
 assert( CGAL::lexicographically_xyz_smaller(p3,p2) );
 assert( !CGAL::lexicographically_xyz_smaller(p4,p3) );

 assert( CGAL::compare_x(p2,p3) == CGAL::EQUAL );
 assert( CGAL::compare_x(p2,p4) == CGAL::SMALLER );
 assert( CGAL::compare_x(p4,p5) == CGAL::LARGER );
 assert( CGAL::compare_y(p2,p1) == CGAL::EQUAL );
 assert( CGAL::compare_y(p3,p2) == CGAL::SMALLER );
 assert( CGAL::compare_y(p4,p5) == CGAL::LARGER );
 assert( CGAL::compare_x(p1,p3) == CGAL::EQUAL );
 assert( CGAL::compare_x(p2,p4) == CGAL::SMALLER );
 assert( CGAL::compare_x(p4,p5) == CGAL::LARGER );

 assert( CGAL::x_equal(p1,p1) );
 assert( CGAL::x_equal(p2,p3) );
 assert( !CGAL::x_equal(p2,p4) );

 assert( CGAL::y_equal(p1,p2) );
 assert( !CGAL::y_equal(p1,p3) );

 assert( CGAL::z_equal(p1,p3) );
 assert( !CGAL::z_equal(p4,p5) );

 std::cout <<'.';

 CGAL::Point_3<R> p6 ( RT(6), RT(4), RT(7) );
 assert( CGAL::coplanar( p1, p2, p3, p6) );
 assert( CGAL::coplanar( p1, p1, p3, p4) );
 assert( CGAL::coplanar( p4, p1, p5, p5 + (p4-p1) ) );
 assert( !CGAL::coplanar( p4, p1, p2, p3 ) );

 assert( !CGAL::collinear( p1, p2, p3 ) );
 assert( CGAL::collinear( p1, p2, p2 + (p2-p1) ) );

 // ordered: arg1 - arg2 - arg3
 assert( CGAL::are_ordered_along_line( p1, p2, p2 + (p2-p1)) );
 assert( CGAL::are_ordered_along_line( p1, p2, p2) );
 assert( !CGAL::are_ordered_along_line( p1, p2 + (p2-p1), p2) );
 assert( !CGAL::are_ordered_along_line( p1, p5, p2 ) );
 assert( CGAL::are_ordered_along_line( p2, p2, p2) );

 assert( CGAL::collinear_are_ordered_along_line( p1, p2, p2 + (p2-p1)) );
 assert( !CGAL::collinear_are_ordered_along_line( p1, p2 + (p2-p1), p2) );

 assert( CGAL::collinear_are_ordered_along_line( p1, p1, p1));
 assert( !CGAL::collinear_are_ordered_along_line( p1, p4, p1));
 assert( !CGAL::collinear_are_ordered_along_line( p1, p3, p1));

 // strictly ordered: ordered && args pairwise distinct
 assert( CGAL::are_strictly_ordered_along_line( p1, p2, p2 + (p2-p1)) );
 assert( !CGAL::are_strictly_ordered_along_line( p1, p2, p2) );
 assert( !CGAL::are_strictly_ordered_along_line( p1, p2 + (p2-p1), p2) );
 assert( !CGAL::are_strictly_ordered_along_line( p1, p5, p2 ) );
 assert( !CGAL::are_strictly_ordered_along_line( p2, p2, p2) );

 assert( CGAL::collinear_are_strictly_ordered_along_line(p1, p2, p2 + (p2-p1)));
 assert( !CGAL::collinear_are_strictly_ordered_along_line(p1, p2, p2));
 assert( !CGAL::collinear_are_strictly_ordered_along_line(p1, p2 + (p2-p1), p2));
 assert( !CGAL::collinear_are_strictly_ordered_along_line(p1, p1, p1));
 assert( !CGAL::collinear_are_strictly_ordered_along_line(p1, p4, p1));

 assert( CGAL::collinear( p3, p2, p3 ) );
 assert( CGAL::collinear( p2, p2, p3 ) );
 assert( CGAL::collinear( p2, p3, p3 ) );

 std::cout << '.';

 CGAL::Point_3<R> pe0( RT(1), RT(0), RT(0) );
 CGAL::Point_3<R> pe1( RT(0), RT(1), RT(0) );
 CGAL::Point_3<R> pe2( RT(0), RT(0), RT(1) );

 assert( CGAL::orientation( CGAL::Point_3<R>(CGAL::ORIGIN), pe0, pe1, pe2 ) \
                           == CGAL::POSITIVE);

 assert( CGAL::orientation( p1, p2, p3, p6 ) == CGAL::ZERO );

 CGAL::Point_3<R> p7( RT(-8), RT(0), RT(0), RT(-2) );
 CGAL::Point_3<R> p8( RT(8), RT(4), RT(0), RT(2) );
 CGAL::Point_3<R> p9( RT(0), RT(12), RT(0), RT(4) );

 assert( CGAL::orientation( p7, p8, p9, p4) == CGAL::POSITIVE );
 assert( CGAL::orientation( p7, p9, p8, p5) == CGAL::POSITIVE );
 assert( CGAL::orientation( p7, p8, p9, p5) == CGAL::NEGATIVE );
 assert( CGAL::orientation( p8, p7, p9, p4) == CGAL::NEGATIVE );

 // normal of three points
 CGAL::Vector_3<R> res( RT(0), RT(0), RT(1) );
 assert( CGAL::normal( CGAL::Point_3<R>(CGAL::ORIGIN), pe0, pe1) \
                           == res);
 std::cout <<'.';

 CGAL::Point_3<R> p10( RT(0), RT(0), RT(16), RT(8) );

// CGAL::side_of_bounded_sphere()
 assert( CGAL::side_of_bounded_sphere(p7,p8,p9,p10,p1) ==CGAL::ON_UNBOUNDED_SIDE);
 assert( CGAL::side_of_bounded_sphere(p7,p9,p8,p10,p1) ==CGAL::ON_UNBOUNDED_SIDE);
 CGAL::Point_3<R> p0(CGAL::ORIGIN);
 assert( CGAL::side_of_bounded_sphere(p7,p8,p9,p10,p0) ==CGAL::ON_BOUNDED_SIDE);
 CGAL::Vector_3<R> v001( RT(0), RT(0), RT(1) );
 CGAL::Vector_3<R> v010( RT(0), RT(1), RT(0) );
 CGAL::Vector_3<R> v100( RT(1), RT(0), RT(0) );
 assert( CGAL::side_of_bounded_sphere(p3 + v001, p3-v001, p3+v010, p3-v100, \
                                      p3 - v010) == CGAL::ON_BOUNDARY );
// CGAL::side_of_bounded_sphere() is further tested in
// _test_fct_points_implicit_sphere(const R& )

 assert( CGAL::compare_distance_to_point(p3, p3 + v001, p3+v010) ==
                                                        CGAL::EQUAL );
 assert( CGAL::compare_distance_to_point(p0, p1, p2) == CGAL::LARGER );
 assert( CGAL::compare_distance_to_point(p0, p3, p1) == CGAL::SMALLER );
 assert( CGAL::compare_distance_to_point(p1, p3, p5) == CGAL::SMALLER );
 assert( CGAL::compare_distance(p3, p3 + v001, p3+v010) ==
                                                        CGAL::EQUAL );
 assert( CGAL::compare_distance(p0, p1, p2) == CGAL::LARGER );
 assert( CGAL::compare_distance(p0, p3, p1) == CGAL::SMALLER );
 assert( CGAL::compare_distance(p1, p3, p5) == CGAL::SMALLER );
 assert( CGAL::has_larger_distance_to_point(p0, p1, p2) );
 assert( CGAL::has_larger_distance_to_point(p3, p5, p1) );
 assert( CGAL::has_smaller_distance_to_point(p0, p2, p1) );
 assert( CGAL::has_smaller_distance_to_point(p3, p1, p5) );

 assert( CGAL::compare_distance(p3, p3 + v001, p3+v010,p3) ==
                                                        CGAL::EQUAL );
 assert( CGAL::compare_distance(p0, p1, p2,p0) == CGAL::LARGER );
 assert( CGAL::compare_distance(p0, p3, p1,p0) == CGAL::SMALLER );
 assert( CGAL::compare_distance(p1, p3, p5,p1) == CGAL::SMALLER );

 assert( CGAL::compare_squared_distance(p3, p3 + v001, CGAL::squared_distance(p3+v010,p3)) ==
                                                        CGAL::EQUAL );
 assert( CGAL::compare_squared_distance(p0, p1, CGAL::squared_distance(p2,p0)) == CGAL::LARGER );
 assert( CGAL::compare_squared_distance(p0, p3, CGAL::squared_distance(p1,p0)) == CGAL::SMALLER );
 assert( CGAL::compare_squared_distance(p1, p3, CGAL::squared_distance(p5,p1)) == CGAL::SMALLER );

 {
   CGAL::Point_3<R> p0(-2,0,0), p1(2,0,0), p2(0,2,0), p3(0,0,2);
   FT four(4);
   assert( CGAL::compare_squared_radius(p0, 0) == CGAL::EQUAL );
   assert( CGAL::compare_squared_radius(p0, -1) == CGAL::POSITIVE );
   assert( CGAL::compare_squared_radius(p0,  1) == CGAL::NEGATIVE );
   assert( CGAL::compare_squared_radius(p0, p1, four) == CGAL::EQUAL );
   assert( CGAL::compare_squared_radius(p0, p1, p2, four) == CGAL::EQUAL );
   assert( CGAL::compare_squared_radius(p0, p1, p2, p3, four) == CGAL::EQUAL );
 }

 {
   CGAL::Point_3<R> p0(0,0,1), p1(0,0,2), p2(1,0,1), p3(1,0,2), p4(1,0,3);
   assert( CGAL::compare_slope(p0, p2, p1, p3) == CGAL::EQUAL );
   assert( CGAL::compare_slope(p0, p2, p0, p3) == CGAL::SMALLER );
   assert( CGAL::compare_slope(p0, p2, p1, p2) == CGAL::LARGER );
   assert( CGAL::compare_slope(p0, p3, p0, p2) == CGAL::LARGER );
   assert( CGAL::compare_slope(p0, p3, p0, p4) == CGAL::SMALLER );
   assert( CGAL::compare_slope(p0, p3, p1, p2) == CGAL::LARGER );
   assert( CGAL::compare_slope(p1, p2, p0, p2) == CGAL::SMALLER );
   assert( CGAL::compare_slope(p1, p2, p0, p3) == CGAL::SMALLER );
   assert( CGAL::compare_slope(p1, p2, p4, p0) == CGAL::LARGER );

 }


 assert(CGAL::l_infinity_distance(p1,p2) == FT(11));
 assert(CGAL::l_infinity_distance(p1,p5) == FT(6));
 // More tests, that require sqrt().
 {
     typedef ::CGAL::Algebraic_structure_traits<FT> AST;
     static const bool has_sqrt =
         ! ::boost::is_same< ::CGAL::Null_functor, typename AST::Sqrt >::value;
     _test_fct_point_sqrt_3(R(), ::CGAL::Boolean_tag<has_sqrt>());
 }
 std::cout << "done" << std::endl;
 return true;
}
#endif // CGAL__TEST_FCT_POINT_3_H
