#ifndef CGAL__TEST_FCT_POINT_D_C
#define CGAL__TEST_FCT_POINT_D_C

#include <CGAL/_test_fct_point_d.h>

template <class R>
bool
_test_fct_point_d(const R& )
{
 std::cout << "Testing functions Point_d" ;

 typedef typename  R::RT    RT;
 typedef typename  R::FT    FT;

 CGAL::Point_d<R> p1(3, RT(18), RT(15), RT(-21), RT(3) ); // 6,5,-7
 CGAL::Point_d<R> p2(3, RT(18), RT(15), RT( 12), RT(3) ); // 6,5,4
 CGAL::Point_d<R> p3(3, RT(18), RT(12), RT(-21), RT(3) ); // 6,4,-7
 CGAL::Point_d<R> p4(3, RT(28), RT(40), RT( 20), RT(4) ); // 7,10,5
 CGAL::Point_d<R> p5(3, RT(12), RT(-4), RT(-20), RT(4) ); // 3,-1,-5

 assert( CGAL::compare_lexicographically_d(p1,p2) == CGAL::SMALLER );
 assert( CGAL::compare_lexicographically_d(p3,p2) == CGAL::SMALLER );
 assert( CGAL::compare_lexicographically_d(p3,p1) == CGAL::SMALLER );
 assert( CGAL::compare_lexicographically_d(p3,p2) == CGAL::SMALLER );
 assert( CGAL::compare_lexicographically_d(p2,p1) == CGAL::LARGER );
 assert( CGAL::compare_lexicographically_d(p2,p3) == CGAL::LARGER );
 assert( CGAL::compare_lexicographically_d(p4,p3) == CGAL::LARGER );
 assert( CGAL::compare_lexicographically_d(p4,p4) == CGAL::EQUAL );

 assert( CGAL::lexicographically_d_smaller_or_equal(p1,p1) );
 assert( CGAL::lexicographically_d_smaller_or_equal(p3,p1) );
 assert( CGAL::lexicographically_d_smaller_or_equal(p3,p2) );
 assert( CGAL::lexicographically_d_smaller_or_equal(p3,p4) );

 assert( !CGAL::lexicographically_d_smaller(p3,p3) );
 assert( CGAL::lexicographically_d_smaller(p3,p2) );
 assert( !CGAL::lexicographically_d_smaller(p4,p3) );

 assert( CGAL::compare_x(p2,p3) == CGAL::EQUAL );
 assert( CGAL::compare_x(p2,p4) == CGAL::SMALLER );
 assert( CGAL::compare_x(p4,p5) == CGAL::LARGER );
 assert( CGAL::compare_x(p2,p1,1) == CGAL::EQUAL );
 assert( CGAL::compare_x(p3,p2,1) == CGAL::SMALLER );
 assert( CGAL::compare_x(p4,p5,1) == CGAL::LARGER );
 assert( CGAL::compare_x(p1,p3) == CGAL::EQUAL );
 assert( CGAL::compare_x(p2,p4) == CGAL::SMALLER );
 assert( CGAL::compare_x(p4,p5) == CGAL::LARGER );

 assert( CGAL::x_equal(p1,p1) );
 assert( CGAL::x_equal(p2,p3) );
 assert( !CGAL::x_equal(p2,p4) );

 assert( CGAL::x_equal(p1,p2,1) );
 assert( !CGAL::x_equal(p1,p3,1) );

 assert( CGAL::x_equal(p1,p3,2) );
 assert( !CGAL::x_equal(p4,p5,2) );

 std::cout <<'.';

 CGAL::Point_d<R> p6 (3, RT(6), RT(4), RT(7) );
 CGAL::Point_d<R> h1[4] = { p1, p2, p3, p6 };
 assert( CGAL::coplanar( h1+0, h1+4 ) );
 CGAL::Point_d<R> h2[4] = { p1, p1, p3, p4 };
 assert( CGAL::coplanar( h2+0, h2+4 ) );
 CGAL::Point_d<R> h3[4] = { p4, p1, p5, p5 + (p4-p1) };
 assert( CGAL::coplanar( h3+0, h3+4 ) );
 CGAL::Point_d<R> h4[4] = { p4, p1, p2, p3 };
 assert( !CGAL::coplanar( h4+0, h4+4 ) );

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

 CGAL::Point_d<R> pe0(3, RT(1), RT(0), RT(0) );
 CGAL::Point_d<R> pe1(3, RT(0), RT(1), RT(0) );
 CGAL::Point_d<R> pe2(3, RT(0), RT(0), RT(1) );

 CGAL::Point_d<R> h5[4] = { CGAL::Point_d<R>(3, CGAL::ORIGIN), pe0, pe1, pe2 };
 assert( CGAL::orientation( h5+0, h5+4 ) \
                           == CGAL::POSITIVE);

 CGAL::Point_d<R> h6[4] = { p1, p2, p3, p6 };
 assert( CGAL::orientation( h6+0, h6+4 ) == CGAL::ZERO );

 CGAL::Point_d<R> p7(3, RT(-8), RT(0), RT(0), RT(-2) );
 CGAL::Point_d<R> p8(3, RT(8), RT(4), RT(0), RT(2) );
 CGAL::Point_d<R> p9(3, RT(0), RT(12), RT(0), RT(4) );

 CGAL::Point_d<R> h7[4] = { p7, p8, p9, p4 };
 assert( CGAL::orientation( h7+0, h7+4 ) == CGAL::POSITIVE );
 CGAL::Point_d<R> h8[4] = { p7, p9, p8, p5 };
 assert( CGAL::orientation( h8+0, h8+4 ) == CGAL::POSITIVE );
 CGAL::Point_d<R> h9[4] = { p7, p8, p9, p5 };
 assert( CGAL::orientation( h9+0, h9+4 ) == CGAL::NEGATIVE );
 CGAL::Point_d<R> h10[4] = { p8, p7, p9, p4 };
 assert( CGAL::orientation( h10+0, h10+4 ) == CGAL::NEGATIVE );

 std::cout <<'.';

 CGAL::Point_d<R> p10(3, RT(0), RT(0), RT(16), RT(8) );

// CGAL::side_of_bounded_sphere()
 CGAL::Point_d<R> h11[4] = { p7,p8,p9,p10 };
 assert( CGAL::side_of_bounded_sphere(h11+0,h11+4,p1)==CGAL::ON_UNBOUNDED_SIDE);
 CGAL::Point_d<R> h12[4] = { p7,p9,p8,p10 };
 assert( CGAL::side_of_bounded_sphere(h12+0,h12+4,p1)==CGAL::ON_UNBOUNDED_SIDE);
 CGAL::Point_d<R> p0(3, CGAL::ORIGIN);
 assert( CGAL::side_of_bounded_sphere(h11+0,h11+4,p0)==CGAL::ON_BOUNDED_SIDE);
 CGAL::Vector_d<R> v001(3, RT(0), RT(0), RT(1) );
 CGAL::Vector_d<R> v010(3, RT(0), RT(1), RT(0) );
 CGAL::Vector_d<R> v100(3, RT(1), RT(0), RT(0) );
 CGAL::Point_d<R> h14[4] = { p3 + v001, p3-v001, p3+v010, p3-v100 };
 assert( CGAL::side_of_bounded_sphere(h14+0,h14+4, \
                                      p3 - v010) == CGAL::ON_BOUNDARY );
// CGAL::side_of_bounded_sphere() is further tested in 
// _test_fct_points_implicit_sphere(const R& )

 std::cout <<'.';
 
 assert( CGAL::cmp_dist_to_point(p3, p3 + v001, p3+v010) == CGAL::EQUAL );
 assert( CGAL::cmp_dist_to_point(p0, p1, p2) == CGAL::LARGER );
 assert( CGAL::cmp_dist_to_point(p0, p3, p1) == CGAL::SMALLER );
 assert( CGAL::cmp_dist_to_point(p1, p3, p5) == CGAL::SMALLER );
 assert( CGAL::has_larger_dist_to_point(p0, p1, p2) );
 assert( CGAL::has_larger_dist_to_point(p3, p5, p1) );
 assert( CGAL::has_smaller_dist_to_point(p0, p2, p1) );
 assert( CGAL::has_smaller_dist_to_point(p3, p1, p5) );
 
 std::cout << "done" << std::endl;
 return true;
}
#endif // CGAL__TEST_FCT_POINT_D_C
