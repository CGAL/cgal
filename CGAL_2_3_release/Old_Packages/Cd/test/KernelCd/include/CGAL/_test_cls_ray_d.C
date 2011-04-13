// revision      : 2.2.2
// revision_date : 28 Sep 1999 
// author(s)     : Stefan Schirra <Stefan.Schirra@mpi-sb.mpg.de>

#ifndef CGAL__TEST_CLS_RAY_D_C
#define CGAL__TEST_CLS_RAY_D_C

#include <CGAL/_test_cls_ray_d.h>

template <class R>
bool
_test_cls_ray_d(const R& )
{
 std::cout << "Testing class Ray_d" ;

 typedef typename  R::RT    RT;
 typedef typename  R::FT    FT;

 typename R::Ray_d ir;
 CGAL::Ray_d<R>  r1(ir);

 RT  n1 =  8;
 RT  n2 = 20;
 RT  n3 =  4;
 RT  n4 = 10;
 RT  n5 =  5;
 RT  n6 = 20;
 RT  n7 = -2;
 RT  n8 =  3;

 CGAL::Point_d<R> p1(3, n1, n2, n3, n3);
 CGAL::Point_d<R> p2(3, n4, n5, n6, n5);
 CGAL::Point_d<R> p3(3, n7, n2, n4, n7);

 CGAL::Ray_d<R> r2( p1, p2 );
 CGAL::Ray_d<R> r3( p2, p1 );
 CGAL::Ray_d<R> r4( r2 );
 r1 = r4;
 CGAL::Direction_d<R> dir( p2 - p1 );
 CGAL::Ray_d<R> r7(p1, dir);

 assert( r1 == r1 );
 assert( r4 == r2 );
 assert( r1 == r4 );
 assert( r1 == r2 );
 assert( r7 == r2 );
 assert( r2 != r3 );

 std::cout <<'.';

 CGAL::Ray_d<R> r5 (p3, p3 + (p1 - p3) );
 assert( r5.has_on( p1 ) );
 assert( r5.has_on( p3 ) );
 assert( r5.has_on( p3 + (p1 - p3) ) );
 assert( r3.has_on( p2 + (p1 - p2) + (p1 - p2) ) );
 assert( r2.has_on( r2.second_point() ));
 assert( r5.has_on( r5.second_point() ));
 assert( r4.has_on( r4.point(1) ));
 assert( r4.has_on( r4.point(3) ));

 std::cout <<'.';

 assert( r5.source() == p3 );
 assert( r2.source() != r3.source() );
 assert( r7.direction() == dir );
 assert( r2.direction() == CGAL::Direction_d<R>(r2.point(2) - r2.point(1) ));
 assert( r2.direction() == r3.opposite().direction() );
 assert( r1.supporting_line() == r2.supporting_line() );
 CGAL::Line_d<R> lin(p1,p2);
 assert( r2.supporting_line() == lin );

 std::cout << '.';

 CGAL::Ray_d<R> r8( p3, dir );
 CGAL::Ray_d<R> r9( p3, -dir );
 assert( r8.opposite() == r9 );
 assert( r9.opposite() == r8 );
 CGAL::Ray_d<R> sdeg(p3,p3);
 assert( sdeg.is_degenerate() );

 std::cout << "done" << std::endl;
 return true;
}

#endif // CGAL__TEST_CLS_RAY_d_C
