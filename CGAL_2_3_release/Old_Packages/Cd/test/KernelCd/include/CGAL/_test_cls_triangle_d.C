// revision      : 2.2.2
// revision_date : 28 Sep 1999 
// author(s)     : Stefan Schirra <Stefan.Schirra@mpi-sb.mpg.de>


#ifndef CGAL__TEST_CLS_TRIANGLE_D_C
#define CGAL__TEST_CLS_TRIANGLE_D_C

#include <CGAL/_test_cls_triangle_d.h>

template <class R>
bool
_test_cls_triangle_d(const R& )
{
 std::cout << "Testing class Triangle_d" ;

 typedef typename  R::RT    RT;
 typedef typename  R::FT    FT;

 typename R::Triangle_d it;
 CGAL::Triangle_d<R> t0(it);

 RT n0 =  0;
 RT n1 = 12;
 RT n2 = 16;
 RT n3 = -4;
 RT n4 =  2;
 RT n5 =  3;
 RT n6 = 30;
 RT n7 =  9;
 RT n8 = 24;
 RT n9 =  8;

 CGAL::Point_d<R> p1(3, n1, n2, n3, n4);
 CGAL::Point_d<R> p2(3, n2, n9, n3,-n3);
 CGAL::Point_d<R> p3(3, n5, n6, n1, n5);
 CGAL::Point_d<R> p4(3, n7, n7, n8, n5);

 CGAL::Point_d<R> ps3(3, n0, n0, n7, n5);
 CGAL::Point_d<R> ps2(3, n0, n7, n0, n5);
 CGAL::Point_d<R> ps1(3, n7, n0, n0, n5);

 CGAL::Triangle_d<R> t1(p1,p2,p3);
 CGAL::Triangle_d<R> t2(p4,p2,p3);
 CGAL::Triangle_d<R> t3(ps1,ps2,ps3);
 CGAL::Triangle_d<R> t4(ps2,ps1,ps3);
 CGAL::Triangle_d<R> t5( t1 );
 t0 = t1;

 assert( t0 == t0 );
 assert( t0 == t1 );
 assert( t5 == t1 );
 assert( t2 != t4 );
 assert( t3 != t4 );

 std::cout <<'.';

 CGAL::Plane_d<R> pl1(3, p1,p2,p3);
 CGAL::Plane_d<R> pl2(3, p4,p2,p3);
 assert( t1.supporting_plane() == pl1 );
 assert( t2.supporting_plane() == pl2 );
 assert( t3.supporting_plane() == t4.supporting_plane().opposite() );

 std::cout <<'.';

 assert( t1.has_on(p3) );
 assert( t1.has_on(p2) );
 assert( t2.has_on(p4) );
 assert( ! t1.has_on(p4) );
 CGAL::Point_d<R> pt(3, n7, n7, n7, n7);
 assert( t3.has_on( pt ) );
 assert( t4.has_on( pt ) );

 assert( t1.vertex(0) == p1 );
 assert( t1.vertex(1) == p2 );
 assert( t1.vertex(2) == p3 );
 assert( t4[0] == ps2 );
 assert( t4[1] == ps1 );
 assert( t4[2] == ps3 );

 std::cout <<'.';

 CGAL::Triangle_d<R> tdeg1( p3,p3,p1);
 CGAL::Triangle_d<R> tdeg2( p3,p3,p3);
 assert( tdeg1.is_degenerate() );
 assert( tdeg2.is_degenerate() );

 std::cout << "done" << std::endl;
 return true;
}

#endif // CGAL__TEST_CLS_TRIANGLE_D_C
