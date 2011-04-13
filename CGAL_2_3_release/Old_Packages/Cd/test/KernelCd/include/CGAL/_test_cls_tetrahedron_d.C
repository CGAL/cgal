// revision      : 2.2.2
// revision_date : 28 Sep 1999 
// author(s)     : Stefan Schirra <Stefan.Schirra@mpi-sb.mpg.de>

#ifndef CGAL__TEST_CLS_TETRAHEDRON_D_C
#define CGAL__TEST_CLS_TETRAHEDRON_D_C
#include <CGAL/_test_cls_tetrahedron_d.h>

template <class R>
bool
_test_cls_tetrahedron_d(const R& )
{
 std::cout << "Testing class Tetrahedron_d" ;

 typedef typename  R::RT    RT;
 typedef typename  R::FT    FT;

 typename R::Tetrahedron_d it;
 CGAL::Tetrahedron_d<R>  t0(it);

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
 CGAL::Point_d<R> ps0(3, CGAL::ORIGIN );

 CGAL::Tetrahedron_d<R> t1(p1,p2,p3,p4);
 CGAL::Tetrahedron_d<R> t2(p2,p1,p3,p4);
 CGAL::Tetrahedron_d<R> t3(ps0,ps1,ps2,ps3); // positive oriented
 CGAL::Tetrahedron_d<R> t4(ps0,ps1,ps3,ps2); // negative oriented
 CGAL::Tetrahedron_d<R> t5(ps0,p1,p3,ps2);
 CGAL::Tetrahedron_d<R> t6(t3);
 CGAL::Tetrahedron_d<R> the(p2,p3,p1,p4);
 CGAL::Tetrahedron_d<R> td1(p2,p3,p3,p4);
 CGAL::Tetrahedron_d<R> td2(p2,p2,p3,p4);
 t0 = t4;

 assert( t4 == t4 );
 assert( t4 == t0 );
 assert( t6 == t3 );
 assert( t4 != t2 );
 assert( t4 != t3 );
 assert( t1 == the );
 assert( td1 == td2 );


 std::cout << '.';

 assert( t5.vertex(0) == ps0 );
 assert( t5.vertex(1) == p1 );
 assert( t5.vertex(2) == p3 );
 assert( t5.vertex(3) == ps2 );
 assert( t5.vertex(4) == ps0 );
 assert( t5.vertex(5) == p1 );
 assert( t1[0] == p1 );
 assert( t1[1] == p2 );
 assert( t1[2] == p3 );
 assert( t1[3] == p4 );
 assert( t1[4] == p1 );
 assert( t1[9] == p2 );

 assert( t3.orientation() == CGAL::POSITIVE );
 assert( t4.orientation() == CGAL::NEGATIVE );
 CGAL::Tetrahedron_d<R> t7( p1,p1,p2,p3);
 CGAL::Tetrahedron_d<R> t8( p2,p2,p2,p2);
 assert( t7.is_degenerate() );
 assert( t8.is_degenerate() );

 std::cout << '.';

 CGAL::Point_d<R> p5(3,n6,n6,n6,n4);
 CGAL::Point_d<R> p6(3, n4,n4,n4,n9);
 CGAL::Point_d<R> p7(3, n7,n7,n7,n7);
 assert( t3.has_on_unbounded_side( p5 ));
 assert( t3.has_on_bounded_side( p6 ));
 assert( t3.has_on_boundary( p7 ));
 assert( t4.has_on_unbounded_side( p5 ));
 assert( t4.has_on_bounded_side( p6 ));
 assert( t4.has_on_boundary( p7 ));
 assert( t2.has_on_unbounded_side( p5 ));
 assert( t4.bounded_side( p5 ) == CGAL::ON_UNBOUNDED_SIDE );
 assert( t4.bounded_side( p6 ) == CGAL::ON_BOUNDED_SIDE );
 assert( t4.bounded_side( p7 ) == CGAL::ON_BOUNDARY );

 std::cout << '.';

 assert( t3.oriented_side( p3 ) != t4.oriented_side( p3 ) );
 assert( t4.has_on_positive_side( p5 ));
 assert( t4.has_on_negative_side( p6 ));
 assert( t3.has_on_positive_side( p6 ));
 assert( t3.has_on_negative_side( p5 ));
 assert( t4.oriented_side( p5 ) == CGAL::ON_POSITIVE_SIDE );
 assert( t4.oriented_side( p6 ) == CGAL::ON_NEGATIVE_SIDE );
 assert( t4.oriented_side( p7 ) == CGAL::ON_ORIENTED_BOUNDARY );
 assert( t3.oriented_side( p7 ) == CGAL::ON_ORIENTED_BOUNDARY );
 assert( t3.oriented_side( p6 ) == CGAL::ON_POSITIVE_SIDE );
 assert( t3.oriented_side( p5 ) == CGAL::ON_NEGATIVE_SIDE );
 assert( t2.has_on_boundary(p2) );
 assert( t2.bounded_side(p3) == CGAL::ON_BOUNDARY );
 assert( t2.oriented_side(p4) == CGAL::ON_ORIENTED_BOUNDARY );
 CGAL::Point_d<R> p8(3, n3, n0, n0, n3);
 CGAL::Point_d<R> p9(3, n0, n3, n0, n3);
 assert( t3.has_on_boundary( p8 ) );
 assert( t3.has_on_boundary( p9 ) );
 assert( t4.has_on_boundary( p8 ) );
 assert( t4.has_on_boundary( p9 ) );
 assert( t3.bounded_side(p8) == CGAL::ON_BOUNDARY );
 assert( t3.oriented_side(p8) == CGAL::ON_ORIENTED_BOUNDARY );
 assert( t4.bounded_side(p9) == CGAL::ON_BOUNDARY );
 assert( t4.oriented_side(p9) == CGAL::ON_ORIENTED_BOUNDARY );

 std::cout << "done" << std::endl;
 return true;
}
#endif // CGAL__TEST_CLS_TETRAHEDRON_D_C
