// revision      : 2.2.2
// revision_date : 28 Sep 1999 
// author(s)     : Stefan Schirra <Stefan.Schirra@mpi-sb.mpg.de>

#ifndef CGAL__TEST_CLS_PLANE_D_C
#define CGAL__TEST_CLS_PLANE_D_C
#include <CGAL/_test_cls_plane_d.h>

template <class gnuR>
bool
_test_cls_plane_d(const gnuR& )
{
 std::cout << "Testing class Plane_d" ;

 typedef typename  gnuR::RT    RT;
 typedef typename  gnuR::FT    FT;

 typename gnuR::Plane_d  ip;
 CGAL::Plane_d<gnuR> pl0(ip);

 RT  x1 = -1;
 RT  x2 =  4;

 CGAL::Point_d<gnuR> p1(3, x1,x1,x1),
                 p2(3,x2,x2,x2),
                 px(3,RT(1),RT(0),RT(0)),
                 py(3,RT(0),RT(1),RT(0)),
                 pz(3,RT(0),RT(0),RT(1));

 CGAL::Point_d<gnuR> p3(p1);
 CGAL::Direction_d<gnuR> d1(3, RT(5), RT(-5), RT(10) );
 CGAL::Vector_d<gnuR>    v1 = d1.to_vector();

 CGAL::Plane_d<gnuR> pld( p1, d1);
 CGAL::Plane_d<gnuR> plv( p1, v1);

 CGAL::Plane_d<gnuR> xy_pl(3,CGAL::ORIGIN,px,py);
 CGAL::Plane_d<gnuR> plc(xy_pl);
 CGAL::Plane_d<gnuR> pl1(3,p1,pz,p2);
 CGAL::Plane_d<gnuR> pla;
 CGAL::Plane_d<gnuR> xy_pl_eq(3, RT(0), RT(0),
                              RT(1), RT(0) );
 CGAL::Plane_d<gnuR> neg_xy_pl_eq(3, RT(0), RT(0),
                                  RT(-1), RT(0) );

 std::cout << '.';

 pla = xy_pl;
 assert( pla == pla );
 assert( plc ==  xy_pl );
 assert( pla == xy_pl );
 assert( plc == pla );
 assert( xy_pl == xy_pl_eq );
 assert( xy_pl != neg_xy_pl_eq);
 assert( pld == plv );

 CGAL::Point_d<gnuR> origo(3, CGAL::ORIGIN );
 CGAL::Vector_d<gnuR> v_zdir(3, RT(0), RT(0), RT(12) );
 CGAL::Plane_d<gnuR> xy_pl_v(origo, v_zdir);

 assert( xy_pl == xy_pl_v );

 assert( xy_pl_eq.a() == RT(0) );
 assert( xy_pl_eq.b() == RT(0) );
 assert( xy_pl_eq.c() == RT(1) );
 assert( xy_pl_eq.d() == RT(0) );
 assert( CGAL::Plane_d<gnuR>(3, pl1.a(), pl1.b(), pl1.c(), pl1.d() ) == pl1 );

 std::cout <<'.';

 // assert( xy_pl_eq.perpendicular_line(px) == xy_pl.perpendicular_line(px) );
 // assert( xy_pl.perpendicular_line(px) == neg_xy_pl_eq.perpendicular_line(px).opposite() );
 // assert( neg_xy_pl_eq.opposite() == xy_pl );
 // assert( plr != plr.opposite() );

 CGAL::Point_d<gnuR> pp0(3, RT(4), RT(6), RT( 0), RT(2) );
 CGAL::Point_d<gnuR> pp1(3, RT(4), RT(6), RT(-8), RT(2) );
 assert( xy_pl_eq.projection( pp0 ) == pp0 );
 assert( xy_pl_eq.has_on( xy_pl_eq.projection( pp0 ) )  );
 assert( xy_pl_eq.has_on( xy_pl_eq.projection( pp1 ) )  );
 assert( pl1.has_on( pl1.projection( pp0 ) )  );
 assert( pl1.has_on( pl1.projection( pp1 ) )  );

 std::cout <<'.';

 assert( plc.has_on_boundary(plc.point()) );
 assert( plc.orthogonal_direction() == pla.orthogonal_direction() );
 assert( plc.perpendicular_line( plc.point() )
  == CGAL::Line_d<gnuR>( plc.point(), plc.orthogonal_direction()) );
 assert( CGAL::Line_d<gnuR>( pl1.point(), pl1.point()+pl1.orthogonal_vector() )
  == CGAL::Line_d<gnuR>( pl1.point(), pl1.orthogonal_direction()) );
 CGAL::Point_d<gnuR>  gnup(3,RT(345),RT(23),RT(0));
 assert( xy_pl.has_on_boundary( gnup ) );

 std::cout <<'.';

 CGAL::Vector_d<gnuR> nov = pl1.orthogonal_vector();
 CGAL::Vector_d<gnuR> vb1 = pl1.base(0);
 CGAL::Vector_d<gnuR> vb2 = pl1.base(1);
 assert( (nov*pl1.base(0)) == FT(0) && (nov*pl1.base(1)) == FT(0) );
 assert( (pl1.base(1)*pl1.base(0)) == FT(0) );
 assert( pl1.has_on_boundary(pl1.point() + pl1.base(0)) );
 assert( pl1.has_on_boundary(pl1.point() + pl1.base(1)) );

 std::cout << '.';

 assert( neg_xy_pl_eq.oriented_side( p1 ) == CGAL::ON_POSITIVE_SIDE );
 assert( xy_pl.oriented_side( p1 ) == CGAL::ON_NEGATIVE_SIDE );
 assert( pl1.oriented_side( p1 ) == CGAL::ON_ORIENTED_BOUNDARY );
 assert( pl1.oriented_side( pz ) == CGAL::ON_ORIENTED_BOUNDARY );
 assert( pl1.oriented_side( p2 ) == CGAL::ON_ORIENTED_BOUNDARY );
 assert( xy_pl.has_on_positive_side( p2 ) );
 assert( xy_pl.has_on_negative_side( p1 ) );

 CGAL::Plane_d<gnuR> pldeg(3, p1, p1, p2);
 assert( pldeg.is_degenerate() );

 std::cout << "done" << std::endl;
 return true;
}

#endif // CGAL__TEST_CLS_PLANE_D_C
