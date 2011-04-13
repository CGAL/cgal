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
// file          : _test_cls_plane_3.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#ifndef CGAL__TEST_CLS_PLANE_3_C
#define CGAL__TEST_CLS_PLANE_3_C

#include <CGAL/_test_cls_plane_3.h>

template <class R>
bool
_test_cls_plane_3(const R& )
{
 std::cout << "Testing class Plane_3" ;

 typedef typename  R::RT    RT;
 typedef typename  R::FT    FT;

 typename R::Plane_3  ip;
 CGAL::Plane_3<R> pl0(ip);

 RT  x1 = -1;
 RT  x2 =  4;

 CGAL::Point_3<R> p1(x1,x1,x1),
                 p2(x2,x2,x2),
                 px(RT(1),RT(0),RT(0)),
                 py(RT(0),RT(1),RT(0)),
                 pz(RT(0),RT(0),RT(1));

 CGAL::Point_3<R> p3(p1);
 CGAL::Direction_3<R> d1( RT(5), RT(-5), RT(10) );
#if (__GNUG__ == 2) && (__GNUC_MINOR__==91)
 CGAL::Vector_3<R>    v1 = d1.to_vector();
#else
 CGAL::Vector_3<R>    v1 = d1.vector();
#endif // egcs 2.91.66

 CGAL::Plane_3<R> pld( p1, d1);
 CGAL::Plane_3<R> plv( p1, v1);

 CGAL::Line_3<R>  l,
                 lp(p1,p2),
                 lxy(px,py);

 CGAL::Segment_3<R>  s12(p1,p2);
 CGAL::Ray_3<R>      r12(p1,p2);

 CGAL::Plane_3<R> xy_pl(CGAL::ORIGIN,px,py);
 CGAL::Plane_3<R> plc(xy_pl);
 CGAL::Plane_3<R> pl1(p1,pz,p2);
 CGAL::Plane_3<R> pls(s12,pz);
 CGAL::Plane_3<R> plr(r12,pz);
 CGAL::Plane_3<R> pla;
 CGAL::Plane_3<R> xy_pl_eq( RT(0), RT(0),
                              RT(1), RT(0) );
 CGAL::Plane_3<R> neg_xy_pl_eq( RT(0), RT(0),
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
 assert( pls == plr );

 CGAL::Point_3<R> origo( CGAL::ORIGIN );
 CGAL::Vector_3<R> v_zdir( RT(0), RT(0), RT(12) );
 CGAL::Plane_3<R> xy_pl_v( origo, v_zdir);

 assert( xy_pl == xy_pl_v );

 assert( xy_pl_eq.a() == RT(0) );
 assert( xy_pl_eq.b() == RT(0) );
 assert( xy_pl_eq.c() == RT(1) );
 assert( xy_pl_eq.d() == RT(0) );
 assert( CGAL::Plane_3<R>(pl1.a(), pl1.b(), pl1.c(), pl1.d() ) == pl1 );

 std::cout <<'.';

 assert( xy_pl_eq.perpendicular_line(px) == xy_pl.perpendicular_line(px) );
 assert( xy_pl.perpendicular_line(px) == \
         neg_xy_pl_eq.perpendicular_line(px).opposite() );
 assert( neg_xy_pl_eq.opposite() == xy_pl );
 assert( plr != plr.opposite() );

 CGAL::Point_3<R> pp0( RT(4), RT(6), RT( 0), RT(2) );
 CGAL::Point_3<R> pp1( RT(4), RT(6), RT(-8), RT(2) );
 assert( xy_pl_eq.projection( pp0 ) == pp0 );
 assert( xy_pl_eq.has_on( xy_pl_eq.projection( pp0 ) )  );
 assert( xy_pl_eq.has_on( xy_pl_eq.projection( pp1 ) )  );
 assert( pl1.has_on( pl1.projection( pp0 ) )  );
 assert( pl1.has_on( pl1.projection( pp1 ) )  );

 assert( plc.has_on(plc.point()) );
 assert( plc.orthogonal_direction() == pla.orthogonal_direction() );
 assert( plc.perpendicular_line( plc.point() ) == \
         CGAL::Line_3<R>( plc.point(), plc.orthogonal_direction()) );
 assert( CGAL::Line_3<R>( pl1.point(), pl1.point()+pl1.orthogonal_vector() )\
         == CGAL::Line_3<R>( pl1.point(), pl1.orthogonal_direction()) );
 CGAL::Point_3<R>  gnup(RT(345),RT(23),RT(0));
 assert( xy_pl.has_on( gnup ) );

 CGAL::Vector_3<R> nov = pl1.orthogonal_vector();
 CGAL::Vector_3<R> vb1 = pl1.base1();
 CGAL::Vector_3<R> vb2 = pl1.base2();
 assert( (nov*pl1.base1()) == FT(0)&&(nov*pl1.base2()) == FT(0) );
 assert( (pl1.base2()*pl1.base1()) == FT(0) );
 assert( pl1.has_on(pl1.point() + pl1.base1()) );
 assert( pl1.has_on(pl1.point() + pl1.base2()) );

 std::cout << '.';

 assert( pl1.has_on( pl1.to_3d( pl1.to_2d( pl1.point() ))) );
 assert( pl1.has_on( pl1.to_3d( pl1.to_2d( pz ))) );

 assert( neg_xy_pl_eq.oriented_side( p1 ) == CGAL::ON_POSITIVE_SIDE );
 assert( xy_pl.oriented_side( p1 ) == CGAL::ON_NEGATIVE_SIDE );
 assert( pl1.oriented_side( p1 ) == CGAL::ON_ORIENTED_BOUNDARY );
 assert( xy_pl.has_on_positive_side( p2 ) );
 assert( xy_pl.has_on_negative_side( p1 ) );

 CGAL::Plane_3<R> pldeg(p1, p1, p2);
 assert( pldeg.is_degenerate() );

 std::cout << "done" << std::endl;
 return true;
}

#endif // CGAL__TEST_CLS_PLANE_3_C
