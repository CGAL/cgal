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
// source        : test_kernel_3.fw
// file          : _test_cls_plane_3.C
// revision      : 2.0.5
// revision_date : 24 Mar 1999 
// author(s)     : Stefan Schirra <Stefan.Schirra@mpi-sb.mpg.de>
//
// coordinator   : MPI, Saarbruecken
// ============================================================================


#ifndef CGAL__TEST_CLS_PLANE_3_C
#define CGAL__TEST_CLS_PLANE_3_C
#ifndef CGAL__TEST_CLS_PLANE_3_H
#include <CGAL/_test_cls_plane_3.h>
#endif // CGAL__TEST_CLS_PLANE_3_H

template <class gnuR>
bool
_test_cls_plane_3(const gnuR& )
{
 cout << "Testing class Plane_3" ;

 typedef typename  gnuR::RT    RT;
 typedef typename  gnuR::FT    FT;

 typename gnuR::Plane_3  ip;
 CGAL::Plane_3<gnuR> pl0(ip);

 RT  x1 = -1;
 RT  x2 =  4;

 CGAL::Point_3<gnuR> p1(x1,x1,x1),
                 p2(x2,x2,x2),
                 px(RT(1),RT(0),RT(0)),
                 py(RT(0),RT(1),RT(0)),
                 pz(RT(0),RT(0),RT(1));

 CGAL::Point_3<gnuR> p3(p1);
 CGAL::Direction_3<gnuR> d1( RT(5), RT(-5), RT(10) );
 CGAL::Vector_3<gnuR>    v1 = d1.to_vector();

 CGAL::Plane_3<gnuR> pld( p1, d1);
 CGAL::Plane_3<gnuR> plv( p1, v1);

 CGAL::Line_3<gnuR>  l,
                 lp(p1,p2),
                 lxy(px,py);

 CGAL::Segment_3<gnuR>  s12(p1,p2);
 CGAL::Ray_3<gnuR>      r12(p1,p2);

 CGAL::Plane_3<gnuR> xy_pl(CGAL::ORIGIN,px,py);
 CGAL::Plane_3<gnuR> plc(xy_pl);
 CGAL::Plane_3<gnuR> pl1(p1,pz,p2);
 CGAL::Plane_3<gnuR> pls(s12,pz);
 CGAL::Plane_3<gnuR> plr(r12,pz);
 CGAL::Plane_3<gnuR> pla;
 CGAL::Plane_3<gnuR> xy_pl_eq( RT(0), RT(0),
                              RT(1), RT(0) );
 CGAL::Plane_3<gnuR> neg_xy_pl_eq( RT(0), RT(0),
                                  RT(-1), RT(0) );

 cout << '.';

 pla = xy_pl;
 assert( pla == pla );
 assert( plc ==  xy_pl );
 assert( pla == xy_pl );
 assert( plc == pla );
 assert( xy_pl == xy_pl_eq );
 assert( xy_pl != neg_xy_pl_eq);
 assert( pld == plv );
 assert( pls == plr );

 CGAL::Point_3<gnuR> origo( CGAL::ORIGIN );
 CGAL::Vector_3<gnuR> v_zdir( RT(0), RT(0), RT(12) );
 CGAL::Plane_3<gnuR> xy_pl_v( origo, v_zdir);

 assert( xy_pl == xy_pl_v );

 assert( xy_pl_eq.a() == RT(0) );
 assert( xy_pl_eq.b() == RT(0) );
 assert( xy_pl_eq.c() == RT(1) );
 assert( xy_pl_eq.d() == RT(0) );
 assert( CGAL::Plane_3<gnuR>(pl1.a(), pl1.b(), pl1.c(), pl1.d() ) == pl1 );

 cout <<'.';

 assert( xy_pl_eq.perpendicular_line(px) == xy_pl.perpendicular_line(px) );
 assert( xy_pl.perpendicular_line(px) == \
         neg_xy_pl_eq.perpendicular_line(px).opposite() );
 assert( neg_xy_pl_eq.opposite() == xy_pl );
 assert( plr != plr.opposite() );

#ifndef CGAL_STRICT_09
 CGAL::Point_3<gnuR> pp0( RT(4), RT(6), RT( 0), RT(2) );
 CGAL::Point_3<gnuR> pp1( RT(4), RT(6), RT(-8), RT(2) );
 assert( xy_pl_eq.projection( pp0 ) == pp0 );
 assert( xy_pl_eq.has_on( xy_pl_eq.projection( pp0 ) )  );
 assert( xy_pl_eq.has_on( xy_pl_eq.projection( pp1 ) )  );
 assert( pl1.has_on( pl1.projection( pp0 ) )  );
 assert( pl1.has_on( pl1.projection( pp1 ) )  );
#endif // CGAL_STRICT_09

 assert( plc.has_on_boundary(plc.point()) );
 assert( plc.orthogonal_direction() == pla.orthogonal_direction() );
 assert( plc.perpendicular_line( plc.point() ) == \
         CGAL::Line_3<gnuR>( plc.point(), plc.orthogonal_direction()) );
 assert( CGAL::Line_3<gnuR>( pl1.point(), pl1.point()+pl1.orthogonal_vector() )\
         == CGAL::Line_3<gnuR>( pl1.point(), pl1.orthogonal_direction()) );
 CGAL::Point_3<gnuR>  gnup(RT(345),RT(23),RT(0));
 assert( xy_pl.has_on_boundary( gnup ) );

 CGAL::Vector_3<gnuR> nov = pl1.orthogonal_vector();
 CGAL::Vector_3<gnuR> vb1 = pl1.base1();
 CGAL::Vector_3<gnuR> vb2 = pl1.base2();
 assert( (nov*pl1.base1()) == FT(0)&&(nov*pl1.base2()) == FT(0) );
 assert( (pl1.base2()*pl1.base1()) == FT(0) );
 assert( pl1.has_on_boundary(pl1.point() + pl1.base1()) );
 assert( pl1.has_on_boundary(pl1.point() + pl1.base2()) );

 cout << '.';

#ifdef CGAL_CARTESIAN_BOTH_23_KERNEL
 assert( pl1.has_on_boundary( pl1.to_3d( pl1.to_2d( pl1.point() ))) );
 assert( pl1.has_on_boundary( pl1.to_3d( pl1.to_2d( pz ))) );
#endif

 assert( neg_xy_pl_eq.oriented_side( p1 ) == CGAL::ON_POSITIVE_SIDE );
 assert( xy_pl.oriented_side( p1 ) == CGAL::ON_NEGATIVE_SIDE );
 assert( pl1.oriented_side( p1 ) == CGAL::ON_ORIENTED_BOUNDARY );
 assert( xy_pl.has_on_positive_side( p2 ) );
 assert( xy_pl.has_on_negative_side( p1 ) );

 CGAL::Plane_3<gnuR> pldeg(p1, p1, p2);
 assert( pldeg.is_degenerate() );

 cout << "done" << endl;
 return true;
}

#endif // CGAL__TEST_CLS_PLANE_3_C
