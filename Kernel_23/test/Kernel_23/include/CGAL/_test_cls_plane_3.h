// Copyright (c) 1999  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Stefan Schirra
 

#ifndef CGAL__TEST_CLS_PLANE_3_H
#define CGAL__TEST_CLS_PLANE_3_H

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
 CGAL::Vector_3<R>    v1 = d1.vector();

 CGAL::Plane_3<R> pld( p1, d1);
 CGAL::Plane_3<R> plv( p1, v1);

 CGAL::Line_3<R>  l,
                 lp(p1,p2),
                 lxy(px,py);

 CGAL::Segment_3<R>  s12(p1,p2);
 CGAL::Ray_3<R>      r12(p1,p2);

 CGAL::Point_3<R> origo( CGAL::ORIGIN );

 CGAL::Plane_3<R> xy_pl(origo,px,py);
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
 CGAL_test_assert( pla == pla );
 CGAL_test_assert( plc ==  xy_pl );
 CGAL_test_assert( pla == xy_pl );
 CGAL_test_assert( plc == pla );
 CGAL_test_assert( xy_pl == xy_pl_eq );
 CGAL_test_assert( xy_pl != neg_xy_pl_eq);
 CGAL_test_assert( pld == plv );
 CGAL_test_assert( pls == plr );

 CGAL::Vector_3<R> v_zdir( RT(0), RT(0), RT(12) );
 CGAL::Plane_3<R> xy_pl_v( origo, v_zdir);

 CGAL_test_assert( xy_pl == xy_pl_v );

 CGAL_test_assert( xy_pl_eq.a() == RT(0) );
 CGAL_test_assert( xy_pl_eq.b() == RT(0) );
 CGAL_test_assert( xy_pl_eq.c() == RT(1) );
 CGAL_test_assert( xy_pl_eq.d() == RT(0) );
 CGAL_test_assert( CGAL::Plane_3<R>(pl1.a(), pl1.b(), pl1.c(), pl1.d() ) == pl1 );

 std::cout <<'.';

 CGAL_test_assert( xy_pl_eq.perpendicular_line(px) == xy_pl.perpendicular_line(px) );
 CGAL_test_assert( xy_pl.perpendicular_line(px) == \
         neg_xy_pl_eq.perpendicular_line(px).opposite() );
 CGAL_test_assert( neg_xy_pl_eq.opposite() == xy_pl );
 CGAL_test_assert( plr != plr.opposite() );

 CGAL::Point_3<R> pp0( RT(4), RT(6), RT( 0), RT(2) );
 CGAL::Point_3<R> pp1( RT(4), RT(6), RT(-8), RT(2) );
 CGAL_test_assert( xy_pl_eq.projection( pp0 ) == pp0 );
 CGAL_test_assert( xy_pl_eq.has_on( xy_pl_eq.projection( pp0 ) )  );
 CGAL_test_assert( xy_pl_eq.has_on( xy_pl_eq.projection( pp1 ) )  );
 CGAL_test_assert( pl1.has_on( pl1.projection( pp0 ) )  );
 CGAL_test_assert( pl1.has_on( pl1.projection( pp1 ) )  );

 CGAL_test_assert( plc.has_on(plc.point()) );
 CGAL_test_assert( plc.orthogonal_direction() == pla.orthogonal_direction() );
 CGAL_test_assert( plc.perpendicular_line( plc.point() ) == \
         CGAL::Line_3<R>( plc.point(), plc.orthogonal_direction()) );
 CGAL_test_assert( CGAL::Line_3<R>( pl1.point(), pl1.point()+pl1.orthogonal_vector() )\
         == CGAL::Line_3<R>( pl1.point(), pl1.orthogonal_direction()) );
 CGAL::Point_3<R>  gnup(RT(345),RT(23),RT(0));
 CGAL_test_assert( xy_pl.has_on( gnup ) );

 CGAL::Vector_3<R> nov = pl1.orthogonal_vector();
 CGAL::Vector_3<R> vb1 = pl1.base1();
 CGAL::Vector_3<R> vb2 = pl1.base2();
 CGAL_test_assert( (nov*pl1.base1()) == FT(0)&&(nov*pl1.base2()) == FT(0) );
 CGAL_test_assert( (pl1.base2()*pl1.base1()) == FT(0) );
 CGAL_test_assert( pl1.has_on(pl1.point() + pl1.base1()) );
 CGAL_test_assert( pl1.has_on(pl1.point() + pl1.base2()) );

 std::cout << '.';

 CGAL_test_assert( pl1.has_on( pl1.to_3d( pl1.to_2d( pl1.point() ))) );
 CGAL_test_assert( pl1.has_on( pl1.to_3d( pl1.to_2d( pz ))) );

 CGAL_test_assert( neg_xy_pl_eq.oriented_side( p1 ) == CGAL::ON_POSITIVE_SIDE );
 CGAL_test_assert( xy_pl.oriented_side( p1 ) == CGAL::ON_NEGATIVE_SIDE );
 CGAL_test_assert( pl1.oriented_side( p1 ) == CGAL::ON_ORIENTED_BOUNDARY );
 CGAL_test_assert( xy_pl.has_on_positive_side( p2 ) );
 CGAL_test_assert( xy_pl.has_on_negative_side( p1 ) );

 CGAL::Plane_3<R> pldeg(p1, p1, p2);
 CGAL_test_assert( pldeg.is_degenerate() );

 std::cout << "done" << std::endl;
 return true;
}

#endif // CGAL__TEST_CLS_PLANE_3_H
