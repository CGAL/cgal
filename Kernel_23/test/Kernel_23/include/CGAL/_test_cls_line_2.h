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
 

#ifndef CGAL__TEST_CLS_LINE_2_H
#define CGAL__TEST_CLS_LINE_2_H

template <class R>
bool
_test_cls_line_2(const R& )
{
 std::cout << "Testing class Line_2";

 typedef typename  R::RT    RT;
 typedef typename  R::FT    gnuFT;

 RT n0 =  0;
 RT n1 =  1;
 RT n2 =  2;
 RT n3 =  3;
 RT n4 =  4;
 RT n5 =  5;
 RT n6 =  6;
 RT n8 =  8;

 CGAL::Point_2<R> p1( n2, n4, n2 );   // ( 1, 2 )
 CGAL::Point_2<R> p2( n6, n8, n2 );   // ( 3, 4 )
 CGAL::Point_2<R> p3(-n6, n6, n3 );   // (-2, 2 )
 CGAL::Point_2<R> p4( n8, n4, n2 );   // ( 4, 2 )

 typename R::Line_2 il;
 CGAL::Line_2<R>  l0(il);
 CGAL::Line_2<R>  l12( p1, p2 );
 CGAL::Line_2<R>  l21( p2, p1 );
 CGAL::Line_2<R>  l34( p3, p4 );
 CGAL::Line_2<R>  l43( p4, p3 );
 CGAL::Line_2<R>  lc ( l12 );
 l0 = l12 ;

 CGAL::Segment_2<R> s12( p1, p2);
 CGAL::Line_2<R>  ls12( s12 );
 CGAL::Ray_2<R>  r21( p2, p1 );
 CGAL::Line_2<R>  lr21( r21 );
 CGAL::Direction_2<R> d12 (p2 - p1);
 CGAL::Direction_2<R> d21 (p1 - p2);
 CGAL::Line_2<R>  ld12(p1, d12);
 CGAL::Line_2<R>  ld21(p2, d21);

 CGAL::Vector_2<R> v12 (p2 - p1);
 CGAL::Vector_2<R> v21 (p1 - p2);
 CGAL::Line_2<R>  ld12v(p1, v12);
 CGAL::Line_2<R>  ld21v(p2, v21);

 CGAL_test_assert(ld12 == ld12v);
 CGAL_test_assert(ld21 == ld21v);

 std::cout << '.';

 CGAL_test_assert( l12 == l12 );
 CGAL_test_assert( l0  == l12 );
 CGAL_test_assert( l12 == lc  );
 CGAL_test_assert( l21 == lr21 );
 CGAL_test_assert( l12 == ls12 );
 CGAL_test_assert( l12 == ld12 );
 CGAL_test_assert( l12 != l21 );
 CGAL_test_assert( l12 != ld21 );
 CGAL_test_assert( lr21 != ls12 );
 CGAL_test_assert( l34 != l43 );
 CGAL_test_assert( l12 != l43 );

 CGAL_test_assert( l34.opposite() == l43);
 CGAL_test_assert( l43.opposite() == l34);
 CGAL_test_assert( l43.opposite().opposite() == l43);
 CGAL_test_assert( ld12 == ld21.opposite() );

 CGAL::Line_2<R> labc( n2, n1, n4);
 CGAL_test_assert( labc.a() == n2 );
 CGAL_test_assert( labc.b() == n1 );
 CGAL_test_assert( labc.c() == n4 );

 CGAL_test_assert( l12.direction() == d12 );
 CGAL_test_assert( l21.direction() == d21 );
 CGAL_test_assert( ld21.direction() == d21 );
 CGAL_test_assert( ld21.direction() ==  - ld12.direction() );
 CGAL_test_assert( labc.direction() == CGAL::Direction_2<R>(labc.b(), - labc.a() ) );

 CGAL_test_assert( ld12v.to_vector().direction() == v12.direction() );
 CGAL_test_assert( ld21v.to_vector().direction() == v21.direction() );

 std::cout << '.';

 CGAL_test_assert( l43.has_on( l43.point(0) ) );
 CGAL_test_assert( lr21.has_on( lr21.point(1) ) );
 CGAL_test_assert( ld21.has_on( ld21.point(-2) ) );
 CGAL_test_assert( lr21.has_on( r21.source() ) );
 CGAL_test_assert( labc.has_on( labc.point(0) ) );

 CGAL_test_assert( l43.is_horizontal() );
 CGAL_test_assert( ! l34.is_vertical() );
 CGAL_test_assert( CGAL::Line_2<R>( n1, n0, n3 ).is_vertical() );
 CGAL_test_assert( CGAL::Line_2<R>( n0, n2, n3 ).is_horizontal() );
 CGAL_test_assert( ! lr21.is_horizontal() );

 CGAL_test_assert( ld12.y_at_x( gnuFT(3) ) == gnuFT( 4) );
 CGAL_test_assert( lr21.y_at_x( gnuFT(3) ) == gnuFT( 4) );
 CGAL_test_assert( ld12.y_at_x( gnuFT(1) ) == gnuFT( 2) );
 CGAL_test_assert( l12.y_at_x( gnuFT(5) ) == gnuFT( 6) );
 CGAL_test_assert( l34.y_at_x( gnuFT(8) ) == gnuFT( 2) );

 CGAL_test_assert( l12.x_at_y( gnuFT(0) ) == gnuFT( -1 ) );
 CGAL_test_assert( ls12.x_at_y( gnuFT(4) ) == gnuFT( 3 ) );
 CGAL_test_assert( l21.x_at_y( gnuFT(6) ) == gnuFT( 5 ) );
 CGAL_test_assert( ld21.x_at_y( gnuFT(2) ) == gnuFT( 1 ) );

 CGAL::Direction_2<R> up( n0, n1 );
 CGAL::Aff_transformation_2<R> rot90(CGAL::ROTATION, up, n1, RT(100) );
 CGAL::Line_2<R> l12perp1( l12.perpendicular( p1 ) );
 CGAL::Line_2<R> l21perp1( l21.perpendicular( p1 ) );
 CGAL::Line_2<R> labcperp( labc.perpendicular( labc.point(1) ) );
 CGAL_test_assert( l12perp1.opposite() == l21perp1 );
 CGAL_test_assert( labcperp.direction() == CGAL::Direction_2<R>( labc.a(), labc.b()) );
 CGAL_test_assert( l12perp1.has_on( p1 ) );
 CGAL_test_assert( l21perp1.has_on( p1 ) );
 CGAL::Line_2<R> l12perp4( l12.perpendicular( p4 ) );
 CGAL_test_assert( l12perp4.has_on( p4 ) );
 CGAL_test_assert( l12.direction().transform( rot90 ) == l12perp4.direction() );

 CGAL_test_assert( CGAL::Line_2<R>( n0, n0, n6 ).is_degenerate() );
 CGAL_test_assert( CGAL::Line_2<R>( p1, p1 ).is_degenerate() );
 CGAL_test_assert( ! CGAL::Line_2<R>( p1, p3 ).is_degenerate() );
 CGAL_test_assert( ! l34.is_degenerate() );

 std::cout << '.';

 CGAL::Point_2<R> p5( n5, n6 );
 CGAL_test_assert( l12.oriented_side(p3) == CGAL::ON_POSITIVE_SIDE );
 CGAL_test_assert( l12.oriented_side(p4) == CGAL::ON_NEGATIVE_SIDE );
 CGAL_test_assert( l12.oriented_side(p2) == CGAL::ON_ORIENTED_BOUNDARY );
 CGAL_test_assert( l12.oriented_side(p5) == CGAL::ON_ORIENTED_BOUNDARY );
 CGAL_test_assert( l21.oriented_side(p3) == CGAL::ON_NEGATIVE_SIDE );
 CGAL_test_assert( l21.oriented_side(p5) == CGAL::ON_ORIENTED_BOUNDARY );

 CGAL_test_assert( l21.has_on_negative_side( p3 ) );
 CGAL_test_assert( l12.has_on_positive_side( p3 ) );
 CGAL_test_assert( l34.has_on_positive_side( p2 ) );
 CGAL_test_assert( l43.has_on( CGAL::Point_2<R>( n8, n2 )) );
 CGAL_test_assert( l43.has_on_boundary( CGAL::Point_2<R>( n8, n2 )) );
 CGAL_test_assert( lr21.has_on( CGAL::Point_2<R>( -n1, n0 )) );

 std::cout << '.';

 CGAL_test_assert( l21.has_on( l21.projection( p3 )) );
 CGAL_test_assert( l21.has_on( l21.projection( p4 )) );
 CGAL_test_assert( l21.has_on( l21.projection( p5 )) );
 CGAL_test_assert( l34.has_on( l34.projection( p3 )) );
 CGAL_test_assert( l34.has_on( l34.projection( p4 )) );
 CGAL_test_assert( l34.has_on( l34.projection( p5 )) );


 std::cout << "done" << std::endl;
 return true;

}
#endif // CGAL__TEST_CLS_LINE_2_H
