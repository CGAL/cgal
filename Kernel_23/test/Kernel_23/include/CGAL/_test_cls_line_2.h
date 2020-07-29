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


#ifndef CGAL__TEST_CLS_LINE_2_H
#define CGAL__TEST_CLS_LINE_2_H

#include <CGAL/use.h>

template <class R>
bool
_test_cls_line_2(const R& )
{
 std::cout << "Testing class Line_2";

 typedef typename  R::RT    RT;
 typedef typename  R::FT    FT;

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

 typename R::Line_2 il0; CGAL_USE(il0); // test default-construction
 typename R::Line_2 il(p1, p2);
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

 assert(ld12 == ld12v);
 assert(ld21 == ld21v);

 std::cout << '.';

 assert( l12 == l12 );
 assert( l0  == l12 );
 assert( l12 == lc  );
 assert( l21 == lr21 );
 assert( l12 == ls12 );
 assert( l12 == ld12 );
 assert( l12 != l21 );
 assert( l12 != ld21 );
 assert( lr21 != ls12 );
 assert( l34 != l43 );
 assert( l12 != l43 );

 assert( l34.opposite() == l43);
 assert( l43.opposite() == l34);
 assert( l43.opposite().opposite() == l43);
 assert( ld12 == ld21.opposite() );

 CGAL::Line_2<R> labc( n2, n1, n4);
 assert( labc.a() == n2 );
 assert( labc.b() == n1 );
 assert( labc.c() == n4 );

 assert( l12.direction() == d12 );
 assert( l21.direction() == d21 );
 assert( ld21.direction() == d21 );
 assert( ld21.direction() ==  - ld12.direction() );
 assert( labc.direction() == CGAL::Direction_2<R>(labc.b(), - labc.a() ) );

 assert( ld12v.to_vector().direction() == v12.direction() );
 assert( ld21v.to_vector().direction() == v21.direction() );

 std::cout << '.';

 assert( l43.has_on( l43.point(0) ) );
 assert( lr21.has_on( lr21.point(1) ) );
 assert( ld21.has_on( ld21.point(-2) ) );
 assert( ld21.has_on( ld21.point(FT(1)/FT(2)) ) );
 assert( lr21.has_on( r21.source() ) );
 assert( labc.has_on( labc.point(0) ) );

 assert( l43.is_horizontal() );
 assert( ! l34.is_vertical() );
 assert( CGAL::Line_2<R>( n1, n0, n3 ).is_vertical() );
 assert( CGAL::Line_2<R>( n0, n2, n3 ).is_horizontal() );
 assert( ! lr21.is_horizontal() );

 assert( ld12.y_at_x( FT(3) ) == FT( 4) );
 assert( lr21.y_at_x( FT(3) ) == FT( 4) );
 assert( ld12.y_at_x( FT(1) ) == FT( 2) );
 assert( l12.y_at_x( FT(5) ) == FT( 6) );
 assert( l34.y_at_x( FT(8) ) == FT( 2) );

 assert( l12.x_at_y( FT(0) ) == FT( -1 ) );
 assert( ls12.x_at_y( FT(4) ) == FT( 3 ) );
 assert( l21.x_at_y( FT(6) ) == FT( 5 ) );
 assert( ld21.x_at_y( FT(2) ) == FT( 1 ) );

 CGAL::Direction_2<R> up( n0, n1 );
 CGAL::Aff_transformation_2<R> rot90(CGAL::ROTATION, up, n1, RT(100) );
 CGAL::Line_2<R> l12perp1( l12.perpendicular( p1 ) );
 CGAL::Line_2<R> l21perp1( l21.perpendicular( p1 ) );
 CGAL::Line_2<R> labcperp( labc.perpendicular( labc.point(1) ) );
 assert( l12perp1.opposite() == l21perp1 );
 assert( labcperp.direction() == CGAL::Direction_2<R>( labc.a(), labc.b()) );
 assert( l12perp1.has_on( p1 ) );
 assert( l21perp1.has_on( p1 ) );
 CGAL::Line_2<R> l12perp4( l12.perpendicular( p4 ) );
 assert( l12perp4.has_on( p4 ) );
 assert( l12.direction().transform( rot90 ) == l12perp4.direction() );

 assert( CGAL::Line_2<R>( n0, n0, n6 ).is_degenerate() );
 assert( CGAL::Line_2<R>( p1, p1 ).is_degenerate() );
 assert( ! CGAL::Line_2<R>( p1, p3 ).is_degenerate() );
 assert( ! l34.is_degenerate() );

 std::cout << '.';

 CGAL::Point_2<R> p5( n5, n6 );
 assert( l12.oriented_side(p3) == CGAL::ON_POSITIVE_SIDE );
 assert( l12.oriented_side(p4) == CGAL::ON_NEGATIVE_SIDE );
 assert( l12.oriented_side(p2) == CGAL::ON_ORIENTED_BOUNDARY );
 assert( l12.oriented_side(p5) == CGAL::ON_ORIENTED_BOUNDARY );
 assert( l21.oriented_side(p3) == CGAL::ON_NEGATIVE_SIDE );
 assert( l21.oriented_side(p5) == CGAL::ON_ORIENTED_BOUNDARY );

 assert( l21.has_on_negative_side( p3 ) );
 assert( l12.has_on_positive_side( p3 ) );
 assert( l34.has_on_positive_side( p2 ) );
 assert( l43.has_on( CGAL::Point_2<R>( n8, n2 )) );
 assert( l43.has_on_boundary( CGAL::Point_2<R>( n8, n2 )) );
 assert( lr21.has_on( CGAL::Point_2<R>( -n1, n0 )) );

 std::cout << '.';

 assert( l21.has_on( l21.projection( p3 )) );
 assert( l21.has_on( l21.projection( p4 )) );
 assert( l21.has_on( l21.projection( p5 )) );
 assert( l34.has_on( l34.projection( p3 )) );
 assert( l34.has_on( l34.projection( p4 )) );
 assert( l34.has_on( l34.projection( p5 )) );


 std::cout << "done" << std::endl;
 return true;

}
#endif // CGAL__TEST_CLS_LINE_2_H
