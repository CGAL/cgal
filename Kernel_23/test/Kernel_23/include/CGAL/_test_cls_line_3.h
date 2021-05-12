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


#ifndef CGAL__TEST_CLS_LINE_3_H
#define CGAL__TEST_CLS_LINE_3_H

#include <CGAL/use.h>
#include <boost/type_traits/is_same.hpp>

template <class R>
bool
_test_cls_line_3(const R& )
{
 std::cout << "Testing class Line_3" ;

 typedef typename  R::RT    RT;
 typedef typename  R::FT    FT;

 const bool nonexact = boost::is_same<RT, double>::value;

 typename R::Line_3 il;
 CGAL::Line_3<R> l0( il ); CGAL_USE(l0);
 CGAL::Line_3<R> l1; CGAL_USE(l1);

 RT n1 =  3;
 RT n2 = 53;
 RT n3 =-14;
 RT n4 = 28;
 RT n5 = 16;
 RT n6 = -6;
 RT n7 = 11;
 RT n8 =-17;
 RT n9 = 30;

 CGAL::Point_3<R> p1( n1, n2, n3);
 CGAL::Point_3<R> p2( n4, n5, n6);
 CGAL::Point_3<R> p3( n7, n8, n9);

 CGAL::Line_3<R> l2(p1,p2);
 CGAL::Line_3<R> l3(l2);

 assert( l2 == l2);
 assert( l2 == l3);
 assert( CGAL::parallel(l2, l3) );

 CGAL::Direction_3<R> dir( n9, n3, n1);
 CGAL::Line_3<R> l4(p3, dir);
 assert( l2 != l4);
 assert( ! CGAL::parallel(l2, l4) );

 CGAL::Vector_3<R> vec( n9, n3, n1);
 CGAL::Line_3<R> l4v(p3, vec);
 assert( l4 == l4v);
 assert( l4.to_vector() == vec);

 CGAL::Segment_3<R> seg(p1,p2);
 CGAL::Ray_3<R>     ray(p2,p1);
 CGAL::Line_3<R>    l5(seg);
 CGAL::Line_3<R>    l6(ray);
 assert( l2 == l5);

 std::cout <<'.';

 assert( l2 == l5 );
 assert( l2.direction() == l5.direction() );
 assert( l5.direction() ==  - l6.direction() );
 assert( l5.has_on( p1 ) );
 assert( l5.has_on( p2 ) );
 assert( l5.has_on( l5.point() ));
 assert( l6.has_on( p1 ) );
 assert( l6.has_on( p2 ) );
 assert( l6.has_on( l5.point() ));
 assert( l5.opposite() == l6 );
 assert( l2.opposite() == l6 );
 assert( l5 != l6 );

 CGAL::Plane_3<R> pl = l6.perpendicular_plane( l6.point() );
 CGAL::Plane_3<R> plstrich(l6.point(), l6.direction() );
 assert( pl == plstrich );
 assert( pl.orthogonal_direction() == l6.direction() );
 CGAL::Plane_3<R> plzweistrich(l6.point(), l5.direction() );
 assert( plzweistrich.opposite() == pl );

 std::cout << '.';

 assert( l4.point(2) - l4.point(1) == l4.point(1) - l4.point(0) );
 assert( (l4.point(FT(1)/FT(3)) - l4.point(0) == l4.point(1+FT(1)/FT(3)) - l4.point(1)) || nonexact );

 CGAL::Point_3<R> p1l4proj = l4.projection(p1);
 assert( l4.has_on( p1l4proj ) || nonexact );
 assert( l4.perpendicular_plane( p1l4proj ).has_on( p1l4proj ) || nonexact );
 assert( l4.perpendicular_plane( p1l4proj ).has_on( p1 ) || nonexact );
 CGAL::Point_3<R> p4 = l4.projection(p2);
 CGAL::Point_3<R> p5 = l4.projection(p3);
 assert(  ( l4.direction() == ( p5 - p4 ).direction() )\
        ||( l4.direction() == ( p4 - p5 ).direction() )  || nonexact );
 assert( l5.direction() == - l6.direction() );

 std::cout <<'.';

 assert( l2.has_on(p1) );
 assert( l2.has_on(p2) );
 assert( l4.has_on(p4) || nonexact );
 assert( l4.has_on(p5) );
 assert( CGAL::Line_3<R>(p1,p1).is_degenerate() );

 {
  CGAL::Point_3<R> p(0, 0, 0);
  CGAL::Point_3<R> q(2, 0, 0);
  CGAL::Point_3<R> r(0, 2, 0);
  CGAL::Line_3<R> l(CGAL::Point_3<R>(1, 1, 0), CGAL::Vector_3<R>(0, 0, 1));
  assert( l == CGAL::equidistant_line(p, q, r) );
 }

 std::cout << "done" << std::endl;
 return true;
}

#endif // CGAL__TEST_CLS_LINE_3_H
