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


#ifndef CGAL__TEST_MF_PLANE_3_TO_2D_H
#define CGAL__TEST_MF_PLANE_3_TO_2D_H

#include <boost/type_traits/is_same.hpp>

template <class R>
bool
_test_mf_plane_3_to_2d(const R& )
{
 std::cout << "Testing member function Plane_3::to_2d" ;

 typedef typename  R::RT    RT;
 typedef CGAL::Plane_3< R>   Plane_3;
 typedef CGAL::Point_3< R>   Point_3;
 typedef CGAL::Point_2< R>   Point_2;

 const bool nonexact = boost::is_same<RT, double>::value;

 RT  n0 =  0;
 RT  n1 =  7;
 RT  n2 = 21;
 RT  n3 = 14;
 RT  n4 =-10;
 RT  n5 =  5;
 RT  n6 = 20;
 RT  n7 = -1;
 RT  n8 =  3;

 Point_3 p1( n1, n2, n3, n1);
 Point_3 p2( n4, n5, n6, n5);
 Point_3 p3( n2, n8, n2, n8);
 Point_3 p4 = p3 + (p2 - p1);

 Plane_3 pl1( p1, p2, p3);
 assert( pl1.has_on( pl1.to_3d( pl1.to_2d( pl1.point() ))) || nonexact );
 assert( pl1.has_on( pl1.to_3d( pl1.to_2d( p4 ))) || nonexact );
 assert( p1 == pl1.to_3d( pl1.to_2d( p1)) || nonexact );
 assert( p2 == pl1.to_3d( pl1.to_2d( p2)) || nonexact );
 assert( p3 == pl1.to_3d( pl1.to_2d( p3)) || nonexact );
 assert( p4 == pl1.to_3d( pl1.to_2d( p4)) || nonexact );

 std::cout << '.';

 Plane_3 pl2( p2, p1, p3);
 assert( pl2.has_on( pl2.to_3d( pl2.to_2d( pl2.point() ))) || nonexact );
 assert( pl2.has_on( pl2.to_3d( pl2.to_2d( p4 ))) || nonexact );
 assert( p1 == pl2.to_3d( pl2.to_2d( p1)) || nonexact );
 assert( p2 == pl2.to_3d( pl2.to_2d( p2)) || nonexact );
 assert( p3 == pl2.to_3d( pl2.to_2d( p3)) || nonexact );
 assert( p4 == pl2.to_3d( pl2.to_2d( p4)) || nonexact );

 Point_3 p5( n2, n8, n0, n7);
 Point_3 p6( n4, n5, n0, n8);
 Plane_3 pl3( p4, p5, p6);
 assert( p4 == pl3.to_3d( pl3.to_2d( p4)) || nonexact );
 assert( p5 == pl3.to_3d( pl3.to_2d( p5)) );
 assert( p6 == pl3.to_3d( pl3.to_2d( p6)) || nonexact );
 Plane_3 pl4( p4, p6, p5);
 assert( p4 == pl4.to_3d( pl4.to_2d( p4)) || nonexact );
 assert( p5 == pl4.to_3d( pl4.to_2d( p5)) );
 assert( p6 == pl4.to_3d( pl4.to_2d( p6)) || nonexact );

 Point_3 p7 = CGAL::midpoint( p1, p2);
 Point_3 p8 = CGAL::midpoint( p3, p3 + (p2-p1) );
 Point_3 p9 = CGAL::midpoint( p1, p3);
 Point_3 p10= CGAL::midpoint( p2, p2 + (p3-p1) );
 assert( pl1.has_on( p7 ));
 assert( pl1.has_on( p8 ));
 assert( pl1.has_on( p9 ));
 assert( pl1.has_on( p10));

 std::cout << '.';

 Point_2 pp7 = pl1.to_2d( p7 );
 Point_2 pp8 = pl1.to_2d( p8 );
 Point_2 pp9 = pl1.to_2d( p9 );
 Point_2 pp10= pl1.to_2d( p10);
 CGAL::Segment_2<R> sp1( pp7, pp8);
 CGAL::Segment_2<R> sp2( pp9, pp10);
 Point_2 pp;
 assert( CGAL::assign( pp, CGAL::intersection( sp1, sp2)) );
 assert( sp1.has_on( pp) || nonexact );
 assert( sp2.has_on( pp) || nonexact );
 Point_3 p = pl1.to_3d( pp);
 assert( pl1.has_on( p ) || nonexact );
 CGAL::Segment_3<R> s1( p7, p8);
 CGAL::Segment_3<R> s2( p9, p10);
 assert( s1.has_on( p) || nonexact );
 assert( s2.has_on( p) || nonexact );

 std::cout << '.' << std::endl;
 return true;
}
#endif // CGAL__TEST_MF_PLANE_3_TO_2D_H
