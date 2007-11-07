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
 

#ifndef CGAL__TEST_CLS_RAY_3_H
#define CGAL__TEST_CLS_RAY_3_H

template <class R>
bool
_test_cls_ray_3(const R& )
{
 std::cout << "Testing class Ray_3" ;

 typedef typename  R::RT    RT;
 typedef typename  R::FT    FT;

 typename R::Ray_3 ir;
 CGAL::Ray_3<R>  r1(ir);

 RT  n1 =  8;
 RT  n2 = 20;
 RT  n3 =  4;
 RT  n4 = 10;
 RT  n5 =  5;
 RT  n6 = 20;
 RT  n7 = -2;

 CGAL::Point_3<R> p1( n1, n2, n3, n3);
 CGAL::Point_3<R> p2( n4, n5, n6, n5);
 CGAL::Point_3<R> p3( n7, n2, n4, n7);

 CGAL::Ray_3<R> r2( p1, p2 );
 CGAL::Ray_3<R> r3( p2, p1 );
 CGAL::Ray_3<R> r4( r2 );
 r1 = r4;
 CGAL::Direction_3<R> dir( p2 - p1 );
 CGAL::Vector_3<R> vec( p2 - p1 );
 CGAL::Line_3<R> l( p1, p2 );
 CGAL::Ray_3<R> r7(p1, dir);

 CGAL_test_assert( r1 == r1 );
 CGAL_test_assert( r4 == r2 );
 CGAL_test_assert( r1 == r4 );
 CGAL_test_assert( r1 == r2 );
 CGAL_test_assert( r7 == r2 );
 CGAL_test_assert( r2 != r3 );

 CGAL_test_assert(   CGAL::parallel(r1, r1) );
 CGAL_test_assert(   CGAL::parallel(r4, r2) );
 CGAL_test_assert(   CGAL::parallel(r2, r3) );

 CGAL::Ray_3<R> r7l(p1, l);
 CGAL::Ray_3<R> r7v(p1, vec);

 CGAL_test_assert(r7 == r7l);
 CGAL_test_assert(r7 == r7v);

 std::cout <<'.';

 CGAL::Ray_3<R> r5 (p3, p3 + (p1 - p3) );
 CGAL_test_assert( r5.has_on( p1 ) );
 CGAL_test_assert( r5.has_on( p3 ) );
 CGAL_test_assert( r5.has_on( p3 + (p1 - p3) ) );
 CGAL_test_assert( r3.has_on( p2 + (p1 - p2) + (p1 - p2) ) );
 CGAL_test_assert( r2.has_on( r2.second_point() ));
 CGAL_test_assert( r5.has_on( r5.second_point() ));
 CGAL_test_assert( r4.has_on( r4.point(1) ));
 CGAL_test_assert( r4.has_on( r4.point(3) ));

 CGAL_test_assert( ! CGAL::parallel(r1, r5) );
 CGAL_test_assert( ! CGAL::parallel(r2, r5) );

 std::cout <<'.';

 CGAL_test_assert( r5.source() == p3 );
 CGAL_test_assert( r2.source() != r3.source() );
 CGAL_test_assert( r7.direction() == dir );
 CGAL_test_assert( r2.direction() == CGAL::Direction_3<R>(r2.point(2) - r2.point(1) ));
 CGAL_test_assert( r2.direction() == r3.opposite().direction() );
 CGAL_test_assert( r2.to_vector().direction() == r3.opposite().to_vector().direction() );
 CGAL_test_assert( r1.supporting_line() == r2.supporting_line() );
 CGAL::Line_3<R> lin(p1,p2);
 CGAL_test_assert( r2.supporting_line() == lin );

 std::cout << '.';

 CGAL::Ray_3<R> r8( p3, dir );
 CGAL::Ray_3<R> r9( p3, -dir );
 CGAL_test_assert( r8.opposite() == r9 );
 CGAL_test_assert( r9.opposite() == r8 );
 CGAL::Ray_3<R> sdeg(p3,p3);
 CGAL_test_assert( sdeg.is_degenerate() );

 std::cout << "done" << std::endl;
 return true;
}

#endif // CGAL__TEST_CLS_RAY_3_H
