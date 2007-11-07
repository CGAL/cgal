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
 

#ifndef CGAL__TEST_FCT_VECTOR_2_H
#define CGAL__TEST_FCT_VECTOR_2_H

template <class R>
bool
_test_fct_vector_2(const R& )
{
 std::cout << "Testing functions Vector_2" ;

 typedef typename  R::RT    RT;
 typedef typename  R::FT    FT;

 RT  n1( 12 );
 RT  n2( -4 );
 RT  n3(  6 );
 RT  n4(  2 );
 RT  n5(  9 );
 RT  n6(-18 );
 RT  n8(  3 );
 RT n10( -8 );
 RT n11( 24 );

 CGAL::Vector_2<R>  v0(CGAL::NULL_VECTOR);   // ( 0, 0)
 CGAL::Vector_2<R>  v1(n1, n2, n4);         // ( 6,-2)
 CGAL::Vector_2<R>  v2(n5, n6, n8);         // ( 3,-6)
 CGAL::Vector_2<R>  v3(n5, n10);            // ( 9,-8)
 CGAL::Vector_2<R>  v4(n8, -n2 );           // ( 3, 4)
 CGAL::Vector_2<R> mv4(-n8, n2 );           // (-3,-4)

 RT  n15( -5 );
 CGAL::Vector_2<R>  v8( n1,  n3, -n4);      // (-6,-3)
 CGAL::Vector_2<R>  v9( RT(0), n15);        // ( 0,-5)
 CGAL::Vector_2<R>  v10( n8,-n2, n4);       // (1.5,2)
 CGAL::Vector_2<R>  v11(-n6, n11, n8);      // ( 6, 8)

 CGAL_test_assert( orientation(v0, v0) == CGAL::COLLINEAR );
 CGAL_test_assert( orientation(v1, v1) == CGAL::COLLINEAR );
 CGAL_test_assert( orientation(v1, v0) == CGAL::COLLINEAR );
 CGAL_test_assert( orientation(v1, v2) == CGAL::RIGHT_TURN );
 CGAL_test_assert( orientation(v2, v1) == CGAL::LEFT_TURN );

 CGAL_test_assert( determinant(v0, v0) == 0 );
 CGAL_test_assert( determinant(v1, v1) == 0 );
 CGAL_test_assert( determinant(v1, v0) == 0 );
 CGAL_test_assert( determinant(v1, v2) == -30 );
 CGAL_test_assert( determinant(v2, v1) == 30 );

 CGAL_test_assert( v1 + v2 == v3 );
 CGAL_test_assert( v1 - v2 == v4 );
 CGAL_test_assert( v3 - v1 == v2 );
 CGAL_test_assert( v3 - v2 == v1 );
 CGAL_test_assert( v4 + v2 == v1 );
 CGAL_test_assert( v4 + v2 == v3 - v2);
 CGAL_test_assert( v9 == v1 + v8 );

 std::cout << '.';

 CGAL_test_assert( (-(- v1)) == v1 );
 CGAL_test_assert( -v4 == mv4);
 CGAL_test_assert( mv4 == v2 - v1);
 CGAL_test_assert( -( v1 - v2) == mv4);
 CGAL_test_assert( v1 + v0 == v1 );
 CGAL_test_assert( v0 - v4 == mv4 );
 CGAL_test_assert( v0 + v4 == v4 );
 CGAL_test_assert( v2 - v0 == v2 );

 std::cout << '.';

 CGAL_test_assert( v1 * v2 == FT(30) );
 CGAL_test_assert( v1 * v0 == FT(0) );
 CGAL_test_assert( v1.squared_length() == FT(40) );
 CGAL_test_assert( CGAL::Vector_2<R>( n1, n2) == v1 * RT(2));
 CGAL_test_assert( CGAL::Vector_2<R>( n5, n6) == v2 * RT(3));
 CGAL_test_assert( CGAL::Vector_2<R>( n1, n2) == RT(2) * v1);
 CGAL_test_assert( CGAL::Vector_2<R>( n5, n6) == RT(3) * v2);
 CGAL_test_assert( CGAL::Vector_2<R>( n1, n2) == v1 * FT(2));
 CGAL_test_assert( CGAL::Vector_2<R>( n5, n6) == v2 * FT(3));
 CGAL_test_assert( CGAL::Vector_2<R>( n1, n2) == FT(2) * v1);
 CGAL_test_assert( CGAL::Vector_2<R>( n5, n6) == FT(3) * v2);
 CGAL_test_assert( v2 / RT(3) == CGAL::Vector_2<R>( RT(1), -n4) );
 CGAL_test_assert( (v2 * RT(3)) / RT(3) == v2 );
 CGAL_test_assert( (v2 / RT(3)) * RT(3) == v2 );

#ifdef CGAL_VECTOR_MULTIPLICATION_FROM_LEFT
 CGAL_test_assert( CGAL::Vector_2<R>( n1, n2) == RT(2) * v1 );
 CGAL_test_assert( (CGAL::Vector_2<R>&)(RT(2) * v1) == CGAL::Vector_2<R>( n1, n2) );
 CGAL_test_assert( CGAL::Vector_2<R>( n5, n6) == RT(3) * v2 );
#endif // CGAL_VECTOR_MULTIPLICATION_FROM_LEFT

 CGAL_test_assert( (v4 / (FT(n1)/FT(n3))) == v10 );
 CGAL_test_assert( (v4 * (FT(n1)/FT(n3))) == v11 );
 CGAL_test_assert( (v4 / (FT(n3)/FT(n1))) == v11 );
 CGAL_test_assert( (v4 * (FT(n3)/FT(n1))) == v10 );

 std::cout << '.';

 CGAL_test_assert( v2.cartesian(0) == v2[0] );
 CGAL_test_assert( v2.cartesian(1) == v2[1] );

 CGAL::Point_2<R> p0(CGAL::ORIGIN);
 CGAL::Point_2<R> p1 = CGAL::ORIGIN + v1;
 CGAL::Point_2<R> p2 = CGAL::ORIGIN + v2;
 CGAL::Point_2<R> p3 = CGAL::ORIGIN + v3;

 CGAL_test_assert( CGAL::ORIGIN + v2 == CGAL::Point_2<R>( n5, n6, n8) );
 CGAL_test_assert( CGAL::ORIGIN - v2 == CGAL::Point_2<R>( -n5, -n6, n8) );
 CGAL_test_assert( p1 - p1 == v0 );
 CGAL_test_assert( p1 - p0 == p1 - CGAL::ORIGIN);
 CGAL_test_assert( p0 - p1 == CGAL::ORIGIN - p1);
 CGAL_test_assert( p1 - p2 == v4 );
 CGAL_test_assert( p2 + v4 == p1 );
 CGAL_test_assert( p3 - v1 == p2 );
 CGAL_test_assert( p3 - p1 == v2 );

 std::cout << "done" << std::endl;
 return true;
}
#endif // CGAL__TEST_FCT_VECTOR_2_H
