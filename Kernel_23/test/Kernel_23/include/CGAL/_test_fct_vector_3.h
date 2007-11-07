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
 

#ifndef CGAL__TEST_FCT_VECTOR_3_H
#define CGAL__TEST_FCT_VECTOR_3_H

template <class R>
bool
_test_fct_vector_3(const R& )
{
 std::cout << "Testing functions Vector_3" ;

 typedef typename  R::RT    RT;
 typedef typename  R::FT    FT;

 RT  n0(  0 );
 RT  n1( 12 );
 RT  n2( -4 );
 RT  n3(  6 );
 RT  n4(  2 );
 RT  n5(  9 );
 RT  n6(-18 );
 RT  n7( 36 );
 RT  n8(  3 );
 RT  n9( 15 );
 RT n10( -8 );
 RT n11( 24 );
 RT n12( 54 );


 CGAL::Vector_3<R>  v0(CGAL::NULL_VECTOR);  // ( 0, 0, 0)
 CGAL::Vector_3<R>  v1(n1, n2, n3, n4);    // ( 6,-2, 3)
 CGAL::Vector_3<R>  v2(n5, n6, n7, n8);    // ( 3,-6,12)
 CGAL::Vector_3<R>  v3(n5, n10, n9);       // ( 9,-8,15)
 CGAL::Vector_3<R>  v4(n8, -n2, -n5);      // ( 3, 4,-9)
 CGAL::Vector_3<R> mv4(-n8, n2, n5);       // (-3,-4, 9)

 CGAL::Vector_3<R>  v10( n8,-n2,-n5, n4);  // (1.5,2,-4.5)
 CGAL::Vector_3<R>  v11(-n6, n11,-n12, n8);// ( 6, 8, -18)
 CGAL::Vector_3<R>  v12(n1, n2, -n3, n4);  // ( 6,-2, -3)

 CGAL_test_assert( orientation(v0, v0, v0) == CGAL::COPLANAR );
 CGAL_test_assert( orientation(v1, v1, v1) == CGAL::COPLANAR );
 CGAL_test_assert( orientation(v1, v1, v2) == CGAL::COPLANAR );
 CGAL_test_assert( orientation(v1, v2, v1) == CGAL::COPLANAR );
 CGAL_test_assert( orientation(v2, v1, v1) == CGAL::COPLANAR );
 CGAL_test_assert( orientation(v1, v2, v3) == CGAL::COPLANAR );
 CGAL_test_assert( orientation(v1, v2, v12) == CGAL::POSITIVE );
 CGAL_test_assert( orientation(v1, v12, v2) == CGAL::NEGATIVE );

 CGAL_test_assert( determinant(v0, v0, v0) == 0 );
 CGAL_test_assert( determinant(v1, v1, v1) == 0 );
 CGAL_test_assert( determinant(v1, v1, v2) == 0 );
 CGAL_test_assert( determinant(v1, v2, v1) == 0 );
 CGAL_test_assert( determinant(v2, v1, v1) == 0 );
 CGAL_test_assert( determinant(v1, v2, v3) == 0 );
 CGAL_test_assert( determinant(v1, v2, v12) == 180 );
 CGAL_test_assert( determinant(v1, v12, v2) == -180 );

 CGAL_test_assert( v1 + v2 == v3 );
 CGAL_test_assert( v1 - v2 == v4 );
 CGAL_test_assert( v3 - v1 == v2 );
 CGAL_test_assert( v3 - v2 == v1 );
 CGAL_test_assert( v4 + v2 == v1 );
 CGAL_test_assert( v4 + v2 == v3 - v2);

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

 CGAL_test_assert( v1.squared_length() == FT(49) );
 CGAL_test_assert( v1 * v2 == FT(66) );
 CGAL_test_assert( v1 * v0 == FT(0) );
 CGAL_test_assert( CGAL::Vector_3<R>( n1, n2, n3) == v1 * RT(2));
 CGAL_test_assert( CGAL::Vector_3<R>( n5, n6, n7) == v2 * RT(3));
 CGAL_test_assert( CGAL::Vector_3<R>( n1, n2, n3) == RT(2) * v1);
 CGAL_test_assert( CGAL::Vector_3<R>( n5, n6, n7) == RT(3) * v2);
 CGAL_test_assert( CGAL::Vector_3<R>( n1, n2, n3) == v1 * FT(2));
 CGAL_test_assert( CGAL::Vector_3<R>( n5, n6, n7) == v2 * FT(3));
 CGAL_test_assert( CGAL::Vector_3<R>( n1, n2, n3) == FT(2) * v1);
 CGAL_test_assert( CGAL::Vector_3<R>( n5, n6, n7) == FT(3) * v2);
 CGAL_test_assert( v2 / RT(3) == CGAL::Vector_3<R>( RT(1), -n4, -n2) );
 CGAL_test_assert( (v2 * RT(3)) / RT(3) == v2 );
 CGAL_test_assert( (v2 / RT(3)) * RT(3) == v2 );

 CGAL_test_assert( (v4 / (FT(n1)/FT(n3))) == v10 );
 CGAL_test_assert( (v4 * (FT(n1)/FT(n3))) == v11 );
 CGAL_test_assert( (v4 / (FT(n3)/FT(n1))) == v11 );
 CGAL_test_assert( (v4 * (FT(n3)/FT(n1))) == v10 );

 std::cout << '.';

 CGAL_test_assert( v2.cartesian(0) == v2[0] );
 CGAL_test_assert( v2.cartesian(1) == v2[1] );
 CGAL_test_assert( v2.cartesian(2) == v2[2] );


 CGAL::Point_3<R> p0(CGAL::ORIGIN);
 CGAL::Point_3<R> p1 = CGAL::ORIGIN + v1;
 CGAL::Point_3<R> p2 = CGAL::ORIGIN + v2;
 CGAL::Point_3<R> p3 = CGAL::ORIGIN + v3;

 CGAL_test_assert( CGAL::ORIGIN + v2 == CGAL::Point_3<R>( n5, n6, n7, n8) );
 CGAL_test_assert( CGAL::ORIGIN - v2 == CGAL::Point_3<R>( -n5, -n6, -n7, n8) );
 CGAL_test_assert( p1 - p1 == v0 );
 CGAL_test_assert( p1 - p0 == p1 - CGAL::ORIGIN);
 CGAL_test_assert( p0 - p1 == CGAL::ORIGIN - p1);
 CGAL_test_assert( p1 - p2 == v4 );
 CGAL_test_assert( p2 + v4 == p1 );
 CGAL_test_assert( p3 - v1 == p2 );
 CGAL_test_assert( p3 - p1 == v2 );

 std::cout << '.';

 CGAL::Vector_3<R> evx( n5, n0, n0, n5);
 CGAL::Vector_3<R> evy( n0, n5, n0, n5);
 CGAL::Vector_3<R> evz( n0, n0, n5, n5);
 CGAL_test_assert( CGAL::cross_product(evx,evy) == evz );
 CGAL_test_assert( CGAL::cross_product(evy,evx) == - evz );
 RT ns1(90);
 RT ns2(189);
 RT ns3(60);
 RT ns4(126);
 CGAL::Vector_3<R> vs1(-n6,ns2,ns1,-n8);   // (-6,-63,-30)
 CGAL_test_assert( CGAL::cross_product(v1,v2) == vs1 );
 CGAL::Vector_3<R> vs2(n1,ns4,ns3,n4);     // (6,63,30)
 CGAL_test_assert( CGAL::cross_product(v1,v4) == vs2 );
 CGAL_test_assert( CGAL::cross_product(v1,mv4) == - vs2 );

 {
  CGAL::Point_3<R> p0(4,2,1,1);
  CGAL::Point_3<R> p1(2,4,1,1);
  CGAL::Segment_3<R> s1(p0, p1);
  CGAL::Ray_3<R> r1(p0, p1);
  CGAL::Line_3<R> l1(p0, p1);
  CGAL::Vector_3<R> v8(p0, p1);
  CGAL::Vector_3<R> v9(s1);
  CGAL::Vector_3<R> v10(r1);
  CGAL::Vector_3<R> v11(l1);

  CGAL_test_assert( v8 == (p1-p0) );
  CGAL_test_assert( v8 == v9 );
  CGAL_test_assert( v10.direction() == v8.direction() );
  CGAL_test_assert( v11.direction() == v8.direction() );
 }

 std::cout << "done" << std::endl;
 return true;
}
#endif // CGAL__TEST_FCT_VECTOR_3_H
