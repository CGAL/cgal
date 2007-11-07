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
 

#ifndef CGAL__TEST_CLS_VECTOR_3_H
#define CGAL__TEST_CLS_VECTOR_3_H

template <class R>
bool
_test_cls_vector_3(const R& )
{
 std::cout << "Testing class Vector_3" ;

 typedef typename  R::RT    RT;
 typedef typename  R::FT    FT;

 typename R::Vector_3       iv;

 CGAL::Vector_3<R>  v1;
 CGAL::Vector_3<R>  v2(iv);
 CGAL::Vector_3<R>  v0(CGAL::NULL_VECTOR);

 RT  n1( 12 );
 RT  n2( -4 );
 RT  n3(  6 );
 RT  n4(  2 );

 CGAL::Vector_3<R>  v3(n1, n2, n3);
 CGAL::Vector_3<R>  v4(n1, n2, n3, n4);
 CGAL::Vector_3<R>  v5(n1, n2, n3, n4);
 CGAL::Vector_3<R>  v6( v5 );
                   v1 = v4;
 CGAL::Vector_3<R>  v7(-n1, -n2, -n3, -n4);

 std::cout << '.';

 CGAL_test_assert( v3 == CGAL::Vector_3<R>(FT(n1), FT(n2), FT(n3)) );
 CGAL_test_assert( v3 == CGAL::Vector_3<R>(12, -4, 6) );

 CGAL_test_assert( v4 == v5 );
 CGAL_test_assert( v5 == v6 );
 CGAL_test_assert( v4 == v6 );
 CGAL_test_assert( v1 == v6 );
 CGAL_test_assert( v0 == CGAL::NULL_VECTOR);
 CGAL_test_assert( v5 == v7 );

 CGAL_test_assert( v3 != v4 );
 CGAL_test_assert( v0 != v1 );
 CGAL_test_assert( v1 != CGAL::NULL_VECTOR);

 CGAL_test_assert( v3.hx() == n1 );   // don't replace v3
 CGAL_test_assert( v3.hy() == n2 );
 CGAL_test_assert( v3.hz() == n3 );

 CGAL_test_assert( FT( v5.hx()) / FT(v5.hw()) == FT( n1) / FT( n4) );
 CGAL_test_assert( FT( v5.hy()) / FT(v5.hw()) == FT( n2) / FT( n4) );
 CGAL_test_assert( FT( v5.hz()) / FT(v5.hw()) == FT( n3) / FT( n4) );

 CGAL_test_assert( v5.x() == FT( n1) / FT( n4) );
 CGAL_test_assert( v5.y() == FT( n2) / FT( n4) );
 CGAL_test_assert( v5.z() == FT( n3) / FT( n4) );

 std::cout << '.';

 CGAL_test_assert( v3.homogeneous(0) == v3.hx() );  // don't replace v3
 CGAL_test_assert( v3.homogeneous(1) == v3.hy() );
 CGAL_test_assert( v3.homogeneous(2) == v3.hz() );
 CGAL_test_assert( v3.homogeneous(3) == v3.hw() );
 CGAL_test_assert( v6.cartesian(0) == v6.x() );
 CGAL_test_assert( v6.cartesian(1) == v6.y() );
 CGAL_test_assert( v6.cartesian(2) == v6.z() );

 std::cout << '.';

 CGAL_test_assert( v0.dimension() == 3 );

 std::cout << "done" << std::endl;
 return true;
}
#endif // CGAL__TEST_CLS_VECTOR_3_H
