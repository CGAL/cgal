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

template <class NT>
bool
_test_cls_quotient(const NT& )
{
  std::cout << "testing class quotient" ;

  typedef CGAL::Quotient<NT>  QT;

  NT n0(0);
  NT n1(1);
  NT n2(2);
  NT n3(3);
  NT n4(4);
  NT n5(5);
  NT n6(6);
  QT q0(0);
  QT q1(1);
  QT q2(2);
  QT q1_2(n1,n2);
  QT q1_3(n1,n3);
  QT q3_4(n3,n4);
  QT q3_6(n3,n6);
  QT q1_4(n1,n4);

  std::cout << '.';

  assert( q0 + q0 == q0 );
  assert( q1 + q0 == q1 );
  assert( q0 + q1 == q1 );
  assert( q1 + n0 == q1 );
  assert( n0 + q1 == q1 );
  assert( n1 + q0 == q1 );
  assert( q0 + n1 == q1 );
  assert( q1 + q1 == q2 );
  assert( q2 - q1 == q1 );
  assert( q2 - n1 == q1 );
  assert( q0 + q0 != q1 );

  assert( q1_2 + q1_2 == n1 );
  assert( q1_2 - q1_4 == q1_4 );
  assert( q1_2 * q1_2 == q1_4 );
  assert( q1_2 == q1 / q2 );
  assert( q1_2 == n1 / q2 );
  assert( q1_3 == q1 / n3 );
  assert( q1_2 + q1_3 == QT(n5,n6) );
  assert( q1_2 == q3_6 );
  assert( q1_2 * n5 == n2 + q1_4 + q1_4 );

  std::cout << '.';

  assert( q1_2 < q3_4 );
  assert( q1_3 < q3_4 );
  assert( n1 < QT(n6,n4) );
  assert( q3_6 < n2 );

  std::cout << q1_2 << "done" << std::endl;
  return true;
}

