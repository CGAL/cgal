// Copyright (c) 1999  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Stefan Schirra
 

#ifndef CGAL__TEST_CLS_POINT_3_H
#define CGAL__TEST_CLS_POINT_3_H

#include <CGAL/Bbox_3.h>
#include <cassert>

template <class R>
bool
_test_cls_point_3(const R& )
{
 std::cout << "Testing class Point_3" ;

 typedef typename  R::RT    RT;
 typedef typename  R::FT    FT;

 typename R::Point_3       ip;
 typedef typename R::Point_3::Cartesian_const_iterator CCI;

 CGAL::Point_3<R>  p1;
 CGAL::Point_3<R>  p2(ip);
 CGAL::Point_3<R>  p0(CGAL::ORIGIN);

 RT  n1(-35 );
 RT  n2( 50 );
 RT  n3(-20 );
 RT  n4(  5 );

 CGAL::Point_3<R>  p3(n1, n2, n3);
 CGAL::Point_3<R>  p4(n1, n2, n3, n4);
 CGAL::Point_3<R>  p5(n1, n2, n3, n4);
 CGAL::Point_3<R>  p6( p5 );
                  p1 = p4;

 std::cout << '.';

 assert( p3 == CGAL::Point_3<R>(FT(n1), FT(n2), FT(n3)) );
 assert( p3 == CGAL::Point_3<R>(-35, 50, -20) );

 assert( p4 == p5 );
 assert( p5 == p6 );
 assert( p4 == p6 );
 assert( p1 == p6 );

 assert( p4 <= p5 );
 assert( p4 >= p5 );
 assert( ! (p4 < p5) );
 assert( ! (p4 > p5) );

 assert( p3 != p4 );
 assert( p0 != p1 );

 assert( p3 < p4);
 assert( p4 > p3);

 assert( p0 == CGAL::ORIGIN);
 assert( p1 != CGAL::ORIGIN);
 // Doesn't work; Point_2::operator== can't be used :(
#ifdef ENHANCED
 assert( CGAL::ORIGIN == p0 );
 assert( CGAL::ORIGIN != p1 );
#endif // ENHANCED

 assert( p3.hx() == n1 );   // don't replace p3
 assert( p3.hy() == n2 );
 assert( p3.hz() == n3 );

 assert( FT(p5.hx()) / FT(p5.hw()) == FT( n1) / FT( n4) );
 assert( FT(p5.hy()) / FT(p5.hw()) == FT( n2) / FT( n4) );
 assert( FT(p5.hz()) / FT(p5.hw()) == FT( n3) / FT( n4) );

 assert( p5.x() == FT( n1) / FT( n4 ) );
 assert( p5.y() == FT( n2) / FT( n4 ) );
 assert( p5.z() == FT( n3) / FT( n4 ) );

 std::cout << '.';

 assert( p3.homogeneous(0) == p3.hx() );  // don't replace p3
 assert( p3.homogeneous(1) == p3.hy() );
 assert( p3.homogeneous(2) == p3.hz() );
 assert( p3.homogeneous(3) == p3.hw() );
 assert( p6.cartesian(0) == p6.x() );
 assert( p6.cartesian(1) == p6.y() );
 assert( p6.cartesian(2) == p6.z() );
 assert( p6[0] == p6.x() );
 assert( p6[1] == p6.y() );
 assert( p6[2] == p6.z() );

 std::cout << '.';

 assert( p0.dimension() == 3 );

 // now we test the Coordinate iterator
  const CGAL::Point_3<R> p(1, 2, 3);

  CCI it = p.cartesian_begin();

  // Default constructor
  CCI itt;

  // Copy constructor
  CCI itc(it);

  assert(itc == it);

  // Assignment
  itt = it;

  // Equality
  assert(itt == it);

  assert(itt - it == 0);

  // Increment
  itt++;

  // Inequality
  assert(itt != it);

  assert(it < itt);
  assert(itt - it == 1);

  itt++;
  assert(itt - it == 2);

  // dereference
  assert(*it == FT(1));

  it++;
  assert(*it == FT(2));
  it++;
  assert(*it == FT(3));

  it++;
  CCI it2 = p.cartesian_end();

  assert(it == it2);
  it--;
  it--;
  it--;
  assert(*it == FT(1));
  it += 1;
  assert(*it == FT(2));
  it -= 1;
  assert(*it == FT(1));

  it2 = it + 1;
  it2--;
  assert(it == it2);
  it++;
  it2 = it - 1;
  it2++;
  assert(it == it2);


 std::cout << '.';

 CGAL::Bbox_3 bb = p3.bbox();
 assert(bb.xmin() <= -35.0);
 assert(bb.xmax() >= -35.0);
 assert(bb.ymin() <= 50.0);
 assert(bb.ymax() >= 50.0);
 assert(bb.zmin() <= -20.0);
 assert(bb.zmax() >= -20.0);

 std::cout << "done" << std::endl;
 return true;
}
#endif // CGAL__TEST_CLS_POINT_3_H
