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
 

#ifndef CGAL__TEST_FCT_COPLANAR_3_H
#define CGAL__TEST_FCT_COPLANAR_3_H

template <class R>
bool
_test_fct_coplanar_3(const R& )
{
  typedef typename R::RT     RT;
  typedef CGAL::Point_3<R>   Point;
    RT RT0(0);
    RT RT1(1);
    RT RT2(2);
    RT RT3(3);
    RT RT4(4);
    RT RT6(6);
    RT RT8(8);
  
  Point p = Point( RT1, RT0, RT1, RT2);
  Point q = Point( RT4, RT1, RT2, RT8);
  Point r = Point( RT3, RT1, RT3, RT6);
  Point s = p + (q - r);
  assert( CGAL::coplanar( p,q,r,s));
  assert( CGAL::coplanar_orientation( p,q,s)   != CGAL::COLLINEAR );
  assert( CGAL::coplanar_orientation( p,q,r)   != CGAL::COLLINEAR );
  assert( CGAL::coplanar_orientation( p,q,r,s) == CGAL::NEGATIVE );
  assert( CGAL::coplanar_orientation( p,q,s,r) == CGAL::NEGATIVE );
  assert( CGAL::coplanar_orientation( p,q,s) !=
          CGAL::coplanar_orientation( p,q,r) );
  assert( CGAL::coplanar_orientation( p,q,r,r) == CGAL::POSITIVE );
  assert( CGAL::coplanar_side_of_bounded_circle( p,q,r,s) ==
          CGAL::ON_UNBOUNDED_SIDE );
  assert( CGAL::coplanar_side_of_bounded_circle( q,p,r,s) ==
          CGAL::ON_UNBOUNDED_SIDE );
  assert( CGAL::coplanar_side_of_bounded_circle( r,q,p,s) ==
          CGAL::ON_UNBOUNDED_SIDE );
  s = p + RT2*( q - p);
  assert( CGAL::coplanar_orientation( p,q,r,s) == CGAL::COLLINEAR );
  assert( CGAL::coplanar_orientation( p,q,s)   == CGAL::COLLINEAR );
  s = p - (q - r);
  assert( CGAL::coplanar_orientation( p,q,r,s) == CGAL::POSITIVE );
  assert( CGAL::coplanar_orientation( p,q,s)   != CGAL::COLLINEAR );
  assert( CGAL::coplanar_orientation( p,q,r) == 
          CGAL::coplanar_orientation( p,q,s) );
  p = Point( RT0, RT1, RT1, RT2);
  q = Point( RT1, RT4, RT2, RT8);
  r = Point( RT1, RT3, RT3, RT6);
  s = p + (q - r);
  assert( CGAL::coplanar_orientation( p,q,s)   != CGAL::COLLINEAR );
  assert( CGAL::coplanar_orientation( p,q,r,s) == CGAL::NEGATIVE );
  assert( CGAL::coplanar_orientation( p,q,s,r) == CGAL::NEGATIVE );
  assert( CGAL::coplanar_orientation( p,q,r,r) == CGAL::POSITIVE );
  assert( CGAL::coplanar_orientation( p,q,s) !=
          CGAL::coplanar_orientation( p,q,r) );
  s = p + RT2*( q - p);
  assert( CGAL::coplanar_orientation( p,q,r,s) == CGAL::COLLINEAR );
  assert( CGAL::coplanar_orientation( p,q,s)   == CGAL::COLLINEAR );
  s = p - (q - r);
  assert( CGAL::coplanar_orientation( p,q,r,s) == CGAL::POSITIVE );
  assert( CGAL::coplanar_orientation( p,q,s) ==
          CGAL::coplanar_orientation( p,q,r) );
  p = Point( RT0, RT1, RT1, RT2);
  q = Point( RT1, RT2, RT4, RT8);
  r = Point( RT1, RT3, RT3, RT6);
  s = p + (q - r);
  assert( CGAL::coplanar_orientation( p,q,r,s) == CGAL::NEGATIVE );
  assert( CGAL::coplanar_orientation( p,q,s,r) == CGAL::NEGATIVE );
  assert( CGAL::coplanar_orientation( p,q,r,r) == CGAL::POSITIVE );
  s = p + RT2*( q - p);
  assert( CGAL::coplanar_orientation( p,q,r,s) == CGAL::COLLINEAR );
  s = p - (q - r);
  assert( CGAL::coplanar_orientation( p,q,r,s) == CGAL::POSITIVE );
  return true;
}

#endif // CGAL__TEST_FCT_COPLANAR_3_H
