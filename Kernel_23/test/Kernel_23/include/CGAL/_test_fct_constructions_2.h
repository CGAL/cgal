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
 

#ifndef CGAL__TEST_FCT_CONSTRUCTIONS_2_H
#define CGAL__TEST_FCT_CONSTRUCTIONS_2_H

template <class R>
bool
_test_fct_constructions_2(const R&)
{
  typedef typename R::RT     RT;
  typedef CGAL::Point_2<R>    Point;
  typedef CGAL::Triangle_2<R> Triangle;
  typedef CGAL::Vector_2<R>   Vector;

  RT RT0(0);
  RT RT1(1);
  RT RT2(2);
  RT RT4(4);
  RT RT8(8);

  Point p( RT4, RT8, RT2);   // ( 2, 4)
  Point pne = p + Vector( RT1, RT1 );
  Point pnw = p + Vector(-RT1, RT1 );
  Point pse = p + Vector( RT1,-RT1 );
  Point psw = p + Vector(-RT1,-RT1 );

  // midpoint
  assert( CGAL::midpoint( pne, psw) == p);
  assert( CGAL::midpoint( pnw, pse) == p);

  // circumcenter
  assert( CGAL::circumcenter( pne, pne ) == pne);
  assert( CGAL::circumcenter( pne, pse ) == CGAL::midpoint(pne, pse) );
  assert( CGAL::circumcenter( pne, pse, pnw) == p);
  assert( CGAL::circumcenter( psw, pse, pnw) == p);
  assert( CGAL::circumcenter( Triangle(pne, pse, pnw)) == p);
  assert( CGAL::circumcenter( Triangle(psw, pse, pnw)) == p);

  // centroid
  Point pe = p + Vector(RT1, RT0);
  assert( CGAL::centroid( pne, pse, pe) == pe);
  assert( CGAL::centroid( pne, psw, pse, pnw) == p);
  assert( CGAL::centroid( Triangle(pne, pse, pe)) == pe);

  // barycenter
  assert( CGAL::barycenter( pne, 1, pe ) == pne );
  assert( CGAL::barycenter( pne, 0, pe ) == pe );
  assert( CGAL::barycenter( pne, 1, pe, 1 ) == CGAL::midpoint(pne, pe) );
  assert( CGAL::barycenter( pne, 1, pse, 0, pe) == pne);
  assert( CGAL::barycenter( pne, 0, pse, 1, pe) == pse);
  assert( CGAL::barycenter( pne, 0, pse, 0, pe) == pe);
  assert( CGAL::barycenter( pne, 1, pse, 1, pe, 1) == pe);
  assert( CGAL::barycenter( pne, 1, psw, 0, pse, 0, pnw) == pne);
  assert( CGAL::barycenter( pne, 0, psw, 1, pse, 0, pnw) == psw);
  assert( CGAL::barycenter( pne, 0, psw, 0, pse, 1, pnw) == pse);
  assert( CGAL::barycenter( pne, 0, psw, 0, pse, 0, pnw) == pnw);
  assert( CGAL::barycenter( pne, 1, psw, 1, pse, 1, pnw, 1) == p);

  // general position intersection point

  return true;
}

#endif // CGAL__TEST_FCT_CONSTRUCTIONS_2_H
