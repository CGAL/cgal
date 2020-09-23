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


#ifndef CGAL__TEST_FCT_CONSTRUCTIONS_3_H
#define CGAL__TEST_FCT_CONSTRUCTIONS_3_H

template <class R>
bool
_test_fct_constructions_3(const R& r)
{
  typedef typename R::RT               RT;
  typedef typename R::Point_3          Point;
  typedef typename R::Weighted_point_3 Weighted_point;
  typedef typename R::Segment_3        Segment;
  typedef typename R::Ray_3        Ray;
  typedef typename R::Plane_3          Plane;
  typedef typename R::Vector_3         Vector;
  typedef typename R::Triangle_3       Triangle;
  typedef typename R::Tetrahedron_3    Tetrahedron;

  RT RT0(0);
  RT RT1(1);
  RT RT2(2);
  RT RT4(4);
  RT RT8(8);

  Point p( RT4, RT8, -RT2, RT2);   // ( 2, 4, -1)
  Point p111 = p + Vector( RT1, RT1, RT1 );
  Point p011 = p + Vector(-RT1, RT1, RT1 );
  Point p101 = p + Vector( RT1,-RT1, RT1 );
  Point p001 = p + Vector(-RT1,-RT1, RT1 );
  Point p000 = p + Vector(-RT1,-RT1,-RT1 );
  Point p100 = p + Vector( RT1,-RT1,-RT1 );
  Point p110 = p + Vector( RT1, RT1,-RT1 );
  Point p010 = p + Vector(-RT1, RT1,-RT1 );

  Point p2   = p + Vector(-RT1, RT0, RT0 );
  Point p3   = p + Vector( RT1, RT0, RT0 );
  Point p4   = p + Vector( RT0, RT1, RT0 );

  Weighted_point wp000(p000);
  Weighted_point wp010(p010);
  Weighted_point wp001(p001);
  Weighted_point wp111(p111);
  Weighted_point wp2(p2);
  Weighted_point wp3(p3);

  Weighted_point wp000_b(p000, 0);
  Weighted_point wp100_b(p100, 4);
  Weighted_point wp010_b(p010, 4);
  Weighted_point wp001_b(p001, 4);

  // midpoint
  assert( CGAL::midpoint( p111, p000) == p);
  assert( CGAL::midpoint( p110, p001) == p);
  assert( CGAL::midpoint( p010, p101) == p);
  assert( CGAL::midpoint( p100, p011) == p);

  // circumcenter
  assert( CGAL::circumcenter( p111, p001, p010, p000) == p);
  assert( CGAL::circumcenter( p101, p001, p010, p100) == p);
  assert( CGAL::circumcenter( p001, p000, p110, p100) == p);
  assert( CGAL::circumcenter( Tetrahedron(p111, p001, p010, p000) ) == p);
  assert( CGAL::circumcenter( Tetrahedron(p101, p001, p010, p100) ) == p);
  assert( CGAL::circumcenter( Tetrahedron(p001, p000, p110, p100) ) == p);

  assert( CGAL::circumcenter( p2, p2 ) == p2);
  assert( CGAL::circumcenter( p2, p3 ) == CGAL::midpoint(p2, p3) );
  assert( CGAL::circumcenter( p2, p3, p4 ) == p);
  assert( CGAL::circumcenter( Triangle(p2, p3, p4) ) == p);

  // centroid
  Point p_11 = p + Vector(RT0, RT1, RT1);
  assert( CGAL::centroid( p111, p010, p101, p000) == p);
  assert( CGAL::centroid( p111, p_11, p011 ) == p_11);
  assert( CGAL::centroid( Tetrahedron(p111, p010, p101, p000)) == p);
  assert( CGAL::centroid( Triangle(p111, p_11, p011) ) == p_11);

  // barycenter
  assert( CGAL::barycenter( p111, 1, p011 ) == p111 );
  assert( CGAL::barycenter( p111, 0, p011 ) == p011 );
  assert( CGAL::barycenter( p111, 1, p011, 1 ) == CGAL::midpoint(p111, p011) );
  assert( CGAL::barycenter( p111, 1, p_11, 0, p011) == p111);
  assert( CGAL::barycenter( p111, 0, p_11, 1, p011) == p_11);
  assert( CGAL::barycenter( p111, 0, p_11, 0, p011) == p011);
  assert( CGAL::barycenter( p111, 1, p_11, 1, p011, 1) == CGAL::centroid(p111, p_11, p011) );
  assert( CGAL::barycenter( p111, 1, p010, 0, p_11, 0, p000) == p111);
  assert( CGAL::barycenter( p111, 0, p010, 1, p_11, 0, p000) == p010);
  assert( CGAL::barycenter( p111, 0, p010, 0, p_11, 1, p000) == p_11);
  assert( CGAL::barycenter( p111, 0, p010, 0, p_11, 0, p000) == p000);
  assert( CGAL::barycenter( p111, 1, p010, 1, p_11, 1, p000, 1) == CGAL::centroid(p111, p010, p_11, p000) );

  // orthogonal_vector
  Point p0(RT0, RT0, RT0), px1(RT1, RT0, RT0), py1(RT0, RT1, RT0);
  Vector vz1(RT0, RT0, RT1);
  Vector orth = orthogonal_vector(p0, px1, py1);
  assert( (vz1 * orth) > 0 );
  assert( parallel(Segment(p0, p0+orth), Segment(p0, p0+vz1)) );

  orth = orthogonal_vector(Plane(p0, px1, py1));
  assert( (vz1 * orth) > 0 );
  assert( parallel(Segment(p0, p0+orth), Segment(p0, p0+vz1)) );

  // weighted circumcenter
  assert( CGAL::weighted_circumcenter( wp2, wp3 ) == CGAL::midpoint(p2, p3) );
  assert( CGAL::weighted_circumcenter( wp000, wp001, wp111) == p);
  assert( CGAL::weighted_circumcenter( wp000, wp001, wp010, wp111) == p);

  assert( CGAL::weighted_circumcenter( wp000_b, wp100_b) == wp000_b);
  assert( CGAL::weighted_circumcenter( wp000_b, wp100_b, wp010_b) == wp000_b);
  assert( CGAL::weighted_circumcenter( wp000_b, wp100_b, wp010_b, wp001_b) == wp000_b);

    // projected point
  Ray ray(Point(0,0,0), Point (1,1,0));
  Segment s(Point(0,0,0), Point (1,1,0));
  assert( r.construct_projected_point_3_object()(ray, Point(-1,0,0)) == Point(0,0,0));
  assert( r.construct_projected_point_3_object()(s, Point(-1,0,0)) == Point(0,0,0));
  assert( r.construct_projected_point_3_object()(s, Point(2,0,0)) == Point(1,1,0));
  return true;
}

#endif // CGAL__TEST_FCT_CONSTRUCTIONS_3_H
