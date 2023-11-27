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


#ifndef CGAL__TEST_FCT_POINTS_IMPLICIT_SPHERE_H
#define CGAL__TEST_FCT_POINTS_IMPLICIT_SPHERE_H

template <class R>
bool
_test_fct_points_implicit_sphere(const R&)
{
  typedef typename R::RT    RT;
  typedef typename R::FT    FT;
  typedef CGAL::Tetrahedron_3<R>  Tetrahedron;

  const bool nonexact = std::is_floating_point<FT>::value;

  const RT RT0(0);
  const RT RT4(4);
  const FT FT1(1);
  const RT RT1(1);
  CGAL::Point_3<R> org( RT0, RT0, RT0);
  CGAL::Point_3<R>   p( RT0, RT1, RT0);
  CGAL::Point_3<R>   q( RT0, RT0, RT1);
  CGAL::Point_3<R>   r( RT1, RT0, RT0);
  CGAL::Point_3<R>   s( RT0, RT1, RT0);
  CGAL::Point_3<R> out( RT4, RT4, RT1);


  RT sin;
  RT cos;
  RT den;
  RT num = RT1;
  RT eps = RT(100);

  CGAL::rational_rotation_approximation( RT(12), RT(23),
                                         sin, cos, den,
                                         num, eps );
  CGAL::Aff_transformation_3<R> rot_z( cos, sin, RT0, RT0,
                                      -sin, cos, RT0, RT0,
                                       RT0, RT0, den, RT0,
                                                      den );

  CGAL::rational_rotation_approximation( RT(17), RT(4),
                                         sin, cos, den,
                                         num, eps );
  CGAL::Aff_transformation_3<R> rot_y( cos, RT0, sin, RT0,
                                       RT0, den, RT0, RT0,
                                      -sin, RT0, cos, RT0,
                                                      den );

  assert( CGAL::squared_distance(   p, org ) == FT1 );
  CGAL::Point_3<R> tpt = p.transform(rot_z);
  assert( CGAL::squared_distance( tpt, org ) == FT1 );
  p = tpt.transform(rot_z);
  assert( CGAL::squared_distance(   p, org ) == FT1 || nonexact );

  CGAL::rational_rotation_approximation( RT(35), RT(-8),
                                         sin, cos, den,
                                         num, eps );
  CGAL::Aff_transformation_3<R> rot_x( den, RT0, RT0, RT0,
                                       RT0, cos, sin, RT0,
                                       RT0,-sin, cos, RT0,
                                                      den );

  assert( CGAL::squared_distance(   q, org ) == FT1 );
  tpt = q.transform(rot_x);
  assert( CGAL::squared_distance( tpt, org ) == FT1 || nonexact );
  q = tpt.transform(rot_y);
  assert( CGAL::squared_distance(   q, org ) == FT1 || nonexact );

  CGAL::rational_rotation_approximation( RT(9), RT(-8),
                                         sin, cos, den,
                                         num, eps );
  rot_y = CGAL::Aff_transformation_3<R>( cos, RT0, sin, RT0,
                                         RT0, den, RT0, RT0,
                                        -sin, RT0, cos, RT0,
                                                        den );

  assert( CGAL::squared_distance(   r, org ) == FT1 );
  tpt = r.transform(rot_z);
  assert( CGAL::squared_distance( tpt, org ) == FT1 || nonexact);
  r = tpt.transform(rot_y);
  assert( CGAL::squared_distance(   r, org ) == FT1 || nonexact );

  CGAL::rational_rotation_approximation( RT(-19), RT(-1),
                                         sin, cos, den,
                                         num, eps );
  rot_z = CGAL::Aff_transformation_3<R>(  cos, sin, RT0, RT0,
                                         -sin, cos, RT0, RT0,
                                          RT0, RT0, den, RT0,
                                                         den );
  assert( CGAL::squared_distance(   s, org ) == FT1 );
  tpt = s.transform(rot_z);
  assert( CGAL::squared_distance( tpt, org ) == FT1 );
  s = tpt.transform(rot_x);
  assert( CGAL::squared_distance(   s, org ) == FT1 );


  assert( CGAL::side_of_oriented_sphere(p,q,r,s,p) \
          == CGAL::ON_ORIENTED_BOUNDARY);
  assert( CGAL::orientation(p,r,q,s) == CGAL::POSITIVE );
  assert( CGAL::side_of_oriented_sphere(p,r,q,s,org) \
          == CGAL::ON_POSITIVE_SIDE);
  assert( CGAL::side_of_oriented_sphere(p,q,r,s,org) \
          == CGAL::ON_NEGATIVE_SIDE);
  assert( CGAL::side_of_bounded_sphere(p,q,r,s,org) \
          == CGAL::ON_BOUNDED_SIDE);
  assert( CGAL::side_of_bounded_sphere(p,r,q,s,org) \
          == CGAL::ON_BOUNDED_SIDE);
  assert( CGAL::side_of_bounded_sphere(p,q,r,s,q) \
          == CGAL::ON_BOUNDARY);
  assert( CGAL::side_of_bounded_sphere(p,q,r,s,out) \
          == CGAL::ON_UNBOUNDED_SIDE);
  assert( CGAL::side_of_bounded_sphere(p,r,q,s,out) \
          == CGAL::ON_UNBOUNDED_SIDE);

  CGAL::Point_3<R> ex( RT1, RT0, RT0);
  CGAL::Point_3<R> ey( RT0, RT1, RT0);
  CGAL::Point_3<R> ez( RT0, RT0, RT1);
  CGAL::Point_3<R> oz( RT0, RT0, -RT1);
  assert( CGAL::circumcenter(ex, ey, ez, oz) == org );
  assert( CGAL::circumcenter(p,q,r,s) == org || nonexact );
  assert( CGAL::circumcenter(p,r,q,s) == org || nonexact );
  assert( CGAL::circumcenter(Tetrahedron(ex, ey, ez, oz)) == org );
  assert( CGAL::circumcenter(Tetrahedron(p,q,r,s)) == org || nonexact );
  assert( CGAL::circumcenter(Tetrahedron(p,r,q,s)) == org || nonexact );

  CGAL::Vector_3<R>  v( RT(12), RT(4), RT(-4), RT(2) );
  CGAL::Point_3<R>   pt = p + v;
  CGAL::Point_3<R>   qt = q + v;
  CGAL::Point_3<R>   rt = r + v;
  CGAL::Point_3<R>   st = s + v;
  CGAL::Point_3<R>   c = org + v;
  CGAL::Point_3<R>   ot= out + v;

  assert( CGAL::side_of_oriented_sphere(pt,qt,rt,st,pt) \
          == CGAL::ON_ORIENTED_BOUNDARY);
  assert( CGAL::orientation(pt,rt,qt,st) == CGAL::POSITIVE );
  assert( CGAL::side_of_oriented_sphere(pt,rt,qt,st,c) \
          == CGAL::ON_POSITIVE_SIDE);
  assert( CGAL::side_of_oriented_sphere(pt,qt,rt,st,c) \
          == CGAL::ON_NEGATIVE_SIDE);
  assert( CGAL::side_of_bounded_sphere(pt,qt,rt,st,c) \
          == CGAL::ON_BOUNDED_SIDE);
  assert( CGAL::side_of_bounded_sphere(pt,rt,qt,st,c) \
          == CGAL::ON_BOUNDED_SIDE);
  assert( CGAL::side_of_bounded_sphere(pt,qt,rt,st,qt) \
          == CGAL::ON_BOUNDARY);
  assert( CGAL::side_of_bounded_sphere(pt,qt,rt,st,ot) \
          == CGAL::ON_UNBOUNDED_SIDE);
  assert( CGAL::side_of_bounded_sphere(pt,rt,qt,st,ot) \
          == CGAL::ON_UNBOUNDED_SIDE);

  assert( CGAL::circumcenter(pt,qt,rt,st) == c || nonexact );
  assert( CGAL::circumcenter(pt,rt,qt,st) == c || nonexact );
  assert( CGAL::circumcenter(Tetrahedron(pt,qt,rt,st)) == c || nonexact );
  assert( CGAL::circumcenter(Tetrahedron(pt,rt,qt,st)) == c || nonexact );

  // Now test side_of_bounded_sphere(p, q, t).

  CGAL::Point_3<R> p1 (RT(100),RT(100),RT(100),RT(10));
  CGAL::Point_3<R> p2 (RT(-100),RT(-100),RT(-100),RT(10));
  CGAL::Point_3<R> p3 (RT(37),RT(42),RT(56),RT(12));
  CGAL::Point_3<R> pt1 (CGAL::ORIGIN);
  CGAL::Point_3<R> pt2 (RT(-100),RT(-100),RT(100),RT(10));
  CGAL::Point_3<R> pt3 (RT(-100),RT(100),RT(100),RT(10));

  assert( CGAL::side_of_bounded_sphere(p1, p2, pt1) == CGAL::ON_BOUNDED_SIDE);
  assert( CGAL::side_of_bounded_sphere(p1, p2, pt2) == CGAL::ON_BOUNDARY);
  assert( CGAL::side_of_bounded_sphere(p1, p2, p3, pt1) == CGAL::ON_BOUNDED_SIDE);
  assert( CGAL::side_of_bounded_sphere(p1, p2, pt2, pt3) == CGAL::ON_BOUNDARY);

  // Now test squared_radius().

  assert( CGAL::squared_radius(p1, p2, pt2, pt3) == FT(300));
  assert( CGAL::squared_radius(p1, p2, pt2) == FT(300));
  assert( CGAL::squared_radius(p1, p2) == FT(300));
  assert( CGAL::squared_radius(p1) == FT(0));

  return true;
}

#endif // CGAL__TEST_FCT_POINTS_IMPLICIT_SPHERE_H
