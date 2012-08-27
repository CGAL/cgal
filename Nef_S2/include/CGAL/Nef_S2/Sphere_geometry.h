// Copyright (c) 1997-2002  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Michael Seel       <seel@mpi-inf.mpg.de>
//                 Peter Hachenberger <hachenberger@mpi-inf.mpg.de>

#ifndef CGAL_SPHERE_GEOMETRY_H
#define CGAL_SPHERE_GEOMETRY_H

#include <CGAL/basic.h>
#include <CGAL/intersection_3.h>
#include <list>

#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 113
#include <CGAL/Nef_2/debug.h>

namespace CGAL {

template <class R> class Sphere_point;
template <class R> class Sphere_segment;
template <class R> class Sphere_triangle;
template <class R> class Sphere_circle;
template <class R> class Sphere_direction;

} //namespace CGAL

#include <CGAL/Nef_S2/Sphere_point.h>
#include <CGAL/Nef_S2/Sphere_circle.h>
#include <CGAL/Nef_S2/Sphere_direction.h>
#include <CGAL/Nef_S2/Sphere_segment.h>
#include <CGAL/Nef_S2/Sphere_triangle.h>
#include <CGAL/Nef_S2/sphere_predicates.h>

namespace CGAL {

template <typename R_>
struct Positive_halfsphere_geometry {

typedef R_                      R;
typedef CGAL::Sphere_point<R>   Point_2;
typedef CGAL::Sphere_segment<R> Segment_2;

int axis;

Positive_halfsphere_geometry() : axis(2) {}
Positive_halfsphere_geometry(int check_sphere) : axis(check_sphere) {}

Point_2 source(const Segment_2& s) const
{ return s.source(); }
Point_2 target(const Segment_2& s) const
{ return s.target(); }
Segment_2 construct_segment(const Point_2& p, const Point_2& q) const
{ return Segment_2(p,q); }

void xz_pi_half_rotate(Point_2& p) const
{ p = Point_2(-p.hz(),p.hy(),p.hx()); }

void zx_pi_half_rotate(Point_2& p) const
{ p = Point_2(-p.hz(),p.hy(),p.hx()); }

void xy_pi_half_rotate(Point_2& p) const
{ p = Point_2(-p.hy(),p.hx(),p.hz()); }

int orientation(const Point_2& p1, const Point_2& p2,
                const Point_2& p3) const { 
  
  int sor = CGAL::spherical_orientation(p1,p2,p3); 
  if (sor) return sor;
  Point_2 pp1(p1), pp2(p2), pp3(p3); 
  switch(axis) {
  case 0:
    if ( !( p1.hx() == 0 && p2.hx() == 0 && p3.hx() == 0) ) return sor;
    if ( p1.hz()<0 ) zx_pi_half_rotate(pp1);
    if ( p2.hz()<0 ) zx_pi_half_rotate(pp2);
    if ( p3.hz()<0 ) zx_pi_half_rotate(pp3);
    return CGAL::spherical_orientation(pp1,pp2,pp3);
    break;   
  case 1:
    if ( !( p1.hy() == 0 && p2.hy() == 0 && p3.hy() == 0) ) return sor;
    if ( p1.hx()>0 ) xy_pi_half_rotate(pp1);
    if ( p2.hx()>0 ) xy_pi_half_rotate(pp2);
    if ( p3.hx()>0 ) xy_pi_half_rotate(pp3);
    return CGAL::spherical_orientation(pp1,pp2,pp3);
    break;
  case 2:
    if ( !( p1.hz() == 0 && p2.hz() == 0 && p3.hz() == 0) ) return sor;
    // sor==0 we perturb any point in the xy-plane with x>0
    // by a negative rotation around the y-axis 
    // our perturbation is big :-) we take PI/2 :
    if ( p1.hx()>0 ) xz_pi_half_rotate(pp1);
    if ( p2.hx()>0 ) xz_pi_half_rotate(pp2);
    if ( p3.hx()>0 ) xz_pi_half_rotate(pp3);
    return CGAL::spherical_orientation(pp1,pp2,pp3);
    break;
  }
  return 0;
}

int orientation(const Segment_2& s, const Point_2& p) const
{ return orientation(s.source(),s.target(),p); }

bool is_degenerate(const Segment_2& s) const
{ return s.is_degenerate(); }

int compare_xy(const Point_2& p1, const Point_2& p2) const {
  CGAL_NEF_TRACEN("compare_xy " << axis << ":" << p1 << " / " << p2);
  return CGAL::spherical_compare(p1,p2,axis,+1);
}

Point_2 intersection(const Segment_2& s1, const Segment_2& s2) const
{ if (s1.sphere_circle() != s2.sphere_circle().opposite()) 
    return s1.intersection(s2); 
  CGAL_assertion(s1.target()==s2.target());
  return s1.target();
}

}; // Positive_halfsphere_geometry<R>

template <typename R>
struct Negative_halfsphere_geometry :
  public Positive_halfsphere_geometry<R> {

typedef Positive_halfsphere_geometry<R> Base;
typedef typename Base::Point_2   Point_2;
typedef typename Base::Segment_2 Segment_2;

  using Base::xz_pi_half_rotate;
  using Base::zx_pi_half_rotate;
  using Base::xy_pi_half_rotate;

Negative_halfsphere_geometry() : Base() {}
Negative_halfsphere_geometry(int check_sphere) : Base(check_sphere) {}
  
int orientation(const Point_2& p1, const Point_2& p2,
                const Point_2& p3) const { 

  int sor = CGAL::spherical_orientation(p1,p2,p3); 
  if (sor) return sor;
  Point_2 pp1(p1), pp2(p2), pp3(p3); 
  switch(((Base*) this)->axis) {
  case 0:
    if ( !( p1.hx() == 0 && p2.hx() == 0 && p3.hx() == 0) ) return sor;
    if ( p1.hz()>0 ) zx_pi_half_rotate(pp1);
    if ( p2.hz()>0 ) zx_pi_half_rotate(pp2);
    if ( p3.hz()>0 ) zx_pi_half_rotate(pp3);
    return CGAL::spherical_orientation(pp1,pp2,pp3);
    break;   
  case 1:
    if ( !( p1.hy() == 0 && p2.hy() == 0 && p3.hy() == 0) ) return sor;
    if ( p1.hx()<0 ) xy_pi_half_rotate(pp1);
    if ( p2.hx()<0 ) xy_pi_half_rotate(pp2);
    if ( p3.hx()<0 ) xy_pi_half_rotate(pp3);
    return CGAL::spherical_orientation(pp1,pp2,pp3);
    break;
  case 2:
    if ( !( p1.hz() == 0 && p2.hz() == 0 && p3.hz() == 0) ) return sor;
    // sor==0 we perturb any point in the xy-plane with x>0
    // by a negative rotation around the y-axis 
    // our perturbation is big :-) we take PI/2 :
    if ( p1.hx()<0 ) xz_pi_half_rotate(pp1);
    if ( p2.hx()<0 ) xz_pi_half_rotate(pp2);
    if ( p3.hx()<0 ) xz_pi_half_rotate(pp3);
    return CGAL::spherical_orientation(pp1,pp2,pp3);
    break;
  }
  return 0;
}

int orientation(const Segment_2& s, const Point_2& p) const
{ return orientation(s.source(),s.target(),p); }

int compare_xy(const Point_2& p1, const Point_2& p2) const
{ return CGAL::spherical_compare(p1,p2,this->axis,-1); }

}; // Negative_halfsphere_geometry<R>

template <typename R_>
struct Sphere_geometry {

typedef R_                          R;
typedef typename R_::RT             RT;
typedef typename R_::FT             FT;
typedef CGAL::Sphere_point<R>       Sphere_point;
typedef CGAL::Sphere_segment<R>     Sphere_segment;
typedef CGAL::Sphere_circle<R>      Sphere_circle;
typedef CGAL::Sphere_direction<R>   Sphere_direction;
typedef CGAL::Sphere_triangle<R>    Sphere_triangle;
typedef typename R::Point_3                  Point_3;
typedef typename R::Plane_3                  Plane_3;
typedef typename R::Aff_transformation_3     Aff_transformation_3;
typedef CGAL::Positive_halfsphere_geometry<R> Positive_halfsphere_geometry;
typedef CGAL::Negative_halfsphere_geometry<R> Negative_halfsphere_geometry;

Sphere_point source(const Sphere_segment& s) const
{ return s.source(); }

Sphere_point target(const Sphere_segment& s) const
{ return s.target(); }

Sphere_segment construct_segment(const Sphere_point& p, 
                                 const Sphere_point& q) const
{ return Sphere_segment(p,q); }

Sphere_segment construct_segment(const Sphere_point& p, 
                                 const Sphere_point& q,
                                 const Plane_3& h) const
{ return Sphere_segment(p,q,Sphere_circle(h)); }

Plane_3 affine_representation(const Plane_3& h, const Point_3& p) const
{ RT wp = p.hw();
  return Plane_3(wp*h.a(),wp*h.b(),wp*h.c(),
                 -(p.hx()*h.a() + p.hy()*h.b() + p.hz()*h.c())); }

Plane_3 linear_representation(const Plane_3& h) const
{ return Plane_3(h.a(),h.b(),h.c(),0); }

/*
Positive_halfsphere_geometry PHG;
const Positive_halfsphere_geometry& 
get_positive_halfsphere_geometry() const
{ return PHG; }

Negative_halfsphere_geometry NHG;
const Negative_halfsphere_geometry& 
get_negative_halfsphere_geometry() const
{ return NHG; }
*/

const Positive_halfsphere_geometry& 
get_positive_halfsphere_geometry(int a) const {
  return Positive_halfsphere_geometry(a); 
}

const Negative_halfsphere_geometry& 
get_negative_halfsphere_geometry(int a) const {
  return Negative_halfsphere_geometry(a);
}

};



} //namespace CGAL
#endif //CGAL_SPHERE_GEOMETRY_H
