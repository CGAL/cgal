// ======================================================================
//
// Copyright (c) 1997-2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 2000, September 20
//
// file          : include/CGAL/Width_default_traits_3.h
// package       : Width_3 (1.6)
// maintainer    : Thomas Herrmann <herrmann@ifor.math.ethz.ch>
// chapter       : Geometric Optimisation
//
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Thomas Herrmann
// coordinator   : ETH Zuerich (Bernd Gaertner <gaertner@inf.ethz.ch>)
//
// implementation: 3D Width of a Point Set
// ======================================================================

#ifndef CGAL_WIDTH_DEFAULT_TRAITS_3_H
#define CGAL_WIDTH_DEFAULT_TRAITS_3_H

#include <CGAL/Point_3.h>
#include <CGAL/Vector_3.h>
#include <CGAL/Plane_3.h>
#include <CGAL/dd_geo/chull_traits_3.h>

CGAL_BEGIN_NAMESPACE

template <class _RepCls>
class Width_default_traits_3 {
 public:
  typedef _RepCls R;
  typedef typename R::RT RT;
  typedef CGAL::Point_3<R> Point;
  typedef CGAL::Vector_3<R> Vector;
  typedef CGAL::Plane_3<R> Plane;
  typedef chull_traits_3<R> ChullTraits;

  Width_default_traits_3(){}
  ~Width_default_traits_3(){}

  RT get_hx(const Point& p) const {
    return p.hx();
  }
  RT get_hy(const Point& p) const {
    return p.hy();
  }
  RT get_hz(const Point& p) const {
    return p.hz();
  }
  RT get_hw(const Point& p) const {
    return p.hw();
  }
  void get_point_coordinates(const Point& p, 
			     RT& px, RT& py, RT& pz, RT& ph) const {
    px=get_hx(p);
    py=get_hy(p);
    pz=get_hz(p);
    ph=get_hw(p);
  }
  RT get_a(const Plane& f) const {
    return f.a();
  }
  RT get_b(const Plane& f) const {
    return f.b();
  }
  RT get_c(const Plane& f) const {
    return f.c();
  }
  RT get_d(const Plane& f) const {
    return f.d();
  }
  void get_plane_coefficients(const Plane& f, 
			      RT& a, RT& b, RT& c, RT& d) const {
    a=get_a(f);
    b=get_b(f);
    c=get_c(f);
    d=get_d(f);
  }

  Point make_point(const RT& hx, const RT& hy, const RT& hz, 
		   const RT& hw) const  {
    return Point(hx,hy,hz,hw);
  }
  Plane make_plane(const RT& a, const RT& b, const RT& c, const RT& d) const {
    return Plane(a,b,c,d);
  }
  Vector make_vector(const RT& a, const RT& b, const RT& c) const {
    return Vector(a,b,c);
  }

  void inverse_normal(Vector& nor) {
    nor = -nor;
  }
  void opposite_plane(Plane& pln) {
    pln=pln.opposite();
  }
  Vector orthogonal_vector(Plane& pln) {
    return pln.orthogonal_vector();
  }
};

CGAL_END_NAMESPACE

#endif //CGAL_WIDTH_DEFAULT_TRAITS_3_H
