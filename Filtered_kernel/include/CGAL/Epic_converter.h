// Copyright (c) 2017  GeometryFactory
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Andreas Fabri, Laurent Rineau

#ifndef CGAL_EPIC_CONVERTER_H
#define CGAL_EPIC_CONVERTER_H


#include <CGAL/config.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

namespace CGAL {

template<typename IK>
class Epic_converter {
  typedef typename Exact_predicates_inexact_constructions_kernel::Point_2 Point_2;
  typedef typename Exact_predicates_inexact_constructions_kernel::Direction_2 Direction_2;
  typedef typename Exact_predicates_inexact_constructions_kernel::Vector_2 Vector_2;
  typedef typename Exact_predicates_inexact_constructions_kernel::Weighted_point_2 Weighted_point_2;
  typedef typename Exact_predicates_inexact_constructions_kernel::Segment_2 Segment_2;
  typedef typename Exact_predicates_inexact_constructions_kernel::Line_2 Line_2;
  typedef typename Exact_predicates_inexact_constructions_kernel::Ray_2 Ray_2;
  typedef typename Exact_predicates_inexact_constructions_kernel::Triangle_2 Triangle_2;
  typedef typename Exact_predicates_inexact_constructions_kernel::Circle_2 Circle_2;
  typedef typename Exact_predicates_inexact_constructions_kernel::Iso_rectangle_2 Iso_rectangle_2;

  typedef typename Exact_predicates_inexact_constructions_kernel::Line_3 Line_3;
  typedef typename Exact_predicates_inexact_constructions_kernel::Plane_3 Plane_3;
  typedef typename Exact_predicates_inexact_constructions_kernel::Triangle_3 Triangle_3;
  typedef typename Exact_predicates_inexact_constructions_kernel::Tetrahedron_3 Tetrahedron_3;
  typedef typename Exact_predicates_inexact_constructions_kernel::Ray_3 Ray_3;
  typedef typename Exact_predicates_inexact_constructions_kernel::Point_3 Point_3;
  typedef typename Exact_predicates_inexact_constructions_kernel::Direction_3 Direction_3;
  typedef typename Exact_predicates_inexact_constructions_kernel::Vector_3 Vector_3;
  typedef typename Exact_predicates_inexact_constructions_kernel::Segment_3 Segment_3;
  typedef typename Exact_predicates_inexact_constructions_kernel::Weighted_point_3 Weighted_point_3;
  typedef typename Exact_predicates_inexact_constructions_kernel::Sphere_3 Sphere_3;
  typedef typename Exact_predicates_inexact_constructions_kernel::Circle_3 Circle_3;
  typedef typename Exact_predicates_inexact_constructions_kernel::Iso_cuboid_3 Iso_cuboid_3;

  typedef typename IK::FT IK_FT;
public:



  std::pair<double,bool> operator()(const typename IK::FT n) const
  {
    double d;
    internal::init_double(d, (IK_FT*)(0));
    if(fit_in_double(n,d)){
      return std::make_pair(d,true);
    }
    return std::make_pair(0,false);
  }

  std::pair<Bbox_2,bool> operator()(const Bbox_2 b) const
  {
    return std::make_pair(b,true);
  }

  std::pair<Bbox_3,bool> operator()(const Bbox_3 b) const
  {
    return std::make_pair(b,true);
  }

  std::pair<Point_2,bool> operator()(const typename IK::Point_2& p) const
  {
    double x, y;
    internal::init_double(x, y, (IK_FT*)(0));
    if(fit_in_double(p.x(),x) && fit_in_double(p.y(),y)){
      return std::make_pair(Point_2(x,y),true);
    }
    return std::make_pair(ORIGIN,false);
  }

  std::pair<Vector_2,bool> operator()(const typename IK::Vector_2& v) const
  {
    double x, y;
    internal::init_double(x, y, (IK_FT*)(0));
    if(fit_in_double(v.x(),x) && fit_in_double(v.y(),y)){
      return std::make_pair(Vector_2(x,y),true);
    }
    return std::make_pair(Vector_2(),false);
  }

  std::pair<Direction_2,bool> operator()(const typename IK::Direction_2& d) const
  {
    double x, y;
    internal::init_double(x, y, (IK_FT*)(0));
    if(fit_in_double(d.dx(),x) && fit_in_double(d.dy(),y)){
      return std::make_pair(Direction_2(x,y),true);
    }
    return std::make_pair(Direction_2(),false);
  }

  std::pair<Weighted_point_2,bool> operator()(const typename IK::Weighted_point_2& wp) const
  {
    std::pair<Point_2,bool> sp = operator()(wp.point());
    std::pair<double,bool> w = operator()(wp.weight());
    if(sp.second && w.second){
      return std::make_pair(Weighted_point_2(sp.first,w.first),true);
    }
    return std::make_pair(Weighted_point_2(),false);
  }

  std::pair<Segment_2,bool> operator()(const typename IK::Segment_2& s) const
  {
    std::pair<Point_2,bool> sp = operator()(s.source());
    if(! sp.second){
      return std::make_pair(Segment_2(),false);
    }
    std::pair<Point_2,bool> tp = operator()(s.target());
    if(! tp.second){
      return std::make_pair(Segment_2(),false);
    }
    return std::make_pair(Segment_2(sp.first,tp.first), true);
  }

  std::pair<Line_2,bool> operator()(const typename IK::Line_2& li) const
  {
    std::pair<double,bool> a = operator()(li.a()),  b = operator()(li.b()) , c = operator()(li.c());
    if(a.second && b.second && c.second){
      return std::make_pair(Line_2(a.first, b.first, c.first),true);
    }
    return std::make_pair(Line_2(), false);
  }

  std::pair<Ray_2,bool> operator()(const typename IK::Ray_2& r) const
  {
    std::pair<Point_2,bool> sp = operator()(r.source());
    if(! sp.second){
      return std::make_pair(Ray_2(),false);
    }
    std::pair<Point_2,bool> tp = operator()(r.second_point());
    if(! tp.second){
      return std::make_pair(Ray_2(),false);
    }
    return std::make_pair(Ray_2(sp.first,tp.first), true);
  }

  std::pair<Triangle_2,bool> operator()(const typename IK::Triangle_2& t) const
  {
    std::pair<Point_2,bool> v0 = operator()(t.vertex(0));
    if(! v0.second){
      return std::make_pair(Triangle_2(),false);
    }
    std::pair<Point_2,bool> v1 = operator()(t.vertex(1));
    if(! v1.second){
      return std::make_pair(Triangle_2(),false);
    }
    std::pair<Point_2,bool> v2 = operator()(t.vertex(2));
    if(! v2.second){
      return std::make_pair(Triangle_2(),false);
    }
    return std::make_pair(Triangle_2(v0.first,v1.first, v2.first), true);
  }

  std::pair<Circle_2,bool> operator()(const typename IK::Circle_2& ci) const
  {
    std::pair<Point_2,bool> c = operator()(ci.center());
    std::pair<double, bool> sr = operator()(ci.squared_radius());
    if(c.second && sr.second){
      return std::make_pair(Circle_2(c.first, sr.first, ci.orientation()),true);
    }
    return std::make_pair(Circle_2(), false);
  }

  std::pair<Iso_rectangle_2,bool> operator()(const typename IK::Iso_rectangle_2& ir) const
  {
    std::pair<Point_2,bool> sp = operator()((ir.min)());
    if(! sp.second){
      return std::make_pair(Iso_rectangle_2(),false);
    }
    std::pair<Point_2,bool> tp = operator()((ir.max)());
    if(! tp.second){
      return std::make_pair(Iso_rectangle_2(),false);
    }
    return std::make_pair(Iso_rectangle_2(sp.first,tp.first), true);
  }


  std::pair<Line_3,bool> operator()(const typename IK::Line_3& li) const
  {
    std::pair<Point_3,bool> sp = operator()(li.point());
    if(! sp.second){
      return std::make_pair(Line_3(),false);
    }
    std::pair<Vector_3,bool> tp = operator()(li.to_vector());
    if(! tp.second){
      return std::make_pair(Line_3(),false);
    }
    return std::make_pair(Line_3(sp.first,tp.first), true);
  }

  std::pair<Plane_3,bool> operator()(const typename IK::Plane_3& pl) const
  {
    std::pair<double,bool> a = operator()(pl.a()),  b = operator()(pl.b()) , c = operator()(pl.c()) , d = operator()(pl.d());
    if(a.second && b.second && c.second && d.second){
      return std::make_pair(Plane_3(a.first, b.first, c.first, d.first),true);
    }
    return std::make_pair(Plane_3(), false);
  }

  std::pair<Triangle_3,bool> operator()(const typename IK::Triangle_3& t) const
  {
    std::pair<Point_3,bool> v0 = operator()(t.vertex(0));
    if(! v0.second){
      return std::make_pair(Triangle_3(),false);
    }
    std::pair<Point_3,bool> v1 = operator()(t.vertex(1));
    if(! v1.second){
      return std::make_pair(Triangle_3(),false);
    }
    std::pair<Point_3,bool> v2 = operator()(t.vertex(2));
    if(! v2.second){
      return std::make_pair(Triangle_3(),false);
    }
    return std::make_pair(Triangle_3(v0.first,v1.first, v2.first), true);
  }

  std::pair<Tetrahedron_3,bool> operator()(const typename IK::Tetrahedron_3& t) const
  {
    std::pair<Point_3,bool> v0 = operator()(t.vertex(0));
    if(! v0.second){
      return std::make_pair(Tetrahedron_3(),false);
    }
    std::pair<Point_3,bool> v1 = operator()(t.vertex(1));
    if(! v1.second){
      return std::make_pair(Tetrahedron_3(),false);
    }
    std::pair<Point_3,bool> v2 = operator()(t.vertex(2));
    if(! v2.second){
      return std::make_pair(Tetrahedron_3(),false);
    }
    std::pair<Point_3,bool> v3 = operator()(t.vertex(3));
    if(! v3.second){
      return std::make_pair(Tetrahedron_3(),false);
    }
    return std::make_pair(Tetrahedron_3(v0.first,v1.first, v2.first, v3.first), true);
  }

  std::pair<Ray_3,bool> operator()(const typename IK::Ray_3& r) const
  {
    std::pair<Point_3,bool> sp = operator()(r.source());
    if(! sp.second){
      return std::make_pair(Ray_3(),false);
    }
    std::pair<Point_3,bool> tp = operator()(r.second_point());
    if(! tp.second){
      return std::make_pair(Ray_3(),false);
    }
    return std::make_pair(Ray_3(sp.first,tp.first), true);
  }

  std::pair<Point_3,bool> operator()(const typename IK::Point_3& p) const
  {
    double x, y, z;
    internal::init_double(x, y, z, (IK_FT*)(0));
    if(fit_in_double(p.x(),x) && fit_in_double(p.y(),y) && fit_in_double(p.z(),z)){
      return std::make_pair(Point_3(x,y,z),true);
    }
    return std::make_pair(ORIGIN,false);
  }

  std::pair<Vector_3,bool> operator()(const typename IK::Vector_3& v) const
  {
    double x, y, z;
    internal::init_double(x, y, z, (IK_FT*)(0));
    if(fit_in_double(v.x(),x) && fit_in_double(v.y(),y) && fit_in_double(v.z(),z)){
      return std::make_pair(Vector_3(x,y,z),true);
    }
    return std::make_pair(Vector_3(),false);
  }

  std::pair<Direction_3,bool> operator()(const typename IK::Direction_3& d) const
  {
    double x, y, z;
    internal::init_double(x, y, z, (IK_FT*)(0));
    if(fit_in_double(d.dx(),x) && fit_in_double(d.dy(),y) && fit_in_double(d.dz(),z)){
      return std::make_pair(Direction_3(x,y,z),true);
    }
    return std::make_pair(Direction_3(),false);
  }

  std::pair<Segment_3,bool> operator()(const typename IK::Segment_3& s) const
  {
    std::pair<Point_3,bool> sp = operator()(s.source());
    if(! sp.second){
      return std::make_pair(Segment_3(),false);
    }
    std::pair<Point_3,bool> tp = operator()(s.target());
    if(! tp.second){
      return std::make_pair(Segment_3(),false);
    }
    return std::make_pair(Segment_3(sp.first,tp.first), true);
  }

  std::pair<Weighted_point_3,bool> operator()(const typename IK::Weighted_point_3& wp) const
  {
    std::pair<Point_3,bool> sp = operator()(wp.point());
    std::pair<double,bool> w = operator()(wp.weight());
    if(sp.second && w.second){
      return std::make_pair(Weighted_point_3(sp.first,w.first),true);
    }
    return std::make_pair(Weighted_point_3(),false);
  }

  std::pair<Sphere_3,bool> operator()(const typename IK::Sphere_3& s) const
  {
    std::pair<Point_3,bool> c = operator()(s.center());
    std::pair<double, bool> sr = operator()(s.squared_radius());
    if(c.second && sr.second){
      return std::make_pair(Sphere_3(c.first, sr.first, s.orientation()),true);
    }
    return std::make_pair(Sphere_3(), false);
  }

  std::pair<Circle_3,bool> operator()(const typename IK::Circle_3& ci) const
  {
    std::pair<Sphere_3, bool> sr = operator()(ci.diametral_sphere());
    std::pair<Plane_3,bool> c = operator()(ci.supporting_plane());
    if(c.second && sr.second){
      return std::make_pair(Circle_3(sr.first, c.first),true);
    }
    return std::make_pair(Circle_3(), false);
  }

  std::pair<Iso_cuboid_3,bool> operator()(const typename IK::Iso_cuboid_3& ic) const
  {
    std::pair<Point_3,bool> sp = operator()((ic.min)());
    if(! sp.second){
      return std::make_pair(Iso_cuboid_3(),false);
    }
    std::pair<Point_3,bool> tp = operator()((ic.max)());
    if(! tp.second){
      return std::make_pair(Iso_cuboid_3(),false);
    }
    return std::make_pair(Iso_cuboid_3(sp.first,tp.first), true);
  }


};

} // CGAL

#endif // CGAL_EPIC_CONVERTER_H
