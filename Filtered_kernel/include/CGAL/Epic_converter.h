// Copyright (c) 2017  GeometryFactory
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
// Author(s)     : Andreas Fabri, Laurent Rineau

#ifndef CGAL_EPIC_CONVERTER_H
#define CGAL_EPIC_CONVERTER_H


#include <CGAL/basic.h>
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
public:



  std::pair<double,bool> operator()(const typename IK::FT n) const
  {
    double d;
    if(fit_in_double(n,d)){
      return std::make_pair(d,true);
    }
    return std::make_pair(0,false);
  }

  
  std::pair<Point_2,bool> operator()(const typename IK::Point_2& p) const
  {
    double x, y;
    if(fit_in_double(p.x(),x) && fit_in_double(p.y(),y)){
      return std::make_pair(Point_2(x,y),true);
    }
    return std::make_pair(ORIGIN,false);
  }
  
  std::pair<Vector_2,bool> operator()(const typename IK::Vector_2& p) const
  {
    double x, y;
    if(fit_in_double(p.x(),x) && fit_in_double(p.y(),y)){
      return std::make_pair(Vector_2(x,y),true);
    }
    return std::make_pair(Vector_2(0,0),false);
  }
  
  std::pair<Direction_2,bool> operator()(const typename IK::Direction_2& p) const
  {
    double x, y;
    if(fit_in_double(p.dx(),x) && fit_in_double(p.dy(),y)){
      return std::make_pair(Direction_2(x,y),true);
    }
    return std::make_pair(Direction_2(0,0),false);
  }
  
  std::pair<Weighted_point_2,bool> operator()(const typename IK::Weighted_point_2& p) const
  {
    std::pair<Point_2,bool> sp = operator()(p.point());
    std::pair<double,bool> w = operator()(p.weight());
    if(sp.second && w.second){
      return std::make_pair(Weighted_point_2(sp.first,w.first),true);
    }
    return std::make_pair(Weighted_point_2(ORIGIN,0),false);
  }
  
  std::pair<Segment_2,bool> operator()(const typename IK::Segment_2& s) const
  {
    std::pair<Point_2,bool> sp = operator()(s.source());
    if(! sp.second){
      return std::make_pair(Segment_2(sp.first,sp.first),false);
    }
    std::pair<Point_2,bool> tp = operator()(s.target());
    if(! tp.second){
      return std::make_pair(Segment_2(tp.first,tp.first),false);
    }
    return std::make_pair(Segment_2(sp.first,tp.first), true);
  }

  std::pair<Line_2,bool> operator()(const typename IK::Line_2& s) const
  {
    std::pair<double,bool> a = operator()(s.a()),  b = operator()(s.b()) , c = operator()(s.c());
    if(a.second && b.second && c.second){
      return std::make_pair(Line_2(a.first, b.first, c.first),true);
    }
    return std::make_pair(Line_2(0,0,0), false);
  }
  
  std::pair<Ray_2,bool> operator()(const typename IK::Ray_2& s) const
  {
    std::pair<Point_2,bool> sp = operator()(s.source());
    if(! sp.second){
      return std::make_pair(Ray_2(ORIGIN,ORIGIN),false);
    }
    std::pair<Point_2,bool> tp = operator()(s.second_point());
    if(! tp.second){
      return std::make_pair(Ray_2(ORIGIN,ORIGIN),false);
    }
    return std::make_pair(Ray_2(sp.first,tp.first), true);
  }
  
  std::pair<Triangle_2,bool> operator()(const typename IK::Triangle_2& s) const
  {
    std::pair<Point_2,bool> v0 = operator()(s.vertex(0));
    if(! v0.second){
      return std::make_pair(Triangle_2(ORIGIN, ORIGIN, ORIGIN),false);
    }
    std::pair<Point_2,bool> v1 = operator()(s.vertex(1));
    if(! v1.second){
      return std::make_pair(Triangle_2(ORIGIN,ORIGIN, ORIGIN),false);
    }
    std::pair<Point_2,bool> v2 = operator()(s.vertex(2));
    if(! v2.second){
      return std::make_pair(Triangle_2(ORIGIN,ORIGIN, ORIGIN),false);
    }
    return std::make_pair(Triangle_2(v0.first,v1.first, v2.first), true);
  }
  
  std::pair<Circle_2,bool> operator()(const typename IK::Circle_2& s) const
  {
    std::pair<Point_2,bool> c = operator()(s.center());
    std::pair<double, bool> sr = operator()(s.squared_radius());
    if(c.second && sr.second){
      return std::make_pair(Circle_2(c.first, sr.first, s.orientation()),true);
    }
    return std::make_pair(Circle_2(ORIGIN, 0, s.orientation()), false);
  }

  std::pair<Iso_rectangle_2,bool> operator()(const typename IK::Iso_rectangle_2& s) const
  {
    std::pair<Point_2,bool> sp = operator()((s.min)());
    if(! sp.second){
      return std::make_pair(Iso_rectangle_2(ORIGIN,ORIGIN),false);
    }
    std::pair<Point_2,bool> tp = operator()((s.max)());
    if(! tp.second){
      return std::make_pair(Iso_rectangle_2(ORIGIN,ORIGIN),false);
    }
    return std::make_pair(Iso_rectangle_2(sp.first,tp.first), true);
  }

  
  std::pair<Line_3,bool> operator()(const typename IK::Line_3& s) const
  {
    std::pair<Point_3,bool> sp = operator()(s.point());
    if(! sp.second){
      return std::make_pair(Line_3(),false);
    }
    std::pair<Vector_3,bool> tp = operator()(s.to_vector());
    if(! tp.second){
      return std::make_pair(Line_3(ORIGIN,NULL_VECTOR),false);
    }
    return std::make_pair(Line_3(sp.first,tp.first), true);
  }

  std::pair<Plane_3,bool> operator()(const typename IK::Plane_3& s) const
  {
    std::pair<double,bool> a = operator()(s.a()),  b = operator()(s.b()) , c = operator()(s.c()) , d = operator()(s.d());
    if(a.second && b.second && c.second && d.second){
      return std::make_pair(Plane_3(a.first, b.first, c.first, d.first),true);
    }
    return std::make_pair(Plane_3(0,0,0,0), false);
  }

  std::pair<Triangle_3,bool> operator()(const typename IK::Triangle_3& s) const
  {
    std::pair<Point_3,bool> v0 = operator()(s.vertex(0));
    if(! v0.second){
      return std::make_pair(Triangle_3(ORIGIN, ORIGIN, ORIGIN),false);
    }
    std::pair<Point_3,bool> v1 = operator()(s.vertex(1));
    if(! v1.second){
      return std::make_pair(Triangle_3(ORIGIN,ORIGIN, ORIGIN),false);
    }
    std::pair<Point_3,bool> v2 = operator()(s.vertex(2));
    if(! v2.second){
      return std::make_pair(Triangle_3(ORIGIN,ORIGIN, ORIGIN),false);
    }
    return std::make_pair(Triangle_3(v0.first,v1.first, v2.first), true);
  }
  
  std::pair<Tetrahedron_3,bool> operator()(const typename IK::Tetrahedron_3& s) const
  {
    std::pair<Point_3,bool> v0 = operator()(s.vertex(0));
    if(! v0.second){
      return std::make_pair(Tetrahedron_3(ORIGIN, ORIGIN, ORIGIN, ORIGIN),false);
    }
    std::pair<Point_3,bool> v1 = operator()(s.vertex(1));
    if(! v1.second){
      return std::make_pair(Tetrahedron_3(ORIGIN,ORIGIN, ORIGIN, ORIGIN),false);
    }
    std::pair<Point_3,bool> v2 = operator()(s.vertex(2));
    if(! v2.second){
      return std::make_pair(Tetrahedron_3(ORIGIN,ORIGIN, ORIGIN, ORIGIN),false);
    }
    std::pair<Point_3,bool> v3 = operator()(s.vertex(3));
    if(! v3.second){
      return std::make_pair(Tetrahedron_3(ORIGIN,ORIGIN, ORIGIN, ORIGIN),false);
    }
    return std::make_pair(Tetrahedron_3(v0.first,v1.first, v2.first, v3.first), true);
  }
  
  std::pair<Ray_3,bool> operator()(const typename IK::Ray_3& s) const
  {
    std::pair<Point_3,bool> sp = operator()(s.source());
    if(! sp.second){
      return std::make_pair(Ray_3(ORIGIN,ORIGIN),false);
    }
    std::pair<Point_3,bool> tp = operator()(s.second_point());
    if(! tp.second){
      return std::make_pair(Ray_3(ORIGIN,ORIGIN),false);
    }
    return std::make_pair(Ray_3(sp.first,tp.first), true);
  }
  
  std::pair<Point_3,bool> operator()(const typename IK::Point_3& p) const
  {
    double x, y, z;
    if(fit_in_double(p.x(),x) && fit_in_double(p.y(),y) && fit_in_double(p.z(),z)){
      return std::make_pair(Point_3(x,y,z),true);
    }
    return std::make_pair(ORIGIN,false);
  }
  
  std::pair<Vector_3,bool> operator()(const typename IK::Vector_3& p) const
  {
    double x, y, z;
    if(fit_in_double(p.x(),x) && fit_in_double(p.y(),y) && fit_in_double(p.z(),z)){
      return std::make_pair(Vector_3(x,y,z),true);
    }
    return std::make_pair(Vector_3(0,0,0),false);
  }
  
  std::pair<Direction_3,bool> operator()(const typename IK::Direction_3& p) const
  {
    double x, y, z;
    if(fit_in_double(p.dx(),x) && fit_in_double(p.dy(),y) && fit_in_double(p.dz(),z)){
      return std::make_pair(Direction_3(x,y,z),true);
    }
    return std::make_pair(Direction_3(0,0,0),false);
  }
  
  std::pair<Segment_3,bool> operator()(const typename IK::Segment_3& s) const
  {
    std::pair<Point_3,bool> sp = operator()(s.source());
    if(! sp.second){
      return std::make_pair(Segment_3(sp.first,sp.first),false);
    }
    std::pair<Point_3,bool> tp = operator()(s.target());
    if(! tp.second){
      return std::make_pair(Segment_3(tp.first,tp.first),false);
    }
    return std::make_pair(Segment_3(sp.first,tp.first), true);
  }

  std::pair<Weighted_point_3,bool> operator()(const typename IK::Weighted_point_3& p) const
  {
    std::pair<Point_3,bool> sp = operator()(p.point());
    std::pair<double,bool> w = operator()(p.weight());
    if(sp.second && w.second){
      return std::make_pair(Weighted_point_3(sp.first,w.first),true);
    }
    return std::make_pair(Weighted_point_3(ORIGIN,0),false);
  }
  
  std::pair<Sphere_3,bool> operator()(const typename IK::Sphere_3& s) const
  {
    std::pair<Point_3,bool> c = operator()(s.center());
    std::pair<double, bool> sr = operator()(s.squared_radius());
    if(c.second && sr.second){
      return std::make_pair(Sphere_3(c.first, sr.first, s.orientation()),true);
    }
    return std::make_pair(Sphere_3(ORIGIN, 0, s.orientation()), false);
  }

    std::pair<Circle_3,bool> operator()(const typename IK::Circle_3& s) const
  {
    std::pair<Sphere_3, bool> sr = operator()(s.diametral_sphere());
    std::pair<Plane_3,bool> c = operator()(s.supporting_plane());
    if(c.second && sr.second){
      return std::make_pair(Circle_3(sr.first, c.first),true);
    }
    return std::make_pair(Circle_3(), false);
  }

  std::pair<Iso_cuboid_3,bool> operator()(const typename IK::Iso_cuboid_3& s) const
  {
    std::pair<Point_3,bool> sp = operator()((s.min)());
    if(! sp.second){
      return std::make_pair(Iso_cuboid_3(ORIGIN,ORIGIN),false);
    }
    std::pair<Point_3,bool> tp = operator()((s.max)());
    if(! tp.second){
      return std::make_pair(Iso_cuboid_3(ORIGIN,ORIGIN),false);
    }
    return std::make_pair(Iso_cuboid_3(sp.first,tp.first), true);
  }


};
  
} // CGAL

#endif // CGAL_EPIC_CONVERTER_H
