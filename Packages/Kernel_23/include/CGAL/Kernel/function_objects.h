// Copyright (c) 1999,2002  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Stefan Schirra, Sylvain Pion

#ifndef CGAL_KERNEL_FUNCTION_OBJECTS_H
#define CGAL_KERNEL_FUNCTION_OBJECTS_H

#include <CGAL/functional_base.h>
#include <CGAL/Origin.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Kernel/Cartesian_coordinate_iterator_2.h>
#include <CGAL/Kernel/Cartesian_coordinate_iterator_3.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/intersection_2.h>
#include <CGAL/intersection_3.h>

CGAL_BEGIN_NAMESPACE

namespace CommonKernelFunctors {

  template <typename K>
  class Are_ordered_along_line_2
  {
    typedef typename K::Point_2     Point_2;
    typedef typename K::Collinear_2 Collinear_2;
    typedef typename K::Collinear_are_ordered_along_line_2 
    Collinear_are_ordered_along_line_2;
  
    Collinear_2 c;
    Collinear_are_ordered_along_line_2 cao;
  public:
    typedef bool             result_type;
    typedef Arity_tag< 3 >   Arity;

    Are_ordered_along_line_2() {}
    Are_ordered_along_line_2(const Collinear_2& c_, 
			     const Collinear_are_ordered_along_line_2& cao_)
      : c(c_), cao(cao_)
    {}

    bool
    operator()(const Point_2& p, const Point_2& q, const Point_2& r) const
    { return c(p, q, r) && cao(p, q, r); }
  };

  template <typename K>
  class Are_ordered_along_line_3
  {
    typedef typename K::Point_3     Point_3;
    typedef typename K::Collinear_3 Collinear_3;
    typedef typename K::Collinear_are_ordered_along_line_3 
    Collinear_are_ordered_along_line_3;
  
    Collinear_3 c;
    Collinear_are_ordered_along_line_3 cao;
  public:
    typedef bool             result_type;
    typedef Arity_tag< 3 >   Arity;

    Are_ordered_along_line_3() {}
    Are_ordered_along_line_3(const Collinear_3& c_, 
			     const Collinear_are_ordered_along_line_3& cao_)
      : c(c_), cao(cao_)
    {}

    bool
    operator()(const Point_3& p, const Point_3& q, const Point_3& r) const
    { return c(p, q, r) && cao(p, q, r); }
  };

  template <typename K>
  class Are_strictly_ordered_along_line_2
  {
    typedef typename K::Point_2 Point_2;
    typedef typename K::Collinear_2 Collinear_2;
    typedef typename K::Collinear_are_strictly_ordered_along_line_2 
    Collinear_are_strictly_ordered_along_line_2;
  
    Collinear_2 c;
    Collinear_are_strictly_ordered_along_line_2 cao;
  public:
    typedef bool             result_type;
    typedef Arity_tag< 3 >   Arity;

    Are_strictly_ordered_along_line_2() {}
    Are_strictly_ordered_along_line_2(
				      const Collinear_2& c_, 
				      const Collinear_are_strictly_ordered_along_line_2& cao_)
      : c(c_), cao(cao_)
    {}

    bool
    operator()(const Point_2& p, const Point_2& q, const Point_2& r) const
    { return c(p, q, r) && cao(p, q, r); }
  };

  template <typename K>
  class Are_strictly_ordered_along_line_3
  {
    typedef typename K::Point_3     Point_3;
    typedef typename K::Collinear_3 Collinear_3;
    typedef typename K::Collinear_are_strictly_ordered_along_line_3 
    Collinear_are_strictly_ordered_along_line_3;
  
    Collinear_3 c;
    Collinear_are_strictly_ordered_along_line_3 cao;
  public:
    typedef bool             result_type;
    typedef Arity_tag< 3 >   Arity;

    Are_strictly_ordered_along_line_3() {}
    Are_strictly_ordered_along_line_3(
				      const Collinear_3& c_, 
				      const Collinear_are_strictly_ordered_along_line_3& cao_)
      : c(c_), cao(cao_)
    {}

    bool
    operator()(const Point_3& p, const Point_3& q, const Point_3& r) const
    { return c(p, q, r) && cao(p, q, r); }
  };

  template <typename K>
  class Assign_2
  {
    typedef typename K::Object_2 Object_2;
  public:
    typedef bool             result_type;
    typedef Arity_tag< 2 >   Arity;

    template <class T>
    bool
    operator()(T& t, const Object_2& o) const
    { return assign(t, o); }
  };

  template <typename K>
  class Assign_3
  {
    typedef typename K::Object_3 Object_3;
  public:
    typedef bool             result_type;
    typedef Arity_tag< 2 >   Arity;

    template <class T>
    bool
    operator()(T& t, const Object_3& o) const
    { return assign(t, o); }
  };

  template <typename K>
  class Bounded_side_2
  {
    typedef typename K::Point_2         Point_2;
    typedef typename K::Circle_2        Circle_2;
    typedef typename K::Triangle_2      Triangle_2;
    typedef typename K::Iso_rectangle_2 Iso_rectangle_2;
  public:
    typedef Bounded_side     result_type;
    typedef Arity_tag< 2 >   Arity;

    Bounded_side
    operator()( const Circle_2& c, const Point_2& p) const
    { return c.bounded_side(p); }

    Bounded_side
    operator()( const Triangle_2& t, const Point_2& p) const
    { return t.bounded_side(p); }

    Bounded_side
    operator()( const Iso_rectangle_2& r, const Point_2& p) const
    { return r.bounded_side(p); }
  };

  template <typename K>
  class Compute_area_3
  {
    typedef typename K::FT                FT;
    typedef typename K::Point_3           Point_3;
    typedef typename K::Triangle_3        Triangle_3;
  public:
    typedef FT               result_type;
    typedef Arity_tag< 1 >   Arity;

    FT
    operator()( const Triangle_3& t ) const
    {
	return CGAL_NTS sqrt(K().compute_squared_area_3_object()(t));
    }

    FT
    operator()( const Point_3& p, const Point_3& q, const Point_3& r ) const
    {
	return CGAL_NTS sqrt(K().compute_squared_area_3_object()(p, q, r));
    }
  };

  template <typename K>
  class Compute_squared_distance_2
  {
    typedef typename K::FT   FT;
  public:
    typedef FT               result_type;
    typedef Arity_tag< 2 >   Arity;

    // There are 25 combinaisons, we use a template.
    template <class T1, class T2>
    FT
    operator()( const T1& t1, const T2& t2) const
    { return CGALi::squared_distance(t1, t2, K()); }
  };

  template <typename K>
  class Compute_squared_distance_3
  {
    typedef typename K::FT   FT;
  public:
    typedef FT               result_type;
    typedef Arity_tag< 2 >   Arity;

    // There are 25 combinaisons, we use a template.
    template <class T1, class T2>
    FT
    operator()( const T1& t1, const T2& t2) const
    { return CGALi::squared_distance(t1, t2, K()); }
  };

  template <typename K>
  class Compute_squared_length_2
  {
    typedef typename K::FT          FT;
    typedef typename K::Segment_2   Segment_2;
    typedef typename K::Vector_2    Vector_2;
  public:
    typedef FT               result_type;
    typedef Arity_tag< 1 >   Arity;

    FT
    operator()( const Vector_2& v) const
    { return v.squared_length(); }

    FT
    operator()( const Segment_2& s) const
    { return s.squared_length(); }
  };

  template <typename K>
  class Compute_squared_length_3
  {
    typedef typename K::FT          FT;
    typedef typename K::Segment_3   Segment_3;
    typedef typename K::Vector_3    Vector_3;
  public:
    typedef FT               result_type;
    typedef Arity_tag< 1 >   Arity;

    FT
    operator()( const Vector_3& v) const
    { return v.squared_length(); }

    FT
    operator()( const Segment_3& s) const
    { return s.squared_length(); }
  };

  template <typename K>
  class Construct_center_2
  {
    typedef typename K::Point_2   Point_2;
    typedef typename K::Circle_2  Circle_2;
  public:
    typedef Point_2          result_type;
    typedef Arity_tag< 1 >   Arity;

    Point_2
    operator()(const Circle_2& c) const
    { return c.center(); }
  };

  template <typename K>
  class Construct_center_3
  {
    typedef typename K::Point_3   Point_3;
    typedef typename K::Sphere_3  Sphere_3;
  public:
    typedef Point_3          result_type;
    typedef Arity_tag< 1 >   Arity;

    Point_3
    operator()(const Sphere_3& s) const
    { return s.center(); }
  };

  template <typename K>
  class Construct_circle_2
  {
    typedef typename K::FT          FT;
    typedef typename K::Point_2     Point_2;
    typedef typename K::Circle_2    Circle_2;
  public:
    typedef Circle_2         result_type;
    typedef Arity_tag< 3 >   Arity;

    Circle_2
    operator()() const
    { return Circle_2(); }

    Circle_2
    operator()( const Point_2& center, const FT& squared_radius,
	        Orientation orientation = COUNTERCLOCKWISE) const
    { return Circle_2(center, squared_radius, orientation); }

    Circle_2
    operator()( const Point_2& p, const Point_2& q, const Point_2& r) const
    { return Circle_2(p, q, r); }

    Circle_2
    operator()( const Point_2& p, const Point_2& q,
	        Orientation orientation = COUNTERCLOCKWISE) const
    { return Circle_2(p, q, orientation); }

    Circle_2
    operator()( const Point_2& center,
	        Orientation orientation = COUNTERCLOCKWISE) const
    { return Circle_2(center, orientation); }
  };

  template <typename K>
  class Construct_direction_2
  {
    typedef typename K::Direction_2     Direction_2;
    typedef typename K::Vector_2        Vector_2;
    typedef typename K::Line_2          Line_2;
    typedef typename K::Ray_2           Ray_2;
    typedef typename K::Segment_2       Segment_2;
    typedef typename K::RT              RT;
  public:
    typedef Direction_2       result_type;
    typedef Arity_tag< 1 >    Arity;

    Direction_2
    operator()() const
    { return Direction_2(); }

#ifndef CGAL_NO_DEPRECATED_CODE
    Direction_2
    operator()(const RT& x, const RT& y) const
    { return Direction_2(x, y); }
#endif // CGAL_NO_DEPRECATED_CODE

    Direction_2
    operator()(const Vector_2& v) const
    { return Direction_2(v); }

    Direction_2
    operator()(const Line_2& l) const
    { return Direction_2(l); }

    Direction_2
    operator()(const Ray_2& r) const
    { return Direction_2(r); }

    Direction_2
    operator()(const Segment_2& s) const
    { return Direction_2(s); }
  };

  template <typename K>
  class Construct_direction_3
  {
    typedef typename K::Direction_3     Direction_3;
    typedef typename K::Vector_3        Vector_3;
    typedef typename K::Line_3          Line_3;
    typedef typename K::Ray_3           Ray_3;
    typedef typename K::Segment_3       Segment_3;
    typedef typename K::RT              RT;
  public:
    typedef Direction_3       result_type;
    typedef Arity_tag< 1 >    Arity;

    Direction_3
    operator()() const
    { return Direction_3(); }

#ifndef CGAL_NO_DEPRECATED_CODE
    Direction_3
    operator()(const RT& x, const RT& y, const RT& z) const
    { return Direction_3(x, y, z); }
#endif // CGAL_NO_DEPRECATED_CODE

    Direction_3
    operator()(const Vector_3& v) const
    { return Direction_3(v); }

    Direction_3
    operator()(const Line_3& l) const
    { return Direction_3(l); }

    Direction_3
    operator()(const Ray_3& r) const
    { return Direction_3(r); }

    Direction_3
    operator()(const Segment_3& s) const
    { return Direction_3(s); }
  };

  template <typename K>
  class Construct_iso_cuboid_3
  {
    typedef typename K::Point_3       Point_3;
    typedef typename K::Iso_cuboid_3  Iso_cuboid_3;
  public:
    typedef Iso_cuboid_3      result_type;
    typedef Arity_tag< 2 >    Arity;

    Iso_cuboid_3
    operator()() const
    { return Iso_cuboid_3(); }

    Iso_cuboid_3
    operator()(const Point_3& p, const Point_3& q) const
    { return Iso_cuboid_3(p, q); }

    Iso_cuboid_3
    operator()(const Point_3 &left,   const Point_3 &right,
               const Point_3 &bottom, const Point_3 &top,
               const Point_3 &far_,   const Point_3 &close) const
    { return Iso_cuboid_3(left, right, bottom, top, far_, close); }
  };

  template <typename K>
  class Construct_iso_rectangle_2
  {
    typedef typename K::Point_2          Point_2;
    typedef typename K::Iso_rectangle_2  Iso_rectangle_2;
  public:
    typedef Iso_rectangle_2   result_type;
    typedef Arity_tag< 2 >    Arity;

    Iso_rectangle_2
    operator()() const
    { return Iso_rectangle_2(); }

    Iso_rectangle_2
    operator()(const Point_2& p, const Point_2& q) const
    { return Iso_rectangle_2(p, q); }

    Iso_rectangle_2
    operator()(const Point_2 &left,   const Point_2 &right,
               const Point_2 &bottom, const Point_2 &top) const
    { return Iso_rectangle_2(left, right, bottom, top); }
  };

  template <typename K>
  class Construct_max_vertex_2
  {
    typedef typename K::Point_2          Point_2;
    typedef typename K::Iso_rectangle_2  Iso_rectangle_2;
  public:
    typedef Point_2   result_type;
    typedef Arity_tag< 1 >    Arity;

    Point_2
    operator()(const Iso_rectangle_2& r) const
    { return r.max(); }
  };

  template <typename K>
  class Construct_min_vertex_2
  {
    typedef typename K::Point_2          Point_2;
    typedef typename K::Iso_rectangle_2  Iso_rectangle_2;
  public:
    typedef Point_2   result_type;
    typedef Arity_tag< 1 >    Arity;

    Point_2
    operator()(const Iso_rectangle_2& r) const
    { return r.min(); }
  };


  template <typename K>
  class Construct_max_vertex_3
  {
    typedef typename K::Point_3          Point_3;
    typedef typename K::Iso_cuboid_3  Iso_cuboid_3;
  public:
    typedef Point_3   result_type;
    typedef Arity_tag< 1 >    Arity;

    Point_3
    operator()(const Iso_cuboid_3& r) const
    { return r.max(); }
  };

  template <typename K>
  class Construct_min_vertex_3
  {
    typedef typename K::Point_3          Point_3;
    typedef typename K::Iso_cuboid_3  Iso_cuboid_3;
  public:
    typedef Point_3   result_type;
    typedef Arity_tag< 1 >    Arity;

    Point_3
    operator()(const Iso_cuboid_3& r) const
    { return r.min(); }
  };


  template <typename K>
  class Construct_object_2
  {
    typedef typename K::Object_2   Object_2;
  public:
    typedef Object_2         result_type;
    typedef Arity_tag< 1 >   Arity;

    template <class Cls>
    Object_2
    operator()( const Cls& c) const
    { return make_object(c); }
  };

  template <typename K>
  class Construct_object_3
  {
    typedef typename K::Object_3   Object_3;
  public:
    typedef Object_3         result_type;
    typedef Arity_tag< 1 >   Arity;

    template <class Cls>
    Object_3
    operator()( const Cls& c) const
    { return make_object(c); }
  };

  template <typename K>
  class Construct_opposite_circle_2
  {
    typedef typename K::Circle_2   Circle_2;
  public:
    typedef Circle_2         result_type;
    typedef Arity_tag< 1 >   Arity;

    Circle_2
    operator()( const Circle_2& c) const
    { return c.opposite(); }
  };

  template <typename K>
  class Construct_opposite_direction_2
  {
    typedef typename K::Direction_2  Direction_2;
  public:
    typedef Direction_2      result_type;
    typedef Arity_tag< 1 >   Arity;

    Direction_2
    operator()( const Direction_2& d) const
    { return -d; }
  };

  template <typename K>
  class Construct_opposite_direction_3
  {
    typedef typename K::Direction_3  Direction_3;
  public:
    typedef Direction_3      result_type;
    typedef Arity_tag< 1 >   Arity;

    Direction_3
    operator()( const Direction_3& d) const
    { return -d; }
  };

  template <typename K>
  class Construct_opposite_line_2
  {
    typedef typename K::Line_2   Line_2;
  public:
    typedef Line_2           result_type;
    typedef Arity_tag< 1 >   Arity;

    Line_2
    operator()( const Line_2& l) const
    { return l.opposite(); }
  };

  template <typename K>
  class Construct_opposite_line_3
  {
    typedef typename K::Line_3   Line_3;
  public:
    typedef Line_3           result_type;
    typedef Arity_tag< 1 >   Arity;

    Line_3
    operator()( const Line_3& l) const
    { return l.opposite(); }
  };

  template <typename K>
  class Construct_opposite_plane_3
  {
    typedef typename K::Plane_3   Plane_3;
  public:
    typedef Plane_3          result_type;
    typedef Arity_tag< 1 >   Arity;

    Plane_3
    operator()( const Plane_3& p) const
    { return p.opposite(); }
  };

  template <typename K>
  class Construct_opposite_ray_2
  {
    typedef typename K::Ray_2   Ray_2;
  public:
    typedef Ray_2            result_type;
    typedef Arity_tag< 1 >   Arity;

    Ray_2
    operator()( const Ray_2& r) const
    { return r.opposite(); }
  };

  template <typename K>
  class Construct_opposite_ray_3
  {
    typedef typename K::Ray_3   Ray_3;
  public:
    typedef Ray_3            result_type;
    typedef Arity_tag< 1 >   Arity;

    Ray_3
    operator()( const Ray_3& r) const
    { return r.opposite(); }
  };

  template <typename K>
  class Construct_opposite_segment_2
  {
    typedef typename K::Segment_2  Segment_2;
  public:
    typedef Segment_2        result_type;
    typedef Arity_tag< 1 >   Arity;

    Segment_2
    operator()( const Segment_2& s) const
    { return s.opposite(); }
  };

  template <typename K>
  class Construct_opposite_segment_3
  {
    typedef typename K::Segment_3  Segment_3;
  public:
    typedef Segment_3        result_type;
    typedef Arity_tag< 1 >   Arity;

    Segment_3
    operator()( const Segment_3& s) const
    { return s.opposite(); }
  };

  template <typename K>
  class Construct_opposite_sphere_3
  {
    typedef typename K::Sphere_3   Sphere_3;
  public:
    typedef Sphere_3         result_type;
    typedef Arity_tag< 1 >   Arity;

    Sphere_3
    operator()( const Sphere_3& s) const
    { return s.opposite(); }
  };

  template <typename K>
  class Construct_opposite_triangle_2
  {
    typedef typename K::Triangle_2  Triangle_2;
  public:
    typedef Triangle_2       result_type;
    typedef Arity_tag< 1 >   Arity;

    Triangle_2
    operator()( const Triangle_2& t) const
    { return t.opposite(); }
  };

  template <typename K>
  class Construct_perpendicular_direction_2
  {
    typedef typename K::Direction_2   Direction_2;
  public:
    typedef Direction_2      result_type;
    typedef Arity_tag< 2 >   Arity;

    Direction_2
    operator()( const Direction_2& d, Orientation o) const
    { return d.perpendicular(o); }
  };

  template <typename K>
  class Construct_perpendicular_line_2
  {
    typedef typename K::Line_2    Line_2;
    typedef typename K::Point_2   Point_2;
  public:
    typedef Line_2           result_type;
    typedef Arity_tag< 2 >   Arity;

    Line_2
    operator()( const Line_2& l, const Point_2& p) const
    { return l.perpendicular(p); }
  };

  template <typename K>
  class Construct_perpendicular_line_3
  {
    typedef typename K::Line_3    Line_3;
    typedef typename K::Point_3   Point_3;
    typedef typename K::Plane_3   Plane_3;
  public:
    typedef Line_3           result_type;
    typedef Arity_tag< 2 >   Arity;

    Line_3
    operator()( const Plane_3& pl, const Point_3& p) const
    { return pl.perpendicular_line(p); }
  };

  template <typename K>
  class Construct_perpendicular_plane_3
  {
    typedef typename K::Line_3    Line_3;
    typedef typename K::Point_3   Point_3;
    typedef typename K::Plane_3   Plane_3;
  public:
    typedef Plane_3          result_type;
    typedef Arity_tag< 2 >   Arity;

    Plane_3
    operator()( const Line_3& l, const Point_3& p) const
    { return l.perpendicular_plane(p); }
  };

  template <typename K>
  class Construct_perpendicular_vector_2
  {
    typedef typename K::Vector_2   Vector_2;
  public:
    typedef Vector_2         result_type;
    typedef Arity_tag< 2 >   Arity;

    Vector_2
    operator()( const Vector_2& v, Orientation o) const
    { return v.perpendicular(o); }
  };

  template <typename K>
  class Construct_plane_3
  {
    typedef typename K::RT           RT;
    typedef typename K::Point_3      Point_3;
    typedef typename K::Vector_3     Vector_3;
    typedef typename K::Direction_3  Direction_3;
    typedef typename K::Line_3       Line_3;
    typedef typename K::Ray_3        Ray_3;
    typedef typename K::Segment_3    Segment_3;
    typedef typename K::Plane_3      Plane_3;
  public:
    typedef Plane_3          result_type;
    typedef Arity_tag< 2 >   Arity;

    Plane_3
    operator()() const
    { return Plane_3(); }

    Plane_3
    operator()(const RT& a, const RT& b, const RT& c, const RT& d) const
    { return Plane_3(a, b, c, d); }

    Plane_3
    operator()(const Point_3& p, const Point_3& q, const Point_3& r) const
    { return Plane_3(p, q, r); }

    Plane_3
    operator()(const Point_3& p, const Direction_3& d) const
    { return Plane_3(p, d); }

    Plane_3
    operator()(const Point_3& p, const Vector_3& v) const
    { return Plane_3(p, v); }

    Plane_3
    operator()(const Line_3& l, const Point_3& p) const
    { return Plane_3(l, p); }

    Plane_3
    operator()(const Ray_3& r, const Point_3& p) const
    { return Plane_3(r, p); }

    Plane_3
    operator()(const Segment_3& s, const Point_3& p) const
    { return Plane_3(s, p); }
  };

  template <typename K>
  class Construct_point_on_2
  {
    typedef typename K::Point_2    Point_2;
    typedef typename K::Segment_2  Segment_2;
    typedef typename K::Line_2     Line_2;
    typedef typename K::Ray_2      Ray_2;
  public:
    typedef Point_2          result_type;
    typedef Arity_tag< 2 >   Arity;

    Point_2
    operator()( const Line_2& l, int i) const
    { return l.point(i); }

    Point_2
    operator()( const Segment_2& s, int i) const
    { return s.point(i); }

    Point_2
    operator()( const Ray_2& r, int i) const
    { return r.point(i); }
  };

  template <typename K>
  class Construct_point_on_3
  {
    typedef typename K::Point_3    Point_3;
    typedef typename K::Segment_3  Segment_3;
    typedef typename K::Line_3     Line_3;
    typedef typename K::Ray_3      Ray_3;
    typedef typename K::Plane_3    Plane_3;
  public:
    typedef Point_3          result_type;
    typedef Arity_tag< 2 >   Arity;

    Point_3
    operator()( const Line_3& l, int i) const
    { return l.point(i); }

    Point_3
    operator()( const Segment_3& s, int i) const
    { return s.point(i); }

    Point_3
    operator()( const Ray_3& r, int i) const
    { return r.point(i); }

    Point_3
    operator()( const Plane_3& p) const
    { return p.point(); }
  };

  template <typename K>
  class Construct_point_2
  {
    typedef typename K::RT         RT;
    typedef typename K::Point_2    Point_2;
  public:
    typedef Point_2          result_type;
    typedef Arity_tag< 1 >   Arity;

    Point_2
    operator()() const
    { return Point_2(); }

    Point_2
    operator()(Origin o) const
    { return Point_2(o); }

    // Reactivated, as some functors in Cartesian/function_objects.h
    // need it for constructions
    //#ifndef CGAL_NO_DEPRECATED_CODE
    Point_2
    operator()(const RT& x, const RT& y) const
    { return Point_2(x, y); }

    Point_2
    operator()(const RT& x, const RT& y, const RT& w) const
    { return Point_2(x, y, w); }
    //#endif // CGAL_NO_DEPRECATED_CODE
  };

  template <typename K>
  class Construct_point_3
  {
    typedef typename K::RT         RT;
    typedef typename K::Point_3    Point_3;
  public:
    typedef Point_3          result_type;
    typedef Arity_tag< 1 >   Arity;

    Point_3
    operator()() const
    { return Point_3(); }

    Point_3
    operator()(Origin o) const
    { return Point_3(o); }


    // Reactivated, as some functors in Cartesian/function_objects.h
    // need it for constructions
    //#ifndef CGAL_NO_DEPRECATED_CODE
    Point_3
    operator()(const RT& x, const RT& y, const RT& z) const
    { return Point_3(x, y, z); }

    Point_3
    operator()(const RT& x, const RT& y, const RT& z, const RT& w) const
    { return Point_3(x, y, z, w); }
    //#endif // CGAL_NO_DEPRECATED_CODE
  };

  template <typename K>
  class Construct_projected_point_2
  {
    typedef typename K::Point_2    Point_2;
    typedef typename K::Line_2     Line_2;
  public:
    typedef Point_2          result_type;
    typedef Arity_tag< 2 >   Arity;

    Point_2
    operator()( const Line_2& l, const Point_2& p ) const
    { return l.projection(p); }
  };

  template <typename K>
  class Construct_projected_xy_point_2
  {
    typedef typename K::Point_2    Point_2;
    typedef typename K::Point_3    Point_3;
    typedef typename K::Plane_3    Plane_3;
  public:
    typedef Point_2          result_type;
    typedef Arity_tag< 2 >   Arity;

    Point_2
    operator()( const Plane_3& h, const Point_3& p) const
    {  return h.to_2d(p); }
  };

  template <typename K>
  class Construct_ray_2
  {
    typedef typename K::Point_2      Point_2;
    typedef typename K::Vector_2     Vector_2;
    typedef typename K::Direction_2  Direction_2;
    typedef typename K::Line_2       Line_2;
    typedef typename K::Ray_2        Ray_2;
  public:
    typedef Ray_2            result_type;
    typedef Arity_tag< 2 >   Arity;

    Ray_2
    operator()() const
    {  return Ray_2(); }

    Ray_2
    operator()(const Point_2& p, const Point_2& q) const
    {  return Ray_2(p, q); }

    Ray_2
    operator()(const Point_2& p, const Vector_2& v) const
    {  return Ray_2(p, v); }

    Ray_2
    operator()(const Point_2& p, const Direction_2& d) const
    {  return Ray_2(p, d); }

    Ray_2
    operator()(const Point_2& p, const Line_2& l) const
    {  return Ray_2(p, l); }
  };

  template <typename K>
  class Construct_ray_3
  {
    typedef typename K::Point_3      Point_3;
    typedef typename K::Vector_3     Vector_3;
    typedef typename K::Direction_3  Direction_3;
    typedef typename K::Line_3       Line_3;
    typedef typename K::Ray_3        Ray_3;
  public:
    typedef Ray_3            result_type;
    typedef Arity_tag< 2 >   Arity;

    Ray_3
    operator()() const
    {  return Ray_3(); }

    Ray_3
    operator()(const Point_3& p, const Point_3& q) const
    {  return Ray_3(p, q); }

    Ray_3
    operator()(const Point_3& p, const Vector_3& v) const
    {  return Ray_3(p, v); }

    Ray_3
    operator()(const Point_3& p, const Direction_3& d) const
    {  return Ray_3(p, d); }

    Ray_3
    operator()(const Point_3& p, const Line_3& l) const
    {  return Ray_3(p, l); }
  };

  template <typename K>
  class Construct_segment_2
  {
    typedef typename K::Segment_2  Segment_2;
    typedef typename K::Point_2    Point_2;
  public:
    typedef Segment_2        result_type;
    typedef Arity_tag< 2 >   Arity;

    Segment_2
    operator()() const
    {  return Segment_2(); }

    Segment_2
    operator()( const Point_2& p, const Point_2& q) const
    {  return Segment_2(p, q); }
  };

  template <typename K>
  class Construct_segment_3
  {
    typedef typename K::Segment_3  Segment_3;
    typedef typename K::Point_3    Point_3;
  public:
    typedef Segment_3        result_type;
    typedef Arity_tag< 2 >   Arity;

    Segment_3
    operator()() const
    {  return Segment_3(); }

    Segment_3
    operator()( const Point_3& p, const Point_3& q) const
    {  return Segment_3(p, q); }
  };

  template <typename K>
  class Construct_sphere_3
  {
    typedef typename K::FT         FT;
    typedef typename K::Point_3    Point_3;
    typedef typename K::Sphere_3   Sphere_3;
  public:
    typedef Sphere_3        result_type;
    typedef Arity_tag< 4 >   Arity;

    Sphere_3
    operator()() const
    {  return Sphere_3(); }

    Sphere_3
    operator()( const Point_3& center, const FT& squared_radius,
	        Orientation orientation = COUNTERCLOCKWISE) const
    {  return Sphere_3(center, squared_radius, orientation); }

    Sphere_3
    operator()( const Point_3& p, const Point_3& q,
	        const Point_3& r, const Point_3& s) const
    {  return Sphere_3(p, q, r, s); }

    Sphere_3
    operator()( const Point_3& p, const Point_3& q, const Point_3& r,
	        Orientation orientation = COUNTERCLOCKWISE) const
    {  return Sphere_3(p, q, r, orientation); }

    Sphere_3
    operator()( const Point_3& p, const Point_3& q,
	        Orientation orientation = COUNTERCLOCKWISE) const
    {  return Sphere_3(p, q, orientation); }

    Sphere_3
    operator()( const Point_3& center,
	        Orientation orientation = COUNTERCLOCKWISE) const
    {  return Sphere_3(center, orientation); }
  };

#ifndef CGAL_NO_DEPRECATED_CODE
  template <typename K>
  class Construct_supporting_line_2
  {
    typedef typename K::Line_2     Line_2;
    typedef typename K::Ray_2      Ray_2;
    typedef typename K::Segment_2  Segment_2;
  public:
    typedef Line_2           result_type;
    typedef Arity_tag< 1 >   Arity;

    Line_2
    operator()( const Ray_2& r) const
    { return r.supporting_line(); }

    Line_2
    operator()( const Segment_2& s) const
    { return s.supporting_line(); }
  };

  template <typename K>
  class Construct_supporting_line_3
  {
    typedef typename K::Line_3     Line_3;
    typedef typename K::Ray_3      Ray_3;
    typedef typename K::Segment_3  Segment_3;
  public:
    typedef Line_3           result_type;
    typedef Arity_tag< 1 >   Arity;

    Line_3
    operator()( const Ray_3& r) const
    { return r.supporting_line(); }

    Line_3
    operator()( const Segment_3& s) const
    { return s.supporting_line(); }
  };
#endif // CGAL_NO_DEPRECATED_CODE

  template <typename K>
  class Construct_supporting_plane_3
  {
    typedef typename K::Triangle_3  Triangle_3;
    typedef typename K::Plane_3     Plane_3;
  public:
    typedef Plane_3          result_type;
    typedef Arity_tag< 1 >   Arity;

    Plane_3
    operator()( const Triangle_3& t) const
    { return t.supporting_plane(); }
  };

  template <typename K>
  class Construct_tetrahedron_3
  {
    typedef typename K::Tetrahedron_3   Tetrahedron_3;
    typedef typename K::Point_3         Point_3;
  public:
    typedef Tetrahedron_3    result_type;
    typedef Arity_tag< 4 >   Arity;

    Tetrahedron_3
    operator()() const
    { return Tetrahedron_3(); }

    Tetrahedron_3
    operator()( const Point_3& p, const Point_3& q,
	        const Point_3& r, const Point_3& s) const
    { return Tetrahedron_3(p, q, r, s); }
  };

  template <typename K>
  class Construct_triangle_2
  {
    typedef typename K::Triangle_2   Triangle_2;
    typedef typename K::Point_2      Point_2;
  public:
    typedef Triangle_2       result_type;
    typedef Arity_tag< 3 >   Arity;

    Triangle_2
    operator()() const
    { return Triangle_2(); }

    Triangle_2
    operator()( const Point_2& p, const Point_2& q, const Point_2& r) const
    { return Triangle_2(p, q, r); }
  };

  template <typename K>
  class Construct_triangle_3
  {
    typedef typename K::Triangle_3   Triangle_3;
    typedef typename K::Point_3      Point_3;
  public:
    typedef Triangle_3       result_type;
    typedef Arity_tag< 3 >   Arity;

    Triangle_3
    operator()() const
    { return Triangle_3(); }

    Triangle_3
    operator()( const Point_3& p, const Point_3& q, const Point_3& r) const
    { return Triangle_3(p, q, r); }
  };

  template <typename K>
  class Construct_vertex_2
  {
    typedef typename K::Point_2          Point_2;
    typedef typename K::Segment_2        Segment_2;
    typedef typename K::Iso_rectangle_2  Iso_rectangle_2;
    typedef typename K::Triangle_2       Triangle_2;
  public:
    typedef Point_2          result_type;
    typedef Arity_tag< 2 >   Arity;

    const Point_2 &
    operator()( const Segment_2& s, int i) const
    { return s.vertex(i); }

    const Point_2 &
    operator()( const Triangle_2& t, int i) const
    { return t.vertex(i); }

    const Point_2
    operator()( const Iso_rectangle_2& r, int i) const
    { return r.vertex(i); }
  };

  template <typename K>
  class Construct_vertex_3
  {
    typedef typename K::Point_3          Point_3;
    typedef typename K::Segment_3        Segment_3;
    typedef typename K::Iso_cuboid_3     Iso_cuboid_3;
    typedef typename K::Triangle_3       Triangle_3;
    typedef typename K::Tetrahedron_3    Tetrahedron_3;
  public:
    typedef Point_3          result_type;
    typedef Arity_tag< 2 >   Arity;

    const Point_3 &
    operator()( const Segment_3& s, int i) const
    { return s.vertex(i); }

    const Point_3 &
    operator()( const Triangle_3& t, int i) const
    { return t.vertex(i); }

    const Point_3
    operator()( const Iso_cuboid_3& r, int i) const
    { return r.vertex(i); }

    const Point_3 &
    operator()( const Tetrahedron_3& t, int i) const
    { return t.vertex(i); }
  };


  template <typename K>
  class Construct_bbox_2
  {
    typedef typename K::Point_2          Point_2;
    typedef typename K::Segment_2        Segment_2;
    typedef typename K::Iso_rectangle_2  Iso_rectangle_2;
    typedef typename K::Triangle_2       Triangle_2;
    typedef typename K::Circle_2         Circle_2;
  public:
    typedef Bbox_2          result_type;
    typedef Arity_tag< 1 >   Arity;

    Bbox_2
    operator()( const Point_2& p) const
    { return p.bbox(); }

    Bbox_2
    operator()( const Segment_2& s) const
    { return s.bbox(); }

    
    Bbox_2
    operator()( const Triangle_2& t) const
    { return t.bbox(); }

    Bbox_2
    operator()( const Iso_rectangle_2& r) const
    { return r.bbox(); }

    Bbox_2
    operator()( const Circle_2& c) const
    { return c.bbox(); }
  };


  template <typename K>
  class Construct_bbox_3
  {
    typedef typename K::Point_3          Point_3;
    typedef typename K::Segment_3        Segment_3;
    typedef typename K::Iso_cuboid_3     Iso_cuboid_3;
    typedef typename K::Triangle_3       Triangle_3;
    typedef typename K::Tetrahedron_3    Tetrahedron_3;
    typedef typename K::Sphere_3         Sphere_3;
  public:
    typedef Bbox_3          result_type;
    typedef Arity_tag< 1 >   Arity;

    Bbox_3
    operator()( const Point_3& p) const
    { return p.bbox(); }

    Bbox_3
    operator()( const Segment_3& s) const
    { return s.bbox(); }

    
    Bbox_3
    operator()( const Triangle_3& t) const
    { return t.bbox(); }

    Bbox_3
    operator()( const Iso_cuboid_3& r) const
    { return r.bbox(); }

    Bbox_3
    operator()( const Tetrahedron_3& t) const
    { return t.bbox(); }

    Bbox_3
    operator()( const Sphere_3& s) const
    { return s.bbox(); }
  };

  template <typename K>
  class Construct_cartesian_const_iterator_2
  {
    typedef typename K::Point_2          Point_2;
    typedef typename K::Cartesian_const_iterator_2
    Cartesian_const_iterator_2;
    
  public:
    typedef Cartesian_const_iterator_2 result_type;
    typedef Arity_tag< 1 >   Arity;

    Cartesian_const_iterator_2
    operator()( const Point_2& p) const
      {
	return p.cartesian_begin();
      }
    
    Cartesian_const_iterator_2
    operator()( const Point_2& p, int) const
    {
      return p.cartesian_end();
    }
  };

  template <typename K>
  class Construct_cartesian_const_iterator_3
  {
    typedef typename K::Point_3          Point_3;
    typedef typename K::Cartesian_const_iterator_3
    Cartesian_const_iterator_3;
    
  public:
    typedef Cartesian_const_iterator_3 result_type;
    typedef Arity_tag< 1 >   Arity;

    Cartesian_const_iterator_3
    operator()( const Point_3& p) const
      {
	return p.cartesian_begin();
      }
    
    Cartesian_const_iterator_3
    operator()( const Point_3& p, int) const
    {
      return p.cartesian_end();
    }
  };

  template <typename K>
  class Coplanar_3
  {
    typedef typename K::Point_3       Point_3;
    typedef typename K::Orientation_3 Orientation_3;
    Orientation_3 o;
  public:
    typedef bool             result_type;
    typedef Arity_tag< 4 >   Arity;

    Coplanar_3() {}
    Coplanar_3(const Orientation_3& o_) : o(o_) {}

    bool
    operator()( const Point_3& p, const Point_3& q,
	        const Point_3& r, const Point_3& s) const
    { 
      return o(p, q, r, s) == COPLANAR;
    }
  };

  template <typename K>
  class Counterclockwise_in_between_2
  {
    typedef typename K::Direction_2  Direction_2;
  public:
    typedef bool             result_type;
    typedef Arity_tag< 3 >   Arity;

    bool
    operator()( const Direction_2& p, const Direction_2& q,
	        const Direction_2& r) const
    {
        if ( q < p)
            return ( p < r )||( r <= q );
        else
            return ( p < r )&&( r <= q );
    }
  };

  template <typename K>
  class Do_intersect_2
  {
  public:
    typedef bool             result_type;
    typedef Arity_tag< 2 >   Arity;

    // There are 36 combinaisons, so I use a template.
    template <class T1, class T2>
    bool
    operator()(const T1& t1, const T2& t2) const
    { return CGALi::do_intersect(t1, t2, K()); }
  };

  template <typename K>
  class Do_intersect_3
  {
  public:
    typedef bool             result_type;
    typedef Arity_tag< 2 >   Arity;

    // There are x combinaisons, so I use a template.
    template <class T1, class T2>
    bool
    operator()(const T1& t1, const T2& t2) const
    { return CGALi::do_intersect(t1, t2, K()); }
  };

  template <typename K>
  class Equal_2
  {
  public:
    typedef bool             result_type;
    typedef Arity_tag< 2 >   Arity;

    // template to replace n different versions
    template <typename T>
    bool
    operator()(const T& p, const T& q) const
    { return p == q; }
  };

  template <typename K>
  class Equal_3
  {
  public:
    typedef bool             result_type;
    typedef Arity_tag< 2 >   Arity;

    // template to replace n different versions
    template <typename T>
    bool
    operator()(const T& p, const T& q) const
    { return p == q; }
  };

  template <typename K>
  class Has_on_boundary_2
  {
    typedef typename K::Point_2          Point_2;
    typedef typename K::Iso_rectangle_2  Iso_rectangle_2;
    typedef typename K::Circle_2         Circle_2;
    typedef typename K::Triangle_2       Triangle_2;
  public:
    typedef bool             result_type;
    typedef Arity_tag< 2 >   Arity;

    bool
    operator()( const Circle_2& c, const Point_2& p) const
    { return c.has_on_boundary(p); }

    bool
    operator()( const Triangle_2& t, const Point_2& p) const
    { return t.has_on_boundary(p); }

    bool
    operator()( const Iso_rectangle_2& r, const Point_2& p) const
    { return r.has_on_boundary(p); }
  };

  template <typename K>
  class Has_on_boundary_3
  {
    typedef typename K::Point_3          Point_3;
    typedef typename K::Iso_cuboid_3     Iso_cuboid_3;
    typedef typename K::Sphere_3         Sphere_3;
    typedef typename K::Tetrahedron_3    Tetrahedron_3;
    typedef typename K::Plane_3          Plane_3;
  public:
    typedef bool             result_type;
    typedef Arity_tag< 2 >   Arity;

    bool
    operator()( const Sphere_3& s, const Point_3& p) const
    { return s.has_on_boundary(p); }

    bool
    operator()( const Tetrahedron_3& t, const Point_3& p) const
    { return t.has_on_boundary(p); }

    bool
    operator()( const Iso_cuboid_3& c, const Point_3& p) const
    { return c.has_on_boundary(p); }
  };

  template <typename K>
  class Has_on_bounded_side_2
  {
    typedef typename K::Point_2          Point_2;
    typedef typename K::Iso_rectangle_2  Iso_rectangle_2;
    typedef typename K::Circle_2         Circle_2;
    typedef typename K::Triangle_2       Triangle_2;
  public:
    typedef bool             result_type;
    typedef Arity_tag< 2 >   Arity;

    bool
    operator()( const Circle_2& c, const Point_2& p) const
    { return c.has_on_bounded_side(p); }

    bool
    operator()( const Triangle_2& t, const Point_2& p) const
    { return t.has_on_bounded_side(p); }

    bool
    operator()( const Iso_rectangle_2& r, const Point_2& p) const
    { return r.has_on_bounded_side(p); }
  };

  template <typename K>
  class Has_on_bounded_side_3
  {
    typedef typename K::Point_3          Point_3;
    typedef typename K::Iso_cuboid_3     Iso_cuboid_3;
    typedef typename K::Sphere_3         Sphere_3;
    typedef typename K::Tetrahedron_3    Tetrahedron_3;
  public:
    typedef bool             result_type;
    typedef Arity_tag< 2 >   Arity;

    bool
    operator()( const Sphere_3& s, const Point_3& p) const
    { return s.has_on_bounded_side(p); }

    bool
    operator()( const Tetrahedron_3& t, const Point_3& p) const
    { return t.has_on_bounded_side(p); }

    bool
    operator()( const Iso_cuboid_3& c, const Point_3& p) const
    { return c.has_on_bounded_side(p); }
  };

  template <typename K>
  class Has_on_negative_side_2
  {
    typedef typename K::Point_2          Point_2;
    typedef typename K::Line_2           Line_2;
    typedef typename K::Circle_2         Circle_2;
    typedef typename K::Triangle_2       Triangle_2;
  public:
    typedef bool             result_type;
    typedef Arity_tag< 2 >   Arity;

    bool
    operator()( const Circle_2& c, const Point_2& p) const
    { return c.has_on_negative_side(p); }

    bool
    operator()( const Triangle_2& t, const Point_2& p) const
    { return t.has_on_negative_side(p); }

    bool
    operator()( const Line_2& l, const Point_2& p) const
    { return l.has_on_negative_side(p); }
  };

  template <typename K>
  class Has_on_negative_side_3
  {
    typedef typename K::Point_3          Point_3;
    typedef typename K::Plane_3          Plane_3;
    typedef typename K::Sphere_3         Sphere_3;
    typedef typename K::Tetrahedron_3    Tetrahedron_3;
  public:
    typedef bool             result_type;
    typedef Arity_tag< 2 >   Arity;

    bool
    operator()( const Sphere_3& s, const Point_3& p) const
    { return s.has_on_negative_side(p); }

    bool
    operator()( const Tetrahedron_3& t, const Point_3& p) const
    { return t.has_on_negative_side(p); }

    bool
    operator()( const Plane_3& pl, const Point_3& p) const
    { return pl.has_on_negative_side(p); }
  };

  template <typename K>
  class Has_on_positive_side_2
  {
    typedef typename K::Point_2          Point_2;
    typedef typename K::Line_2           Line_2;
    typedef typename K::Circle_2         Circle_2;
    typedef typename K::Triangle_2       Triangle_2;
  public:
    typedef bool             result_type;
    typedef Arity_tag< 2 >   Arity;

    bool
    operator()( const Circle_2& c, const Point_2& p) const
    { return c.has_on_positive_side(p); }

    bool
    operator()( const Triangle_2& t, const Point_2& p) const
    { return t.has_on_positive_side(p); }

    bool
    operator()( const Line_2& l, const Point_2& p) const
    { return l.has_on_positive_side(p); }
  };

  template <typename K>
  class Has_on_positive_side_3
  {
    typedef typename K::Point_3          Point_3;
    typedef typename K::Plane_3          Plane_3;
    typedef typename K::Sphere_3         Sphere_3;
    typedef typename K::Tetrahedron_3    Tetrahedron_3;
  public:
    typedef bool             result_type;
    typedef Arity_tag< 2 >   Arity;

    bool
    operator()( const Sphere_3& s, const Point_3& p) const
    { return s.has_on_positive_side(p); }

    bool
    operator()( const Tetrahedron_3& t, const Point_3& p) const
    { return t.has_on_positive_side(p); }

    bool
    operator()( const Plane_3& pl, const Point_3& p) const
    { return pl.has_on_positive_side(p); }
  };

  template <typename K>
  class Has_on_unbounded_side_2
  {
    typedef typename K::Point_2          Point_2;
    typedef typename K::Iso_rectangle_2  Iso_rectangle_2;
    typedef typename K::Circle_2         Circle_2;
    typedef typename K::Triangle_2       Triangle_2;
  public:
    typedef bool             result_type;
    typedef Arity_tag< 2 >   Arity;

    bool
    operator()( const Circle_2& c, const Point_2& p) const
    { return c.has_on_unbounded_side(p); }

    bool
    operator()( const Triangle_2& t, const Point_2& p) const
    { return t.has_on_unbounded_side(p); }

    bool
    operator()( const Iso_rectangle_2& r, const Point_2& p) const
    { return r.has_on_unbounded_side(p); }
  };

  template <typename K>
  class Has_on_unbounded_side_3
  {
    typedef typename K::Point_3          Point_3;
    typedef typename K::Iso_cuboid_3     Iso_cuboid_3;
    typedef typename K::Sphere_3         Sphere_3;
    typedef typename K::Tetrahedron_3    Tetrahedron_3;
  public:
    typedef bool             result_type;
    typedef Arity_tag< 2 >   Arity;

    bool
    operator()( const Sphere_3& s, const Point_3& p) const
    { return s.has_on_unbounded_side(p); }

    bool
    operator()( const Tetrahedron_3& t, const Point_3& p) const
    { return t.has_on_unbounded_side(p); }

    bool
    operator()( const Iso_cuboid_3& c, const Point_3& p) const
    { return c.has_on_unbounded_side(p); }
  };

  template <typename K>
  class Has_on_2
  {
    typedef typename K::Point_2          Point_2;
    typedef typename K::Line_2           Line_2;
    typedef typename K::Ray_2            Ray_2;
    typedef typename K::Segment_2        Segment_2;
  public:
    typedef bool             result_type;
    typedef Arity_tag< 2 >   Arity;

    bool
    operator()( const Line_2& l, const Point_2& p) const
    { return l.has_on(p); }

    bool
    operator()( const Ray_2& r, const Point_2& p) const
    { return r.has_on(p); }

    bool
    operator()( const Segment_2& s, const Point_2& p) const
    { return s.has_on(p); }
  };

  template <typename K>
  class Intersect_2
  {
    typedef typename K::Object_2    Object_2;
  public:
    typedef Object_2         result_type;
    typedef Arity_tag< 2 >   Arity;

    // 25 possibilities, so I keep the template.
    template <class T1, class T2>
    Object_2
    operator()(const T1& t1, const T2& t2) const
    { return CGALi::intersection(t1, t2, K()); }
  };

  template <typename K>
  class Intersect_3
  {
    typedef typename K::Object_3    Object_3;
    typedef typename K::Plane_3     Plane_3;
  public:
    typedef Object_3         result_type;
    typedef Arity_tag< 2 >   Arity;

    // n possibilities, so I keep the template.
    template <class T1, class T2>
    Object_3
    operator()(const T1& t1, const T2& t2) const
    { return CGALi::intersection(t1, t2, K() ); }

    Object_3
    operator()(const Plane_3& pl1, const Plane_3& pl2, const Plane_3& pl3)const
    { return CGALi::intersection(pl1, pl2, pl3, K() ); }
  };

  template <typename K>
  class Is_degenerate_2
  {
    typedef typename K::Circle_2          Circle_2;
    typedef typename K::Iso_rectangle_2   Iso_rectangle_2;
    typedef typename K::Line_2            Line_2;
    typedef typename K::Ray_2             Ray_2;
    typedef typename K::Segment_2         Segment_2;
    typedef typename K::Triangle_2        Triangle_2;
  public:
    typedef bool             result_type;
    typedef Arity_tag< 1 >   Arity;

    bool
    operator()( const Circle_2& c) const
    { return c.is_degenerate(); }

    bool
    operator()( const Iso_rectangle_2& r) const
    { return r.is_degenerate(); }

    bool
    operator()( const Line_2& l) const
    { return l.is_degenerate(); }

    bool
    operator()( const Ray_2& r) const
    { return r.is_degenerate(); }

    bool
    operator()( const Segment_2& s) const
    { return s.is_degenerate(); }

    bool
    operator()( const Triangle_2& t) const
    { return t.is_degenerate(); }
  };

  template <typename K>
  class Is_degenerate_3
  {
    typedef typename K::Iso_cuboid_3      Iso_cuboid_3;
    typedef typename K::Line_3            Line_3;
    typedef typename K::Plane_3           Plane_3;
    typedef typename K::Ray_3             Ray_3;
    typedef typename K::Segment_3         Segment_3;
    typedef typename K::Sphere_3          Sphere_3;
    typedef typename K::Triangle_3        Triangle_3;
    typedef typename K::Tetrahedron_3     Tetrahedron_3;
  public:
    typedef bool             result_type;
    typedef Arity_tag< 1 >   Arity;

    bool
    operator()( const Iso_cuboid_3& c) const
    { return c.is_degenerate(); }

    bool
    operator()( const Line_3& l) const
    { return l.is_degenerate(); }

    bool
    operator()( const Plane_3& pl) const
    { return pl.is_degenerate(); }

    bool
    operator()( const Ray_3& r) const
    { return r.is_degenerate(); }

    bool
    operator()( const Segment_3& s) const
    { return s.is_degenerate(); }

    bool
    operator()( const Sphere_3& s) const
    { return s.is_degenerate(); }

    bool
    operator()( const Triangle_3& t) const
    { return t.is_degenerate(); }

    bool
    operator()( const Tetrahedron_3& t) const
    { return t.is_degenerate(); }
  };

  template <typename K>
  class Is_horizontal_2
  {
    typedef typename K::Line_2    Line_2;
    typedef typename K::Segment_2 Segment_2;
    typedef typename K::Ray_2     Ray_2;
  public:
    typedef bool             result_type;
    typedef Arity_tag< 1 >   Arity;

    bool
    operator()( const Line_2& l) const
    { return l.is_horizontal(); }

    bool
    operator()( const Segment_2& s) const
    { return s.is_horizontal(); }

    bool
    operator()( const Ray_2& r) const
    { return r.is_horizontal(); }
  };

  template <typename K>
  class Is_vertical_2
  {
    typedef typename K::Line_2    Line_2;
    typedef typename K::Segment_2 Segment_2;
    typedef typename K::Ray_2     Ray_2;
  public:
    typedef bool             result_type;
    typedef Arity_tag< 1 >   Arity;

    bool
    operator()( const Line_2& l) const
    { return l.is_vertical(); }

    bool
    operator()( const Segment_2& s) const
    { return s.is_vertical(); }

    bool
    operator()( const Ray_2& r) const
    { return r.is_vertical(); }
  };

  template <typename K>
  class Left_turn_2
  {
    typedef typename K::Point_2        Point_2;
    typedef typename K::Orientation_2  Orientation_2;
    Orientation_2 o;
  public:
    typedef bool             result_type;
    typedef Arity_tag< 3 >   Arity;

    Left_turn_2() {}
    Left_turn_2(const Orientation_2& o_) : o(o_) {}
  
    bool
    operator()(const Point_2& p, const Point_2& q, const Point_2& r) const
    { return o(p, q, r) == LEFT_TURN; }
  };

  template <typename K>
  class Less_rotate_ccw_2
  {
    typedef typename K::Point_2        Point_2;
    typedef typename K::Orientation_2  Orientation_2;
    typedef typename K::Collinear_are_ordered_along_line_2 
    Collinear_are_ordered_along_line_2;
    Orientation_2 o;
    Collinear_are_ordered_along_line_2 co;
  public:
    typedef bool             result_type;
    typedef Arity_tag< 3 >   Arity;

    Less_rotate_ccw_2() {}
    Less_rotate_ccw_2(const Orientation_2& o_, 
		      const Collinear_are_ordered_along_line_2& co_) 
      : o(o_), co(co_)
    {}

    bool
    operator()(const Point_2& r, const Point_2& p, const Point_2& q) const
    {
      Orientation ori = o(r, p, q);
      if ( ori == LEFT_TURN )
	return true;
      else if ( ori == RIGHT_TURN )
	return false;
      else
	{
	  if (p == r) return false;
	  if (q == r) return true;
	  if (p == q) return false;
	  return co( r, q, p);
	}
    }
  };

  template <typename K>
  class Oriented_side_2
  {
    typedef typename K::Point_2     Point_2;
    typedef typename K::Circle_2    Circle_2;
    typedef typename K::Line_2      Line_2;
    typedef typename K::Triangle_2  Triangle_2;
  public:
    typedef Oriented_side    result_type;
    typedef Arity_tag< 2 >   Arity;

    Oriented_side
    operator()( const Circle_2& c, const Point_2& p) const
    { return c.oriented_side(p); }

    Oriented_side
    operator()( const Line_2& l, const Point_2& p) const
    { return l.oriented_side(p); }

    Oriented_side
    operator()( const Triangle_2& t, const Point_2& p) const
    { return t.oriented_side(p); }
  };

  template <typename K>
  class Oriented_side_3
  {
    typedef typename K::Point_3        Point_3;
    typedef typename K::Tetrahedron_3  Tetrahedron_3;
    typedef typename K::Plane_3        Plane_3;
    typedef typename K::Sphere_3       Sphere_3;
  public:
    typedef Oriented_side    result_type;
    typedef Arity_tag< 2 >   Arity;

    Oriented_side
    operator()( const Sphere_3& s, const Point_3& p) const
    { return s.oriented_side(p); }

    Oriented_side
    operator()( const Plane_3& pl, const Point_3& p) const
    { return pl.oriented_side(p); }

    Oriented_side
    operator()( const Tetrahedron_3& t, const Point_3& p) const
    { return t.oriented_side(p); }
  };

} // namespace CommonKernelFunctors

CGAL_END_NAMESPACE

#endif // CGAL_KERNEL_FUNCTION_OBJECTS_H
