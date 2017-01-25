// Copyright (c) 1999-2003  ETH Zurich (Switzerland).
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
// Author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch> and
//                 Emo Welzl <emo@inf.ethz.ch>

#ifndef CGAL_MIN_QUADRILATERAL_TRAITS_2_H
#define CGAL_MIN_QUADRILATERAL_TRAITS_2_H 1

#include <CGAL/license/Bounding_volumes.h>


#include <CGAL/basic.h>
#include <CGAL/Optimisation/assertions.h>
#include <CGAL/Point_2.h>
#include <CGAL/Direction_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/intersections.h>
#include <CGAL/utility.h>
#include <functional>
#include <vector>

namespace CGAL {

namespace Optimisation {
  template < class Kernel_ >
  struct Min_rectangle_2 {
    typedef Kernel_                       Kernel;
    typedef typename Kernel::Point_2      Point_2;
    typedef typename Kernel::Direction_2  Direction_2;

    Min_rectangle_2(const Point_2& q1,
                    const Direction_2& e,
                    const Point_2& q2,
                    const Point_2& q3,
                    const Point_2& q4)
    : p1(q1), p2(q2), p3(q3), p4(q4), d(e)
    {}

    Point_2 p1, p2, p3, p4;
    Direction_2 d;
  };

  template < class Kernel_ >
  struct Min_parallelogram_2 {
    typedef Kernel_                       Kernel;
    typedef typename Kernel::Point_2      Point_2;
    typedef typename Kernel::Direction_2  Direction_2;

    Min_parallelogram_2(const Point_2& q1, const Direction_2& e1,
                        const Point_2& q2, const Direction_2& e2,
                        const Point_2& q3, const Point_2& q4)
    : p1(q1), p2(q2), p3(q3), p4(q4), d1(e1), d2(e2)
    {}

    Point_2 p1, p2, p3, p4;
    Direction_2 d1, d2;
  };

  template < class Kernel >
  std::ostream& operator<<(std::ostream& o,
                           const Min_parallelogram_2<Kernel>& p)
  {
    if (p.d1 == p.d2)
      o << "Dirs equal!\n";
    return o << "Para(\t" << p.p1 << ",\n"
             << p.d1 << ",\n"
             << p.p2 << ",\n"
             << p.d2 << ",\n"
             << p.p3 << ",\n"
             << p.p4 << ")";
  }
} // namespace Optimisation 

template < class K_ >
struct Min_quadrilateral_default_traits_2 {
  // types inherited from Kernel
  typedef K_                                         Kernel;
  typedef typename Kernel::RT                        RT;
  typedef typename Kernel::Point_2                   Point_2;
  typedef typename Kernel::Vector_2                  Vector_2;
  typedef typename Kernel::Direction_2               Direction_2;
  typedef typename Kernel::Line_2                    Line_2;

  // predicates and constructions inherited from Kernel
  typedef typename Kernel::Equal_2                   Equal_2;
  typedef typename Kernel::Less_xy_2                 Less_xy_2;
  typedef typename Kernel::Less_yx_2                 Less_yx_2;
  typedef typename Kernel::Has_on_negative_side_2    Has_on_negative_side_2;
  typedef typename Kernel::Compare_angle_with_x_axis_2
    Compare_angle_with_x_axis_2;
  typedef typename Kernel::Construct_vector_2        Construct_vector_2;
  typedef typename Kernel::Construct_direction_2     Construct_direction_2;
  typedef typename Kernel::Construct_line_2          Construct_line_2;
  typedef typename Kernel::Construct_perpendicular_vector_2
    Construct_perpendicular_vector_2;
  typedef typename Kernel::Construct_opposite_direction_2
    Construct_opposite_direction_2;

  // used for expensive precondition checks only:
  typedef typename Kernel::Orientation_2             Orientation_2;

protected:
  // used internally
  Construct_line_2              line;
  typename Kernel::Intersect_2  isec;

public:
  // new types
  typedef Optimisation::Min_rectangle_2<Kernel>      Rectangle_2;
  typedef Optimisation::Min_parallelogram_2<Kernel>  Parallelogram_2;
  typedef Triple<Point_2,Direction_2,Point_2>        Strip_2;

  // new predicates
  struct Area_less_rectangle_2
  : public std::binary_function< Rectangle_2, Rectangle_2, bool >
  {
    RT
    area_numerator(const Rectangle_2& r, Cartesian_tag) const
    {
      return
      (r.d.dx() * (r.p3.y() - r.p1.y()) + r.d.dy() * (r.p1.x() - r.p3.x())) *
      (r.d.dy() * (r.p2.y() - r.p4.y()) + r.d.dx() * (r.p2.x() - r.p4.x()));
    }
    
    RT
    area_denominator(const Rectangle_2& r, Cartesian_tag) const
    { return CGAL_NTS square(r.d.dx()) + CGAL_NTS square(r.d.dy()); }
    
    RT
    area_numerator(const Rectangle_2& r, Homogeneous_tag) const
    {
      return
      (r.d.dx() * (r.p3.hy() * r.p1.hw() - r.p1.hy() * r.p3.hw()) +
       r.d.dy() * (r.p1.hx() * r.p3.hw() - r.p3.hx() * r.p1.hw())) *
      (r.d.dy() * (r.p2.hy() * r.p4.hw() - r.p4.hy() * r.p2.hw()) +
       r.d.dx() * (r.p2.hx() * r.p4.hw() - r.p4.hx() * r.p2.hw()));
    }
    
    RT
    area_denominator(const Rectangle_2& r, Homogeneous_tag) const
    {
      return r.p1.hw() * r.p2.hw() * r.p3.hw() * r.p4.hw() *
        (CGAL_NTS square(r.d.dx()) + CGAL_NTS square(r.d.dy()));
    }
  
    bool
    operator()(const Rectangle_2& p, const Rectangle_2& q) const
    {
      typedef typename Kernel::Rep_tag Rep_tag;
      Rep_tag tag  CGAL_SUNPRO_INITIALIZE(= Rep_tag());
      return area_numerator(p, tag) * area_denominator(q, tag) <
             area_denominator(p, tag) * area_numerator(q, tag);
    }
  };
  struct Area_less_parallelogram_2
  : public std::binary_function< Parallelogram_2,
                                      Parallelogram_2,
                                      bool >
  {
    RT
    area_numerator(const Parallelogram_2& r, Cartesian_tag) const
    {
      return
      (r.d1.dx() * (r.p3.y() - r.p1.y()) -
       r.d1.dy() * (r.p3.x() - r.p1.x())) *
      (r.d2.dx() * (r.p4.y() - r.p2.y()) -
       r.d2.dy() * (r.p4.x() - r.p2.x()));
    }
    
    RT
    area_denominator(const Parallelogram_2& r, Cartesian_tag) const
    { return r.d1.dx() * r.d2.dy() - r.d1.dy() * r.d2.dx(); }
    
    RT
    area_numerator(const Parallelogram_2& r, Homogeneous_tag) const
    {
      return
      (r.d1.dx() * (r.p3.hy() * r.p1.hw() - r.p1.hy() * r.p3.hw()) -
       r.d1.dy() * (r.p3.hx() * r.p1.hw() - r.p1.hx() * r.p3.hw())) *
      (r.d2.dx() * (r.p4.hy() * r.p2.hw() - r.p2.hy() * r.p4.hw()) -
       r.d2.dy() * (r.p4.hx() * r.p2.hw() - r.p2.hx() * r.p4.hw()));
    }
    
    RT
    area_denominator(const Parallelogram_2& r, Homogeneous_tag) const
    {
      return r.p1.hw() * r.p2.hw() * r.p3.hw() * r.p4.hw() *
        (r.d1.dx() * r.d2.dy() - r.d1.dy() * r.d2.dx());
    }
  
    bool
    operator()(const Parallelogram_2& p, const Parallelogram_2& q) const
    {
      typedef typename Kernel::Rep_tag Rep_tag;
      Rep_tag tag  CGAL_SUNPRO_INITIALIZE(= Rep_tag());
      return area_numerator(p, tag) * area_denominator(q, tag) <
             area_denominator(p, tag) * area_numerator(q, tag);
    }
  };
  struct Width_less_strip_2
  : public std::binary_function< Strip_2, Strip_2, bool >
  {
    RT
    width_numerator(const Strip_2& r, Cartesian_tag) const
    {
      return CGAL_NTS square(
        r.second.dx() * (r.third.y() - r.first.y()) +
        r.second.dy() * (r.first.x() - r.third.x()));
    }
    
    RT
    width_denominator(const Strip_2& r, Cartesian_tag) const
    { return CGAL_NTS square(r.second.dx()) + CGAL_NTS square(r.second.dy()); }
    
    RT
    width_numerator(const Strip_2& r, Homogeneous_tag) const
    {
      return CGAL_NTS square(
        r.second.dx() *
          (r.third.hy() * r.first.hw() - r.first.hy() * r.third.hw()) +
        r.second.dy() *
          (r.first.hx() * r.third.hw() - r.third.hx() * r.first.hw()));
    }
    
    RT
    width_denominator(const Strip_2& r, Homogeneous_tag) const {
      return CGAL_NTS square(r.first.hw() * r.third.hw()) *
        (CGAL_NTS square(r.second.dx()) + CGAL_NTS square(r.second.dy()));
    }
  
    bool
    operator()(const Strip_2& p, const Strip_2& q) const
    {
      typedef typename Kernel::Rep_tag Rep_tag;
      Rep_tag tag  CGAL_SUNPRO_INITIALIZE(= Rep_tag());
      return width_numerator(p, tag) * width_denominator(q, tag) <
             width_denominator(p, tag) * width_numerator(q, tag);
    }
  };

  // new constructions
  struct Construct_vector_from_direction_2
  : public std::unary_function<Direction_2,Vector_2>
  {
    Vector_2 operator()(const Direction_2& d) const { return d.vector(); }
  };
  struct Construct_rectangle_2
  {
    Rectangle_2
    operator()(const Point_2& p1,
               const Direction_2& d1,
               const Point_2& p2,
               const Point_2& p3,
               const Point_2& p4) const
    { return Rectangle_2(p1, d1, p2, p3, p4); }
  };
  template < class OutputIterator >
  OutputIterator
  copy_rectangle_vertices_2(const Rectangle_2& r, OutputIterator o) const
  {
    return copy_parallelogram_vertices_2(
      construct_parallelogram_2_object()(
        r.p1,
        r.d,
        r.p2,
        construct_direction_2_object()(
          construct_perpendicular_vector_2_object()(
            construct_vector_from_direction_2_object()(r.d),
            CLOCKWISE)),
        r.p3,
        r.p4),
        o);
  }
  struct Construct_parallelogram_2
  {
    Parallelogram_2
    operator()(const Point_2& p1,
               const Direction_2& d1,
               const Point_2& p2,
               const Direction_2& d2,
               const Point_2& p3,
               const Point_2& p4) const
    { return Parallelogram_2(p1, d1, p2, d2, p3, p4); }
  };
  template < class OutputIterator >
  OutputIterator
  copy_parallelogram_vertices_2(
    const Parallelogram_2& r, OutputIterator o) const
  {
    Point_2 tmp;
    Line_2  tmpl;
    Object  tmpo;
  
    tmpo = isec(line(r.p1, r.d1), line(r.p2, r.d2));
    if (assign(tmp, tmpo)) {
      *o++ = tmp;
    } else {
      CGAL_optimisation_assertion_code(bool test1 =)
      assign(tmpl, tmpo);
      CGAL_optimisation_assertion(test1);
      *o++ = r.p1;
    }
    tmpo = isec(line(r.p3, r.d1), line(r.p2, r.d2));
    if (assign(tmp, tmpo)) {
      *o++ = tmp;
    } else {
      CGAL_optimisation_assertion_code(bool test1 =)
      assign(tmpl, tmpo);
      CGAL_optimisation_assertion(test1);
      *o++ = r.p2;
    }
    tmpo = isec(line(r.p3, r.d1), line(r.p4, r.d2));
    if (assign(tmp, tmpo)) {
      *o++ = tmp;
    } else {
      CGAL_optimisation_assertion_code(bool test1 =)
      assign(tmpl, tmpo);
      CGAL_optimisation_assertion(test1);
      *o++ = r.p3;
    }
    tmpo = isec(line(r.p1, r.d1), line(r.p4, r.d2));
    if (assign(tmp, tmpo)) {
      *o++ = tmp;
    } else {
      CGAL_optimisation_assertion_code(bool test1 =)
      assign(tmpl, tmpo);
      CGAL_optimisation_assertion(test1);
      *o++ = r.p3;
    }
    return o;
  }
  struct Construct_strip_2
  {
    Strip_2
    operator()(const Point_2& p1,
               const Direction_2& d1,
               const Point_2& p2) const
    { return Strip_2(p1, d1, p2); }
  };
  template < class OutputIterator >
  OutputIterator
  copy_strip_lines_2(const Strip_2& r, OutputIterator o) const
  {
    *o++ = line(r.first, r.second);
    *o++ = line(r.third, r.second);
    return o;
  } 

  Equal_2      equal_2_object()   const { return Equal_2(); }
  Less_xy_2    less_xy_2_object() const { return Less_xy_2(); }
  Less_yx_2    less_yx_2_object() const { return Less_yx_2(); }
  
  Has_on_negative_side_2 has_on_negative_side_2_object() const
  { return Has_on_negative_side_2(); }
  
  Compare_angle_with_x_axis_2 compare_angle_with_x_axis_2_object() const
  { return Compare_angle_with_x_axis_2(); }
  
  Area_less_rectangle_2 area_less_rectangle_2_object() const
  { return Area_less_rectangle_2(); }
  
  Area_less_parallelogram_2 area_less_parallelogram_2_object() const
  { return Area_less_parallelogram_2(); }
  
  Width_less_strip_2 width_less_strip_2_object() const
  { return Width_less_strip_2(); }
  
  Construct_vector_2 construct_vector_2_object() const
  { return Construct_vector_2(); }
  
  Construct_vector_from_direction_2
  construct_vector_from_direction_2_object() const
  { return Construct_vector_from_direction_2(); }
  
  Construct_perpendicular_vector_2
  construct_perpendicular_vector_2_object() const
  { return Construct_perpendicular_vector_2(); }
  
  Construct_direction_2 construct_direction_2_object() const
  { return Construct_direction_2(); }
  
  Construct_opposite_direction_2
  construct_opposite_direction_2_object() const
  { return Construct_opposite_direction_2(); }
  
  Construct_line_2 construct_line_2_object() const
  { return Construct_line_2(); }
  
  Construct_rectangle_2 construct_rectangle_2_object() const
  { return Construct_rectangle_2(); }
  
  Construct_parallelogram_2 construct_parallelogram_2_object() const
  { return Construct_parallelogram_2(); }
  
  Construct_strip_2 construct_strip_2_object() const
  { return Construct_strip_2(); }
  
  Orientation_2 orientation_2_object() const { return Orientation_2(); }
  

};

} //namespace CGAL

#endif // ! (CGAL_MIN_QUADRILATERAL_TRAITS_2_H)

// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------
