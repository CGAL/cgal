// Copyright (c) 1997-2013 INRIA Sophia-Antipolis (France).
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
// Author(s)     : Nico Kruithof <Nico@nghk.nl>

#ifndef CGAL_PERIODIC_2_TRIANGULATION_TRAITS_2_H
#define CGAL_PERIODIC_2_TRIANGULATION_TRAITS_2_H

#include <CGAL/license/Periodic_2_triangulation_2.h>


#include <CGAL/Point_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/Line_2.h>
#include <CGAL/Ray_2.h>
#include <CGAL/predicates_on_points_2.h>
#include <CGAL/basic_constructions_2.h>
#include <CGAL/distance_predicates_2.h>

#include <CGAL/triangulation_assertions.h>
#include <CGAL/Segment_2_Segment_2_intersection.h>

#include <CGAL/Periodic_2_offset_2.h>

namespace CGAL
{

template < class K, class Predicate_ >
class Traits_with_offsets_2_adaptor
{
  typedef K Kernel;
  typedef Predicate_ Predicate;

  typedef typename Kernel::Point_2       Point;
  typedef typename Kernel::Offset_2      Offset;

  // Use the construct_point_2 predicate from the kernel to convert the periodic points to Euclidean points
  typedef typename Kernel::Construct_point_2 Construct_point_2;
public:
  typedef typename Kernel::Iso_rectangle_2  Iso_rectangle_2;

public:
  typedef typename Predicate::result_type result_type;

  Traits_with_offsets_2_adaptor(const Iso_rectangle_2 * dom) : _domain(dom) { }

  result_type operator()(const Point& p0, const Point& p1,
                         const Offset& o0, const Offset& o1) const
  {
    return Predicate()(pp(p0, o0), pp(p1, o1));
  }
  result_type operator()(const Point& p0, const Point& p1, const Point& p2,
                         const Offset& o0, const Offset& o1, const Offset& o2) const
  {
    return Predicate()(pp(p0, o0), pp(p1, o1), pp(p2, o2));
  }
  result_type operator()(const Point& p0, const Point& p1,
                         const Point& p2, const Point& p3,
                         const Offset& o0, const Offset& o1,
                         const Offset& o2, const Offset& o3) const
  {
    return Predicate()(pp(p0, o0), pp(p1, o1), pp(p2, o2), pp(p3, o3));
  }

  result_type operator()(const Point& p0, const Point& p1) const
  {
    return Predicate()(p0, p1);
  }
  result_type operator()(const Point& p0, const Point& p1,
                         const Point& p2) const
  {
    return Predicate()(p0, p1, p2);
  }
  result_type operator()(const Point& p0, const Point& p1,
                         const Point& p2, const Point& p3) const
  {
    return Predicate()(p0, p1, p2, p3);
  }

private:
  Point pp(const Point &p, const Offset &o) const
  {
    return Point(p.x() + (_domain->xmax() - _domain->xmin()) * o.x(),
                 p.y() + (_domain->ymax() - _domain->ymin()) * o.y());
  }
public:
  const Iso_rectangle_2* _domain;
};

template < typename K, typename Construct_point_2_base>
class Periodic_2_construct_point_2 : public Construct_point_2_base
{
  typedef Construct_point_2_base           Base;
  typedef K                                Kernel;

public:
  typedef typename Kernel::Point_2         Point;
  typedef typename Kernel::Offset_2        Offset;
  typedef typename Kernel::Iso_rectangle_2 Iso_rectangle_2;

  typedef Point       result_type;

  Periodic_2_construct_point_2(const Iso_rectangle_2 & dom) : _dom(dom) { }

  using Base::operator();

  Point operator() ( const Point& p, const Offset& o ) const
  {
    return Point(p.x() + (_dom.xmax() - _dom.xmin()) * o.x(),
                 p.y() + (_dom.ymax() - _dom.ymin()) * o.y());
  }

private:
  const Iso_rectangle_2& _dom;
};


template < class Kernel, class Off = typename CGAL::Periodic_2_offset_2 >
class Periodic_2_triangulation_traits_base_2 : public Kernel
{
public: // Undocumented
  typedef Kernel                                                   K;
  typedef Kernel                                                   Kernel_base;
  typedef Off                                                      Offset_2;
  typedef Periodic_2_triangulation_traits_base_2<Kernel, Off>      Self;

  typedef typename K::FT                   FT;
  typedef typename K::RT                   RT;
  typedef typename K::Line_2               Line_2;
  typedef typename K::Ray_2                Ray_2;
  typedef typename K::Circle_2             Circle_2;
  typedef typename K::Aff_transformation_2 Aff_transformation_2;

public:
  typedef typename K::Point_2              Point_2;
  typedef typename K::Segment_2            Segment_2;
  typedef typename K::Vector_2             Vector_2;
  typedef typename K::Triangle_2           Triangle_2;
  typedef typename K::Iso_rectangle_2      Iso_rectangle_2;
  typedef Offset_2                         Periodic_2_offset_2;

  typedef Traits_with_offsets_2_adaptor<Self, typename K::Less_x_2>                   Less_x_2;
  typedef Traits_with_offsets_2_adaptor<Self, typename K::Less_y_2>                   Less_y_2;
  typedef Traits_with_offsets_2_adaptor<Self, typename K::Compare_x_2>                Compare_x_2;
  typedef Traits_with_offsets_2_adaptor<Self, typename K::Compare_y_2>                Compare_y_2;
  typedef Traits_with_offsets_2_adaptor<Self, typename K::Compare_xy_2>                Compare_xy_2;
  typedef Traits_with_offsets_2_adaptor<Self, typename K::Orientation_2>              Orientation_2;
  typedef Traits_with_offsets_2_adaptor<Self, typename K::Side_of_oriented_circle_2>  Side_of_oriented_circle_2;
  typedef Traits_with_offsets_2_adaptor<Self, typename K::Compute_squared_radius_2>   Compute_squared_radius_2;
  typedef Traits_with_offsets_2_adaptor<Self, typename K::Construct_center_2>         Construct_center_2;
  typedef Traits_with_offsets_2_adaptor<Self, typename K::Construct_circumcenter_2>   Construct_circumcenter_2;
  typedef Traits_with_offsets_2_adaptor<Self, typename K::Construct_bisector_2>       Construct_bisector_2;
  typedef Traits_with_offsets_2_adaptor<Self, typename K::Compare_distance_2>         Compare_distance_2;
  typedef Periodic_2_construct_point_2<Self, typename K::Construct_point_2>           Construct_point_2;
  typedef Traits_with_offsets_2_adaptor<Self, typename K::Construct_segment_2>        Construct_segment_2;
  typedef Traits_with_offsets_2_adaptor<Self, typename K::Construct_triangle_2>       Construct_triangle_2;
  typedef Traits_with_offsets_2_adaptor<Self, typename K::Construct_direction_2>      Construct_direction_2;
  typedef Traits_with_offsets_2_adaptor<Self, typename K::Construct_ray_2>            Construct_ray_2;

  Periodic_2_triangulation_traits_base_2() : _domain(Iso_rectangle_2(0, 0, 1, 1))
  {
  }
  Periodic_2_triangulation_traits_base_2(const Periodic_2_triangulation_traits_base_2 &) {}
  Periodic_2_triangulation_traits_base_2 &operator=
  (const Periodic_2_triangulation_traits_base_2 &)
  {
    return *this;
  }

  // Access
  void set_domain(const Iso_rectangle_2& domain)
  {
    _domain = domain;
  }

  Iso_rectangle_2 get_domain() const
  {
    return _domain;
  }

  // Predicates&_domain

  Compare_x_2
  compare_x_2_object() const
  {
    return Compare_x_2(&_domain);
  }

  Compare_y_2
  compare_y_2_object() const
  {
    return Compare_y_2(&_domain);
  }

  Compare_xy_2
  compare_xy_2_object() const
  {
    return Compare_xy_2(&_domain);
  }

  Orientation_2
  orientation_2_object() const
  {
    return Orientation_2(&_domain);
  }

  Side_of_oriented_circle_2
  side_of_oriented_circle_2_object() const
  {
    return Side_of_oriented_circle_2(&_domain);
  }

  Construct_circumcenter_2
  construct_circumcenter_2_object() const
  {
    return Construct_circumcenter_2(&_domain);
  }

  Compare_distance_2
  compare_distance_2_object() const
  {
    return Compare_distance_2(&_domain);
  }

  Construct_point_2  construct_point_2_object() const
  {
    return Construct_point_2(_domain);
  }

  Construct_segment_2  construct_segment_2_object() const
  {
    return Construct_segment_2(&_domain);
  }

  Construct_triangle_2  construct_triangle_2_object() const
  {
    return Construct_triangle_2(&_domain);
  }

protected:
  Iso_rectangle_2 _domain;
};


// Forward declaration for the filtered traits
template < typename K, typename Off = CGAL::Periodic_2_offset_2, bool Has_filtered_predicates = K::Has_filtered_predicates  >
class Periodic_2_triangulation_traits_2;

} //namespace CGAL


// Partial specialization for Filtered_kernel<CK>.
#include <CGAL/Periodic_2_triangulation_filtered_traits_2.h>
#include <CGAL/Filtered_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Periodic_2_triangulation_filtered_traits_2.h>

namespace CGAL
{

// This declaration is needed to break the cyclic dependency.
template < typename K, typename Off >
class Periodic_2_triangulation_filtered_traits_base_2;

template < class K, class Off, bool b>
class Periodic_2_triangulation_traits_2
  : public Periodic_2_triangulation_traits_base_2<K, Off>
{
};

template < typename K, typename Off >
class Periodic_2_triangulation_traits_2 < K, Off, true>
  : public Periodic_2_triangulation_filtered_traits_2 <K, Off >
{
public:
  typedef K  Kernel;
};

} //namespace CGAL

#endif // CGAL_PERIODIC_2_TRIANGULATION_TRAITS_2_H
