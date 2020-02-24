// Copyright (c) 2009  GeometryFactory (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent Rineau


#ifndef CGAL_INTERNAL_TRIANGULATION_2_PROJECTION_TRAITS_BASE_3_H
#define CGAL_INTERNAL_TRIANGULATION_2_PROJECTION_TRAITS_BASE_3_H

#include <CGAL/license/Triangulation_2.h>


#include <CGAL/Profile_timer.h>
#include <CGAL/intersections.h>
#include <CGAL/predicates/sign_of_determinant.h>

namespace CGAL {

namespace TriangulationProjectionTraitsCartesianFunctors {

template <class Traits>
class Projected_orientation_with_normal_3
{
  // private members
  typename Traits::Vector_3 normal;

  // private type aliases
  typedef typename Traits::K K;
  typedef typename Traits::Point_2 Point;
  typedef typename Traits::Vector_3 Vector_3;
public:
  typedef typename K::Orientation Orientation;
  typedef Orientation result_type;

  Projected_orientation_with_normal_3(const Vector_3& normal_)
    : normal(normal_)
  {
    CGAL_PROFILER("Construct Projected_orientation_with_normal_3.")
    CGAL_TIME_PROFILER("Construct Projected_orientation_with_normal_3")
  }

  Orientation operator()(const Point& p,
			 const Point& q,
			 const Point& r) const
  {
    CGAL_PROFILER("Projected_orientation_with_normal_3::operator()");
    CGAL_TIME_PROFILER("Projected_orientation_with_normal_3::operator()");
    return orientation(q-p, r-p, normal);
  }
}; // end class Projected_orientation_with_normal_3<Traits>

template <class Traits>
class Projected_side_of_oriented_circle_with_normal_3
{
  // private members
  typename Traits::Vector_3 normal;

  // private types aliases
  typedef typename Traits::K K;
  typedef typename Traits::Point_2 Point;
  typedef typename Traits::Vector_3 Vector_3;
  typedef typename Traits::FT FT;

  typedef Projected_side_of_oriented_circle_with_normal_3<Traits> Self;

public:
  typedef typename K::Oriented_side Oriented_side;
  typedef Oriented_side result_type;

  Projected_side_of_oriented_circle_with_normal_3(const Vector_3& normal_)
    : normal(normal_)
  {
    CGAL_PROFILER("Construct Projected_side_of_oriented_circle_with_normal_3.")
    CGAL_TIME_PROFILER("Construct Projected_side_of_oriented_circle_with_normal_3.")
  }

  Oriented_side operator()(const Point& p,
			   const Point& q,
			   const Point& r,
			   const Point& t) const
  {
    CGAL_PROFILER("Projected_side_of_oriented_circle_with_normal_3::operator()")
    CGAL_TIME_PROFILER("Projected_side_of_oriented_circle_with_normal_3::operator()")
    const Vector_3& u = normal;

    const Vector_3 tp = p - t;
    const Vector_3 tq = q - t;
    const Vector_3 tr = r - t;

    const FT tp2 = tp * tp;
    const FT tq2 = tq * tq;
    const FT tr2 = tr * tr;
    const FT u2  =  u * u;

    const FT k_p = tp * u;
    const FT k_q = tq * u;
    const FT k_r = tr * u;

    return sign_of_determinant<FT>(
	tp.x(), tp.y(), tp.z(), (tp2 + k_p) * u2 - k_p * k_p,
	tr.x(), tr.y(), tr.z(), (tr2 + k_r) * u2 - k_r * k_r,
	tq.x(), tq.y(), tq.z(), (tq2 + k_q) * u2 - k_q * k_q,
	 u.x(),  u.y(),  u.z(), u2 * u2);
    // Note that q and r have been swapped in the determinant above, to
    // inverse its sign.
  }
}; // end class Projected_side_of_oriented_circle_with_normal_3

template <class Traits>
class Projected_squared_distance_with_normal_3
{
  // private members
  typename Traits::Vector_3 normal;

  // private types aliases
  typedef typename Traits::K K;
  typedef typename Traits::Point_2 Point;
  typedef typename Traits::Line_2 Line;
  typedef typename Traits::Vector_3 Vector_3;
  typedef typename Traits::FT FT;

public:
  Projected_squared_distance_with_normal_3(const Vector_3& normal_)
    : normal(normal_)
  {
    CGAL_PROFILER("Construct Projected_squared_distance_with_normal_3.")
    CGAL_TIME_PROFILER("Construct Projected_squared_distance_with_normal_3")
  }

  FT operator()(const Point& p, const Point& q)
  {
    return squared_distance(p, q);
  }

  FT operator()(const Line& line, const Point& p)
  {
    CGAL_PROFILER("Projected_squared_distance_with_normal_3::operator()")
    CGAL_TIME_PROFILER("Projected_squared_distance_with_normal_3::operator()")
    const Vector_3& vl = line.to_vector();
    const Point& pl = line.point();
    const Vector_3 v = cross_product(normal,
                                     vl);
    if(v == NULL_VECTOR) {
      // den == 0 if the line is vertical
      // In that case, the distance is the distance to the line
      const Vector_3 w = cross_product(pl - p,
                                       vl);
      return (w * w) / (vl * vl);
    }
    else {
      const FT det = determinant(normal,
                                 vl,
                                 pl - p);
      return (det * det) / ( v * v );
    }
  }
}; // end class Projected_squared_distance_with_normal_3

template <class Traits>
class Projected_intersect_3
{
  // private members
  typename Traits::Vector_3 normal;

  // private types aliases
  typedef typename Traits::K K;
  typedef typename Traits::Point_2 Point;
  typedef typename Traits::Line_2 Line;
  typedef typename Traits::Segment_2 Segment;
  typedef typename K::Plane_3 Plane_3;
  typedef typename Traits::Vector_3 Vector_3;
  typedef typename Traits::FT FT;
public:
  Projected_intersect_3(const Vector_3& normal_)
    : normal(normal_)
  {
    CGAL_PROFILER("Construct Projected_intersect_3")
    CGAL_TIME_PROFILER("Construct Projected_intersect_3")
  }

  Object operator()(const Segment& s1, const Segment& s2)
  {
    CGAL_PROFILER("Projected_intersect_3::operator()")
    CGAL_TIME_PROFILER("Projected_intersect_3::operator()")
    const Vector_3 u1 = cross_product(s1.to_vector(), normal);
    if(u1 == NULL_VECTOR)
      return K().intersect_3_object()(s1.supporting_line(), s2);

    const Vector_3 u2 = cross_product(s2.to_vector(), normal);
    if(u2 == NULL_VECTOR)
      return K().intersect_3_object()(s1, s2.supporting_line());

    const Plane_3 plane_1(s1.source(), u1);
    const Plane_3 plane_2(s2.source(), u2);

    Object planes_intersection = intersection(plane_1, plane_2);
    if(planes_intersection.empty()) {
      std::cerr << "planes_intersection is empty\n";
      return planes_intersection;
    }
    if(const Line* line = object_cast<Line>(&planes_intersection))
    {
      const Point& pi = line->point(0);
      if(cross_product(normal, pi - s1.source())
         * cross_product(normal, pi - s1.target()) > FT(0)
         ||
         cross_product(normal, pi - s2.source())
         * cross_product(normal, pi - s2.target()) > FT(0) )
      {
        // the intersection of the lines is not inside the segments
        std::cerr << "intersection not inside\n";
        return Object();
      }
      else
      {
        // Let the plane passing through s1.source() and with normal
        // the cross product of s1.to_vector() and s2.to_vector(). That
        // plane should intersect *l, now.
        return intersection(*line, Plane_3(s1.source(),
                                           cross_product(s1.to_vector(),
                                                         s2.to_vector())));
      }
    }
    if(object_cast<Plane_3>(&planes_intersection))
    {
      std::cerr << "coplanar lines\n";
      CGAL_error();
      return Object();
    }
    return Object();
  }
}; // end class Projected_intersect_3


template <class Traits>
class Less_along_axis
{
  // private members
  typedef typename Traits::Vector_3 Vector_3;
  typedef typename Traits::Point_2 Point;
  Vector_3 base;
public:
  Less_along_axis(const Vector_3& base) : base(base)
  {
    CGAL_PROFILER("Construct Less_along_axis")
    CGAL_TIME_PROFILER("Construct Less_along_axis")
  }

  typedef bool result_type;

  bool operator() (const Point &p, const Point &q) const {
    return base * (p - q) < 0;
  }
}; // end class Less_along_axis

template <class Traits>
class Compare_along_axis
{
  // private members
  typedef typename Traits::Vector_3 Vector_3;
  typedef typename Traits::Point_2 Point;
  Vector_3 base;
public:
  Compare_along_axis(const Vector_3& base) : base(base)
  {
    CGAL_PROFILER("Construct Compare_along_axis")
    CGAL_TIME_PROFILER("Construct Compare_along_axis")
  }

  typedef Comparison_result result_type;

  Comparison_result operator() (const Point &p, const Point &q) const {
    return compare(base * (p - q), 0);
  }
}; // end class Compare_along_axis

} // end namespace TriangulationProjectionTraitsCartesianFunctors


template < class Kernel >
class Triangulation_2_projection_traits_base_3
{
  typedef Triangulation_2_projection_traits_base_3<Kernel> Self;

  typename Kernel::Vector_3 n, b1, b2;

public:
  typedef typename Kernel::Vector_3 Vector_3;


  explicit Triangulation_2_projection_traits_base_3(const Vector_3& n_)
    : n(n_)
  {
    typedef typename Kernel::FT FT;
    typedef typename Kernel::Vector_3 Vector_3;

    const FT& nx = n.x();
    const FT& ny = n.y();
    const FT& nz = n.z();
    if(CGAL::abs(nz) >= CGAL::abs(ny)) {
      b1 = Vector_3(nz, 0, -nx);
    }
    else {
      b1 = Vector_3(ny, -nx, 0);
    }
    b2 = cross_product(n, b1);
  }

  const Vector_3& normal() const
  {
    return n;
  }

  const Vector_3& base1() const{
    return b1;
  }

  const Vector_3& base2() const{
    return b2;
  }

  typedef Kernel K;
  typedef typename K::FT          FT;
  typedef typename K::Point_3     Point_2;
  typedef typename K::Segment_3   Segment_2;
  typedef typename K::Vector_3    Vector_2;
  typedef typename K::Triangle_3  Triangle_2;
  typedef typename K::Line_3      Line_2;

  typedef typename K::Angle_3                                Angle_2;

  typedef TriangulationProjectionTraitsCartesianFunctors::
    Compare_along_axis<Self>                                 Compare_x_2;
  typedef TriangulationProjectionTraitsCartesianFunctors::
    Compare_along_axis<Self>                                 Compare_y_2;

  typedef TriangulationProjectionTraitsCartesianFunctors::
    Less_along_axis<Self>                                    Less_x_2;
  typedef TriangulationProjectionTraitsCartesianFunctors::
    Less_along_axis<Self>                                    Less_y_2;

  typedef TriangulationProjectionTraitsCartesianFunctors::
    Projected_orientation_with_normal_3<Self>                Orientation_2;

  typedef TriangulationProjectionTraitsCartesianFunctors::
    Projected_side_of_oriented_circle_with_normal_3<Self>    Side_of_oriented_circle_2;

  typedef TriangulationProjectionTraitsCartesianFunctors::
  Projected_squared_distance_with_normal_3<Self>             Compute_squared_distance_2;

  typedef TriangulationProjectionTraitsCartesianFunctors::
  Projected_intersect_3<Self>                                Intersect_2;

  typedef typename K::Construct_point_3   Construct_point_2;
  typedef typename K::Construct_weighted_point_3  Construct_weighted_point_2;
  typedef typename K::Construct_segment_3  Construct_segment_2;
  typedef typename K::Construct_vector_3   Construct_vector_2;
  typedef typename K::Construct_line_3     Construct_line_2;
  typedef typename K::Construct_triangle_3 Construct_triangle_2;

  typedef typename K::Construct_scaled_vector_3     Construct_scaled_vector_2;
  typedef typename K::Construct_translated_point_3  Construct_translated_point_2;
  typedef typename K::Construct_midpoint_3          Construct_midpoint_2;
  typedef typename K::Construct_circumcenter_3      Construct_circumcenter_2;

  typedef typename K::Compute_area_3                Compute_area_2;
  typedef typename K::Construct_bbox_3              Construct_bbox_2;

  Less_x_2
  less_x_2_object() const
  {
    return Less_x_2(this->base1());
  }

  Less_y_2
  less_y_2_object() const
  {
    return Less_y_2(this->base2());
  }

  Compare_x_2
  compare_x_2_object() const
  {
    return Compare_x_2(this->base1());
  }

  Compare_y_2
  compare_y_2_object() const
  {
    return Compare_y_2(this->base2());
  }

  Orientation_2
  orientation_2_object() const
  {
    return Orientation_2(this->normal());
  }

  Side_of_oriented_circle_2
  side_of_oriented_circle_2_object() const
  {
    return Side_of_oriented_circle_2(this->normal());
  }

  Compute_squared_distance_2
  compute_squared_distance_2_object() const
  {
    return Compute_squared_distance_2(this->normal());
  }

  Intersect_2
  intersect_2_object () const
  {
    return Intersect_2(this->normal());
  }

  Angle_2  angle_2_object() const
    {return Angle_2();}

  Construct_point_2  construct_point_2_object() const
    {return Construct_point_2();}

  Construct_weighted_point_2  construct_weighted_point_2_object() const
    {return Construct_weighted_point_2();}

  Construct_segment_2  construct_segment_2_object() const
    {return Construct_segment_2();}

  Construct_vector_2  construct_vector_2_object() const
    {return Construct_vector_2();}

  Construct_scaled_vector_2  construct_scaled_vector_2_object() const
    {return Construct_scaled_vector_2();}

  Construct_midpoint_2  construct_midpoint_2_object() const
    {return Construct_midpoint_2();}

  Construct_circumcenter_2  construct_circumcenter_2_object() const
    {return Construct_circumcenter_2();}

  Construct_translated_point_2  construct_translated_point_2_object() const
    {return Construct_translated_point_2();}

  Construct_line_2  construct_line_2_object() const
    {return Construct_line_2();}

  Construct_triangle_2  construct_triangle_2_object() const
    {return Construct_triangle_2();}

  Compute_area_2 compute_area_2_object() const
  {return Compute_area_2();}


  Construct_bbox_2  construct_bbox_2_object() const
    {return Construct_bbox_2();}


  // Special functor, not in the Kernel concept
  class Projection_to_plan {
    // Remeber: Point_2 is K::Point_3
    const Point_2& plane_point;
    const Vector_3& normal;
  public:
    // Return the projection of a point to a plane passing through
    // the point 'plane_point' and with orthogonal vector normal().
    Projection_to_plan(const Point_2& plane_point_, const Self& self)
      : plane_point(plane_point_),
        normal(self.normal())
    {}

    Point_2 operator()(const Point_2& point) const
    {
      return point +
        ( ( (plane_point - point) * normal ) / (normal * normal) ) * normal;
    }
  }; // end Projection_to_plan

  Projection_to_plan projection_to_plan_object(const Point_2& plane_point) const
  {
    return Projection_to_plan(plane_point, *this);
  }

}; // end class Triangulation_2_projection_traits_base_3<Kernel>

} // end namespace CGAL

#endif // CGAL_INTERNAL_TRIANGULATION_2_PROJECTION_TRAITS_BASE_3_H
