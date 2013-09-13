// Copyright (c) 2003  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Julia Floetotto

#ifndef CGAL_VORONOI_INTERSECTION_2_TRAITS_3_H
#define CGAL_VORONOI_INTERSECTION_2_TRAITS_3_H

#include <CGAL/Origin.h>
#include <CGAL/tags.h>
#include <CGAL/number_utils.h>
#include <CGAL/representation_tags.h>
#include <CGAL/predicates/predicates_for_voronoi_intersection_cartesian_2_3.h>
#include <CGAL/constructions/constructions_for_voronoi_intersection_cartesian_2_3.h>
#include <CGAL/predicates/predicates_for_voronoi_intersection_cartesian_2_3.h>
#include <CGAL/constructions/constructions_for_voronoi_intersection_cartesian_2_3.h>
#include <CGAL/function_objects.h>

namespace CGAL {

template <class Traits>
class Orientation_with_normal_plane_2_3
{
public:
  typedef typename Traits::Point_3     Point;
  typedef typename Traits::Vector_3    Vector;

  Orientation_with_normal_plane_2_3(const Vector& _normal,
				    const Traits& _traits)
    : normal(_normal), traits(_traits) {}

  Orientation operator()(const Point &p, const Point &q, const Point &r) const
    {
      return traits.orientation_3_object()(p,q,q+normal,r);
    }

private:
  Vector  normal;
  const Traits&  traits;
};

template < typename K >
class Side_of_plane_centered_sphere_2_3
{
public:
  typedef Oriented_side  result_type;

  typedef typename K::Point_3      Point;
  typedef typename K::Vector_3     Vector;
  typedef typename K::Plane_3      Plane;
  typedef typename K::Direction_3  Direction;

  Side_of_plane_centered_sphere_2_3(const Point& _p,
				    const Vector& _normal)
    : a(_p), normal(_normal) {}

  Oriented_side operator()(const Point &p, const Point &q,
                           const Point &r, const Point &t) const
    {
      return side_of_plane_centered_sphere(a,normal,p,q,r,t);
    }

  Oriented_side operator()(const Point &p, const Point &q,
                           const Point &r) const
    {
      return side_of_plane_centered_sphere(a,normal,p,q,r);
    }

  Oriented_side operator() (const Point &p, const Point &q) const
    {
      return side_of_plane_centered_sphere(a,normal,p,q);
    }

private:
  Point   a;
  Vector  normal;
};

template < typename K >
class Construct_plane_centered_circumcenter_3
{
public:
  typedef typename K::Point_3      Point;
  typedef typename K::Vector_3     Vector;

  Construct_plane_centered_circumcenter_3(const Point& _p,
					  const Vector& _normal)
    : a(_p), normal(_normal) {}

  Point operator()(const Point &p, const Point &q, const Point &r) const
    {
      return plane_centered_circumcenter_3(a,normal,p,q,r);
    }

private:
  Point   a;
  Vector  normal;
};


template < typename K >
class Construct_plane_intersected_bisector_3
{
public:
  typedef typename K::Point_3      Point;
  typedef typename K::Vector_3     Vector;
  typedef typename K::Line_3       Line;

  Construct_plane_intersected_bisector_3(const Point& _a,
					 const Vector& _normal)
    : a(_a), normal(_normal) {}

  Line operator() (const Point &p, const Point &q) const
    {
      return plane_intersected_bisector_3(a, normal, p, q);
    }

private:
  Point   a;
  Vector  normal;
};


template < typename K >
class Compare_first_projection_3
{
  //compares the projection of two points onto a second (non-trivial)
  // vector in the projection plane:
public:
  typedef typename K::Point_3      Point;
  typedef typename K::Vector_3     Vector;
  typedef typename K::FT           Coord_type;

  Compare_first_projection_3(const Vector& _normal)
    : normal(_normal) {}

  Comparison_result operator() (const Point &p, const Point &q) const
    {
      if(normal.x()!=Coord_type(0))
	return (Comparison_result) CGAL_NTS
	  sign(Vector(normal.y(),-normal.x(),Coord_type(0))*(p-q));
      if(normal.y()!= Coord_type(0))
	return (Comparison_result) CGAL_NTS
	  sign(Vector(-normal.y(),normal.x(),Coord_type(0))*(p-q));

      CGAL_assertion(normal.z()!= Coord_type(0));
      return (Comparison_result) CGAL_NTS
	sign(Vector(-normal.z(),Coord_type(0),normal.x())*(p-q));
    }

private:
  Vector normal;
};

template < typename K >
class Compare_second_projection_3
{
  //compares the projection of two points onto a second (non-trivial)
  // vector in the projection plane:
public:
  typedef typename K::FT           Coord_type;
  typedef typename K::Vector_3     Vector;
  typedef typename K::Point_3      Point;

  Compare_second_projection_3(const Vector& _normal)
    : normal(_normal) {}

  Comparison_result operator()(const Point &p, const Point &q) const
    {
      if(normal.x()!=Coord_type(0))
	return (Comparison_result) CGAL_NTS
	  sign(Vector(normal.z(),Coord_type(0),-normal.x())*(p-q));
      if(normal.y()!= Coord_type(0))
	return (Comparison_result) CGAL_NTS
	  sign(Vector(Coord_type(0),normal.z(),-normal.y())*(p-q));

      CGAL_assertion(normal.z()!= Coord_type(0));
      return (Comparison_result) CGAL_NTS
	sign(Vector(Coord_type(0),-normal.z(),normal.y())*(p-q));
     }

private:
  Vector normal;
};

namespace Interpolation {

template < typename K >
class Compute_area_3
{
  //squareroot of compute_squared_area_3:
  //-> if no sqrt is supported, cast to double
public:
  typedef typename K::FT                     FT;
  typedef typename K::Point_3                Point;
  typedef typename K::Compute_squared_area_3 Compute_squared_area_3;

  FT operator()(const Point& p, const Point& q, const Point& r) const
  {
      typedef typename CGAL::Algebraic_structure_traits<FT> AST;
    FT squared_area = Compute_squared_area_3()(p,q,r);
    return cast_sqrt_to_double(squared_area, typename AST::Algebraic_category() );
  }

private:
  FT cast_sqrt_to_double(const FT& squared_area, Field_with_sqrt_tag) const
    {
      return CGAL_NTS sqrt(squared_area);
    }

  FT cast_sqrt_to_double(const FT& squared_area, Integral_domain_without_division_tag) const
    {
      double approx = CGAL_NTS to_double(squared_area);
      return CGAL_NTS sqrt(approx);
    }
};

} // namespace Interpolation

template < class K >
class Voronoi_intersection_2_traits_3
{
public:
  typedef K Rep;

  typedef typename Rep::RT                          Weight;
  typedef typename Rep::FT                          FT;

  //other types needed:
  typedef typename Rep::Point_3                     Point_2;
  typedef typename Rep::Segment_3                   Segment_2;
  typedef typename Rep::Triangle_3                  Triangle_2;
  typedef typename Rep::Line_3                      Line_2;
  typedef typename Rep::Ray_3                       Ray_2;
  typedef typename Rep::Vector_3                    Vector_2;
  typedef typename Rep::Construct_ray_3             Construct_ray_2;
  typedef typename Rep::Construct_segment_3         Construct_segment_2;
  typedef typename Rep::Construct_triangle_3        Construct_triangle_2;

  typedef typename Rep::Compare_distance_3          Compare_distance_2;
  //if no sqrt is supported, it casts to double:
  typedef Interpolation::Compute_area_3<Rep>                       Compute_area_2;

  //the regular triangulation traits model:
  //Traits::Point_2 is a 3D point!!
  typedef Point_2                        Weighted_point_2;
  typedef Point_2                        Bare_point;


  //specific tests:
  typedef Orientation_with_normal_plane_2_3<Rep>           Orientation_2;
  typedef Side_of_plane_centered_sphere_2_3<Rep>           Power_test_2;

  typedef Construct_plane_centered_circumcenter_3<Rep>
          Construct_weighted_circumcenter_2;

  typedef Construct_plane_intersected_bisector_3<Rep> Construct_radical_axis_2;

  typedef Compare_first_projection_3<Rep>                  Compare_x_2;
  typedef Compare_second_projection_3<Rep>                 Compare_y_2;

  typedef Compare_to_less<Compare_x_2>                     Less_x_2;
  typedef Compare_to_less<Compare_y_2>                     Less_y_2;

  //for certificated coordinate/neighbor computation:
  typedef typename Rep::Less_distance_to_point_3    Less_distance_to_point_2;
  typedef typename Rep::Compute_squared_distance_3  Compute_squared_distance_2;


  //instantiations and creation of functors:
  //for the triangulation:
  Orientation_2
  orientation_2_object() const
  { return Orientation_2(normal, Rep()); }

  Power_test_2
  power_test_2_object() const
  { return Power_test_2(a,normal); }

  Compare_distance_2 compare_distance_2_object() const
  { return Compare_distance_2(); }

  Compare_x_2
  compare_x_2_object() const
  { return Compare_x_2(normal); }

  Compare_y_2
  compare_y_2_object() const
  { return Compare_y_2(normal); }

  Less_x_2
  less_x_2_object() const
  { return compare_to_less(compare_x_2_object());; }

  Less_y_2
  less_y_2_object() const
  { return compare_to_less(compare_y_2_object());; }

  //for the coordinate computation:
  Compute_area_2 compute_area_2_object() const
  { return Compute_area_2(); }

  //for constructions of dual:
  Construct_weighted_circumcenter_2
  construct_weighted_circumcenter_2_object() const
  { return Construct_weighted_circumcenter_2(a,normal); }

  Construct_radical_axis_2  construct_radical_axis_2_object() const
  { return Construct_radical_axis_2(a,normal); }

  Construct_ray_2  construct_ray_2_object() const
  { return Construct_ray_2(); }

  Construct_segment_2  construct_segment_2_object() const
  { return Construct_segment_2(); }

  Construct_triangle_2  construct_triangle_2_object() const
  { return Construct_triangle_2(); }


  //for certification of coordinate/neighbor computation:
  Less_distance_to_point_2 less_distance_to_point_2_object() const
  { return Less_distance_to_point_2(); }

  Compute_squared_distance_2 compute_squared_distance_2_object() const
  { return Compute_squared_distance_2(); }

  Voronoi_intersection_2_traits_3<K>(const Point_2& _p = Point_2(),
				     const Vector_2& _normal = NULL_VECTOR)
      : a(_p), normal(_normal) {}

  const Vector_2& get_normal() const { return normal; }
  const Point_2&  get_point()  const { return a; }

private:
  //defining the intersection plane:
  Point_2 a;
  Vector_2 normal;
};


//put the homogeneous or cartesian tag
template < class Point, class Vector >
inline
Oriented_side
side_of_plane_centered_sphere(const Point &a, const Vector& n,
/*defines the plane*/
	   const Point &p, const Point &q, const Point &r, const Point &t)
{
  typedef typename Point::R::Rep_tag Tag;
  return side_of_plane_centered_sphere(a,n,p,q,r,t, Tag());
}

template < class Point, class Vector >
inline
Oriented_side
side_of_plane_centered_sphere(const Point &a, const Vector& n,
/*defines the plane*/
	   const Point &p, const Point &q, const Point &r)
{
  typedef typename Point::R::Rep_tag Tag;
  return side_of_plane_centered_sphere(a,n,p,q,r,Tag());
}

template < class Point, class Vector >
inline
Oriented_side
side_of_plane_centered_sphere(const Point &a, const Vector& n,
/*defines the plane*/
	   const Point &p, const Point &q)
{
  typedef typename Point::R::RT  RT;
  Comparison_result r = Compare<RT>()(-CGAL_NTS square(Vector(a,p)*n),
				      -CGAL_NTS square(Vector(a,q)*n));

  if(r == LARGER)    return ON_NEGATIVE_SIDE;
  else if (r == SMALLER) return ON_POSITIVE_SIDE;
  return ON_ORIENTED_BOUNDARY;
}

template < class Point, class Vector >
inline
Point
plane_centered_circumcenter_3(const Point &a, const Vector& n,
			      /*defines the plane*/
			      const Point &p, const Point &q, const Point &r)
{
  typedef typename Point::R::Rep_tag Tag;
  return plane_centered_circumcenter_3(a,n,p,q,r, Tag());
}

template < class Point, class Vector >
inline
typename Point::R::Line_3
plane_intersected_bisector_3(const Point &a, const Vector& n,
			     /*defines the plane*/
			     const Point &p, const Point &q)
{
  typedef typename Point::R::Rep_tag Tag;
  return plane_intersected_bisector_3(a,n,p,q, Tag());
}

///-----------------------------------------------------------
// Cartesian variants:
//
template < class Point, class Vector>
inline
Oriented_side
side_of_plane_centered_sphere(const Point &a, const Vector& n,
/*defines the plane*/
	   const Point &p, const Point &q, const Point &r, const Point &t,
	   Cartesian_tag)
{
  return side_of_plane_centered_sphereC3(a.x(), a.y(), a.z(),
                                         n.x(), n.y(), n.z(),
                                         p.x(), p.y(), p.z(),
                                         q.x(), q.y(), q.z(),
                                         r.x(), r.y(), r.z(),
                                         t.x(), t.y(), t.z());
}

template < class Point, class Vector >
inline
Oriented_side
side_of_plane_centered_sphere(const Point &a, const Vector& n,
/*defines the plane*/
	   const Point &p, const Point &q, const Point &r, Cartesian_tag)
{
  return side_of_plane_centered_sphereC3(a.x(), a.y(), a.z(),
                                         n.x(), n.y(), n.z(),
                                         p.x(), p.y(), p.z(),
                                         q.x(), q.y(), q.z(),
                                         r.x(), r.y(), r.z());
}

template < class Point, class Vector >
inline
Point
plane_centered_circumcenter_3(const Point &a, const Vector& n,
			      /*defines the plane*/
			      const Point &p, const Point &q, const Point &r,
			      Cartesian_tag)
{
  typename Point::R::RT x,y,z;
  plane_centered_circumcenterC3(a.x(), a.y(), a.z(),
                                n.x(), n.y(), n.z(),
                                p.x(), p.y(), p.z(),
                                q.x(), q.y(), q.z(),
                                r.x(), r.y(), r.z(),x,y,z);
  return Point(x,y,z);
}

template < class Point,  class Vector>
inline
typename Point::R::Line_3
plane_intersected_bisector_3(const Point &a, const Vector& n,
			     /*defines the plane*/
			     const Point &p, const Point &q, Cartesian_tag)
{
  typename Point::R::RT x1,y1,z1, x2,y2,z2;
  typedef typename Point::R::Line_3 Line;
  bisector_plane_intersectionC3(a.x(), a.y(), a.z(),
                                n.x(), n.y(), n.z(),
                                p.x(), p.y(), p.z(),
                                q.x(), q.y(), q.z(), x1,y1,z1,x2,y2,z2);

  return Line(Point(x1,y1,z1), Point(x2,y2,z2));
}


// Homogeneous variants.

// The 3 following call the cartesian version over FT, because an
// homogeneous special version has not yet been written.

template <class Point,  class Vector  >
inline
Oriented_side
side_of_plane_centered_sphere(const Point &a, const Vector& n,
/*defines the plane*/
	   const Point &p, const Point &q, const Point &r, const Point &t,
	   Homogeneous_tag)
{
    return side_of_plane_centered_sphereC3(a.x(), a.y(), a.z(),
                                           n.x(), n.y(), n.z(),
                                           p.x(), p.y(), p.z(),
                                           q.x(), q.y(), q.z(),
                                           r.x(), r.y(), r.z(),
                                           t.x(), t.y(), t.z());
}

template < class Point, class Vector >
inline
Oriented_side
side_of_plane_centered_sphere(const Point &a, const Vector& n,
/*defines the plane*/
	   const Point &p, const Point &q, const Point &r, Homogeneous_tag)
{
  return side_of_plane_centered_sphereC3(a.x(), a.y(), a.z(),
                                         n.x(), n.y(), n.z(),
                                         p.x(), p.y(), p.z(),
                                         q.x(), q.y(), q.z(),
                                         r.x(), r.y(), r.z());
}

template < class Point, class Vector >
inline
Point
plane_centered_circumcenter_3(const Point &a, const Vector& n,
			      /*defines the plane*/
			      const Point &p, const Point &q, const Point &r,
			      Homogeneous_tag)
{
  typename Point::R::RT x,y,z;
  plane_centered_circumcenterC3(a.x(), a.y(), a.z(),
                                n.x(), n.y(), n.z(),
                                p.x(), p.y(), p.z(),
                                q.x(), q.y(), q.z(),
                                r.x(), r.y(), r.z(),x,y,z);
  return Point(x,y,z);
}

template < class Point, class Vector >
inline
typename Point::R::Line_3
plane_intersected_bisector_3(const Point &a, const Vector& n,
			     /*defines the plane*/
			     const Point &p, const Point &q, Homogeneous_tag)
{
  typename Point::R::RT x1,y1,z1, x2,y2,z2;
  typedef typename Point::R::Line_3 Line;
  bisector_plane_intersectionC3(a.x(), a.y(), a.z(),
                                n.x(), n.y(), n.z(),
                                p.x(), p.y(), p.z(),
                                q.x(), q.y(), q.z(),x1,y1,z1,x2,y2,z2);

  return Line(Point(x1,y1,z1), Point(x2,y2,z2));
}

} //namespace CGAL

#endif // CGAL_VORONOI_INTERSECTION_2_TRAITS_3_H
