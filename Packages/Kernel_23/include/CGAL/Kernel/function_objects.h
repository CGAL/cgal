// ======================================================================
//
// Copyright (c) 1999,2002 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// release       : 
// release_date  : 
//
// file          : include/CGAL/Kernel/function_objects.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra, Sylvain Pion
//
// coordinator   : MPI, Saarbruecken
// ======================================================================

#ifndef CGAL_KERNEL_FUNCTION_OBJECTS_H
#define CGAL_KERNEL_FUNCTION_OBJECTS_H

#include <CGAL/functional_base.h>
#include <CGAL/Quotient.h>
#include <CGAL/Origin.h>

CGAL_BEGIN_NAMESPACE
namespace CGALi {

template <typename K>
class Angle_2
{
    typedef typename K::Point_2 Point_2;
public:
    typedef Angle            result_type;
    typedef Arity_tag< 3 >   Arity;

    Angle
    operator()(const Point_2& p, const Point_2& q, const Point_2& r) const
    { return angle(p, q, r); }
};

template <typename K>
class Angle_3
{
    typedef typename K::Point_3 Point_3;
public:
    typedef Angle            result_type;
    typedef Arity_tag< 3 >   Arity;

    Angle
    operator()(const Point_3& p, const Point_3& q, const Point_3& r) const
    { return angle(p, q, r); }
};

template <typename K>
class Are_ordered_along_line_2
{
    typedef typename K::Point_2 Point_2;
public:
    typedef bool             result_type;
    typedef Arity_tag< 3 >   Arity;

    bool
    operator()(const Point_2& p, const Point_2& q, const Point_2& r) const
    { return are_ordered_along_line(p,q,r); }
};

template <typename K>
class Are_ordered_along_line_3
{
    typedef typename K::Point_3 Point_3;
public:
    typedef bool             result_type;
    typedef Arity_tag< 3 >   Arity;

    bool
    operator()(const Point_3& p, const Point_3& q, const Point_3& r) const
    { return are_ordered_along_line(p,q,r); }
};

template <typename K>
class Are_strictly_ordered_along_line_2
{
    typedef typename K::Point_2 Point_2;
public:
    typedef bool             result_type;
    typedef Arity_tag< 3 >   Arity;

    bool
    operator()(const Point_2& p, const Point_2& q, const Point_2& r) const
    { return are_strictly_ordered_along_line(p,q,r); }
};

template <typename K>
class Are_strictly_ordered_along_line_3
{
    typedef typename K::Point_3 Point_3;
public:
    typedef bool             result_type;
    typedef Arity_tag< 3 >   Arity;

    bool
    operator()(const Point_3& p, const Point_3& q, const Point_3& r) const
    { return are_strictly_ordered_along_line(p, q, r); }
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
class Bounded_side_3
{
    typedef typename K::Point_3         Point_3;
    typedef typename K::Sphere_3        Sphere_3;
    typedef typename K::Tetrahedron_3   Tetrahedron_3;
    typedef typename K::Iso_cuboid_3    Iso_cuboid_3;
public:
    typedef Bounded_side     result_type;
    typedef Arity_tag< 2 >   Arity;

    Bounded_side
    operator()( const Sphere_3& s, const Point_3& p) const
    { return s.bounded_side(p); }

    Bounded_side
    operator()( const Tetrahedron_3& t, const Point_3& p) const
    { return t.bounded_side(p); }

    Bounded_side
    operator()( const Iso_cuboid_3& c, const Point_3& p) const
    { return c.bounded_side(p); }
};

template <typename K>
class Collinear_are_ordered_along_line_2
{
    typedef typename K::Point_2         Point_2;
public:
    typedef bool             result_type;
    typedef Arity_tag< 3 >   Arity;

    bool
    operator()(const Point_2& p, const Point_2& q, const Point_2& r) const
    { return collinear_are_ordered_along_line(p, q, r); }
};

template <typename K>
class Collinear_are_ordered_along_line_3
{
    typedef typename K::Point_3         Point_3;
public:
    typedef bool             result_type;
    typedef Arity_tag< 3 >   Arity;

    bool
    operator()(const Point_3& p, const Point_3& q, const Point_3& r) const
    { return collinear_are_ordered_along_line(p, q, r); }
};

template <typename K>
class Collinear_are_strictly_ordered_along_line_2
{
    typedef typename K::Point_2   Point_2;
public:
    typedef bool             result_type;
    typedef Arity_tag< 3 >   Arity;

    bool
    operator()(const Point_2& p, const Point_2& q, const Point_2& r) const
    { return collinear_are_strictly_ordered_along_line(p, q, r); }
};

template <typename K>
class Collinear_are_strictly_ordered_along_line_3
{
    typedef typename K::Point_3   Point_3;
public:
    typedef bool             result_type;
    typedef Arity_tag< 3 >   Arity;

    bool
    operator()(const Point_3& p, const Point_3& q, const Point_3& r) const
    { return collinear_are_strictly_ordered_along_line(p, q, r); }
};

template <typename K>
class Collinear_has_on_2
{
    typedef typename K::Point_2    Point_2;
    typedef typename K::Ray_2      Ray_2;
    typedef typename K::Segment_2  Segment_2;
public:
    typedef bool             result_type;
    typedef Arity_tag< 2 >   Arity;

    bool
    operator()( const Ray_2& r, const Point_2& p) const
    { return r.collinear_has_on(p); }

    bool
    operator()( const Segment_2& s, const Point_2& p) const
    { return s.collinear_has_on(p); }
};

template <typename K>
class Collinear_2
{
    typedef typename K::Point_2    Point_2;
public:
    typedef bool             result_type;
    typedef Arity_tag< 3 >   Arity;

    bool
    operator()(const Point_2& p, const Point_2& q, const Point_2& r) const
    { return collinear(p, q, r); }
};

template <typename K>
class Collinear_3
{
    typedef typename K::Point_3    Point_3;
public:
    typedef bool             result_type;
    typedef Arity_tag< 3 >   Arity;

    bool
    operator()(const Point_3& p, const Point_3& q, const Point_3& r) const
    { return collinear(p, q, r); }
};

template <typename K>
class Compare_angle_with_x_axis_2
{
    typedef typename K::Direction_2  Direction_2;
public:
    typedef Comparison_result        result_type;
    typedef Arity_tag< 2 >           Arity;

    Comparison_result
    operator()(const Direction_2& p, const Direction_2& q) const
    { return compare_angle_with_x_axis(p, q); }
};

template <typename K>
class Compare_distance_2
{
    typedef typename K::Point_2   Point_2;
public:
    typedef Comparison_result     result_type;
    typedef Arity_tag< 3 >        Arity;

    Comparison_result
    operator()(const Point_2& p, const Point_2& q, const Point_2& r) const
    { return compare_distance_to_point(p, q, r); }
};

template <typename K>
class Compare_distance_3
{
    typedef typename K::Point_3   Point_3;
public:
    typedef Comparison_result     result_type;
    typedef Arity_tag< 3 >        Arity;

    Comparison_result
    operator()(const Point_3& p, const Point_3& q, const Point_3& r) const
    { return compare_distance_to_point(p, q, r); }
};

template <typename K>
class Compare_slope_2
{
    typedef typename K::Line_2     Line_2;
    typedef typename K::Segment_2  Segment_2;
public:
    typedef Comparison_result      result_type;
    typedef Arity_tag< 2 >         Arity;

    Comparison_result
    operator()(const Line_2& l1, const Line_2& l2) const
    { return compare_slopes(l1, l2); }

    Comparison_result
    operator()(const Segment_2& s1, const Segment_2& s2) const
    { return compare_slopes(s1, s2); }
};

template <typename K>
class Compare_x_at_y_2
{
    typedef typename K::Point_2    Point_2;
    typedef typename K::Line_2     Line_2;
public:
    typedef Comparison_result      result_type;

    Comparison_result
    operator()( const Point_2& p, const Line_2& h) const
    { return compare_x_at_y(p, h); }

    Comparison_result
    operator()( const Point_2& p, const Line_2& h1, const Line_2& h2) const
    { return compare_x_at_y(p, h1, h2); }

    Comparison_result
    operator()( const Line_2& l1, const Line_2& l2, const Line_2& h) const
    { return compare_x_at_y(l1, l2, h); }

    Comparison_result
    operator()( const Line_2& l1, const Line_2& l2,
	        const Line_2& h1, const Line_2& h2) const
    { return compare_x_at_y(l1, l2, h1, h2); }
};

template <typename K>
class Compare_xyz_3
{
    typedef typename K::Point_3    Point_3;
public:
    typedef Comparison_result  result_type;
    typedef Arity_tag< 2 >     Arity;

    Comparison_result
    operator()( const Point_3& p, const Point_3& q) const
    { return compare_xyz(p, q); }
};

template <typename K>
class Compare_xy_2
{
    typedef typename K::Point_2    Point_2;
public:
    typedef Comparison_result  result_type;
    typedef Arity_tag< 2 >     Arity;

    Comparison_result
    operator()( const Point_2& p, const Point_2& q) const
    { return compare_xy(p, q); }
};

template <typename K>
class Compare_xy_3
{
    typedef typename K::Point_3    Point_3;
public:
    typedef Comparison_result  result_type;
    typedef Arity_tag< 2 >     Arity;

    Comparison_result
    operator()( const Point_3& p, const Point_3& q) const
    { return compare_xy(p, q); }
};

template <typename K>
class Compare_x_2
{
    typedef typename K::Point_2    Point_2;
    typedef typename K::Line_2     Line_2;
public:
    typedef Comparison_result      result_type;

    Comparison_result
    operator()( const Point_2& p, const Point_2& q) const
    { return compare_x(p, q); }

    Comparison_result
    operator()( const Point_2& p, const Line_2& l1, const Line_2& l2) const
    { return compare_x(p, l1, l2); }

    Comparison_result
    operator()( const Line_2& l, const Line_2& h1, const Line_2& h2) const
    { return compare_x(l, h1, h2); }

    Comparison_result
    operator()( const Line_2& l1, const Line_2& l2,
	        const Line_2& h1, const Line_2& h2) const
    { return compare_x(l1, l2, h1, h2); }
};

template <typename K>
class Compare_x_3
{
    typedef typename K::Point_3    Point_3;
public:
    typedef Comparison_result      result_type;
    typedef Arity_tag< 2 >         Arity;

    Comparison_result
    operator()( const Point_3& p, const Point_3& q) const
    { return compare_x(p, q); }
};

template <typename K>
class Compare_y_at_x_2
{
    typedef typename K::Point_2    Point_2;
    typedef typename K::Line_2     Line_2;
    typedef typename K::Segment_2  Segment_2;
public:
    typedef Comparison_result      result_type;

    Comparison_result
    operator()( const Point_2& p, const Line_2& h) const
    { return compare_y_at_x(p, h); }

    Comparison_result
    operator()( const Point_2& p, const Line_2& h1, const Line_2& h2) const
    { return compare_y_at_x(p, h1, h2); }

    Comparison_result
    operator()( const Line_2& l1, const Line_2& l2, const Line_2& h) const
    { return compare_y_at_x(l1, l2, h); }

    Comparison_result
    operator()( const Line_2& l1, const Line_2& l2,
	        const Line_2& h1, const Line_2& h2) const
    { return compare_y_at_x(l1, l2, h1, h2); }

    Comparison_result
    operator()( const Point_2& p, const Segment_2& s) const
    { return compare_y_at_x(p, s); }

    Comparison_result
    operator()( const Point_2& p,
	        const Segment_2& s1, const Segment_2& s2) const
    { return compare_y_at_x(p, s1, s2); }
};

template <typename K>
class Compare_y_2
{
    typedef typename K::Point_2   Point_2;
    typedef typename K::Line_2    Line_2;
public:
    typedef Comparison_result     result_type;

    Comparison_result
    operator()( const Point_2& p, const Point_2& q) const
    { return compare_y(p, q); }

    Comparison_result
    operator()( const Point_2& p, const Line_2& l1, const Line_2& l2) const
    { return compare_y(p, l1, l2); }

    Comparison_result
    operator()( const Line_2& l, const Line_2& h1, const Line_2& h2) const
    { return compare_y(l, h1, h2); }

    Comparison_result
    operator()( const Line_2& l1, const Line_2& l2,
	        const Line_2& h1, const Line_2& h2) const
    { return compare_y(l1, l2, h1, h2); }
};

template <typename K>
class Compare_y_3
{
    typedef typename K::Point_3   Point_3;
public:
    typedef Comparison_result     result_type;
    typedef Arity_tag< 2 >        Arity;

    Comparison_result
    operator()( const Point_3& p, const Point_3& q) const
    { return compare_y(p, q); }
};

template <typename K>
class Compare_z_3
{
    typedef typename K::Point_3   Point_3;
public:
    typedef Comparison_result     result_type;
    typedef Arity_tag< 2 >        Arity;

    Comparison_result
    operator()( const Point_3& p, const Point_3& q) const
    { return compare_z(p, q); }
};

template <typename K>
class Compute_area_2
{
    typedef typename K::FT                FT;
    typedef typename K::Iso_rectangle_2   Iso_rectangle_2;
    typedef typename K::Triangle_2        Triangle_2;
public:
    typedef FT               result_type;
    typedef Arity_tag< 1 >   Arity;

    FT
    operator()( const Iso_rectangle_2& r ) const
    { return r.area(); }

    FT
    operator()( const Triangle_2& t ) const
    { return t.area(); }
};

template <typename K>
class Compute_squared_area_3
{
    typedef typename K::FT                FT;
    typedef typename K::Triangle_3        Triangle_3;
public:
    typedef FT               result_type;
    typedef Arity_tag< 1 >   Arity;

    FT
    operator()( const Triangle_3& t ) const
    { return t.squared_area(); }
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
    { return squared_distance(t1, t2); }
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
    { return squared_distance(t1, t2); }
};

template <typename K>
class Compute_squared_length_2
{
    typedef typename K::FT          FT;
    typedef typename K::Segment_2   Segment_2;
  public:
    typedef FT               result_type;
    typedef Arity_tag< 1 >   Arity;

    FT
    operator()( const Segment_2& s) const
    { return s.squared_length(); }
};

template <typename K>
class Compute_squared_length_3
{
    typedef typename K::FT          FT;
    typedef typename K::Segment_3   Segment_3;
  public:
    typedef FT               result_type;
    typedef Arity_tag< 1 >   Arity;

    FT
    operator()( const Segment_3& s) const
    { return s.squared_length(); }
};

template <typename K>
class Compute_squared_radius_2
{
    typedef typename K::FT          FT;
    typedef typename K::Point_2     Point_2;
    typedef typename K::Circle_2    Circle_2;
public:
    typedef FT         result_type;

    FT
    operator()( const Circle_2& c) const
    { return c.squared_radius(); }

    FT
    operator()( const Point_2& p, const Point_2& q, const Point_2& r) const
    { return squared_radius(p, q, r); }
};

template <typename K>
class Compute_squared_radius_3
{
    typedef typename K::FT          FT;
    typedef typename K::Point_3     Point_3;
    typedef typename K::Sphere_3    Sphere_3;
public:
    typedef FT         result_type;

    FT
    operator()( const Sphere_3& s) const
    { return s.squared_radius(); }

    FT
    operator()( const Point_3& p, const Point_3& q, const Point_3& r) const
    { return squared_radius(p, q, r); }

    FT
    operator()( const Point_3& p, const Point_3& q,
	        const Point_3& r, const Point_3& s) const
    { return squared_radius(p, q, r, s); }
};

template <typename K>
class Compute_volume_3
{
    typedef typename K::FT             FT;
    typedef typename K::Tetrahedron_3  Tetrahedron_3;
    typedef typename K::Iso_cuboid_3   Iso_cuboid_3;
public:
    typedef FT               result_type;
    typedef Arity_tag< 1 >   Arity;

    FT
    operator()( const Tetrahedron_3& t ) const
    { return t.volume(); }

    FT
    operator()( const Iso_cuboid_3& c ) const
    { return c.volume(); }
};

template <typename K>
class Construct_base_vector_3
{
    typedef typename K::Vector_3   Vector_3;
    typedef typename K::Plane_3    Plane_3;
public:
     typedef Vector_3         result_type;
     typedef Arity_tag< 2 >   Arity;

     Vector_3
     operator()( const Plane_3& h, int index ) const
     {
       if (index == 1)
         return h.base1();
       else 
         return h.base2();
     }
};

template <typename K>
class Construct_bisector_2
{
    typedef typename K::Point_2 Point_2;
    typedef typename K::Line_2  Line_2;
public:
    typedef Line_2           result_type;
    typedef Arity_tag< 2 >   Arity;

    Line_2
    operator()(const Point_2& p, const Point_2& q) const
    { return bisector(p,q); }
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
class Construct_centroid_2
{
    typedef typename K::Point_2  Point_2;
public:
    typedef Point_2          result_type;

    Point_2
    operator()(const Point_2& p, const Point_2& q, const Point_2& r) const
    { return centroid(p, q, r); }

    Point_2
    operator()(const Point_2& p, const Point_2& q, 
               const Point_2& r, const Point_2& s) const
    { return centroid(p, q, r, s); }
};

template <typename K>
class Construct_centroid_3
{
    typedef typename K::Point_3  Point_3;
public:
    typedef Point_3          result_type;

    Point_3
    operator()(const Point_3& p, const Point_3& q, const Point_3& r) const
    { return centroid(p, q, r); }

    Point_3
    operator()(const Point_3& p, const Point_3& q, 
               const Point_3& r, const Point_3& s) const
    { return centroid(p, q, r, s); }
};

template <typename K>
class Construct_circle_2
{
    typedef typename K::FT          FT;
    typedef typename K::Point_2     Point_2;
    typedef typename K::Circle_2    Circle_2;
public:
    typedef Circle_2       result_type;

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
class Construct_circumcenter_2
{
    typedef typename K::Point_2  Point_2;
public:
    typedef Point_2          result_type;

    Point_2
    operator()(const Point_2& p, const Point_2& q, const Point_2& r) const
    { return circumcenter(p, q, r); }
};

template <typename K>
class Construct_circumcenter_3
{
    typedef typename K::Point_3  Point_3;
public:
    typedef Point_3          result_type;

    Point_3
    operator()(const Point_3& p, const Point_3& q, const Point_3& r) const
    { return circumcenter(p, q, r); }

    Point_3
    operator()(const Point_3& p, const Point_3& q,
	       const Point_3& r, const Point_3& s) const
    { return circumcenter(p, q, r, s); }
};

template <typename K>
class Construct_cross_product_vector_3
{
    typedef typename K::Vector_3  Vector_3;
public:
    typedef Vector_3         result_type;
    typedef Arity_tag< 2 >   Arity;

    Vector_3
    operator()(const Vector_3& v, const Vector_3& w) const
    { return cross_product(v, w); }
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
    typedef Iso_cuboid_3   result_type;

    Iso_cuboid_3
    operator()() const
    { return Iso_cuboid_3(); }

    Iso_cuboid_3
    operator()(const Point_3& p, const Point_3& q) const
    { return Iso_cuboid_3(p, q); }
};

template <typename K>
class Construct_iso_rectangle_2
{
    typedef typename K::Point_2          Point_2;
    typedef typename K::Iso_rectangle_2  Iso_rectangle_2;
public:
    typedef Iso_rectangle_2   result_type;

    Iso_rectangle_2
    operator()() const
    { return Iso_rectangle_2(); }

    Iso_rectangle_2
    operator()(const Point_2& p, const Point_2& q) const
    { return Iso_rectangle_2(p, q); }
};

template <typename K>
class Construct_lifted_point_3
{
    typedef typename K::Point_2          Point_2;
    typedef typename K::Point_3          Point_3;
    typedef typename K::Plane_3          Plane_3;
public:
    typedef Point_3          result_type;
    typedef Arity_tag< 2 >   Arity;

    Point_3
    operator()(const Plane_3& h, const Point_2& p) const
    {  return h.to_3d(p); }
};

template <typename K>
class Construct_line_2
{
    typedef typename K::RT           RT;
    typedef typename K::Point_2      Point_2;
    typedef typename K::Direction_2  Direction_2;
    typedef typename K::Segment_2    Segment_2;
    typedef typename K::Ray_2        Ray_2;
    typedef typename K::Line_2       Line_2;
public:
    typedef Line_2  result_type;

    Line_2
    operator()() const
    { return Line_2(); }

#ifndef CGAL_NO_DEPRECATED_CODE
    Line_2
    operator()(const RT& a, const RT& b, const RT& c) const
    { return Line_2(a, b, c); }
#endif // CGAL_NO_DEPRECATED_CODE

    Line_2
    operator()(const Point_2& p, const Point_2& q) const
    { return Line_2(p, q); }

    Line_2
    operator()(const Point_2& p, const Direction_2& d) const
    { return Line_2(p, d); }

    Line_2
    operator()(const Segment_2& s) const
    { return Line_2(s); }

    Line_2
    operator()(const Ray_2& r) const
    { return Line_2(r); }
};

template <typename K>
class Construct_line_3
{
    typedef typename K::Point_3      Point_3;
    typedef typename K::Direction_3  Direction_3;
    typedef typename K::Segment_3    Segment_3;
    typedef typename K::Ray_3        Ray_3;
    typedef typename K::Line_3       Line_3;
public:
    typedef Line_3  result_type;

    Line_3
    operator()() const
    { return Line_3(); }

    Line_3
    operator()(const Point_3& p, const Point_3& q) const
    { return Line_3(p, q); }

    Line_3
    operator()(const Point_3& p, const Direction_3& d) const
    { return Line_3(p, d); }

    Line_3
    operator()(const Segment_3& s) const
    { return Line_3(s); }

    Line_3
    operator()(const Ray_3& r) const
    { return Line_3(r); }
};

template <typename K>
class Construct_midpoint_2
{
    typedef typename K::Point_2   Point_2;
public:
    typedef Point_2          result_type;
    typedef Arity_tag< 2 >   Arity;

    Point_2
    operator()(const Point_2& p, const Point_2& q) const
    { return midpoint(p, q); }
};

template <typename K>
class Construct_midpoint_3
{
    typedef typename K::Point_3   Point_3;
public:
    typedef Point_3          result_type;
    typedef Arity_tag< 2 >   Arity;

    Point_3
    operator()(const Point_3& p, const Point_3& q) const
    { return midpoint(p, q); }
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
class Construct_opposite_vector_2
{
    typedef typename K::Vector_2    Vector_2;
public:
    typedef Vector_2         result_type;
    typedef Arity_tag< 1 >   Arity;

    Vector_2
    operator()( const Vector_2& v) const
    { return -v; }
};

template <typename K>
class Construct_opposite_vector_3
{
    typedef typename K::Vector_3    Vector_3;
public:
    typedef Vector_3         result_type;
    typedef Arity_tag< 1 >   Arity;

    Vector_3
    operator()( const Vector_3& v) const
    { return -v; }
};

template <typename K>
class Construct_orthogonal_vector_3
{
    typedef typename K::Vector_3    Vector_3;
    typedef typename K::Plane_3     Plane_3;
public:
    typedef Vector_3         result_type;
    typedef Arity_tag< 1 >   Arity;

    Vector_3
    operator()( const Plane_3& p ) const
    { return p.orthogonal_vector(); }
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
    typedef typename K::Direction_3  Direction_3;
    typedef typename K::Line_3       Line_3;
    typedef typename K::Ray_3        Ray_3;
    typedef typename K::Segment_3    Segment_3;
    typedef typename K::Plane_3      Plane_3;
public:
    typedef Plane_3    result_type;

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
    typedef Point_2     result_type;

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
    typedef Point_3     result_type;

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
    typedef Point_2    result_type;

    Point_2
    operator()() const
    { return Point_2(); }

    Point_2
    operator()(Origin o) const
    { return Point_2(o); }

#ifndef CGAL_NO_DEPRECATED_CODE
    Point_2
    operator()(const RT& x, const RT& y) const
    { return Point_2(x, y); }

    Point_2
    operator()(const RT& x, const RT& y, const RT& w) const
    { return Point_2(x, y, w); }
#endif // CGAL_NO_DEPRECATED_CODE
};

template <typename K>
class Construct_point_3
{
    typedef typename K::RT         RT;
    typedef typename K::Point_3    Point_3;
public:
    typedef Point_3    result_type;

    Point_3
    operator()() const
    { return Point_3(); }

    Point_3
    operator()(Origin o) const
    { return Point_3(o); }

#ifndef CGAL_NO_DEPRECATED_CODE
    Point_3
    operator()(const RT& x, const RT& y, const RT& z) const
    { return Point_3(x, y, z); }

    Point_3
    operator()(const RT& x, const RT& y, const RT& z, const RT& w) const
    { return Point_3(x, y, z, w); }
#endif // CGAL_NO_DEPRECATED_CODE
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
class Construct_projected_point_3
{
    typedef typename K::Point_3    Point_3;
    typedef typename K::Plane_3    Plane_3;
    typedef typename K::Line_3     Line_3;
public:
    typedef Point_3          result_type;
    typedef Arity_tag< 2 >   Arity;

    Point_3
    operator()( const Line_3& l, const Point_3& p ) const
    { return l.projection(p); }

    Point_3
    operator()( const Plane_3& h, const Point_3& p ) const
    { return h.projection(p); }
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
    typedef typename K::Direction_2  Direction_2;
    typedef typename K::Ray_2        Ray_2;
public:
     typedef Ray_2            result_type;

     Ray_2
     operator()() const
     {  return Ray_2(); }

     Ray_2
     operator()(const Point_2& p, const Point_2& q) const
     {  return Ray_2(p, q); }

     Ray_2
     operator()(const Point_2& p, const Direction_2& d) const
     {  return Ray_2(p, d); }
};

template <typename K>
class Construct_ray_3
{
    typedef typename K::Point_3      Point_3;
    typedef typename K::Direction_3  Direction_3;
    typedef typename K::Ray_3        Ray_3;
public:
     typedef Ray_3            result_type;

     Ray_3
     operator()() const
     {  return Ray_3(); }

     Ray_3
     operator()(const Point_3& p, const Point_3& q) const
     {  return Ray_3(p, q); }

     Ray_3
     operator()(const Point_3& p, const Direction_3& d) const
     {  return Ray_3(p, d); }
};

template <typename K>
class Construct_scaled_vector_2
{
    typedef typename K::RT         RT;
    typedef typename K::Vector_2   Vector_2;
public:
    typedef Vector_2         result_type;
    typedef Arity_tag< 2 >   Arity;

    Vector_2
    operator()( const Vector_2& v, const RT& scale) const
    {  return v * scale; }

    Vector_2
    operator()( const Vector_2& v, const Quotient<RT>& scale) const
    {  return v * scale; }
};

template <typename K>
class Construct_scaled_vector_3
{
    typedef typename K::RT         RT;
    typedef typename K::Vector_3   Vector_3;
public:
    typedef Vector_3         result_type;
    typedef Arity_tag< 2 >   Arity;

    Vector_3
    operator()( const Vector_3& v, const RT& scale) const
    {  return v * scale; }

    Vector_3
    operator()( const Vector_3& v, const Quotient<RT>& scale) const
    {  return v * scale; }
};

template <typename K>
class Construct_segment_2
{
    typedef typename K::Segment_2  Segment_2;
    typedef typename K::Point_2    Point_2;
public:
    typedef Segment_2        result_type;

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

    Tetrahedron_3
    operator()() const
    { return Tetrahedron_3(); }

    Tetrahedron_3
    operator()( const Point_3& p, const Point_3& q,
	        const Point_3& r, const Point_3& s) const
    { return Tetrahedron_3(p, q, r, s); }
};

template <typename K>
class Construct_translated_point_2
{
    typedef typename K::Point_2   Point_2;
    typedef typename K::Vector_2  Vector_2;
public:
    typedef Point_2          result_type;
    typedef Arity_tag< 2 >   Arity;

    Point_2
    operator()( const Point_2& p, const Vector_2& v) const
    {  return p + v; }
};

template <typename K>
class Construct_translated_point_3
{
    typedef typename K::Point_3   Point_3;
    typedef typename K::Vector_3  Vector_3;
public:
    typedef Point_3          result_type;
    typedef Arity_tag< 2 >   Arity;

    Point_3
    operator()( const Point_3& p, const Vector_3& v) const
    {  return p + v; }
};

template <typename K>
class Construct_triangle_2
{
    typedef typename K::Triangle_2   Triangle_2;
    typedef typename K::Point_2      Point_2;
public:
    typedef Triangle_2   result_type;

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
    typedef Triangle_3   result_type;

    Triangle_3
    operator()() const
    { return Triangle_3(); }

    Triangle_3
    operator()( const Point_3& p, const Point_3& q, const Point_3& r) const
    { return Triangle_3(p, q, r); }
};

template <typename K>
class Construct_vector_2
{
    typedef typename K::RT           RT;
    typedef typename K::Vector_2     Vector_2;
    typedef typename K::Point_2      Point_2;
public:
    typedef Vector_2   result_type;

    Vector_2
    operator()() const
    { return Vector_2(); }

    Vector_2
    operator()( const Point_2& p, const Point_2& q) const
    { return Vector_2(p, q); }

    Vector_2
    operator()( Null_vector n) const
    { return Vector_2(n); }

#ifndef CGAL_NO_DEPRECATED_CODE
    Vector_2
    operator()( const RT& x, const RT& y) const
    { return Vector_2(x, y); }

    Vector_2
    operator()( const RT& x, const RT& y, const RT& w) const
    { return Vector_2(x, y, w); }
#endif // CGAL_NO_DEPRECATED_CODE
};

template <typename K>
class Construct_vector_3
{
    typedef typename K::RT           RT;
    typedef typename K::Vector_3     Vector_3;
    typedef typename K::Point_3      Point_3;
public:
    typedef Vector_3   result_type;

    Vector_3
    operator()() const
    { return Vector_3(); }

    Vector_3
    operator()( const Point_3& p, const Point_3& q) const
    { return Vector_3(p, q); }

    Vector_3
    operator()( Null_vector n) const
    { return Vector_3(n); }

#ifndef CGAL_NO_DEPRECATED_CODE
    Vector_3
    operator()( const RT& x, const RT& y, const RT& z) const
    { return Vector_3(x, y, z); }

    Vector_3
    operator()( const RT& x, const RT& y, const RT& z, const RT& w) const
    { return Vector_3(x, y, z, w); }
#endif // CGAL_NO_DEPRECATED_CODE
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

    Point_2
    operator()( const Segment_2& s, int i) const
    { return s.vertex(i); }

    Point_2
    operator()( const Triangle_2& t, int i) const
    { return t.vertex(i); }

    Point_2
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

    Point_3
    operator()( const Segment_3& s, int i) const
    { return s.vertex(i); }

    Point_3
    operator()( const Triangle_3& t, int i) const
    { return t.vertex(i); }

    Point_3
    operator()( const Iso_cuboid_3& r, int i) const
    { return r.vertex(i); }

    Point_3
    operator()( const Tetrahedron_3& t, int i) const
    { return t.vertex(i); }
};

template <typename K>
class Coplanar_orientation_3
{
    typedef typename K::Point_3   Point_3;
public:
    typedef Orientation  result_type;

    Orientation
    operator()(const Point_3& p, const Point_3& q, const Point_3& r) const
    { return coplanar_orientation(p, q, r); }

    Orientation
    operator()( const Point_3& p, const Point_3& q,
	        const Point_3& r, const Point_3& t) const
    { return coplanar_orientation(p, q, r, t); }
};

template <typename K>
class Coplanar_side_of_bounded_circle_3
{
    typedef typename K::Point_3   Point_3;
public:
    typedef Bounded_side     result_type;
    typedef Arity_tag< 4 >   Arity;

    Bounded_side
    operator()( const Point_3& p, const Point_3& q,
	        const Point_3& r, const Point_3& t) const
    { return coplanar_side_of_bounded_circle(p, q, r, t); }
};

template <typename K>
class Coplanar_3
{
    typedef typename K::Point_3   Point_3;
public:
    typedef bool             result_type;
    typedef Arity_tag< 4 >   Arity;

    bool
    operator()( const Point_3& p, const Point_3& q,
	        const Point_3& r, const Point_3& s) const
    { return coplanar(p, q, r, s); }
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
    { return p.counterclockwise_in_between(q, r); }
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
    { return do_intersect(t1, t2); }
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
    { return do_intersect(t1, t2); }
};

template <typename K>
class Equal_xy_3
{
    typedef typename K::Point_3    Point_3;
public:
    typedef bool             result_type;
    typedef Arity_tag< 2 >   Arity;

    bool
    operator()( const Point_3& p, const Point_3& q) const
    { return equal_xy(p, q); }
};

template <typename K>
class Equal_x_2
{
    typedef typename K::Point_2    Point_2;
public:
    typedef bool             result_type;
    typedef Arity_tag< 2 >   Arity;

    bool
    operator()( const Point_2& p, const Point_2& q) const
    { return x_equal(p, q); }
};

template <typename K>
class Equal_x_3
{
    typedef typename K::Point_3    Point_3;
public:
    typedef bool             result_type;
    typedef Arity_tag< 2 >   Arity;

    bool
    operator()( const Point_3& p, const Point_3& q) const
    { return x_equal(p, q); }
};

template <typename K>
class Equal_y_2
{
    typedef typename K::Point_2    Point_2;
public:
    typedef bool             result_type;
    typedef Arity_tag< 2 >   Arity;

    bool
    operator()( const Point_2& p, const Point_2& q) const
    { return y_equal(p, q); }
};

template <typename K>
class Equal_y_3
{
    typedef typename K::Point_3    Point_3;
public:
    typedef bool             result_type;
    typedef Arity_tag< 2 >   Arity;

    bool
    operator()( const Point_3& p, const Point_3& q) const
    { return y_equal(p, q); }
};

template <typename K>
class Equal_z_3
{
    typedef typename K::Point_3    Point_3;
public:
    typedef bool             result_type;
    typedef Arity_tag< 2 >   Arity;

    bool
    operator()( const Point_3& p, const Point_3& q) const
    { return z_equal(p, q); }
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

#ifndef CGAL_NO_DEPRECATED_CODE
    bool
    operator()( const Plane_3& pl, const Point_3& p) const
    { return pl.has_on_boundary(p); }
#endif // CGAL_NO_DEPRECATED_CODE
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
class Has_on_3
{
    typedef typename K::Point_3          Point_3;
    typedef typename K::Line_3           Line_3;
    typedef typename K::Ray_3            Ray_3;
    typedef typename K::Segment_3        Segment_3;
    typedef typename K::Plane_3          Plane_3;
    typedef typename K::Triangle_3       Triangle_3;
public:
    typedef bool             result_type;
    typedef Arity_tag< 2 >   Arity;

    bool
    operator()( const Line_3& l, const Point_3& p) const
    { return l.has_on(p); }

    bool
    operator()( const Ray_3& r, const Point_3& p) const
    { return r.has_on(p); }

    bool
    operator()( const Segment_3& s, const Point_3& p) const
    { return s.has_on(p); }

    bool
    operator()( const Plane_3& pl, const Point_3& p) const
    { return pl.has_on(p); }

    bool
    operator()( const Triangle_3& t, const Point_3& p) const
    { return t.has_on(p); }
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
    { return intersection(t1, t2); }
};

template <typename K>
class Intersect_3
{
    typedef typename K::Object_3    Object_3;
public:
    typedef Object_3         result_type;
    typedef Arity_tag< 2 >   Arity;

    // n possibilities, so I keep the template.
    template <class T1, class T2>
    Object_3
    operator()(const T1& t1, const T2& t2) const
    { return intersection(t1, t2); }
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

// We qualify the call to the global function with the CGAL namespace because
// otherwise ambigiuties arise when Point_2 is a leda_point (or leda_rat_point)
// and using LEDA versions earlier than 4.4 (where the leda namespace was
// introduced).  Compilers cannot decide between the left_turn from LEDA in
// global namespace and the left_turn in CGAL namespace defined in 
// (rat_)leda_for_CGAL_2.h.
template <typename K>
class Left_turn_2
{
    typedef typename K::Point_2   Point_2;
public:
    typedef bool             result_type;
    typedef Arity_tag< 3 >   Arity;

    bool
    operator()(const Point_2& p, const Point_2& q, const Point_2& r) const
    { return CGAL::left_turn(p, q, r); }
};

template <typename K>
class Less_distance_to_point_2
{
    typedef typename K::Point_2   Point_2;
public:
    typedef bool             result_type;
    typedef Arity_tag< 3 >   Arity;

    bool
    operator()(const Point_2& p, const Point_2& q, const Point_2& r) const
    { return has_smaller_distance_to_point(p, q, r); }
};

template <typename K>
class Less_distance_to_point_3
{
    typedef typename K::Point_3   Point_3;
public:
    typedef bool             result_type;
    typedef Arity_tag< 3 >   Arity;

    bool
    operator()(const Point_3& p, const Point_3& q, const Point_3& r) const
    { return has_smaller_distance_to_point(p, q, r); }
};

template <typename K>
class Less_rotate_ccw_2
{
    typedef typename K::Point_2   Point_2;
public:
    typedef bool             result_type;
    typedef Arity_tag< 3 >   Arity;

    bool
    operator()(const Point_2& r, const Point_2& p, const Point_2& q) const
    {
          Orientation ori = orientation(r, p, q);
          if ( ori == LEFT_TURN )
              return true;
          else if ( ori == RIGHT_TURN )
              return false;
          else
          {
              if (p == r) return false;
              if (q == r) return true;
              if (p == q) return false;
              return  collinear_are_ordered_along_line( r, q, p);
          }
    }
};

template <typename K>
class Less_signed_distance_to_line_2
{
    typedef typename K::Point_2   Point_2;
public:
    typedef bool             result_type;
    typedef Arity_tag< 4 >   Arity;

    bool
    operator()(const Point_2& a, const Point_2& b,
               const Point_2& c, const Point_2& d) const
    {
        Comparison_result res = compare_signed_distance_to_line(a, b, c, d);
          if ( res == LARGER )
              return false;
          else if ( res == SMALLER )
              return true;
          else
              return lexicographically_xy_smaller( c, d );
    }
};

template <typename K>
class Less_signed_distance_to_plane_3
{
    typedef typename K::Point_3 Point_3;
    typedef typename K::Plane_3 Plane_3;
public:
    typedef bool             result_type;
    typedef Arity_tag< 3 >   Arity;

    bool
    operator()( const Plane_3& p, const Point_3& q, const Point_3& r) const
    { return has_smaller_signed_distance_to_plane(p, q, r); }
};

template <typename K>
class Less_xyz_3
{
    typedef typename K::Point_3 Point_3;
public:
    typedef bool             result_type;
    typedef Arity_tag< 2 >   Arity;

    bool
    operator()( const Point_3& p, const Point_3& q) const
    { return lexicographically_xyz_smaller(p, q); }
};

template <typename K>
class Less_xy_2
{
    typedef typename K::Point_2 Point_2;
public:
    typedef bool             result_type;
    typedef Arity_tag< 2 >   Arity;

    bool
    operator()( const Point_2& p, const Point_2& q) const
    { return lexicographically_xy_smaller( p, q); }
};

template <typename K>
class Less_xy_3
{
    typedef typename K::Point_3 Point_3;
public:
    typedef bool             result_type;
    typedef Arity_tag< 2 >   Arity;

    bool
    operator()( const Point_3& p, const Point_3& q) const
    { return lexicographically_xy_smaller( p, q); }
};

template <typename K>
class Less_x_2
{
    typedef typename K::Point_2 Point_2;
public:
    typedef bool             result_type;
    typedef Arity_tag< 2 >   Arity;

    bool
    operator()( const Point_2& p, const Point_2& q) const
    { return less_x(p, q); }
};

template <typename K>
class Less_x_3
{
    typedef typename K::Point_3 Point_3;
public:
    typedef bool             result_type;
    typedef Arity_tag< 2 >   Arity;

    bool
    operator()( const Point_3& p, const Point_3& q) const
    { return less_x(p, q); }
};

template <typename K>
class Less_yx_2
{
    typedef typename K::Point_2 Point_2;
public:
    typedef bool             result_type;
    typedef Arity_tag< 2 >   Arity;

    bool
    operator()( const Point_2& p, const Point_2& q) const
    { return lexicographically_yx_smaller(p, q); }
};

template <typename K>
class Less_y_2
{
    typedef typename K::Point_2 Point_2;
public:
    typedef bool             result_type;
    typedef Arity_tag< 2 >   Arity;

    bool
    operator()( const Point_2& p, const Point_2& q) const
    { return less_y(p, q); }
};

template <typename K>
class Less_y_3
{
    typedef typename K::Point_3 Point_3;
public:
    typedef bool             result_type;
    typedef Arity_tag< 2 >   Arity;

    bool
    operator()( const Point_3& p, const Point_3& q) const
    { return less_y(p, q); }
};

template <typename K>
class Less_z_3
{
    typedef typename K::Point_3 Point_3;
public:
    typedef bool             result_type;
    typedef Arity_tag< 2 >   Arity;

    bool
    operator()( const Point_3& p, const Point_3& q) const
    { return less_z(p, q); }
};

template <typename K>
class Orientation_2
{
    typedef typename K::Point_2 Point_2;
public:
    typedef Orientation      result_type;
    typedef Arity_tag< 3 >   Arity;

    Orientation
    operator()(const Point_2& p, const Point_2& q, const Point_2& r) const
    { return orientation(p, q, r); }
};

template <typename K>
class Orientation_3
{
    typedef typename K::Point_3 Point_3;
public:
    typedef Orientation      result_type;
    typedef Arity_tag< 4 >   Arity;

    Orientation
    operator()( const Point_3& p, const Point_3& q,
	        const Point_3& r, const Point_3& s) const
    { return orientation(p, q, r, s); }
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

template <typename K>
class Side_of_bounded_circle_2
{
    typedef typename K::Point_2        Point_2;
public:
    typedef Bounded_side   result_type;

    Bounded_side
    operator()( const Point_2& p, const Point_2& q, const Point_2& t) const
    { return side_of_bounded_circle(p, q, t); }

    Bounded_side
    operator()( const Point_2& p, const Point_2& q,
	        const Point_2& r, const Point_2& t) const
    { return side_of_bounded_circle(p, q, r, t); }
};

template <typename K>
class Side_of_bounded_sphere_3
{
    typedef typename K::Point_3        Point_3;
public:
    typedef Bounded_side   result_type;

    Bounded_side
    operator()( const Point_3& p, const Point_3& q, const Point_3& t) const
    { return side_of_bounded_sphere(p, q, t); }

    Bounded_side
    operator()( const Point_3& p, const Point_3& q,
	        const Point_3& r, const Point_3& t) const
    { return side_of_bounded_sphere(p, q, r, t); }

    Bounded_side
    operator()( const Point_3& p, const Point_3& q, const Point_3& r,
	        const Point_3& s, const Point_3& t) const
    { return side_of_bounded_sphere(p, q, r, s, t); }
};

template <typename K>
class Side_of_oriented_circle_2
{
    typedef typename K::Point_2        Point_2;
public:
    typedef Oriented_side    result_type;
    typedef Arity_tag< 4 >   Arity;

    Oriented_side
    operator()( const Point_2& p, const Point_2& q,
	        const Point_2& r, const Point_2& t) const
    { return side_of_oriented_circle(p, q, r, t); }
};

template <typename K>
class Side_of_oriented_sphere_3
{
    typedef typename K::Point_3        Point_3;
public:
    typedef Oriented_side    result_type;
    typedef Arity_tag< 5 >   Arity;

    Oriented_side
    operator()( const Point_3& p, const Point_3& q, const Point_3& r,
	        const Point_3& s, const Point_3& t) const
    { return side_of_oriented_sphere(p, q, r, s, t); }
};

#ifndef CGAL_NO_DEPRECATED_CODE

template <typename K>
class Equal_xy_2
{
    typedef typename K::Point_2    Point_2;
public:
    typedef bool             result_type;
    typedef Arity_tag< 2 >   Arity;

    bool
    operator()( const Point_2& p, const Point_2& q) const
    { return equal_xy(p, q); }
};

template <typename K>
class Leftturn_2
{
    typedef typename K::Point_2   Point_2;
public:
    typedef bool             result_type;
    typedef Arity_tag< 3 >   Arity;

    bool
    operator()(const Point_2& p, const Point_2& q, const Point_2& r) const
    { return CGAL::left_turn(p, q, r); }
};

template <typename K>
class Equal_xyz_3
{
    typedef typename K::Point_3  Point_3;
public:
    typedef bool             result_type;
    typedef Arity_tag< 2 >   Arity;

    bool
    operator()( const Point_3& p, const Point_3& q) const
    { return equal_xyz(p, q); }
};

template <typename K>
class Construct_direction_of_line_2
{
    typedef typename K::Direction_2   Direction_2;
    typedef typename K::Line_2        Line_2;
public:
    typedef Direction_2      result_type;
    typedef Arity_tag< 1 >   Arity;

    Direction_2
    operator()( const Line_2& l) const
    { return l.direction(); }
};

template <typename K>
class Construct_direction_of_line_3
{
    typedef typename K::Direction_3   Direction_3;
    typedef typename K::Line_3        Line_3;
public:
    typedef Direction_3      result_type;
    typedef Arity_tag< 1 >   Arity;

    Direction_3
    operator()( const Line_3& l) const
    { return l.direction(); }
};

template <typename K>
class Construct_direction_of_ray_2
{
    typedef typename K::Direction_2   Direction_2;
    typedef typename K::Ray_2         Ray_2;
public:
    typedef Direction_2      result_type;
    typedef Arity_tag< 1 >   Arity;

    Direction_2
    operator()( const Ray_2& r) const
    { return r.direction(); }
};

template <typename K>
class Construct_direction_of_ray_3
{
    typedef typename K::Direction_3   Direction_3;
    typedef typename K::Ray_3         Ray_3;
public:
    typedef Direction_3      result_type;
    typedef Arity_tag< 1 >   Arity;

    Direction_3
    operator()( const Ray_3& r) const
    { return r.direction(); }
};

template <typename K>
class Construct_max_point_2
{
    typedef typename K::Point_2  Point_2;
public:
    typedef Point_2          result_type;
    typedef Arity_tag< 1 >   Arity;

    template <class Cls>
    Point_2
    operator()( const Cls& c) const
    { return c.max(); }
};

template <typename K>
class Construct_max_point_3
{
    typedef typename K::Point_3  Point_3;
public:
    typedef Point_3          result_type;
    typedef Arity_tag< 1 >   Arity;

    template <class Cls>
    Point_3
    operator()( const Cls& c) const
    { return c.max(); }
};

template <typename K>
class Construct_min_point_2
{
    typedef typename K::Point_2  Point_2;
public:
    typedef Point_2          result_type;
    typedef Arity_tag< 1 >   Arity;

    template <class Cls>
    Point_2
    operator()( const Cls& c) const
    { return c.min(); }
};

template <typename K>
class Construct_min_point_3
{
    typedef typename K::Point_3  Point_3;
public:
    typedef Point_3          result_type;
    typedef Arity_tag< 1 >   Arity;

    template <class Cls>
    Point_3
    operator()( const Cls& c) const
    { return c.min(); }
};

template <typename K>
class Construct_source_point_2
{
    typedef typename K::Point_2  Point_2;
public:
    typedef Point_2          result_type;
    typedef Arity_tag< 1 >   Arity;

    template <class Cls>
    Point_2
    operator()( const Cls& c) const
    { return c.source(); }
};

template <typename K>
class Construct_source_point_3
{
    typedef typename K::Point_3  Point_3;
public:
    typedef Point_3          result_type;
    typedef Arity_tag< 1 >   Arity;

    template <class Cls>
    Point_3
    operator()( const Cls& c) const
    { return c.source(); }
};

template <typename K>
class Construct_target_point_2
{
    typedef typename K::Point_2  Point_2;
public:
    typedef Point_2          result_type;
    typedef Arity_tag< 1 >   Arity;

    template <class Cls>
    Point_2
    operator()( const Cls& c) const
    { return c.target(); }
};

template <typename K>
class Construct_target_point_3
{
    typedef typename K::Point_3  Point_3;
public:
    typedef Point_3          result_type;
    typedef Arity_tag< 1 >   Arity;

    template <class Cls>
    Point_3
    operator()( const Cls& c) const
    { return c.target(); }
};

template <typename K>
class Construct_second_point_on_2
{
    typedef typename K::Point_2  Point_2;
public:
    typedef Point_2          result_type;
    typedef Arity_tag< 1 >   Arity;

    template <class Cls>
    Point_2
    operator()( const Cls& c) const
    { return c.second_point(); }
};

template <typename K>
class Construct_second_point_on_3
{
    typedef typename K::Point_3  Point_3;
public:
    typedef Point_3          result_type;
    typedef Arity_tag< 1 >   Arity;

    template <class Cls>
    Point_3
    operator()( const Cls& c) const
    { return c.second_point(); }
};

template <typename K>
class Transform_2
{
public:
    template <class Transformation, class ArgumentType>
    ArgumentType
    operator()( const ArgumentType& a, const Transformation& t) const
    { return a.transform(t); }
};

template <typename K>
class Transform_3
{
public:
    template <class Transformation, class ArgumentType>
    ArgumentType
    operator()( const ArgumentType& a, const Transformation& t) const
    { return a.transform(t); }
};

// This one is 100% unused.
template <class ReturnType>
class Call_x_at_y_to_get
{
  public:
    typedef ReturnType     result_type;
    typedef Arity_tag< 2 >   Arity;

    template <class Cls>
    ReturnType
    operator()( const Cls& c, const ReturnType& x) const
    { return c.x_at_y(x); }
};

template <typename K>
class Compute_y_at_x_2
{
    typedef typename K::FT        FT;
    typedef typename K::Line_2    Line_2;
public:
    typedef FT               result_type;
    typedef Arity_tag< 2 >   Arity;

    FT
    operator()( const Line_2& l, const FT& x) const
    { return l.y_at_x(x); }
};

template <typename K>
class Construct_aff_transformation_2
{
    typedef typename K::Aff_transformation_2 ToBeConstructed;
public:
    typedef ToBeConstructed  result_type;

    ToBeConstructed
    operator()() const
    { return ToBeConstructed(); }

    template <class A1> 
    ToBeConstructed
    operator()( const A1& a1) const
    { return ToBeConstructed(a1); }

    template <class A1, class A2> 
    ToBeConstructed
    operator()( const A1& a1, const A2& a2) const
    { return ToBeConstructed(a1,a2); }

    template <class A1, class A2, class A3> 
    ToBeConstructed
    operator()( const A1& a1, const A2& a2, const A3& a3) const
    { return ToBeConstructed(a1,a2,a3); }

    template <class A1, class A2, class A3, class A4> 
    ToBeConstructed
    operator()( const A1& a1, const A2& a2, const A3& a3, const A4& a4) const
    { return ToBeConstructed(a1,a2,a3,a4); }

    template <class A1, class A2, class A3, class A4, class A5> 
    ToBeConstructed
    operator()( const A1& a1, const A2& a2, const A3& a3, const A4& a4,
	    const A5& a5) const
    { return ToBeConstructed(a1,a2,a3,a4,a5); }

    template <class A> 
    ToBeConstructed
    operator()( const A& a1, const A& a2, const A& a3,
                const A& a4, const A& a5, const A& a6 ) const
    { return ToBeConstructed(a1,a2,a3,a4,a5,a6); }

    template <class A> 
    ToBeConstructed
    operator()( const A& a1, const A& a2, const A& a3,
                const A& a4, const A& a5, const A& a6,
                const A& a7 ) const
    { return ToBeConstructed(a1,a2,a3,a4,a5,a6,a7); }

    template <class A> 
    ToBeConstructed
    operator()( const A& a1, const A& a2, const A& a3,
                const A& a4, const A& a5, const A& a6,
                const A& a7, const A& a8, const A& a9) const
    { return ToBeConstructed(a1,a2,a3,a4,a5,a6,a7,a8,a9); }

    template <class A> 
    ToBeConstructed
    operator()( const A& a1, const A& a2, const A& a3,
                const A& a4, const A& a5, const A& a6,
                const A& a7, const A& a8, const A& a9,
                const A& a10) const
    { return ToBeConstructed(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10); }

    template <class A> 
    ToBeConstructed
    operator()( const A& a1, const A& a2, const A& a3,
                const A& a4, const A& a5, const A& a6,
                const A& a7, const A& a8, const A& a9,
                const A& a10,const A& a11,const A& a12) const
    { return ToBeConstructed(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12); }

    template <class A> 
    ToBeConstructed
    operator()( const A& a1, const A& a2, const A& a3,
                const A& a4, const A& a5, const A& a6,
                const A& a7, const A& a8, const A& a9,
                const A& a10,const A& a11,const A& a12,
                const A& a13) const
    { return ToBeConstructed(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13); }

};

template <typename K>
class Construct_aff_transformation_3
{
    typedef typename K::Aff_transformation_3 ToBeConstructed;
public:
    typedef ToBeConstructed  result_type;

    ToBeConstructed
    operator()() const
    { return ToBeConstructed(); }

    template <class A1> 
    ToBeConstructed
    operator()( const A1& a1) const
    { return ToBeConstructed(a1); }

    template <class A1, class A2> 
    ToBeConstructed
    operator()( const A1& a1, const A2& a2) const
    { return ToBeConstructed(a1,a2); }

    template <class A1, class A2, class A3> 
    ToBeConstructed
    operator()( const A1& a1, const A2& a2, const A3& a3) const
    { return ToBeConstructed(a1,a2,a3); }

    template <class A1, class A2, class A3, class A4> 
    ToBeConstructed
    operator()( const A1& a1, const A2& a2, const A3& a3, const A4& a4) const
    { return ToBeConstructed(a1,a2,a3,a4); }

    template <class A1, class A2, class A3, class A4, class A5> 
    ToBeConstructed
    operator()( const A1& a1, const A2& a2, const A3& a3, const A4& a4,
	    const A5& a5) const
    { return ToBeConstructed(a1,a2,a3,a4,a5); }

    template <class A> 
    ToBeConstructed
    operator()( const A& a1, const A& a2, const A& a3,
                const A& a4, const A& a5, const A& a6 ) const
    { return ToBeConstructed(a1,a2,a3,a4,a5,a6); }

    template <class A> 
    ToBeConstructed
    operator()( const A& a1, const A& a2, const A& a3,
                const A& a4, const A& a5, const A& a6,
                const A& a7 ) const
    { return ToBeConstructed(a1,a2,a3,a4,a5,a6,a7); }

    template <class A> 
    ToBeConstructed
    operator()( const A& a1, const A& a2, const A& a3,
                const A& a4, const A& a5, const A& a6,
                const A& a7, const A& a8, const A& a9) const
    { return ToBeConstructed(a1,a2,a3,a4,a5,a6,a7,a8,a9); }

    template <class A> 
    ToBeConstructed
    operator()( const A& a1, const A& a2, const A& a3,
                const A& a4, const A& a5, const A& a6,
                const A& a7, const A& a8, const A& a9,
                const A& a10) const
    { return ToBeConstructed(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10); }

    template <class A> 
    ToBeConstructed
    operator()( const A& a1, const A& a2, const A& a3,
                const A& a4, const A& a5, const A& a6,
                const A& a7, const A& a8, const A& a9,
                const A& a10,const A& a11,const A& a12) const
    { return ToBeConstructed(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12); }

    template <class A> 
    ToBeConstructed
    operator()( const A& a1, const A& a2, const A& a3,
                const A& a4, const A& a5, const A& a6,
                const A& a7, const A& a8, const A& a9,
                const A& a10,const A& a11,const A& a12,
                const A& a13) const
    { return ToBeConstructed(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13); }

};

#endif // CGAL_NO_DEPRECATED_CODE

} // namespace CGALi
CGAL_END_NAMESPACE

#endif // CGAL_KERNEL_FUNCTION_OBJECTS_H
