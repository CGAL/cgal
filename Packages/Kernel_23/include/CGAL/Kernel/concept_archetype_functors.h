// ======================================================================
//
// Copyright (c) 1999,2003 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// release       : 
// release_date  : 
//
// file          : concept_archetype_functors.h
// package       : Kernel_23
// maintainer    : 
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Matthias Baesken
//
// coordinator   : MPI, Saarbruecken
// ======================================================================

#ifndef CGAL_CONCEPT_ARCHETYPE_FUNCTION_OBJECTS_H
#define CGAL_CONCEPT_ARCHETYPE_FUNCTION_OBJECTS_H

#include <CGAL/functional_base.h>
#include <CGAL/Quotient.h>
#include <CGAL/Origin.h>

// see kernel functors for original version
// all deprecated stuff was removed ...
// all functors were put into namespace CGALca 
// (to avoid ambiguities with the original functors)

CGAL_BEGIN_NAMESPACE

namespace CGALca {

template <typename K>
class Angle_2
{
    typedef typename K::Point_2 Point_2;
public:
    typedef Angle            result_type;
    typedef Arity_tag< 3 >   Arity;

    Angle operator()(const Point_2& p, const Point_2& q, 
                     const Point_2& r) const
    { CGAL::Angle a = CGAL::RIGHT;
      return a; 
    }
};

template <typename K>
class Angle_3
{
    typedef typename K::Point_3 Point_3;
public:
    typedef Angle            result_type;
    typedef Arity_tag< 3 >   Arity;

    Angle operator()(const Point_3& p, const Point_3& q, 
                     const Point_3& r) const
    { CGAL::Angle a = CGAL::RIGHT;
      return a; 
    }
};

template <typename K>
class Are_ordered_along_line_2
{
    typedef typename K::Point_2 Point_2;
public:
    typedef bool             result_type;
    typedef Arity_tag< 3 >   Arity;

    bool operator()(const Point_2& p, const Point_2& q, const Point_2& r) const
    { return true; }
};

template <typename K>
class Are_ordered_along_line_3
{
    typedef typename K::Point_3 Point_3;
public:
    typedef bool             result_type;
    typedef Arity_tag< 3 >   Arity;

    bool operator()(const Point_3& p, const Point_3& q, const Point_3& r) const
    { return true; }
};

template <typename K>
class Are_strictly_ordered_along_line_2
{
    typedef typename K::Point_2 Point_2;
public:
    typedef bool             result_type;
    typedef Arity_tag< 3 >   Arity;

    bool operator()(const Point_2& p, const Point_2& q, const Point_2& r) const
    { return true; }
};

template <typename K>
class Are_strictly_ordered_along_line_3
{
    typedef typename K::Point_3 Point_3;
public:
    typedef bool             result_type;
    typedef Arity_tag< 3 >   Arity;

    bool operator()(const Point_3& p, const Point_3& q, const Point_3& r) const
    { return true; }
};

template <typename K>
class Assign_2
{
    typedef typename K::Object_2 Object_2;
public:
    typedef bool             result_type;
    typedef Arity_tag< 2 >   Arity;

    template <class T>
    bool operator()(T& t, const Object_2& o) const
    { return true; }
};

template <typename K>
class Assign_3
{
    typedef typename K::Object_3 Object_3;
public:
    typedef bool             result_type;
    typedef Arity_tag< 2 >   Arity;

    template <class T>
    bool operator()(T& t, const Object_3& o) const
    { return true; }
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

    Bounded_side operator()( const Circle_2& c, const Point_2& p) const
    { return CGAL::ON_BOUNDARY; }

    Bounded_side operator()( const Triangle_2& t, const Point_2& p) const
    { return CGAL::ON_BOUNDARY; }

    Bounded_side operator()( const Iso_rectangle_2& r, const Point_2& p) const
    { return CGAL::ON_BOUNDARY; }
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
    { return CGAL::ON_BOUNDARY; }

    Bounded_side
    operator()( const Tetrahedron_3& t, const Point_3& p) const
    { return CGAL::ON_BOUNDARY; }

    Bounded_side
    operator()( const Iso_cuboid_3& c, const Point_3& p) const
    { return CGAL::ON_BOUNDARY; }
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
    { return true; }
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
    { return true; }
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
    { return true; }
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
    { return true; }
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
    { return true; }

    bool
    operator()( const Segment_2& s, const Point_2& p) const
    { return true; }
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
    { return true; }
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
    { return true; }
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
    { return CGAL::LARGER; }
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
    { return CGAL::LARGER; }
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
    { return CGAL::LARGER; }
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
    { return CGAL::LARGER; }

    Comparison_result
    operator()(const Segment_2& s1, const Segment_2& s2) const
    { return CGAL::LARGER; }
};

template <typename K>
class Compare_x_at_y_2
{
    typedef typename K::Point_2    Point_2;
    typedef typename K::Line_2     Line_2;
public:
    typedef Comparison_result      result_type;
    typedef Arity_tag< 3 >         Arity;

    Comparison_result
    operator()( const Point_2& p, const Line_2& h) const
    { return CGAL::LARGER; }

    Comparison_result
    operator()( const Point_2& p, const Line_2& h1, const Line_2& h2) const
    { return CGAL::LARGER; }

    Comparison_result
    operator()( const Line_2& l1, const Line_2& l2, const Line_2& h) const
    { return CGAL::LARGER; }

    Comparison_result
    operator()( const Line_2& l1, const Line_2& l2,
	        const Line_2& h1, const Line_2& h2) const
    { return CGAL::LARGER; }
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
    { return CGAL::LARGER; }
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
    { return CGAL::LARGER; }
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
    { return CGAL::LARGER; }
};

template <typename K>
class Compare_x_2
{
    typedef typename K::Point_2    Point_2;
    typedef typename K::Line_2     Line_2;
public:
    typedef Comparison_result      result_type;
    typedef Arity_tag< 2 >     Arity;

    Comparison_result
    operator()( const Point_2& p, const Point_2& q) const
    { return CGAL::LARGER; }

    Comparison_result
    operator()( const Point_2& p, const Line_2& l1, const Line_2& l2) const
    { return CGAL::LARGER; }

    Comparison_result
    operator()( const Line_2& l, const Line_2& h1, const Line_2& h2) const
    { return CGAL::LARGER; }

    Comparison_result
    operator()( const Line_2& l1, const Line_2& l2,
	        const Line_2& h1, const Line_2& h2) const
    { return CGAL::LARGER; }
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
    { return CGAL::LARGER; }
};

template <typename K>
class Compare_y_at_x_2
{
    typedef typename K::Point_2    Point_2;
    typedef typename K::Line_2     Line_2;
    typedef typename K::Segment_2  Segment_2;
public:
    typedef Comparison_result      result_type;
    typedef Arity_tag< 3 >         Arity;

    Comparison_result
    operator()( const Point_2& p, const Line_2& h) const
    { return CGAL::LARGER; }

    Comparison_result
    operator()( const Point_2& p, const Line_2& h1, const Line_2& h2) const
    { return CGAL::LARGER; }

    Comparison_result
    operator()( const Line_2& l1, const Line_2& l2, const Line_2& h) const
    { return CGAL::LARGER; }

    Comparison_result
    operator()( const Line_2& l1, const Line_2& l2,
	        const Line_2& h1, const Line_2& h2) const
    { return CGAL::LARGER; }

    Comparison_result
    operator()( const Point_2& p, const Segment_2& s) const
    { return CGAL::LARGER; }

    Comparison_result
    operator()( const Point_2& p,
	        const Segment_2& s1, const Segment_2& s2) const
    { return CGAL::LARGER; }
};

template <typename K>
class Compare_y_2
{
    typedef typename K::Point_2   Point_2;
    typedef typename K::Line_2    Line_2;
public:
    typedef Comparison_result     result_type;
    typedef Arity_tag< 2 >         Arity;

    Comparison_result
    operator()( const Point_2& p, const Point_2& q) const
    { return CGAL::LARGER; }

    Comparison_result
    operator()( const Point_2& p, const Line_2& l1, const Line_2& l2) const
    { return CGAL::LARGER; }

    Comparison_result
    operator()( const Line_2& l, const Line_2& h1, const Line_2& h2) const
    { return CGAL::LARGER; }

    Comparison_result
    operator()( const Line_2& l1, const Line_2& l2,
	        const Line_2& h1, const Line_2& h2) const
    { return CGAL::LARGER; }
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
    { return CGAL::LARGER; }
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
    { return CGAL::LARGER; }
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
    { FT val = 0; return val; }

    FT
    operator()( const Triangle_2& t ) const
    { FT val = 0; return val; }
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
    { FT val = 0; return val; }
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
    { FT val = 0; return val; }
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
    { FT val = 0; return val; }
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
    { FT val = 0; return val; }
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
    { FT val = 0; return val; }
};

template <typename K>
class Compute_squared_radius_2
{
    typedef typename K::FT          FT;
    typedef typename K::Point_2     Point_2;
    typedef typename K::Circle_2    Circle_2;
public:
    typedef FT               result_type;
    typedef Arity_tag< 1 >   Arity;

    FT
    operator()( const Circle_2& c) const
    { FT val = 0; return val; }

    FT
    operator()( const Point_2& p, const Point_2& q, const Point_2& r) const
    { FT val = 0; return val; }
};

template <typename K>
class Compute_squared_radius_3
{
    typedef typename K::FT          FT;
    typedef typename K::Point_3     Point_3;
    typedef typename K::Sphere_3    Sphere_3;
public:
    typedef FT               result_type;
    typedef Arity_tag< 1 >   Arity;

    FT
    operator()( const Sphere_3& s) const
    { FT val = 0; return val; }

    FT
    operator()( const Point_3& p, const Point_3& q, const Point_3& r) const
    { FT val = 0; return val; }

    FT
    operator()( const Point_3& p, const Point_3& q,
	        const Point_3& r, const Point_3& s) const
    { FT val = 0; return val; }
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
    { FT val = 0; return val; }

    FT
    operator()( const Iso_cuboid_3& c ) const
    { FT val = 0; return val; }
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
        Vector_3 v;
	return v;
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
    { Line_2 obj; return obj; }
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
    { Point_2 obj; return obj; }
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
    { Point_3 obj; return obj; }
};

template <typename K>
class Construct_centroid_2
{
    typedef typename K::Point_2  Point_2;
public:
    typedef Point_2          result_type;
    typedef Arity_tag< 3 >   Arity;

    Point_2
    operator()(const Point_2& p, const Point_2& q, const Point_2& r) const
    { Point_2 obj; return obj; }

    Point_2
    operator()(const Point_2& p, const Point_2& q, 
               const Point_2& r, const Point_2& s) const
    { Point_2 obj; return obj; }
};

template <typename K>
class Construct_centroid_3
{
    typedef typename K::Point_3  Point_3;
public:
    typedef Point_3          result_type;
    typedef Arity_tag< 3 >   Arity;

    Point_3
    operator()(const Point_3& p, const Point_3& q, const Point_3& r) const
    { Point_3 obj; return obj; }

    Point_3
    operator()(const Point_3& p, const Point_3& q, 
               const Point_3& r, const Point_3& s) const
    { Point_3 obj; return obj; }
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
    { return Circle_2(); }

    Circle_2
    operator()( const Point_2& p, const Point_2& q, const Point_2& r) const
    { return Circle_2(); }

    Circle_2
    operator()( const Point_2& p, const Point_2& q,
	        Orientation orientation = COUNTERCLOCKWISE) const
    { return Circle_2(); }

    Circle_2
    operator()( const Point_2& center,
	        Orientation orientation = COUNTERCLOCKWISE) const
    { return Circle_2(); }
};

template <typename K>
class Construct_circumcenter_2
{
    typedef typename K::Point_2  Point_2;
public:
    typedef Point_2          result_type;
    typedef Arity_tag< 3 >   Arity;

    Point_2
    operator()(const Point_2& p, const Point_2& q, const Point_2& r) const
    { return Point_2(); }
};

template <typename K>
class Construct_circumcenter_3
{
    typedef typename K::Point_3  Point_3;
public:
    typedef Point_3          result_type;
    typedef Arity_tag< 4 >   Arity;

    Point_3
    operator()(const Point_3& p, const Point_3& q, const Point_3& r) const
    { return Point_3(); }

    Point_3
    operator()(const Point_3& p, const Point_3& q,
	       const Point_3& r, const Point_3& s) const
    { return Point_3(); }
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
    { return Vector_3(); }
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

    Direction_2
    operator()(const Vector_2& v) const
    { return Direction_2(); }

    Direction_2
    operator()(const Line_2& l) const
    { return Direction_2(); }

    Direction_2
    operator()(const Ray_2& r) const
    { return Direction_2(); }

    Direction_2
    operator()(const Segment_2& s) const
    { return Direction_2(); }
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

    Direction_3
    operator()(const Vector_3& v) const
    { return Direction_3(); }

    Direction_3
    operator()(const Line_3& l) const
    { return Direction_3(); }

    Direction_3
    operator()(const Ray_3& r) const
    { return Direction_3(); }

    Direction_3
    operator()(const Segment_3& s) const
    { return Direction_3(); }
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
    { return Iso_cuboid_3(); }
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
    { return Iso_rectangle_2(); }
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
    {  return Point_3(); }
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
    typedef Line_2            result_type;
    typedef Arity_tag< 2 >    Arity;

    Line_2
    operator()() const
    { return Line_2(); }

    Line_2
    operator()(const Point_2& p, const Point_2& q) const
    { return Line_2(); }

    Line_2
    operator()(const Point_2& p, const Direction_2& d) const
    { return Line_2(); }

    Line_2
    operator()(const Segment_2& s) const
    { return Line_2(); }

    Line_2
    operator()(const Ray_2& r) const
    { return Line_2(); }
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
    typedef Line_3            result_type;
    typedef Arity_tag< 2 >    Arity;

    Line_3
    operator()() const
    { return Line_3(); }

    Line_3
    operator()(const Point_3& p, const Point_3& q) const
    { return Line_3(); }

    Line_3
    operator()(const Point_3& p, const Direction_3& d) const
    { return Line_3(); }

    Line_3
    operator()(const Segment_3& s) const
    { return Line_3(); }

    Line_3
    operator()(const Ray_3& r) const
    { return Line_3(); }
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
    { return Point_2(); }
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
    { return Point_3(); }
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
    { Object_2 obj; return obj; }
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
    { Object_3 obj; return obj; }
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
    { return Circle_2(); }
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
    { return Direction_2(); }
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
    { return Direction_3(); }
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
    { return Line_2(); }
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
    { return Line_3(); }
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
    { return Plane_3(); }
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
    { return Ray_2(); }
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
    { return Ray_3(); }
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
    { return Segment_2(); }
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
    { return Segment_3(); }
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
    { return Sphere_3(); }
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
    { return Triangle_2(); }
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
    { return Vector_2(); }
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
    { return Vector_3(); }
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
    { return Vector_3(); }
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
    { return Direction_2(); }
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
    { return Line_2(); }
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
    { return Line_3(); }
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
    { return Plane_3(); }
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
    { return Vector_2(); }
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
    typedef Plane_3          result_type;
    typedef Arity_tag< 2 >   Arity;

    Plane_3
    operator()() const
    { return Plane_3(); }

    Plane_3
    operator()(const RT& a, const RT& b, const RT& c, const RT& d) const
    { return Plane_3(); }

    Plane_3
    operator()(const Point_3& p, const Point_3& q, const Point_3& r) const
    { return Plane_3(); }

    Plane_3
    operator()(const Point_3& p, const Direction_3& d) const
    { return Plane_3(); }

    Plane_3
    operator()(const Line_3& l, const Point_3& p) const
    { return Plane_3(); }

    Plane_3
    operator()(const Ray_3& r, const Point_3& p) const
    { return Plane_3(); }

    Plane_3
    operator()(const Segment_3& s, const Point_3& p) const
    { return Plane_3(); }
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
    { return Point_2(); }

    Point_2
    operator()( const Segment_2& s, int i) const
    { return Point_2(); }

    Point_2
    operator()( const Ray_2& r, int i) const
    { return Point_2(); }
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
    { return Point_3(); }

    Point_3
    operator()( const Segment_3& s, int i) const
    { return Point_3(); }

    Point_3
    operator()( const Ray_3& r, int i) const
    { return Point_3(); }

    Point_3
    operator()( const Plane_3& p) const
    { return Point_3(); }
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
    { return Point_2(); }
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
    { return Point_3(); }
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
    { return Point_2(); }
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
    { return Point_3(); }

    Point_3
    operator()( const Plane_3& h, const Point_3& p ) const
    { return Point_3(); }
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
     {  return Point_2(); }
};

template <typename K>
class Construct_ray_2
{
    typedef typename K::Point_2      Point_2;
    typedef typename K::Direction_2  Direction_2;
    typedef typename K::Ray_2        Ray_2;
public:
     typedef Ray_2            result_type;
     typedef Arity_tag< 2 >   Arity;

     Ray_2
     operator()() const
     {  return Ray_2(); }

     Ray_2
     operator()(const Point_2& p, const Point_2& q) const
     {  return Ray_2(); }

     Ray_2
     operator()(const Point_2& p, const Direction_2& d) const
     {  return Ray_2(); }
};

template <typename K>
class Construct_ray_3
{
    typedef typename K::Point_3      Point_3;
    typedef typename K::Direction_3  Direction_3;
    typedef typename K::Ray_3        Ray_3;
public:
     typedef Ray_3            result_type;
     typedef Arity_tag< 2 >   Arity;

     Ray_3
     operator()() const
     {  return Ray_3(); }

     Ray_3
     operator()(const Point_3& p, const Point_3& q) const
     {  return Ray_3(); }

     Ray_3
     operator()(const Point_3& p, const Direction_3& d) const
     {  return Ray_3(); }
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
    {  return Vector_2(); }

    Vector_2
    operator()( const Vector_2& v, const Quotient<RT>& scale) const
    {  return Vector_2(); }
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
    {  return Vector_3(); }

    Vector_3
    operator()( const Vector_3& v, const Quotient<RT>& scale) const
    {  return Vector_3(); }
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
    {  return Segment_2(); }
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
    {  return Segment_3(); }
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
    {  return Sphere_3(); }

    Sphere_3
    operator()( const Point_3& p, const Point_3& q,
	        const Point_3& r, const Point_3& s) const
    {  return Sphere_3(); }

    Sphere_3
    operator()( const Point_3& p, const Point_3& q, const Point_3& r,
	        Orientation orientation = COUNTERCLOCKWISE) const
    {  return Sphere_3(); }

    Sphere_3
    operator()( const Point_3& p, const Point_3& q,
	        Orientation orientation = COUNTERCLOCKWISE) const
    {  return Sphere_3(); }

    Sphere_3
    operator()( const Point_3& center,
	        Orientation orientation = COUNTERCLOCKWISE) const
    {  return Sphere_3(); }
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
    { return Line_2(); }

    Line_2
    operator()( const Segment_2& s) const
    { return Line_2(); }
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
    { return Line_3(); }

    Line_3
    operator()( const Segment_3& s) const
    { return Line_3(); }
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
    { return Plane_3(); }
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
    { return Tetrahedron_3(); }
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
    {  return Point_2(); }
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
    {  return Point_3(); }
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
    { return Triangle_2(); }
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
    { return Triangle_3(); }
};

template <typename K>
class Construct_vector_2
{
    typedef typename K::RT           RT;
    typedef typename K::Vector_2     Vector_2;
    typedef typename K::Point_2      Point_2;
public:
    typedef Vector_2         result_type;
    typedef Arity_tag< 2 >   Arity;

    Vector_2
    operator()() const
    { return Vector_2(); }

    Vector_2
    operator()( const Point_2& p, const Point_2& q) const
    { return Vector_2(); }

    Vector_2
    operator()( Null_vector n) const
    { return Vector_2(); }
};

template <typename K>
class Construct_vector_3
{
    typedef typename K::RT           RT;
    typedef typename K::Vector_3     Vector_3;
    typedef typename K::Point_3      Point_3;
public:
    typedef Vector_3         result_type;
    typedef Arity_tag< 2 >   Arity;

    Vector_3
    operator()() const
    { return Vector_3(); }

    Vector_3
    operator()( const Point_3& p, const Point_3& q) const
    { return Vector_3(); }

    Vector_3
    operator()( Null_vector n) const
    { return Vector_3(); }
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
    { return Point_2(); }

    Point_2
    operator()( const Triangle_2& t, int i) const
    { return Point_2(); }

    Point_2
    operator()( const Iso_rectangle_2& r, int i) const
    { return Point_2(); }
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
    { return Point_3(); }

    Point_3
    operator()( const Triangle_3& t, int i) const
    { return Point_3(); }

    Point_3
    operator()( const Iso_cuboid_3& r, int i) const
    { return Point_3(); }

    Point_3
    operator()( const Tetrahedron_3& t, int i) const
    { return Point_3(); }
};

template <typename K>
class Coplanar_orientation_3
{
    typedef typename K::Point_3   Point_3;
public:
    typedef Orientation  result_type;
    typedef Arity_tag< 4 >   Arity;

    Orientation
    operator()(const Point_3& p, const Point_3& q, const Point_3& r) const
    { return CGAL::COLLINEAR; }

    Orientation
    operator()( const Point_3& p, const Point_3& q,
	        const Point_3& r, const Point_3& t) const
    { return CGAL::COLLINEAR; }
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
    { return CGAL::ON_BOUNDARY; }
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
    { return true; }
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
    { return true; }
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
    { return true; }
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
    { return true; }
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
    { return true; }
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
    { return true; }
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
    { return true; }
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
    { return true; }
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
    { return true; }
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
    { return true; }
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
    { return true; }
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
    { return true; }
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
    { return true; }

    bool
    operator()( const Triangle_2& t, const Point_2& p) const
    { return true; }

    bool
    operator()( const Iso_rectangle_2& r, const Point_2& p) const
    { return true; }
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
    { return true; }

    bool
    operator()( const Tetrahedron_3& t, const Point_3& p) const
    { return true; }

    bool
    operator()( const Iso_cuboid_3& c, const Point_3& p) const
    { return true; }
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
    { return true; }

    bool
    operator()( const Triangle_2& t, const Point_2& p) const
    { return true; }

    bool
    operator()( const Iso_rectangle_2& r, const Point_2& p) const
    { return true; }
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
    { return true; }

    bool
    operator()( const Tetrahedron_3& t, const Point_3& p) const
    { return true; }

    bool
    operator()( const Iso_cuboid_3& c, const Point_3& p) const
    { return true; }
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
    { return true; }

    bool
    operator()( const Triangle_2& t, const Point_2& p) const
    { return true; }

    bool
    operator()( const Line_2& l, const Point_2& p) const
    { return true; }
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
    { return true; }

    bool
    operator()( const Tetrahedron_3& t, const Point_3& p) const
    { return true; }

    bool
    operator()( const Plane_3& pl, const Point_3& p) const
    { return true; }
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
    { return true; }

    bool
    operator()( const Triangle_2& t, const Point_2& p) const
    { return true; }

    bool
    operator()( const Line_2& l, const Point_2& p) const
    { return true; }
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
    { return true; }

    bool
    operator()( const Tetrahedron_3& t, const Point_3& p) const
    { return true; }

    bool
    operator()( const Plane_3& pl, const Point_3& p) const
    { return true; }
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
    { return true; }

    bool
    operator()( const Triangle_2& t, const Point_2& p) const
    { return true; }

    bool
    operator()( const Iso_rectangle_2& r, const Point_2& p) const
    { return true; }
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
    { return true; }

    bool
    operator()( const Tetrahedron_3& t, const Point_3& p) const
    { return true; }

    bool
    operator()( const Iso_cuboid_3& c, const Point_3& p) const
    { return true; }
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
    { return true; }

    bool
    operator()( const Ray_2& r, const Point_2& p) const
    { return true; }

    bool
    operator()( const Segment_2& s, const Point_2& p) const
    { return true; }
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
    { return true; }

    bool
    operator()( const Ray_3& r, const Point_3& p) const
    { return true; }

    bool
    operator()( const Segment_3& s, const Point_3& p) const
    { return true; }

    bool
    operator()( const Plane_3& pl, const Point_3& p) const
    { return true; }

    bool
    operator()( const Triangle_3& t, const Point_3& p) const
    { return true; }
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
    { return Object_2(); }
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
    { return Object_3(); }
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
    { return true; }

    bool
    operator()( const Iso_rectangle_2& r) const
    { return true; }

    bool
    operator()( const Line_2& l) const
    { return true; }

    bool
    operator()( const Ray_2& r) const
    { return true; }

    bool
    operator()( const Segment_2& s) const
    { return true; }

    bool
    operator()( const Triangle_2& t) const
    { return true; }
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
    { return true; }

    bool
    operator()( const Line_3& l) const
    { return true; }

    bool
    operator()( const Plane_3& pl) const
    { return true; }

    bool
    operator()( const Ray_3& r) const
    { return true; }

    bool
    operator()( const Segment_3& s) const
    { return true; }

    bool
    operator()( const Sphere_3& s) const
    { return true; }

    bool
    operator()( const Triangle_3& t) const
    { return true; }

    bool
    operator()( const Tetrahedron_3& t) const
    { return true; }
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
    { return true; }

    bool
    operator()( const Segment_2& s) const
    { return true; }

    bool
    operator()( const Ray_2& r) const
    { return true; }
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
    { return true; }

    bool
    operator()( const Segment_2& s) const
    { return true; }

    bool
    operator()( const Ray_2& r) const
    { return true; }
};


template <typename K>
class Left_turn_2
{
    typedef typename K::Point_2   Point_2;
public:
    typedef bool             result_type;
    typedef Arity_tag< 3 >   Arity;

    bool
    operator()(const Point_2& p, const Point_2& q, const Point_2& r) const
    { return true; }
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
    { return true; }
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
    { return true; }
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
    { return true; }
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
    { return true; }
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
    { return true; }
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
    { return true; }
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
    { return true; }
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
    { return true; }
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
    { return true; }
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
    { return true; }
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
    { return true; }
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
    { return true; }
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
    { return true; }
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
    { return true; }
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
    { return CGAL::COLLINEAR; }
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
    { return CGAL::COLLINEAR; }
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
    { return CGAL::ON_POSITIVE_SIDE; }

    Oriented_side
    operator()( const Line_2& l, const Point_2& p) const
    { return CGAL::ON_POSITIVE_SIDE; }

    Oriented_side
    operator()( const Triangle_2& t, const Point_2& p) const
    { return CGAL::ON_POSITIVE_SIDE; }
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
    { return CGAL::ON_POSITIVE_SIDE; }

    Oriented_side
    operator()( const Plane_3& pl, const Point_3& p) const
    { return CGAL::ON_POSITIVE_SIDE; }

    Oriented_side
    operator()( const Tetrahedron_3& t, const Point_3& p) const
    { return CGAL::ON_POSITIVE_SIDE; }
};

template <typename K>
class Side_of_bounded_circle_2
{
    typedef typename K::Point_2        Point_2;
public:
    typedef Bounded_side     result_type;
    typedef Arity_tag< 4 >   Arity;

    Bounded_side
    operator()( const Point_2& p, const Point_2& q, const Point_2& t) const
    { return CGAL::ON_BOUNDARY; }

    Bounded_side
    operator()( const Point_2& p, const Point_2& q,
	        const Point_2& r, const Point_2& t) const
    { return CGAL::ON_BOUNDARY; }
};

template <typename K>
class Side_of_bounded_sphere_3
{
    typedef typename K::Point_3        Point_3;
public:
    typedef Bounded_side   result_type;
    typedef Arity_tag< 5 >   Arity;

    Bounded_side
    operator()( const Point_3& p, const Point_3& q, const Point_3& t) const
    { return CGAL::ON_BOUNDARY; }

    Bounded_side
    operator()( const Point_3& p, const Point_3& q,
	        const Point_3& r, const Point_3& t) const
    { return CGAL::ON_BOUNDARY; }

    Bounded_side
    operator()( const Point_3& p, const Point_3& q, const Point_3& r,
	        const Point_3& s, const Point_3& t) const
    { return CGAL::ON_BOUNDARY; }
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
    { return CGAL::ON_POSITIVE_SIDE; }
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
    { return CGAL::ON_POSITIVE_SIDE; }
};


} // namespace CGALca

CGAL_END_NAMESPACE

#endif 
