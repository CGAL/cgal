// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
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
// author(s)     : Stefan Schirra <Stefan.Schirra@@mpi-sb.mpg.de>
//
// coordinator   : MPI, Saarbruecken
// ======================================================================

#ifndef CGAL_KERNEL_FUNCTION_OBJECTS_H
#define CGAL_KERNEL_FUNCTION_OBJECTS_H

#include <CGAL/functional_base.h>

CGAL_BEGIN_NAMESPACE
namespace CGALi {

template <class ToBeConstructed>
class Construct
{
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

template <class ToBeConstructed>
class Call_make_object_to_get
{
  public:
    typedef  ToBeConstructed  result_type;
    typedef  Arity_tag< 1 >   Arity;

    template <class Cls>
    ToBeConstructed
    operator()( const Cls& c) const
    { return make_object(c); }
};

template <class Vector>
class Construct_scaled_vector
{
   public:
     typedef Vector    result_type;
     typedef  Arity_tag< 2 >   Arity;

     template <class NT>
     Vector
     operator()( const Vector& v, const NT& scale) const
     {  return v * scale; }
};

template <class Point>
class Construct_translated_point
{
   public:
     typedef Point    result_type;
     typedef  Arity_tag< 2 >   Arity;

     template <class Vector>
     Point
     operator()( const Point& p, const Vector& v) const
     {  return p + v; }
};

template <class Point_2>
class Construct_projected_xy_point
{
   public:
     typedef Point_2    result_type;
     typedef  Arity_tag< 2 >   Arity;

     template <class Plane_3, class Point_3>
     Point_2
     operator()( const Plane_3& h, const Point_3& p) const
     {  return h.to_2d(p); }
};

template <class ReturnType>
class Call_point_to_get
{
  public:
    typedef ReturnType     result_type;

    template <class Cls>
    ReturnType
    operator()( const Cls& c) const
    { return c.point(); }

    template <class Cls>
    ReturnType
    operator()( const Cls& c, int i) const
    { return c.point(i); }
};

template <class ReturnType>
class Call_second_point_to_get
{
  public:
    typedef ReturnType     result_type;
    typedef  Arity_tag< 1 >   Arity;

    template <class Cls>
    ReturnType
    operator()( const Cls& c) const
    { return c.second_point(); }
};

template <class ReturnType>
class Call_perpendicular_to_get
{
  public:
    typedef ReturnType     result_type;
    typedef  Arity_tag< 2 >   Arity;
  
  // This is only called for Vector_2, Direction_2, and Line_2
  // and in all cases only the two parameter variant is documented
  
    template <class Cls>
    ReturnType
    operator()( const Cls& c) const
    { return c.perpendicular(); }

    template <class Cls, class A1>
    ReturnType
    operator()( const Cls& c, const A1& a1) const
    { return c.perpendicular(a1); }
};

template <class ReturnType>
class Call_perpendicular_plane_to_get
{
  public:
    typedef ReturnType     result_type;
    typedef  Arity_tag< 2 >   Arity;

    template <class Cls, class A1>
    ReturnType
    operator()( const Cls& c, const A1& a1) const
    { return c.perpendicular_plane(a1); }
};

template <class ReturnType>
class Call_perpendicular_line_to_get
{
  public:
    typedef ReturnType     result_type;
    typedef  Arity_tag< 2 >   Arity;

    template <class Cls, class A1>
    ReturnType
    operator()( const Cls& c, const A1& a1) const
    { return c.perpendicular_line(a1); }
};

template <class Vector>
class Call_orthogonal_vector_to_get
{
  public:
    typedef Vector     result_type;
    typedef  Arity_tag< 1 >   Arity;

    template <class Cls>
    Vector
    operator()( const Cls& c ) const
    { return c.orthogonal_vector(); }
};

template <class Point>
class Call_projection_to_get
{
  public:
    typedef Point     result_type;
    typedef  Arity_tag< 2 >   Arity;

    template <class Cls>
    Point
    operator()( const Cls& c, const Point& p ) const
    { return c.projection(p); }
};

template <class ReturnType>
class Call_squared_area_to_get
{
  public:
    typedef ReturnType     result_type;
    typedef  Arity_tag< 1 >   Arity;

    template <class Cls>
    ReturnType
    operator()( const Cls& c ) const
    { return c.squared_area(); }
};

template <class ReturnType>
class Call_area_to_get
{
  public:
    typedef ReturnType     result_type;
    typedef Arity_tag< 1 >   Arity;

    template <class Cls>
    ReturnType
    operator()( const Cls& c ) const
    { return c.area(); }
};

template <class ReturnType>
class Call_volume_to_get
{
  public:
    typedef ReturnType     result_type;
    typedef Arity_tag< 1 >   Arity;

    template <class Cls>
    ReturnType
    operator()( const Cls& c ) const
    { return c.volume(); }
};

template <class Point>
class p_Midpoint
{
  public:
    typedef Point          result_type;
    typedef Arity_tag< 2 >   Arity;

    Point
    operator()(const Point& p, const Point& q) const
    { return midpoint(p,q); }
};

template <class Point>
class p_Center
{
  public:
    typedef Point          result_type;
    typedef Arity_tag< 1 >   Arity;

    template <class C>
    Point
    operator()(const C&c) const
    { return c.center(); }
};

template <class Point>
class p_Circumcenter
{
  public:
    typedef Point          result_type;

    Point
    operator()(const Point& p, const Point& q, const Point& r) const
    { return circumcenter(p,q,r); }

    Point
    operator()(const Point& p, const Point& q, 
               const Point& r, const Point& s) const
    { return circumcenter(p,q,r,s); }
};

template <class Point>
class p_Centroid
{
  public:
    typedef Point          result_type;

    Point
    operator()(const Point& p, const Point& q, const Point& r) const
    { return centroid(p,q,r); }

    Point
    operator()(const Point& p, const Point& q, 
               const Point& r, const Point& s) const
    { return centroid(p,q,r,s); }
};

template <class Point, class Line>
class pl_Bisector
{
  public:
    typedef Line           result_type;
    typedef Arity_tag< 2 >   Arity;

    Line
    operator()(const Point& p, const Point& q) const
    { return bisector(p,q); }
};

template <class Vector>
class v_Opposite
{
   public: 
     typedef Vector        result_type;
     typedef Arity_tag< 1 >   Arity;
   
   Vector
   operator()(const Vector& v) const
   { return -v; }
};

template <class Vector>
class v_Cross_product
{
  public:
    typedef Vector          result_type;
    typedef Arity_tag< 2 >   Arity;

    Vector
    operator()(const Vector& v, const Vector& w) const
    { return cross_product(v, w); }
};

template <class Vector>
class v_Base
{
   public:
     typedef Vector        result_type;
     typedef Arity_tag< 2 >   Arity;

     template <class Plane>
     Vector
     operator()( const Plane& pl, int index )
     {
       if (index == 1)
         return pl.base1();
       else 
         return pl.base2();
     }
};

class Intersect
{
  public:
    typedef CGAL::Object   result_type;
    typedef Arity_tag< 2 >   Arity;

    template <class T1, class T2>
    CGAL::Object
    operator()(const T1& t1, const T2& t2) const
    { return intersection( t1, t2); }
};

class Do_intersect
{
  public:
    typedef bool   result_type;
    typedef Arity_tag< 2 >   Arity;

    template <class T1, class T2>
    bool
    operator()(const T1& t1, const T2& t2) const
    { return do_intersect( t1, t2); }
};

class Assign
{
  public:
    typedef bool           result_type;
    typedef Arity_tag< 2 >   Arity;

    template <class T1>
    bool
    operator()(T1& t1, const CGAL::Object& o) const
    { return assign( t1, o); }
};

template <class ReturnType>
class Call_y_at_x_to_get
{
  public:
    typedef ReturnType     result_type;
    typedef Arity_tag< 2 >   Arity;

    template <class Cls>
    ReturnType
    operator()( const Cls& c, const ReturnType& x) const
    { return c.y_at_x(x); }
};

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

template <class ReturnType>
class Call_squared_length_to_get
{
  public:
    typedef ReturnType     result_type;
    typedef Arity_tag< 1 >   Arity;

    template <class Cls>
    ReturnType
    operator()( const Cls& c) const
    { return c.squared_length(); }
};

template <class ReturnType>
class Call_squared_distance
{
  public:
    typedef ReturnType     result_type;
    typedef Arity_tag< 2 >   Arity;

    template <class T1, class T2>
    ReturnType
    operator()( const T1& t1, const T2& t2) const
    { return squared_distance(t1, t2); }
};

template <class ReturnType>
class Call_squared_radius
{
  public:
    typedef ReturnType     result_type;

    template <class T1>
    ReturnType
    operator()( const T1& t1) const
    { return t1.squared_radius(); }

    template <class T1>
    ReturnType
    operator()( const T1& t1, const T1& t2, const T1& t3) const
    { return squared_radius(t1, t2, t3); }

    template <class T1>
    ReturnType
    operator()( const T1& t1, const T1& t2, const T1& t3, const T1& t4) const
    { return squared_radius(t1, t2, t3, t4); }
};

class p_Angle
{
  public:
    typedef Angle           result_type;
    typedef Arity_tag< 3 >   Arity;

    template <class T>
    Angle
    operator()(const T& p, const T& q, const T& r) const
    { return angle(p, q, r); }
};

template <class Point_3>
class p_Lifted
{
   public:
     typedef Point_3        result_type;
     typedef Arity_tag< 2 >   Arity;

     template <class Plane, class Point_2>
     Point_3
     operator()(const Plane& p, const Point_2& pt) const
     {  return p.to_3d(pt); }
};

class Counterclockwise_in_between
{
  public:
    typedef bool           result_type;
    typedef Arity_tag< 3 >   Arity;

    template <class T>
    bool
    operator()(const T& p, const T& q, const T& r) const
    { return p.counterclockwise_in_between(q,r); }
};

class Collinear
{
  public:
    typedef bool           result_type;
    typedef Arity_tag< 3 >   Arity;

    template <class T>
    bool
    operator()(const T& p, const T& q, const T& r) const
    { return collinear(p,q,r); }
};

class Coplanar
{
  public:
    typedef bool  result_type;
    typedef Arity_tag< 4 >   Arity;

    template <class T>
    bool
    operator()(const T& p, const T& q, const T& r, const T& s) const
    { return coplanar(p,q,r,s); }
};

class Coplanar_orientation
{
  public:
    typedef Orientation  result_type;

    template <class P>
    Orientation
    operator()(const P& p, const P& q, const P& r) const
    { return coplanar_orientation(p,q,r); }

    template <class P>
    Orientation
    operator()(const P& p, const P& q, const P& r, const P& t) const
    { return coplanar_orientation(p,q,r,t); }
};

class Coplanar_side_of_bounded_circle
{
  public:
    typedef Bounded_side  result_type;
    typedef Arity_tag< 4 >   Arity;

    template <class P>
    Bounded_side
    operator()(const P& p, const P& q, const P& r, const P& t) const
    { return coplanar_side_of_bounded_circle(p,q,r,t); }
};

class Side_of_oriented_circle
{
  public:
    typedef Oriented_side  result_type;
    typedef Arity_tag< 4 >   Arity;

    template <class T>
    Oriented_side
    operator()(const T& p, const T& q, const T& r, const T& t) const
    { return side_of_oriented_circle(p,q,r,t); }
};

class Side_of_bounded_circle
{
  public:
    typedef Bounded_side   result_type;

    template <class T>
    Bounded_side
    operator()(const T& p, const T& q, const T& t) const
    { return side_of_bounded_circle(p,q,t); }

    template <class T>
    Bounded_side
    operator()(const T& p, const T& q, const T& r, const T& t) const
    { return side_of_bounded_circle(p,q,r,t); }
};

class Side_of_oriented_sphere
{
  public:
    typedef Oriented_side  result_type;
    typedef Arity_tag< 5 >   Arity;

    template <class T>
    Oriented_side
    operator()(const T& p, const T& q, const T& r, const T& s,
	    const T& t) const
    { return side_of_oriented_sphere(p,q,r,s,t); }
};

class Side_of_bounded_sphere
{
  public:
    typedef Bounded_side   result_type;

    template <class T>
    Bounded_side
    operator()(const T& p, const T& q, const T& t) const
    { return side_of_bounded_sphere(p,q,t); }

    template <class T>
    Bounded_side
    operator()(const T& p, const T& q, const T& r, const T& t) const
    { return side_of_bounded_sphere(p,q,r,t); }

    template <class T>
    Bounded_side
    operator()(const T& p, const T& q, const T& r, const T& s,
	    const T& t) const
    { return side_of_bounded_sphere(p,q,r,s,t); }
};

class Call_is_horizontal
{
  public:
    typedef bool           result_type;
    typedef Arity_tag< 1 >   Arity;

    template <class Cls>
    bool
    operator()( const Cls& c) const
    { return c.is_horizontal(); }
};

class Call_is_vertical
{
  public:
    typedef bool           result_type;
    typedef Arity_tag< 1 >   Arity;

    template <class Cls>
    bool
    operator()( const Cls& c) const
    { return c.is_vertical(); }
};

class Call_is_degenerate
{
  public:
    typedef bool           result_type;
    typedef Arity_tag< 1 >   Arity;

    template <class Cls>
    bool
    operator()( const Cls& c) const
    { return c.is_degenerate(); }
};

class Call_has_on_bounded_side
{
  public:
    typedef bool           result_type;
    typedef Arity_tag< 2 >   Arity;

    template <class Cls, class Arg>
    bool
    operator()( const Cls& c, const Arg& a) const
    { return c.has_on_bounded_side(a); }
};

class Call_has_on_unbounded_side
{
  public:
    typedef bool           result_type;
    typedef Arity_tag< 2 >   Arity;

    template <class Cls, class Arg>
    bool
    operator()( const Cls& c, const Arg& a) const
    { return c.has_on_unbounded_side(a); }
};

class Call_has_on_boundary
{
  public:
    typedef bool           result_type;
    typedef Arity_tag< 2 >   Arity;

    template <class Cls, class Arg>
    bool
    operator()( const Cls& c, const Arg& a) const
    { return c.has_on_boundary(a); }
};

class Call_has_on_positive_side
{
  public:
    typedef bool           result_type;
    typedef Arity_tag< 2 >   Arity;

    template <class Cls, class Arg>
    bool
    operator()( const Cls& c, const Arg& a) const
    { return c.has_on_positive_side(a); }
};

class Call_has_on_negative_side
{
  public:
    typedef bool           result_type;
    typedef Arity_tag< 2 >   Arity;

    template <class Cls, class Arg>
    bool
    operator()( const Cls& c, const Arg& a) const
    { return c.has_on_negative_side(a); }
};

class Call_oriented_side
{
  public:
    typedef Oriented_side   result_type;
    typedef Arity_tag< 2 >   Arity;

    template <class Cls, class Arg>
    Oriented_side
    operator()( const Cls& c, const Arg& a) const
    { return c.oriented_side(a); }
};

class Call_bounded_side
{
  public:
    typedef Bounded_side    result_type;
    typedef Arity_tag< 2 >   Arity;

    template <class Cls, class Arg>
    Bounded_side
    operator()( const Cls& c, const Arg& a) const
    { return c.bounded_side(a); }
};

class Less_x
{
  public:
    typedef bool           result_type;
    typedef Arity_tag< 2 >   Arity;

    template <class T1, class T2>
    bool
    operator()( const T1& a1, const T2& a2) const
    { return less_x(a1,a2); }
};

class Less_y
{
  public:
    typedef bool           result_type;
    typedef Arity_tag< 2 >   Arity;

    template <class T1, class T2>
    bool
    operator()( const T1& a1, const T2& a2) const
    { return less_y(a1,a2); }
};

class Less_z
{
  public:
    typedef bool           result_type;
    typedef Arity_tag< 2 >   Arity;

    template <class T1, class T2>
    bool
    operator()( const T1& a1, const T2& a2) const
    { return less_y(a1,a2); }
};

class Less_xy
{
  public:
    typedef bool           result_type;
    typedef Arity_tag< 2 >   Arity;

    template <class T1, class T2>
    bool
    operator()( const T1& a1, const T2& a2) const
    { return lexicographically_xy_smaller(a1,a2); }
};

class Less_yx
{
  public:
    typedef bool           result_type;
    typedef Arity_tag< 2 >   Arity;

    template <class T1, class T2>
    bool
    operator()( const T1& a1, const T2& a2) const
    { return lexicographically_yx_smaller(a1,a2); }
};

class Less_xyz
{
  public:
    typedef bool           result_type;
    typedef Arity_tag< 2 >   Arity;

    template <class T1, class T2>
    bool
    operator()( const T1& a1, const T2& a2) const
    { return lexicographically_xyz_smaller(a1,a2); }
};

class Equal
{
  public:
    typedef bool           result_type;
    typedef Arity_tag< 2 >   Arity;

    template <class T1, class T2>
    bool
    operator()(const T1& p, const T2& q) const
    { return p == q; }
};

class Equal_x
{
  public:
    typedef bool           result_type;
    typedef Arity_tag< 2 >   Arity;

    template <class T1, class T2>
    bool
    operator()( const T1& a1, const T2& a2) const
    { return x_equal(a1,a2); }
};

class Equal_y
{
  public:
    typedef bool           result_type;
    typedef Arity_tag< 2 >   Arity;

    template <class T1, class T2>
    bool
    operator()( const T1& a1, const T2& a2) const
    { return y_equal(a1,a2); }
};

class Equal_z
{
  public:
    typedef bool           result_type;
    typedef Arity_tag< 2 >   Arity;

    template <class T1, class T2>
    bool
    operator()( const T1& a1, const T2& a2) const
    { return z_equal(a1,a2); }
};

class Equal_xy
{
  public:
    typedef bool           result_type;
    typedef Arity_tag< 2 >   Arity;

    template <class T1, class T2>
    bool
    operator()( const T1& a1, const T2& a2) const
    { return equal_xy(a1,a2); }
};

class Equal_xyz
{
  public:
    typedef bool           result_type;
    typedef Arity_tag< 2 >   Arity;

    template <class T1, class T2>
    bool
    operator()( const T1& a1, const T2& a2) const
    { return equal_xyz(a1,a2); }
};

class Compare_x
{
  public:
    typedef Comparison_result result_type;

    template <class T1, class T2>
    Comparison_result
    operator()( const T1& a1, const T2& a2) const
    { return compare_x(a1,a2); }

    template <class T1, class T2, class T3>
    Comparison_result
    operator()( const T1& a1, const T2& a2, const T3& a3) const
    { return compare_x(a1,a2,a3); }

    template <class T1, class T2, class T3, class T4>
    Comparison_result
    operator()( const T1& a1, const T2& a2, const T3& a3, const T4& a4) const
    { return compare_x(a1,a2,a3,a4); }

};

class Compare_y
{
  public:
    typedef Comparison_result result_type;

    template <class T1, class T2>
    Comparison_result
    operator()( const T1& a1, const T2& a2) const
    { return compare_y(a1,a2); }

    template <class T1, class T2, class T3>
    Comparison_result
    operator()( const T1& a1, const T2& a2, const T3& a3) const
    { return compare_y(a1,a2,a3); }

    template <class T1, class T2, class T3, class T4>
    Comparison_result
    operator()( const T1& a1, const T2& a2, const T3& a3, const T4& a4) const
    { return compare_y(a1,a2,a3,a4); }


};

class Compare_z
{
  public:
    typedef Comparison_result result_type;
    typedef Arity_tag< 2 >   Arity;

    template <class T1, class T2>
    Comparison_result
    operator()( const T1& a1, const T2& a2) const
    { return compare_z(a1,a2); }
};

class Compare_xy
{
  public:
    typedef Comparison_result result_type;
    typedef Arity_tag< 2 >   Arity;

    template <class T1, class T2>
    Comparison_result
    operator()( const T1& a1, const T2& a2) const
    { return compare_xy(a1,a2); }
};

class Compare_xyz
{
  public:
    typedef Comparison_result result_type;
    typedef Arity_tag< 2 >   Arity;

    template <class T1, class T2>
    Comparison_result
    operator()( const T1& a1, const T2& a2) const
    { return compare_xyz(a1,a2); }
};

class Compare_y_at_x
{
  public:
    typedef Comparison_result result_type;

    template <class T1, class T2>
    Comparison_result
    operator()( const T1& a1, const T2& a2) const
    { return compare_y_at_x(a1,a2); }

    template <class T1, class T2, class T3>
    Comparison_result
    operator()( const T1& a1, const T2& a2, const T3& a3) const
    { return compare_y_at_x(a1,a2,a3); }
    
    template <class T1, class T2, class T3, class T4>
    Comparison_result
    operator()( const T1& a1, const T2& a2, const T3& a3, const T4& a4) const
    { return compare_y_at_x(a1,a2,a3,a4); }
};

class Compare_x_at_y
{
  public:
    typedef Comparison_result result_type;

    template <class T1, class T2>
    Comparison_result
    operator()( const T1& a1, const T2& a2) const
    { return compare_x_at_y(a1,a2); }

    template <class T1, class T2, class T3>
    Comparison_result
    operator()( const T1& a1, const T2& a2, const T3& a3) const
    { return compare_x_at_y(a1,a2,a3); }
    
    template <class T1, class T2, class T3, class T4>
    Comparison_result
    operator()( const T1& a1, const T2& a2, const T3& a3, const T4& a4) const
    { return compare_x_at_y(a1,a2,a3,a4); }
};

template <class T>
class Compare_distance
{
  public:
    typedef Comparison_result           result_type;
    typedef Arity_tag< 3 >   Arity;

    Comparison_result
    operator()(const T& p, const T& q, const T& r) const
    { return cmp_dist_to_point(p,q,r); }
};

template <class T>
class Compare_angle_with_x_axis
{
  public:
    typedef Comparison_result           result_type;
    typedef Arity_tag< 2 >   Arity;

    Comparison_result
    operator()(const T& p, const T& q) const
    { return compare_angle_with_x_axis(p,q); }
};

template <class Plane, class Point>
class Less_signed_distance_to_plane
{
  public:
    typedef bool           result_type;
    typedef Arity_tag< 3 >   Arity;

    bool operator()( const Plane & p, const Point& q, const Point& r) const
    { return has_smaller_signed_dist_to_plane(p,q,r); }
};

class Are_ordered_along_line
{
  public:
    typedef bool           result_type;
    typedef Arity_tag< 3 >   Arity;

    template <class T>
    bool
    operator()(const T& p, const T& q, const T& r) const
    { return are_ordered_along_line(p,q,r); }
};

class Are_strictly_ordered_along_line
{
  public:
    typedef bool           result_type;
    typedef Arity_tag< 3 >   Arity;

    template <class T>
    bool
    operator()(const T& p, const T& q, const T& r) const
    { return are_strictly_ordered_along_line(p,q,r); }
};

class Collinear_are_ordered_along_line
{
  public:
    typedef bool           result_type;
    typedef Arity_tag< 3 >   Arity;

    template <class T>
    bool
    operator()(const T& p, const T& q, const T& r) const
    { return collinear_are_ordered_along_line(p,q,r); }
};

class Collinear_are_strictly_ordered_along_line
{
  public:
    typedef bool           result_type;
    typedef Arity_tag< 3 >   Arity;

    template <class T>
    bool
    operator()(const T& p, const T& q, const T& r) const
    { return collinear_are_strictly_ordered_along_line(p,q,r); }
};

#ifndef CGAL_NO_DEPRECATED_CODE
class Call_transform
{
 public:
    template <class Transformation, class ArgumentType>
    ArgumentType
    operator()( const ArgumentType& a, const Transformation& t) const
    { return a.transform(t); }
};
#endif // CGAL_NO_DEPRECATED_CODE

template <class ReturnType>
class Call_source_to_get
{
  public:
    typedef ReturnType     result_type;
    typedef Arity_tag< 1 >   Arity;

    template <class Cls>
    ReturnType
    operator()( const Cls& c) const
    { return c.source(); }
};

template <class ReturnType>
class Call_target_to_get
{
  public:
    typedef ReturnType     result_type;
    typedef Arity_tag< 1 >   Arity;

    template <class Cls>
    ReturnType
    operator()( const Cls& c) const
    { return c.target(); }
};

template <class ReturnType>
class Call_min_to_get
{
  public:
    typedef ReturnType     result_type;
    typedef Arity_tag< 1 >   Arity;

    template <class Cls>
    ReturnType
    operator()( const Cls& c) const
    { return c.min(); }
};

template <class ReturnType>
class Call_max_to_get
{
  public:
    typedef ReturnType     result_type;
    typedef Arity_tag< 1 >   Arity;

    template <class Cls>
    ReturnType
    operator()( const Cls& c) const
    { return c.max(); }
};

template <class ReturnType>
class Call_vertex_to_get
{
  public:
    typedef ReturnType     result_type;
    typedef Arity_tag< 2 >   Arity;

    template <class Cls>
    ReturnType
    operator()( const Cls& c, int i) const
    { return c.vertex(i); }
};

template <class ReturnType>
class Call_direction_to_get
{
  public:
    typedef ReturnType     result_type;
    typedef Arity_tag< 1 >   Arity;

    template <class Cls>
    ReturnType
    operator()( const Cls& c) const
    { return c.direction(); }
};

template <class ReturnType>
class Call_supporting_line_to_get
{
  public:
    typedef ReturnType     result_type;
    typedef Arity_tag< 1 >   Arity;

    template <class Cls>
    ReturnType
    operator()( const Cls& c) const
    { return c.supporting_line(); }
};

template <class ReturnType>
class Call_supporting_plane_to_get
{
  public:
    typedef ReturnType     result_type;
    typedef Arity_tag< 1 >   Arity;

    template <class Cls>
    ReturnType
    operator()( const Cls& c) const
    { return c.supporting_plane(); }
};

template <class ReturnType>
class Call_opposite_to_get
{
  public:
    typedef ReturnType     result_type;
    typedef Arity_tag< 1 >   Arity;

    template <class Cls>
    ReturnType
    operator()( const Cls& c) const
    { return c.opposite(); }
};


class Call_has_on
{
  public:
    typedef bool           result_type;
    typedef Arity_tag< 2 >   Arity;

    template <class Cls, class A1>
    bool
    operator()( const Cls& c, const A1& a1) const
    { return c.has_on(a1); }
};

class Call_collinear_has_on
{
  public:
    typedef bool           result_type;
    typedef Arity_tag< 2 >   Arity;

    template <class Cls, class A1>
    bool
    operator()( const Cls& c, const A1& a1) const
    { return c.collinear_has_on(a1); }
};


} // end namespace CGALi
CGAL_END_NAMESPACE

#endif // CGAL_KERNEL_FUNCTION_OBJECTS_H
