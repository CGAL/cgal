#line 709 "bits.lw"
#ifndef CGAL_I_FUNCTION_OBJECTS_H
#define CGAL_I_FUNCTION_OBJECTS_H

#line 756 "bits.lw"
CGAL_BEGIN_NAMESPACE
#line 762 "bits.lw"
namespace CGALi {

#line 122 "bits.lw"
template <class ToBeConstructed>
class Construct
{
  public:
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
    operator()( const A1& a1, const A2& a2, const A3& a3, const A4) const
    { return ToBeConstructed(a1,a2,a3,a4); }

    template <class A1, class A2, class A3, class A4, class A5> 
    ToBeConstructed
    operator()( const A1& a1, const A2& a2, const A3& a3, const A4, const A5) const
    { return ToBeConstructed(a1,a2,a3,a4,a5); }
};

#line 273 "bits.lw"
template <class ReturnType>
class Call_point_to_get
{
  public:
    template <class Cls>
    ReturnType
    operator()( const Cls& c) const
    { return c.point(); }

    template <class Cls>
    ReturnType
    operator()( const Cls& c, int i) const
    { return c.point(i); }
};
#line 290 "bits.lw"
template <class ReturnType>
class Call_second_point_to_get
{
  public:
    template <class Cls>
    ReturnType
    operator()( const Cls& c) const
    { return c.second_point(); }
};
#line 302 "bits.lw"
template <class ReturnType>
class Call_perpendicular_to_get
{
  public:
    template <class Cls>
    ReturnType
    operator()( const Cls& c) const
    { return c.perpendicular(); }

    template <class Cls, class A1>
    ReturnType
    operator()( const Cls& c, const A1& a1) const
    { return c.perpendicular(a1); }
};
#line 319 "bits.lw"
template <class Point>
class p_Midpoint
{
  public:
    Point
    operator()(const Point& p, const Point& q) const { return midpoint(p,q); }
};
#line 339 "bits.lw"
template <class Point>
class p_Circumcenter
{
  public:
    Point
    operator()(const Point& p, const Point& q, const Point& r) const
    { return circumcenter(p,q,r); }
};
#line 329 "bits.lw"
template <class Point, class Line>
class pl_Bisector
{
  public:
    Line
    operator()(const Point& p, const Point& q) const { return bisector(p,q); }
};
#line 365 "bits.lw"
class Intersect
{
  public:
    template <class T1, class T2>
    ::CGAL::Object
    operator()(const T1& t1, const T2& t2) const
    { return intersection( t1, t2); }
};
#line 376 "bits.lw"
template <class ReturnType>
class Call_y_at_x_to_get
{
  public:
    template <class Cls>
    ReturnType
    operator()( const Cls& c, class FT& x) const
    { return c.y_at_x(x); }
};
#line 388 "bits.lw"
template <class ReturnType>
class Call_squared_length_to_get
{
  public:
    template <class Cls>
    ReturnType
    operator()( const Cls& c) const
    { return c.squared_length(); }
};
#line 512 "bits.lw"
class Collinear
{
  public:
    template <class T>
    bool
    operator()(const T& p, const T& q, const T& r) const
    { return collinear(p,q,r); }
};
#line 523 "bits.lw"
class Side_of_oriented_circle
{
  public:
    template <class T>
    bool
    operator()(const T& p, const T& q, const T& r, const T& t) const
    { return side_of_oriented_circle(p,q,r,t); }
};

class Side_of_bounded_circle
{
  public:
    template <class T>
    bool
    operator()(const T& p, const T& q, const T& r, const T& t) const
    { return side_of_bounded_circle(p,q,r,t); }
};
#line 544 "bits.lw"
class Call_is_horizontal
{
  public:
    template <class Cls>
    bool
    operator()( const Cls& c) const
    { return c.is_horizontal(); }
};

class Call_is_vertical
{
  public:
    template <class Cls>
    bool
    operator()( const Cls& c) const
    { return c.is_vertical(); }
};
#line 564 "bits.lw"
class Call_is_degenerate
{
  public:
    template <class Cls>
    bool
    operator()( const Cls& c) const
    { return c.is_degenerate(); }
};
#line 575 "bits.lw"
class Call_has_on_bounded_side
{
  public:
    template <class Cls, class Arg>
    bool
    operator()( const Cls& c, const Arg& a) const
    { return c.has_on_bounded_side(a); }
};

class Call_has_on_unbounded_side
{
  public:
    template <class Cls, class Arg>
    bool
    operator()( const Cls& c, const Arg& a) const
    { return c.has_on_unbounded_side(a); }
};

class Call_has_on_boundary
{
  public:
    template <class Cls, class Arg>
    bool
    operator()( const Cls& c, const Arg& a) const
    { return c.has_on_boundary(a); }
};

class Call_has_on_positive_side
{
  public:
    template <class Cls, class Arg>
    bool
    operator()( const Cls& c, const Arg& a) const
    { return c.has_on_positive_side(a); }
};

class Call_has_on_negative_side
{
  public:
    template <class Cls, class Arg>
    bool
    operator()( const Cls& c, const Arg& a) const
    { return c.has_on_negative_side(a); }
};

class Call_oriented_side
{
  public:
    template <class Cls, class Arg>
    Oriented_side
    operator()( const Cls& c, const Arg& a) const
    { return c.oriented_side(a); }
};
#line 631 "bits.lw"
class Compare_x
{
  public:
    template <class T1, class T2>
    Comparison_result
    operator()( const T1& a1, const T2& a2) const
    { return compare_x(a1,a2); }
};

class Compare_y
{
  public:
    template <class T1, class T2>
    Comparison_result
    operator()( const T1& a1, const T2& a2) const
    { return compare_y(a1,a2); }
};

class Compare_y_at_x
{
  public:
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

#line 671 "bits.lw"
class Are_ordered_along_line
{
  public:
    template <class T>
    bool
    operator()(const T& p, const T& q, const T& r) const
    { return are_ordered_along_line(p,q,r); }
};

class Are_strictly_ordered_along_line
{
  public:
    template <class T>
    bool
    operator()(const T& p, const T& q, const T& r) const
    { return are_strictly_ordered_along_line(p,q,r); }
};

class Collinear_are_ordered_along_line
{
  public:
    template <class T>
    bool
    operator()(const T& p, const T& q, const T& r) const
    { return collinear_are_ordered_along_line(p,q,r); }
};

class Collinear_are_strictly_ordered_along_line
{
  public:
    template <class T>
    bool
    operator()(const T& p, const T& q, const T& r) const
    { return collinear_are_strictly_ordered_along_line(p,q,r); }
};

#line 765 "bits.lw"
} // end namespace CGALi
#line 759 "bits.lw"
CGAL_END_NAMESPACE

#line 736 "bits.lw"
#endif // CGAL_I_FUNCTION_OBJECTS_H
