// ============================================================================
//
// Copyright (c) 1999, 2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision $
// release_date  : $CGAL_Date $
//
// file          : Min_quadrilateral_traits_2.h
// chapter       : $CGAL_Chapter: Geometric Optimisation $
// package       : $CGAL_Package: Min_quadrilaterals $
// source        : oops.aw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch> and
//                 Emo Welzl <emo@inf.ethz.ch>
//
// maintainer    : Michael Hoffmann <hoffmann@inf.ethz.ch>
// coordinator   : ETH
//
// Default traits for computing minimum enclosing four-gons
// ============================================================================

#if ! (CGAL_MIN_QUADRILATERAL_TRAITS_2_H)
#define CGAL_MIN_QUADRILATERAL_TRAITS_2_H 1

#include <CGAL/basic.h>
#include <CGAL/Optimisation/assertions.h>
#include <CGAL/Point_2.h>
#include <CGAL/Direction_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/intersections.h>
#include <vector>

CGAL_BEGIN_NAMESPACE

template < class R >
struct Min_rectangle_2 {
  typedef typename R::Point_2     Point_2;
  typedef typename R::Direction_2 Direction_2;

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

template < class R >
struct Min_parallelogram_2 {
  typedef typename R::Point_2     Point_2;
  typedef typename R::Direction_2 Direction_2;

  Min_parallelogram_2(const Point_2& q1, const Direction_2& e1,
                      const Point_2& q2, const Direction_2& e2,
                      const Point_2& q3, const Point_2& q4)
  : p1(q1), p2(q2), p3(q3), p4(q4), d1(e1), d2(e2)
  {}

  Point_2 p1, p2, p3, p4;
  Direction_2 d1, d2;
};

template < class R >
std::ostream& operator<<(std::ostream& o, const Min_parallelogram_2< R >& p)
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

template < class R >
struct Min_strip_2 {
  typedef typename R::Point_2     Point_2;
  typedef typename R::Direction_2 Direction_2;

  Min_strip_2(const Point_2& q1, const Direction_2& e, const Point_2& q2)
  : p1(q1), p2(q2), d(e)
  {}

  Point_2 p1, p2;
  Direction_2 d;
};

template < class R_ >
struct Min_quadrilateral_default_traits_2 {
  // Types
  typedef R_                                 R;
  typedef typename R::RT                     RT;
  typedef typename R::Point_2                Point_2;
  typedef typename R::Direction_2            Direction_2;
  typedef typename R::Line_2                 Line_2;
  typedef Min_rectangle_2< R >               Rectangle_2;
  typedef Min_parallelogram_2< R >           Parallelogram_2;
  typedef Min_strip_2< R >                   Strip_2;

  // Predicates
  typedef typename R::Equal_2                Equal_2;
  // std::equal_to< Point_2 >                Equal_2;
  typedef typename R::Less_x_2               Less_x_2;
  typedef typename R::Less_y_2               Less_y_2;
      struct Greater_x_2
      : public CGAL_STD::binary_function< Point_2, Point_2, bool >
      {
        bool
        operator()(const Point_2& p, const Point_2& q) const
        { return p.x() > q.x(); }
      };
  struct Greater_y_2
  : public CGAL_STD::binary_function< Point_2, Point_2, bool >
  {
    bool
    operator()(const Point_2& p, const Point_2& q) const
    { return p.y() > q.y(); }
  };
  struct Right_of_implicit_line_2
  {
    bool
    operator()(const Point_2& p,
               const Point_2& q, const Direction_2& d) const
    {
      typedef CGAL::Line_2< R >  Line_2;
      return Line_2(q, d).has_on_negative_side(p);
    }
  };
  struct Less_rotate_ccw_2
  : public CGAL_STD::binary_function< Direction_2, Direction_2, bool >
  {
    bool
    operator()(const Direction_2& d, const Direction_2& e) const
    { return d.dy() * e.dx() < d.dx() * e.dy(); }
  };
  struct Area_less_rectangle_2
  : public CGAL_STD::binary_function< Rectangle_2, Rectangle_2, bool >
  {
    RT
    area_numerator(const Rectangle_2& r, Cartesian_tag) const
    {
      return
      (r.d.dx() * (r.p3.y() - r.p1.y()) + r.d.dy() * (r.p1.x() - r.p3.x())) *
        (-r.d.dy() * (r.p4.y() - r.p2.y()) + r.d.dx() * (r.p2.x() - r.p4.x()));
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
      (-r.d.dy() * (r.p4.hy() * r.p2.hw() - r.p2.hy() * r.p4.hw()) +
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
      typename R::Rep_tag tag;
      return area_numerator(p, tag) * area_denominator(q, tag) <
        area_denominator(p, tag) * area_numerator(q, tag);
    }
  };
  struct Area_less_parallelogram_2
  : public CGAL_STD::binary_function< Parallelogram_2,
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
      typename R::Rep_tag tag;
      return area_numerator(p, tag) * area_denominator(q, tag) <
        area_denominator(p, tag) * area_numerator(q, tag);
    }
  };
  struct Width_less_strip_2
  : public CGAL_STD::binary_function< Strip_2, Strip_2, bool >
  {
    RT
    width_numerator(const Strip_2& r, Cartesian_tag) const
    {
      return
        r.d.dx() * (r.p2.y() - r.p1.y()) +
        r.d.dy() * (r.p1.x() - r.p2.x());
    }
    
    RT
    width_denominator(const Strip_2& r, Cartesian_tag) const
    { return CGAL_NTS square(r.d.dx()) + CGAL_NTS square(r.d.dy()); }
    
    RT
    width_numerator(const Strip_2& r, Homogeneous_tag) const
    {
      return
        r.d.dx() * (r.p2.hy() * r.p1.hw() - r.p1.hy() * r.p2.hw()) +
        r.d.dy() * (r.p1.hx() * r.p2.hw() - r.p2.hx() * r.p1.hw());
    }
    
    RT
    width_denominator(const Strip_2& r, Homogeneous_tag) const {
      return r.p1.hw() * r.p2.hw() *
        (CGAL_NTS square(r.d.dx()) + CGAL_NTS square(r.d.dy()));
    }
  
    bool
    operator()(const Strip_2& p, const Strip_2& q) const
    {
      typename R::Rep_tag tag;
      return width_numerator(p, tag) * width_denominator(q, tag) <
        width_denominator(p, tag) * width_numerator(q, tag);
    }
  };

  // Constructions
  struct Construct_direction_2
  : public CGAL_STD::binary_function< Point_2, Point_2, Direction_2 >
  {
    Direction_2
    operator()(const Point_2& p, const Point_2& q) const
    { return (q - p).direction(); }
  };
  struct Rotate_direction_by_multiple_of_pi_2
  : public CGAL_STD::binary_function< Direction_2, int, Direction_2 >
  {
    Direction_2
    operator()(const Direction_2& d, int i) const
    {
      typedef typename R::Vector_2 Vector_2;
      CGAL_precondition(i >= 0 && i < 4);
      if (i == 0)
        return d;
      if (i == 1)
        return Direction_2(Vector_2(d.to_vector().hy(),
                                    -d.to_vector().hx(),
                                    d.to_vector().hw()));
      if (i == 2)
        return -d;
      return Direction_2(Vector_2(-d.to_vector().hy(),
                                  d.to_vector().hx(),
                                  d.to_vector().hw()));
    }
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
        rotate_direction_by_multiple_of_pi_2_object()(r.d, 1),
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
    typedef typename R::Line_2  Line_2;
    Point_2 tmp;
    Line_2  tmpl;
    Object  tmpo;
  
    tmpo = intersection(Line_2(r.p1, r.d1), Line_2(r.p2, r.d2));
    if (assign(tmp, tmpo)) {
      *o++ = tmp;
    } else {
      CGAL_optimisation_assertion_code(bool test1 =)
      assign(tmpl, tmpo);
      CGAL_optimisation_assertion(test1);
      *o++ = r.p1;
    }
    tmpo = intersection(Line_2(r.p3, r.d1), Line_2(r.p2, r.d2));
    if (assign(tmp, tmpo)) {
      *o++ = tmp;
    } else {
      CGAL_optimisation_assertion_code(bool test1 =)
      assign(tmpl, tmpo);
      CGAL_optimisation_assertion(test1);
      *o++ = r.p2;
    }
    tmpo = intersection(Line_2(r.p3, r.d1), Line_2(r.p4, r.d2));
    if (assign(tmp, tmpo)) {
      *o++ = tmp;
    } else {
      CGAL_optimisation_assertion_code(bool test1 =)
      assign(tmpl, tmpo);
      CGAL_optimisation_assertion(test1);
      *o++ = r.p3;
    }
    tmpo = intersection(Line_2(r.p1, r.d1), Line_2(r.p4, r.d2));
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
    *o++ = Line_2(r.p1, r.d);
    *o++ = Line_2(r.p2, r.d);
    return o;
  } 

  Equal_2     equal_2_object()     const { return Equal_2(); }
  Less_x_2    less_x_2_object()    const { return Less_x_2(); }
  Less_y_2    less_y_2_object()    const { return Less_y_2(); }
  Greater_x_2 greater_x_2_object() const { return Greater_x_2(); }
  Greater_y_2 greater_y_2_object() const { return Greater_y_2(); }
  
  Right_of_implicit_line_2 right_of_implicit_line_2_object() const
  { return Right_of_implicit_line_2(); }
  
  Less_rotate_ccw_2 less_rotate_ccw_2_object() const
  { return Less_rotate_ccw_2(); }
  
  Area_less_rectangle_2 area_less_rectangle_2_object() const
  { return Area_less_rectangle_2(); }
  
  Area_less_parallelogram_2 area_less_parallelogram_2_object() const
  { return Area_less_parallelogram_2(); }
  
  Width_less_strip_2 width_less_strip_2_object() const
  { return Width_less_strip_2(); }
  
  Construct_direction_2 construct_direction_2_object() const
  { return Construct_direction_2(); }
  
  Rotate_direction_by_multiple_of_pi_2
  rotate_direction_by_multiple_of_pi_2_object() const
  { return Rotate_direction_by_multiple_of_pi_2(); }
  
  Construct_rectangle_2 construct_rectangle_2_object() const
  { return Construct_rectangle_2(); }
  
  Construct_parallelogram_2 construct_parallelogram_2_object() const
  { return Construct_parallelogram_2(); }
  
  Construct_strip_2 construct_strip_2_object() const
  { return Construct_strip_2(); }
};


CGAL_END_NAMESPACE

#endif // ! (CGAL_MIN_QUADRILATERAL_TRAITS_2_H)

// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

