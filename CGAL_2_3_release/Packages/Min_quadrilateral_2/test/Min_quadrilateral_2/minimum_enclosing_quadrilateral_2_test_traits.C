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
// file          : minimum_enclosing_quadrilateral_2_test_traits.C
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
// Test Program: Computing minimum enclosing quadrilaterals
// ============================================================================

#include <CGAL/min_quadrilateral_2.h>
#include <CGAL/predicates/kernel_ftC2.h>
#include <functional>
#include <vector>
#include <iterator>

struct MyTraits {
  struct Point_2      { double xc, yc; };
  struct Direction_2  { double xd, yd; };
  struct Line_2       { Point_2 p; Direction_2 d; };
  struct Rectangle_2  { Point_2 pp1, pp2, pp3, pp4; Direction_2 dd; };
  struct Parallelogram_2 {
    Point_2 pp1, pp2, pp3, pp4;
    Direction_2 dd1, dd2;
  };
  struct Strip_2      { Point_2 pp1, pp2; Direction_2 dd; };
private:
  struct Equal_2 : public std::binary_function< Point_2, Point_2, bool > {
    bool operator()(const Point_2& p, const Point_2& q) const
    { return p.xc == q.xc && p.yc == q.yc; }
  };
  struct Less_x_2
  : public std::binary_function< Point_2, Point_2, bool >
  {
    bool operator()(const Point_2& p, const Point_2& q) const
    { return p.xc < q.xc; }
  };
  struct Less_y_2
  : public std::binary_function< Point_2, Point_2, bool >
  {
    bool operator()(const Point_2& p, const Point_2& q) const
    { return p.yc < q.yc; }
  };
  struct Greater_x_2
  : public std::binary_function< Point_2, Point_2, bool >
  {
    bool operator()(const Point_2& p, const Point_2& q) const
    { return p.xc > q.xc; }
  };
  struct Greater_y_2
  : public std::binary_function< Point_2, Point_2, bool >
  {
    bool operator()(const Point_2& p, const Point_2& q) const
    { return p.yc > q.yc; }
  };
  struct Right_of_implicit_line_2 {
    bool operator()(const Point_2& p,
                    const Point_2& q, const Direction_2& d) const
    { return d.xd * (q.yc - p.yc) < d.yd * (q.xc - p.xc); }
  };
  struct Less_rotate_ccw_2 {
    bool operator()(const Direction_2& d, const Direction_2& e) const
    { return CGAL::SMALLER ==
      CGAL::compare_angle_with_x_axisC2(d.xd, d.yd, e.xd, e.yd); }
  };
  struct Area_less_rectangle_2
  : public std::binary_function< Rectangle_2, Rectangle_2, bool >
  {
    bool operator()(const Rectangle_2& d, const Rectangle_2& e) const
    {
      return
      (d.dd.xd * (d.pp3.yc - d.pp1.yc) +
       d.dd.yd * (d.pp1.xc - d.pp3.xc)) *
      (d.dd.yd * (d.pp2.yc - d.pp4.yc) +
       d.dd.xd * (d.pp2.xc - d.pp4.xc)) *
      (CGAL_NTS square(e.dd.xd) + CGAL_NTS square(e.dd.yd)) <
      (e.dd.xd * (e.pp3.yc - e.pp1.yc) +
       e.dd.yd * (e.pp1.xc - e.pp3.xc)) *
      (e.dd.yd * (e.pp2.yc - e.pp4.yc) +
       e.dd.xd * (e.pp2.xc - e.pp4.xc)) *
      (CGAL_NTS square(d.dd.xd) + CGAL_NTS square(d.dd.yd));
    }
  };
  struct Area_less_parallelogram_2
  : public std::binary_function< Parallelogram_2, Parallelogram_2, bool >
  {
    bool operator()(const Parallelogram_2& d,
                    const Parallelogram_2& e) const
    {
      return
      (d.dd1.xd * (d.pp3.yc - d.pp1.yc) +
       d.dd1.yd * (d.pp1.xc - d.pp3.xc)) *
      (d.dd2.yd * (d.pp2.yc - d.pp4.yc) +
       d.dd2.xd * (d.pp4.xc - d.pp2.xc)) *
      (e.dd1.xd * e.dd2.yd - e.dd1.yd * e.dd2.xd) <
      (e.dd1.xd * (e.pp3.yc - e.pp1.yc) +
       e.dd1.yd * (e.pp1.xc - e.pp3.xc)) *
      (e.dd2.yd * (e.pp2.yc - e.pp4.yc) +
       e.dd2.xd * (e.pp4.xc - e.pp2.xc)) *
      (d.dd1.xd * d.dd2.yd - d.dd1.yd * d.dd2.xd);
    }
  };
  struct Width_less_strip_2
  : public std::binary_function< Strip_2, Strip_2, bool >
  {
    bool operator()(const Strip_2& d, const Strip_2& e) const
    {
      return
      (d.dd.xd * (d.pp2.yc - d.pp1.yc) + d.dd.yd * (d.pp1.xc - d.pp2.xc)) *
        (CGAL_NTS square(e.dd.xd) + CGAL_NTS square(e.dd.yd)) <
      (e.dd.xd * (e.pp2.yc - e.pp1.yc) + e.dd.yd * (e.pp1.xc - e.pp2.xc)) *
        (CGAL_NTS square(d.dd.xd) + CGAL_NTS square(d.dd.yd));
    }
  };
  struct Rotate_direction_by_multiple_of_pi_2
  : public std::binary_function< Direction_2, int, Direction_2 >
  {
    Direction_2 operator()(const Direction_2& d, int i) const
    {
      Direction_2 n;
      if (i == 0)
        return d;
      if (i == 1) {
        n.xd = d.yd;
        n.yd = -d.xd;
        return n;
      }
      if (i == 2) {
        n.xd = -d.xd;
        n.yd = -d.yd;
        return n;
      }
      n.xd = -d.yd;
      n.yd = d.xd;
      return n;
    }
  };
  struct Construct_direction_2
  : public std::binary_function< Point_2, Point_2, Direction_2 >
  {
    Direction_2 operator()(const Point_2& p, const Point_2& q) const
    {
      Direction_2 n;
      n.xd = q.xc - p.xc;
      n.yd = q.yc - p.yc;
      return n;
    }
  };
  struct Construct_rectangle_2 {
    Rectangle_2 operator()(const Point_2& p1,
                           const Direction_2& d,
                           const Point_2& p2,
                           const Point_2& p3,
                           const Point_2& p4) const
    {
      Rectangle_2 n;
      n.pp1 = p1, n.pp2 = p2, n.pp3 = p3, n.pp4 = p4, n.dd = d;
      return n;
    }
  };
  struct Construct_parallelogram_2 {
    Parallelogram_2 operator()(const Point_2& p1,
                               const Direction_2& d1,
                               const Point_2& p2,
                               const Direction_2& d2,
                               const Point_2& p3,
                               const Point_2& p4) const
    {
      Parallelogram_2 n;
      n.pp1 = p1, n.pp2 = p2, n.pp3 = p3, n.pp4 = p4, n.dd1 = d1, n.dd2 = d2;
      return n;
    }
  };
  struct Construct_strip_2 {
    Strip_2 operator()(const Point_2& p1,
                       const Direction_2& d,
                       const Point_2& p2) const
    {
      Strip_2 n;
      n.pp1 = p1, n.pp2 = p2, n.dd = d;
      return n;
    }
  };
public:
  template < class OutputIterator >
  OutputIterator
  copy_rectangle_vertices_2(const Rectangle_2&, OutputIterator o) const
  {
    // no output, this is just a test
    return o;
  }
  template < class OutputIterator >
  OutputIterator
  copy_parallelogram_vertices_2(const Parallelogram_2&,
                                OutputIterator o) const
  {
    // no output, this is just a test
    return o;
  }
  template < class OutputIterator >
  OutputIterator
  copy_strip_lines_2(const Strip_2&, OutputIterator o) const
  {
    // no output, this is just a test
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
  friend class Data;
};

struct Data {
  typedef MyTraits::Point_2 P;
  Data() {
    P p;
    p.xc = 0, p.yc = 2;
    i.push_back(p);
    p.xc = 2;
    i.push_back(p);
    p.yc = 4;
    i.push_back(p);
  }
  std::vector< P > i, o;
};

int main()
{
  MyTraits t;
  Data d;
  CGAL::min_rectangle_2(
    d.i.begin(), d.i.end(), std::back_inserter(d.o), t);
  CGAL::min_parallelogram_2(
    d.i.begin(), d.i.end(), std::back_inserter(d.o), t);
  CGAL::min_strip_2(
    d.i.begin(), d.i.end(), std::back_inserter(d.o), t);
  return 0;
}
// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

