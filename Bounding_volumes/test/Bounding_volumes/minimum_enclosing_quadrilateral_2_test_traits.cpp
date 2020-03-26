// Copyright (c) 1999-2003  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch> and
//                 Emo Welzl <emo@inf.ethz.ch>

#define CGAL_OPTIMISATION_NO_PRECONDITIONS
#include <CGAL/min_quadrilateral_2.h>
#include <CGAL/predicates/kernel_ftC2.h>
#include <CGAL/constructions/kernel_ftC2.h>
#include <functional>
#include <vector>
#include <iterator>

struct MyTraits {
  struct Point_2      { double xc, yc; };
  struct Vector_2     { double xx, yy; };
  struct Direction_2  { double xd, yd; };
  struct Line_2       { double aa, bb, cc; };
  struct Rectangle_2  { Point_2 pp1, pp2, pp3, pp4; Direction_2 dd; };
  struct Parallelogram_2 {
    Point_2 pp1, pp2, pp3, pp4;
    Direction_2 dd1, dd2;
  };
  struct Strip_2      { Point_2 pp1, pp2; Direction_2 dd; };
  struct Equal_2
  : public CGAL::cpp98::binary_function<Point_2,Point_2,bool>
  {
    bool operator()(const Point_2& p, const Point_2& q) const
    { return p.xc == q.xc && p.yc == q.yc; }
  };
  struct Less_xy_2
  : public CGAL::cpp98::binary_function<Point_2,Point_2,bool>
  {
    bool operator()(const Point_2& p, const Point_2& q) const
    { return p.xc < q.xc || (p.xc == q.xc && p.yc < q.yc); }
  };
  struct Less_yx_2
  : public CGAL::cpp98::binary_function<Point_2,Point_2,bool>
  {
    bool operator()(const Point_2& p, const Point_2& q) const
    { return p.yc < q.yc || (p.yc == q.yc && p.xc < q.xc); }
  };
  struct Orientation_2 {
    typedef CGAL::Orientation result_type;

    result_type
    operator()(const Point_2& p, const Point_2& q, const Point_2& r) const {
      return CGAL::orientationC2(p.xc, p.yc, q.xc, q.yc, r.xc, r.yc);
    }
  };
  struct Has_on_negative_side_2
  : public CGAL::cpp98::binary_function<Line_2,Point_2,bool>
  {
    bool operator()(const Line_2& l, const Point_2& p) const {
      return
      CGAL::side_of_oriented_lineC2(l.aa, l.bb, l.cc, p.xc, p.yc)
      ==
      CGAL::ON_NEGATIVE_SIDE;
    }
  };
  struct Compare_angle_with_x_axis_2
  : public CGAL::cpp98::binary_function<Direction_2,Direction_2,CGAL::Comparison_result>
  {
    CGAL::Comparison_result
    operator()(const Direction_2& d, const Direction_2& e) const
    { return CGAL::compare_angle_with_x_axisC2(d.xd, d.yd, e.xd, e.yd); }
  };
  struct Area_less_rectangle_2
  : public CGAL::cpp98::binary_function<Rectangle_2,Rectangle_2,bool>
  {
    bool operator()(const Rectangle_2& d, const Rectangle_2& e) const
    {
      using CGAL::square;
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
  : public CGAL::cpp98::binary_function<Parallelogram_2,Parallelogram_2,bool>
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
  : public CGAL::cpp98::binary_function<Strip_2,Strip_2,bool>
  {
    bool operator()(const Strip_2& d, const Strip_2& e) const
    {
      using CGAL::square;
      return
      (d.dd.xd * (d.pp2.yc - d.pp1.yc) + d.dd.yd * (d.pp1.xc - d.pp2.xc)) *
        (CGAL_NTS square(e.dd.xd) + CGAL_NTS square(e.dd.yd)) <
      (e.dd.xd * (e.pp2.yc - e.pp1.yc) + e.dd.yd * (e.pp1.xc - e.pp2.xc)) *
        (CGAL_NTS square(d.dd.xd) + CGAL_NTS square(d.dd.yd));
    }
  };
  struct Construct_vector_2
  : public CGAL::cpp98::binary_function<Point_2,Point_2,Vector_2>
  {
    Vector_2 operator()(const Point_2& p, const Point_2& q) const
    {
      Vector_2 v;
      v.xx = q.xc - p.xc;
      v.yy = q.yc - p.yc;
      return v;
    }
  };
  struct Construct_vector_from_direction_2
  : public CGAL::cpp98::unary_function<Direction_2,Vector_2>
  {
    Vector_2 operator()(const Direction_2& d) const
    {
      Vector_2 v;
      v.xx = d.xd;
      v.yy = d.yd;
      return v;
    }
  };
  struct Construct_perpendicular_vector_2
  : public CGAL::cpp98::binary_function<Vector_2,CGAL::Orientation,Vector_2>
  {
    Vector_2 operator()(const Vector_2& v, CGAL::Orientation o) const
    {
      Vector_2 vm;
      if (o == CGAL::CLOCKWISE) {
        vm.xx = v.yy;
        vm.yy = -v.xx;
      } else {
        vm.xx = -v.yy;
        vm.yy = v.xx;
      }
      return vm;
    }
  };
  struct Construct_direction_2
  : public CGAL::cpp98::unary_function<Vector_2,Direction_2>
  {
    Direction_2 operator()(const Vector_2& v) const
    {
      Direction_2 d;
      d.xd = v.xx;
      d.yd = v.yy;
      return d;
    }
  };
  struct Construct_opposite_direction_2
  : public CGAL::cpp98::unary_function<Direction_2,Direction_2>
  {
    Direction_2 operator()(const Direction_2& d) const
    {
      Direction_2 dm;
      dm.xd = -d.xd;
      dm.yd = -d.yd;
      return dm;
    }
  };
  struct Construct_line_2
  : public CGAL::cpp98::binary_function<Point_2,Direction_2,Line_2>
  {
    Line_2 operator()(const Point_2& p, const Direction_2& d) const
    {
      Line_2 l;
      CGAL::line_from_point_directionC2(p.xc, p.yc, d.xd, d.yd,
                                        l.aa, l.bb, l.cc);
      return l;
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


  friend struct Data;
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

