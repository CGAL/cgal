// ======================================================================
//
// Copyright (c) The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-2.3-I-24 $
// release_date  : $CGAL_Date: 2000/12/29 $
//
// file          : include/CGAL/Snap_rounding_2_traits.h
// package       : arr (1.73)
// maintainer    : Eli Packer <elip@post.tau.ac.il>
// author(s)     : Eli Packer
// coordinator   : Tel-Aviv University (Dan Halperin <danha@post.tau.ac.il>)
//
// ======================================================================

#ifndef CGAL_SR_2_TRAITS_H
#define CGAL_SR_2_TRAITS_H

#include <CGAL/basic.h>
#include <map>

#include <CGAL/leda_integer.h> 
#include "../../include/CGAL/Snap_rounding_2_utility.h"

CGAL_BEGIN_NAMESPACE

template<class base_rep>
class Snap_rounding_traits : public base_rep {

typedef typename base_rep::FT                    NT;
typedef typename base_rep::Point_2               Point_2;
typedef typename base_rep::Segment_2             Segment_2;
typedef typename base_rep::Iso_rectangle_2       Iso_rectangle_2;

public:

/*! Functor
 */
class Snap_2 {
 public:
  void operator()(Point_2 p,NT pixel_size,NT &x,NT &y)
  {
    NT x_tmp = p.x() / pixel_size;
    NT y_tmp = p.y() / pixel_size;

    x = floor(x_tmp.to_double()) * pixel_size + pixel_size / 2.0;
    y = floor(y_tmp.to_double()) * pixel_size + pixel_size / 2.0;
  }
};

Snap_2 snap_2_object() const { return Snap_2(); }

/*! Functor
 */
class Integer_grid_point_2 {
 public:
  Point_2 operator()(Point_2 p,NT pixel_size)
  {
    NT x = (p.x() - pixel_size / 2.0) / pixel_size;
    NT y = (p.y() - pixel_size / 2.0) / pixel_size;
    Point_2 out_p(x,y);

    return(out_p);
  }
};

Integer_grid_point_2 integer_grid_point_2_object() const
    { return Integer_grid_point_2(); }

/*! Functor
 */
class Segment_direction_2 {
 public:
  double operator()(Segment_2 s)
  {
    double x1 = s.source().x().to_double();
    double y1 = s.source().y().to_double();
    double x2 = s.target().x().to_double();
    double y2 = s.target().y().to_double();

    return(atan((y2 - y1)/(x2 - x1)));
  }
};

Segment_direction_2 segment_direction_2_object() const
    {return Segment_direction_2(); }

/*! Functor
 */

class Bounding_box_of_minkowski_sum_2 {
 private:
  const Snap_rounding_traits<base_rep>* _gt;
  Bounding_box_of_minkowski_sum_2(
      const Snap_rounding_traits<base_rep>* gt) : _gt(gt) {}

  Point_2 small_x_point(Point_2 p1,Point_2 p2)
  {
    Comparison_result c = _gt->compare_x_2_object()(p1,p2);
    if(c == SMALLER)
      return(p1);
    else
      return(p2);
  }

  Point_2 big_x_point(Point_2 p1,Point_2 p2)
  {
    Comparison_result c = _gt->compare_x_2_object()(p1,p2);
    if(c == SMALLER)
      return(p2);
    else
      return(p1);
  }

  Point_2 small_y_point(Point_2 p1,Point_2 p2)
  {
    Comparison_result c = _gt->compare_y_2_object()(p1,p2);
    if(c == SMALLER)
      return(p1);
    else
      return(p2);
  }

  Point_2 big_y_point(Point_2 p1,Point_2 p2)
  {
    Comparison_result c = _gt->compare_y_2_object()(p1,p2);
    if(c == SMALLER)
      return(p2);
    else
      return(p1);
  }

 public:
  Iso_rectangle_2 operator()(Segment_2 s,
                             NT unit_squere,
		             NT angle)
  {
    Point_2 ms1,ms2,ms3,ms4,ms5,ms6;// minkowski sum points

    Comparison_result cx =  _gt->compare_x_2_object()(s.source(),s.target());
    NT x1 = s.source().x(),y1 = s.source().y(),x2 =
       s.target().x(),y2 = s.target().y();

    if(cx == SMALLER) {
      // we use unit_squere instead of unit_squere / 2 in order to
      // find tangency points which are not supported by kd-tree
      ms1 = Point_2(x1 - 0.6 * unit_squere,y1 - 0.6 * unit_squere);
      ms2 = Point_2(x1 - 0.6 * unit_squere,y1 + 0.6 * unit_squere);
      ms3 = Point_2(x1 + 0.6 * unit_squere,y1 - 0.6 * unit_squere);
      ms4 = Point_2(x2 + 0.6 * unit_squere,y2 - 0.6 * unit_squere);
      ms5 = Point_2(x2 + 0.6 * unit_squere,y2 + 0.6 * unit_squere);
      ms6 = Point_2(x2 - 0.6 * unit_squere,y2 + 0.6 * unit_squere);
    } else {
      // we use unit_squere instead of unit_squere / 2 in order to
      // find tangency points which are not supported by kd-tree
      ms1 = Point_2(x1 + 0.6 * unit_squere,y1 - 0.6 * unit_squere);
      ms2 = Point_2(x1 - 0.6 * unit_squere,y1 - 0.6 * unit_squere);
      ms3 = Point_2(x1 + 0.6 * unit_squere,y1 + 0.6 * unit_squere);
      ms4 = Point_2(x2 + 0.6 * unit_squere,y2 + 0.6 * unit_squere);
      ms5 = Point_2(x2 - 0.6 * unit_squere,y2 + 0.6 * unit_squere);
      ms6 = Point_2(x2 - 0.6 * unit_squere,y2 - 0.6 * unit_squere);
    }

    static Snap_rounding_rotation<base_rep> r;
    r(ms1,angle);
    r(ms2,angle);
    r(ms3,angle);
    r(ms4,angle);
    r(ms5,angle);
    r(ms6,angle);

    // query
    Point_2 point_left,point_right,point_bot,point_top;

    point_left = small_x_point(ms1,ms2);
    point_left = small_x_point(point_left,ms3);
    point_left = small_x_point(point_left,ms4);
    point_left = small_x_point(point_left,ms5);
    point_left = small_x_point(point_left,ms6);

    point_right = big_x_point(ms1,ms2);
    point_right = big_x_point(point_right,ms3);
    point_right = big_x_point(point_right,ms4);
    point_right = big_x_point(point_right,ms5);
    point_right = big_x_point(point_right,ms6);

    point_bot = small_y_point(ms1,ms2);
    point_bot = small_y_point(point_bot,ms3);
    point_bot = small_y_point(point_bot,ms4);
    point_bot = small_y_point(point_bot,ms5);
    point_bot = small_y_point(point_bot,ms6);

    point_top = big_y_point(ms1,ms2);
    point_top = big_y_point(point_top,ms3);
    point_top = big_y_point(point_top,ms4);
    point_top = big_y_point(point_top,ms5);
    point_top = big_y_point(point_top,ms6);

    Iso_rectangle_2 rec(point_left,point_right,point_bot,point_top);

    return(rec);
  }

  friend class Snap_rounding_traits<base_rep>;
};

Bounding_box_of_minkowski_sum_2 bounding_box_of_minkowski_sum_2_object() const
    {return Bounding_box_of_minkowski_sum_2(this); }

};

CGAL_END_NAMESPACE

#endif // CGAL_ISR_2_TRAITS_H
