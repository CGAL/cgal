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
#include <list>

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

  void minkowski_sum(std::list<Point_2>& points_list,Segment_2 s,NT unit_squere)
  {
     Comparison_result cx =  _gt->compare_x_2_object()(s.source(),s.target());
     NT x1 = s.source().x(),y1 = s.source().y(),x2 =
       s.target().x(),y2 = s.target().y();
     Point_2 ms1,ms2,ms3,ms4,ms5,ms6;// minkowski sum points

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

     points_list.push_back(ms1);
     points_list.push_back(ms2);
     points_list.push_back(ms3);
     points_list.push_back(ms4);
     points_list.push_back(ms5);
     points_list.push_back(ms6);
  }

 public:
  Iso_rectangle_2 operator()(Segment_2 s,
                             NT unit_squere,
		             NT angle)
  {
    std::list<Point_2> points_list;

    minkowski_sum(points_list,s,unit_squere);

    static Snap_rounding_rotation<base_rep> r;

    typename std::list<Point_2>::iterator iter;

    for(iter = points_list.begin();iter != points_list.end();++iter)
      r(*iter,angle);

    // query
    iter = points_list.begin();
    Point_2 point_left,point_right,point_bot,point_top;
    point_left = point_right = point_bot = point_top = *iter;
    for(++iter;iter != points_list.end();++iter) {
      point_left = small_x_point(point_left,*iter);
      point_right = big_x_point(point_right,*iter);
      point_bot = small_y_point(point_bot,*iter);
      point_top = big_y_point(point_top,*iter);
    }

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
