// Copyright (c) 1997  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$
// $Name$
//
// author(s)     : Eli Packer <elip@post.tau.ac.il>
#ifndef CGAL_SR_2_TRAITS_H
#define CGAL_SR_2_TRAITS_H

#include <CGAL/basic.h>
#include <map>
#include <list>

#include <CGAL/Arr_segment_cached_traits_2.h>

CGAL_BEGIN_NAMESPACE

template<class base_rep>
class Snap_rounding_traits_2 :
    public CGAL::Arr_segment_cached_traits_2<base_rep> {

public: // otherwise Segment_data cannot access the types
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

    x = NT(floor(to_double(x_tmp))) * pixel_size + pixel_size / NT(2.0);
    y = NT(floor(to_double(y_tmp))) * pixel_size + pixel_size / NT(2.0);
    //    x = floor(x_tmp.to_double()) * pixel_size + pixel_size / 2.0;
    //    y = floor(y_tmp.to_double()) * pixel_size + pixel_size / 2.0;
  }
};

Snap_2 snap_2_object() const { return Snap_2(); }

/*! Functor
 */
class Integer_grid_point_2 {
 public:
  Point_2 operator()(Point_2 p,NT pixel_size)
  {
    NT x = (p.x() - pixel_size / NT(2.0)) / pixel_size;
    NT y = (p.y() - pixel_size / NT(2.0)) / pixel_size;
    Point_2 out_p(x,y);

    return(out_p);
  }
};

Integer_grid_point_2 integer_grid_point_2_object() const
    { return Integer_grid_point_2(); }

/*! Functor
 */

class Minkowski_sum_with_pixel_2 {
 private:
  const Snap_rounding_traits_2<base_rep>* _gt;
  Minkowski_sum_with_pixel_2(
      const Snap_rounding_traits_2<base_rep>* gt) : _gt(gt) {}

 public:
  void operator()(std::list<Point_2>& points_list,
                             Segment_2 s,
                             NT unit_square)
  {
     Comparison_result cx =  _gt->compare_x_2_object()(s.source(),s.target());
     Comparison_result cy =  _gt->compare_y_2_object()(s.source(),s.target());
     NT x1 = s.source().x(),y1 = s.source().y(),x2 =
       s.target().x(),y2 = s.target().y();
     Point_2 ms1,ms2,ms3,ms4,ms5,ms6;// minkowski sum points

     if(cx == SMALLER) {
       if(cy == SMALLER) {
         // we use unit_square instead of unit_square / 2 in order to
         // find tangency points which are not supported by kd-tree
         ms1 = Point_2(x1 - unit_square,y1 - unit_square);
         ms2 = Point_2(x1 - unit_square,y1 + unit_square);
         ms3 = Point_2(x1 + unit_square,y1 - unit_square);
         ms4 = Point_2(x2 + unit_square,y2 - unit_square);
         ms5 = Point_2(x2 + unit_square,y2 + unit_square);
         ms6 = Point_2(x2 - unit_square,y2 + unit_square);
       } else {
         ms1 = Point_2(x1 - unit_square,y1 - unit_square);
         ms2 = Point_2(x1 - unit_square,y1 + unit_square);
         ms3 = Point_2(x1 + unit_square,y1 + unit_square);
         ms4 = Point_2(x2 + unit_square,y2 - unit_square);
         ms5 = Point_2(x2 + unit_square,y2 + unit_square);
         ms6 = Point_2(x2 - unit_square,y2 - unit_square);
       }
     } else {
       if(cy == SMALLER) {
         ms1 = Point_2(x1 + unit_square,y1 - unit_square);
         ms2 = Point_2(x1 + unit_square,y1 + unit_square);
         ms3 = Point_2(x1 - unit_square,y1 - unit_square);
         ms4 = Point_2(x2 + unit_square,y2 + unit_square);
         ms5 = Point_2(x2 - unit_square,y2 + unit_square);
         ms6 = Point_2(x2 - unit_square,y2 - unit_square);
       } else {
         ms1 = Point_2(x1 + unit_square,y1 - unit_square);
         ms2 = Point_2(x1 + unit_square,y1 + unit_square);
         ms3 = Point_2(x1 - unit_square,y1 + unit_square);
         ms4 = Point_2(x2 + unit_square,y2 - unit_square);
         ms5 = Point_2(x2 - unit_square,y2 - unit_square);
         ms6 = Point_2(x2 - unit_square,y2 + unit_square);
       }
     }

     points_list.push_back(ms1);
     points_list.push_back(ms2);
     points_list.push_back(ms3);
     points_list.push_back(ms4);
     points_list.push_back(ms5);
     points_list.push_back(ms6);
  }

  friend class Snap_rounding_traits_2<base_rep>;
};

Minkowski_sum_with_pixel_2 minkowski_sum_with_pixel_2_object() const
    {return Minkowski_sum_with_pixel_2(this); }

};

CGAL_END_NAMESPACE

#endif // CGAL_ISR_2_TRAITS_H
