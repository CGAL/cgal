// Copyright (c) 2001,2009,2014  Tel-Aviv University (Israel), Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// author(s)     : Eli Packer <elip@post.tau.ac.il>,
//                 Waqar Khan <wkhan@mpi-inf.mpg.de>

#ifndef CGAL_SNAP_ROUNDING_2_TRAITS_H
#define CGAL_SNAP_ROUNDING_2_TRAITS_H

#include <CGAL/basic.h>
#include <map>
#include <list>

#include <CGAL/Arr_segment_traits_2.h>

namespace CGAL {

template<class Base_kernel>
class Snap_rounding_traits_2 :
  public CGAL::Arr_segment_traits_2<Base_kernel> {

public: // otherwise Segment_data cannot access the types
  typedef typename Base_kernel::FT                          NT;
  typedef typename Base_kernel::FT                          FT;
  typedef typename Base_kernel::Point_2                     Point_2;
  typedef typename Base_kernel::Segment_2                   Segment_2;
  typedef typename Base_kernel::Iso_rectangle_2             Iso_rectangle_2;
  typedef typename Base_kernel::Vector_2                    Vector_2;
  typedef typename Base_kernel::Line_2                      Line_2;
  typedef typename Base_kernel::Aff_transformation_2        Aff_transformation_2;
  typedef typename Base_kernel::Direction_2                 Direction_2;
  typedef typename Base_kernel::Construct_vertex_2          Construct_vertex_2 ;
  typedef typename Base_kernel::Construct_segment_2         Construct_segment_2 ;
  typedef typename Base_kernel::Construct_iso_rectangle_2   Construct_iso_rectangle_2;
  typedef typename Base_kernel::Compare_y_2                 Compare_y_2;

  typedef typename Base_kernel::Construct_min_vertex_2                    Construct_min_vertex_2;
  typedef typename Base_kernel::Construct_max_vertex_2                    Construct_max_vertex_2;
  typedef typename Base_kernel::Cartesian_const_iterator_2                Cartesian_const_iterator_2;
  typedef typename Base_kernel::Construct_cartesian_const_iterator_2      Construct_cartesian_const_iterator_2;

  typedef CGAL::Arr_segment_traits_2<Base_kernel>                         Base_traits;
  typedef typename Base_traits::Compare_x_2                               Compare_x_2;
  typedef CGAL::To_double<NT>                                             To_double;

public:
  /*! Functor */
  class Snap_2 {
  public:
    void operator()(const Point_2& p, NT pixel_size, NT &x, NT &y)
    {
      NT x_tmp = p.x() / pixel_size;
      NT y_tmp = p.y() / pixel_size;

      double x_floor = std::floor(CGAL::to_double(x_tmp));
      double y_floor = std::floor(CGAL::to_double(y_tmp));
      x = NT(x_floor) * pixel_size + pixel_size / NT(2.0);
      y = NT(y_floor) * pixel_size + pixel_size / NT(2.0);
    }
  };

  Snap_2 snap_2_object() const { return Snap_2(); }

  /*! Functor */
  class Integer_grid_point_2 {
  public:
    Point_2 operator()(const Point_2& p, NT pixel_size)
    {
      NT x = (p.x() - pixel_size / NT(2.0)) / pixel_size;
      NT y = (p.y() - pixel_size / NT(2.0)) / pixel_size;
      Point_2 out_p(x,y);

      return(out_p);
    }
  };

  /*! */
  Integer_grid_point_2 integer_grid_point_2_object() const
  {
    return Integer_grid_point_2();
  }

  /*! Functor */
  class Minkowski_sum_with_pixel_2 {
  private:
    typedef Snap_rounding_traits_2<Base_kernel>         Traits;
    typedef std::list<Point_2>                          Point_list;


    const Traits * m_gt;

    Minkowski_sum_with_pixel_2(const Traits * gt) : m_gt(gt) {}

  public:
    void operator()(Point_list & points_list, const Segment_2& s, NT unit_square)
    {

      Construct_vertex_2 construct_ver = m_gt->construct_vertex_2_object();
      Compare_y_2 compare_y = m_gt->compare_y_2_object();
      Compare_x_2 compare_x = m_gt->compare_x_2_object();

      Point_2 src = construct_ver(s, 0);
      Point_2 trg = construct_ver(s, 1);
      Comparison_result cx = compare_x(src, trg);
      Comparison_result cy = compare_y(src, trg);
      NT x1 = src.x();
      NT y1 = src.y();
      NT x2 = trg.x();
      NT y2 = trg.y();
      Point_2 ms1, ms2, ms3, ms4, ms5, ms6;// minkowski sum points

      if (cx == SMALLER) {
        if (cy == SMALLER) {
          // we use unit_square instead of unit_square / 2 in order to
          // find tangency points which are not supported by kd-tree
          ms1 = Point_2(x1 - unit_square, y1 - unit_square);
          ms2 = Point_2(x1 - unit_square, y1 + unit_square);
          ms3 = Point_2(x1 + unit_square, y1 - unit_square);
          ms4 = Point_2(x2 + unit_square, y2 - unit_square);
          ms5 = Point_2(x2 + unit_square, y2 + unit_square);
          ms6 = Point_2(x2 - unit_square, y2 + unit_square);
        } else {
          ms1 = Point_2(x1 - unit_square, y1 - unit_square);
          ms2 = Point_2(x1 - unit_square, y1 + unit_square);
          ms3 = Point_2(x1 + unit_square, y1 + unit_square);
          ms4 = Point_2(x2 + unit_square, y2 - unit_square);
          ms5 = Point_2(x2 + unit_square, y2 + unit_square);
          ms6 = Point_2(x2 - unit_square, y2 - unit_square);
        }
      } else {
        if(cy == SMALLER) {
          ms1 = Point_2(x1 + unit_square, y1 - unit_square);
          ms2 = Point_2(x1 + unit_square, y1 + unit_square);
          ms3 = Point_2(x1 - unit_square, y1 - unit_square);
          ms4 = Point_2(x2 + unit_square, y2 + unit_square);
          ms5 = Point_2(x2 - unit_square, y2 + unit_square);
          ms6 = Point_2(x2 - unit_square, y2 - unit_square);
        } else {
          ms1 = Point_2(x1 + unit_square, y1 - unit_square);
          ms2 = Point_2(x1 + unit_square, y1 + unit_square);
          ms3 = Point_2(x1 - unit_square, y1 + unit_square);
          ms4 = Point_2(x2 + unit_square, y2 - unit_square);
          ms5 = Point_2(x2 - unit_square, y2 - unit_square);
          ms6 = Point_2(x2 - unit_square, y2 + unit_square);
        }
      }

      points_list.push_back(ms1);
      points_list.push_back(ms2);
      points_list.push_back(ms3);
      points_list.push_back(ms4);
      points_list.push_back(ms5);
      points_list.push_back(ms6);
    }

    friend class Snap_rounding_traits_2<Base_kernel>;
  };

  /*! */
  Minkowski_sum_with_pixel_2 minkowski_sum_with_pixel_2_object() const
  {
    return Minkowski_sum_with_pixel_2(this);
  }

  Construct_segment_2  construct_segment_2_object() const
  {
    Base_kernel k;
    return k.construct_segment_2_object();
  }

  Construct_vertex_2 construct_vertex_2_object() const
  {
    Base_kernel k;
    return k.construct_vertex_2_object();
  }

  Compare_y_2 compare_y_2_object() const
  {
    Base_kernel k;
    return k.compare_y_2_object();
  }

  Construct_iso_rectangle_2 construct_iso_rectangle_2_object() const
  {
    Base_kernel k;
    return k.construct_iso_rectangle_2_object();
  }

};

} //namespace CGAL

#endif // CGAL_ISR_2_TRAITS_H
