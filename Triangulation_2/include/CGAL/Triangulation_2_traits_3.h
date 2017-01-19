// Copyright (c) 1997   INRIA Sophia-Antipolis (France).
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
// Author(s)     : Mariette Yvinec

#ifndef CGAL_TRIANGULATION_2_TRAITS_3_H
#define CGAL_TRIANGULATION_2_TRAITS_3_H

#include <CGAL/license/Triangulation_2.h>



#include <CGAL/Point_3.h>
#include <CGAL/Segment_3.h>
#include <CGAL/Triangle_3.h>
#include <CGAL/Kernel/global_functions_3.h>

#include <CGAL/triangulation_assertions.h>

namespace CGAL { 

template<class R>
class Compare_yz_3
{
public:
  typedef typename  R::Point_3     Point;

  Comparison_result operator() (Point p, Point q){
    Comparison_result r;
    r = CGAL::compare_y(p,q);
    if (r == EQUAL) r = CGAL::compare_z(p,q);
    return r;
   }
};
    

template <class R>
class Side_of_oriented_circle_2_3
{
  // 2d triangulation needs a side_of_oriented_circle
   // that in fact is a side_of_bounded_circle
  // meaning that
  // bounded side of circle = positive side
public:
  typedef typename  R::Point_3                   Point;
  typedef typename  R::Coplanar_side_of_bounded_circle_3
                                                 Side_of_bounded_circle_2_3;
  Oriented_side operator() (const Point& p, 
			    const Point& q, 
			    const Point& r,
			    const Point& s) {
    Side_of_bounded_circle_2_3  side;
    Bounded_side bs = side(p,q,r,s);
    return ( bs == ON_UNBOUNDED_SIDE) ? ON_NEGATIVE_SIDE :
      (bs == ON_BOUNDED_SIDE ) ? ON_POSITIVE_SIDE :
      ON_ORIENTED_BOUNDARY;
  }   
};



template < class R >
class Triangulation_2_traits_3 
{
public:
  typedef R Rep;
  typedef typename Rep::Point_3    Point_2;
  typedef typename Rep::Segment_3  Segment_2;
  typedef typename Rep::Triangle_3 Triangle_2;
 
  typedef typename Rep::Compare_x_3               Compare_x_2;
  typedef Compare_yz_3<Rep>                       Compare_y_2;
  typedef typename Rep::Coplanar_orientation_3    Orientation_2;
  typedef Side_of_oriented_circle_2_3<Rep>        Side_of_oriented_circle_2;  
  typedef typename Rep::Construct_segment_3       Construct_segment_2;
  typedef typename Rep::Construct_triangle_3      Construct_triangle_2;

  // for compatibility with previous versions
  typedef Point_2      Point;
  typedef Segment_2    Segment;
  typedef Triangle_2   Triangle;

  Compare_x_2
  compare_x_2_object() const
    { return Compare_x_2();}

  Compare_y_2
  compare_y_2_object() const
    { return Compare_y_2();}
  
  Orientation_2
  orientation_2_object() const
    { return Orientation_2();}

  Side_of_oriented_circle_2
  side_of_oriented_circle_2_object() const
  {return Side_of_oriented_circle_2();}

  Construct_segment_2  construct_segment_2_object() const
    {return Construct_segment_2();}

  Construct_triangle_2  construct_triangle_2_object() const
    {return Construct_triangle_2();}

};

} //namespace CGAL 
#endif // CGAL_TRIANGULATION_2_TRAITS_3_H
