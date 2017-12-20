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
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Mariette Yvinec

#ifndef CGAL_TRIANGULATION_EUCLIDEAN_TRAITS_2_H
#define CGAL_TRIANGULATION_EUCLIDEAN_TRAITS_2_H

#include <CGAL/license/Triangulation_2.h>


#define CGAL_DEPRECATED_HEADER "<CGAL/Triangulation_euclidean_traits_2.h>"
#include <CGAL/internal/deprecation_warning.h>

#include <CGAL/Point_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/Line_2.h>
#include <CGAL/Ray_2.h>
#include <CGAL/predicates_on_points_2.h>
#include <CGAL/basic_constructions_2.h>
#include <CGAL/distance_predicates_2.h>

#include <CGAL/triangulation_assertions.h>
#include <CGAL/Segment_2_Segment_2_intersection.h>

namespace CGAL { 

template < class R >
class Triangulation_euclidean_traits_2 {
public:
  typedef R Rep;
  typedef typename R::Point_2     Point_2;
  typedef typename R::Segment_2   Segment_2;
  typedef typename R::Triangle_2  Triangle_2;
  typedef typename R::Line_2      Line_2;
  typedef typename R::Ray_2       Ray_2;

  typedef typename R::Less_x_2                   Less_x_2;
  typedef typename R::Less_y_2                   Less_y_2;
  typedef typename R::Compare_x_2                Compare_x_2;
  typedef typename R::Compare_y_2                Compare_y_2;
  typedef typename R::Orientation_2              Orientation_2;
  typedef typename R::Side_of_oriented_circle_2  Side_of_oriented_circle_2;
  typedef typename R::Construct_circumcenter_2   Construct_circumcenter_2;
  typedef typename R::Construct_bisector_2       Construct_bisector_2;
  typedef typename R::Compare_distance_2         Compare_distance_2;
  typedef typename R::Construct_segment_2        Construct_segment_2;
  typedef typename R::Construct_triangle_2       Construct_triangle_2;
  typedef typename R::Construct_direction_2      Construct_direction_2;
  typedef typename R::Construct_ray_2            Construct_ray_2;
  
  //for natural_neighbor_coordinates_2
  typedef typename R::FT                         FT;
  typedef typename R::Equal_x_2                  Equal_x_2;
  typedef typename R::Compute_area_2             Compute_area_2;
  Compute_area_2 compute_area_2_object () const {return Compute_area_2();}
  
  // for compatibility with previous versions
  typedef Point_2      Point;
  typedef Segment_2    Segment;
  typedef Triangle_2   Triangle;
  typedef Ray_2        Ray;
  typedef Line_2       Line;

  Triangulation_euclidean_traits_2() {}
  Triangulation_euclidean_traits_2(const Triangulation_euclidean_traits_2 &) {}
  Triangulation_euclidean_traits_2 &operator=
      (const Triangulation_euclidean_traits_2 &)
  {return *this;}
 
  Less_x_2
  less_x_2_object() const
    { return Less_x_2();}

  Less_y_2
  less_y_2_object() const
    { return Less_y_2();}
  
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
 
  Construct_circumcenter_2
  construct_circumcenter_2_object() const
    { return Construct_circumcenter_2();}

  Construct_bisector_2
  construct_bisector_2_object() const
    {return Construct_bisector_2();}
  
  Compare_distance_2
  compare_distance_2_object() const
    {return Compare_distance_2();}

  Construct_segment_2  construct_segment_2_object() const
    {return Construct_segment_2();}

  Construct_triangle_2  construct_triangle_2_object() const
    {return Construct_triangle_2();}

  Construct_direction_2  construct_direction_2_object() const
    {return Construct_direction_2();}

  Construct_ray_2  construct_ray_2_object() const
    {return Construct_ray_2();}

};

} //namespace CGAL 

#endif // CGAL_TRIANGULATION_EUCLIDEAN_TRAITS_2_H
