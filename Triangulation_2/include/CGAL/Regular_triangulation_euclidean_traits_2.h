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
// Author(s)     : Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>
//                 Sylvain Pion

#ifndef CGAL_REGULAR_TRIANGULATION_EUCLIDEAN_TRAITS_2_H
#define CGAL_REGULAR_TRIANGULATION_EUCLIDEAN_TRAITS_2_H

#include <CGAL/basic.h>


namespace CGAL { 

template < class R, class W = typename R::RT>
class Regular_triangulation_euclidean_traits_2
  : public R
{
public:
  typedef R                                     Kernel;
  typedef R                                     Rep;
  typedef typename R::FT                        Weight;
  typedef R                                     Traits;
  typedef typename Traits::Point_2              Bare_point;
  typedef typename Traits::Point_2              Point_2;
  typedef typename R::Weighted_point_2          Weighted_point_2;
  // This is required for the point() function of vertex base class
  // to be correctly return a weighted_point;
  // patch 27/11/00
  // 18/03/03 I put now the same typedef in Regulat_triangulation_2
  // for the need of hierarchy
  // don't know if this is definitive
  //typedef Weighted_point                        Point_2;

  typedef Regular_triangulation_euclidean_traits_2<R, W>   Self;

  typedef typename R::Power_side_of_oriented_power_circle_2 Power_side_of_oriented_power_circle_2;
  typedef typename R::Compare_power_distance_2  Compare_power_distance_2;

  // construction objects
  typedef typename R::Construct_weighted_circumcenter_2
                                            Construct_weighted_circumcenter_2;
  typedef typename R::Construct_radical_axis_2 Construct_radical_axis_2;
  
  Power_side_of_oriented_power_circle_2 
  power_side_of_power_oriented_circle_2_object() const
    {  return Power_side_of_oriented_power_circle_2();}

  Compare_power_distance_2
  compare_power_distance_2_object() const {
    return Compare_power_distance_2();
  }

  //constructions for dual:
  Construct_weighted_circumcenter_2
  construct_weighted_circumcenter_2_object() const
    {return Construct_weighted_circumcenter_2();}
  
  Construct_radical_axis_2
  construct_radical_axis_2_object() const
    {return Construct_radical_axis_2();}

};
 

} //namespace CGAL

#endif // CGAL_REGULAR_TRIANGULATION_EUCLIDEAN_TRAITS_2_H
