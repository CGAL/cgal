// Copyright (c) 2006-2009   INRIA Sophia-Antipolis (France).
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
// Author(s)     : Nico Kruithof <Nico.Kruithof@sophia.inria.fr>
//                 Manuel Caroli <Manuel.Caroli@sophia.inria.fr>

#ifndef CGAL_PERIODIC_4_HYPERBOLIC_TRIANGULATION_TRAITS_2_H
#define CGAL_PERIODIC_4_HYPERBOLIC_TRIANGULATION_TRAITS_2_H


#include <CGAL/Traits_with_hyperbolic_offsets_adaptor.h>
#include <CGAL/Triangulation_hyperbolic_traits_2.h>
#include <CGAL/Square_root_2_field.h>
#include <CGAL/Hyperbolic_octagon_group.h>
#include <CGAL/Periodic_4_construct_hyperbolic_point_2.h>


namespace CGAL { 

template < class THT2 >
class Periodic_4_hyperbolic_triangulation_traits_base_2
  : public THT2
{
public:
  typedef THT2                                                   K;
  typedef HyperbolicOctagonGroup                                 Offset;
  typedef Periodic_4_hyperbolic_triangulation_traits_base_2< K >    Self;  

  typedef typename K::RT                RT;
  typedef typename K::FT                FT;
  typedef typename K::Point_2           Point_2;
  //typedef typename K::Vector_3          Vector_3;
  typedef Offset                        Periodic_4_hyperbolic_offset_2;
  typedef typename K::Circle_2          Circle_2;

  // The next typedef is there for backward compatibility
  // Some users take their point type from the traits class.
  // Before this type was Point
  typedef Point_2 Point;


  // Triangulation predicates
  typedef Traits_with_hyperbolic_offsets_adaptor<Self, typename K::Compare_xy_2>
      Compare_xy_2;
  typedef Traits_with_hyperbolic_offsets_adaptor<Self, typename K::Orientation_2>
      Orientation_2;
  
  // Triangulation constructions
  typedef Periodic_4_construct_hyperbolic_point_2<Self, typename K::Construct_point_2>
      Construct_hyperbolic_point_2;
  typedef Traits_with_hyperbolic_offsets_adaptor<Self, typename K::Construct_segment_2>
      Construct_segment_2;
  typedef Traits_with_hyperbolic_offsets_adaptor<Self, typename K::Construct_triangle_2>
      Construct_hyperbolic_triangle_2;

  // Access
  void set_domain(const Circle_2& domain) {
    _domain = domain;
  }
  
  Circle_2 get_domain() const {
    return _domain;
  }

  // Operations
  Compare_xy_2
  compare_xy_2_object() const {
    return Compare_xy_2(&_domain);
  }
  Orientation_2
  orientation_2_object() const {
    return Orientation_2(&_domain);
  }

  Construct_hyperbolic_point_2
  construct_hyperbolic_point_2_object() const {
    return Construct_point_2(_domain);
  }
  Construct_segment_2
  construct_segment_2_object() const {
    return Construct_segment_2(&_domain);
  }


protected:
  Circle_2 _domain;
};

namespace future_release
{
template < typename K >
class Periodic_4_hyperbolic_triangulation_traits_2;
}

} //namespace CGAL


#include <CGAL/Periodic_3_Delaunay_triangulation_traits_3.h>


namespace CGAL
{
template < typename K >
class Periodic_4_hyperbolic_Delaunay_triangulation_traits_2;

}

#endif // CGAL_PERIODIC_4_HYPERBOLIC_TRIANGULATION_TRAITS_2_H
