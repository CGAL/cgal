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

#ifndef CGAL_PERIODIC_4_HYPERBOLIC_DELAUNAY_TRIANGULATION_TRAITS_2_H
#define CGAL_PERIODIC_4_HYPERBOLIC_DELAUNAY_TRIANGULATION_TRAITS_2_H

//#include <CGAL/basic.h>
#include <CGAL/Traits_with_hyperbolic_offsets_adaptor.h>
#include <CGAL/Periodic_4_construct_hyperbolic_point_2.h>
#include <CGAL/Periodic_4_hyperbolic_triangulation_traits_2.h>
#include <CGAL/Hyperbolic_octagon_group.h>

namespace CGAL { 

template < class THT2 >
class Periodic_4_hyperbolic_Delaunay_triangulation_traits_base_2
  : public Periodic_3_triangulation_traits_base_3<THT2>
{
public:
  typedef THT2                                                               K;
  typedef HyperbolicOctagonGroup                                             Offset;
  typedef Periodic_4_hyperbolic_triangulation_traits_base_2< K >             Base;
  typedef Periodic_4_hyperbolic_Delaunay_triangulation_traits_base_2< K >    Self;

  typedef typename Base::RT                   RT;
  typedef typename Base::FT                   FT;
  typedef typename Base::Point_2              Point_2;
  typedef typename Base::Vector_2             Vector_2;
  //typedef typename Base::Periodic_4_offset_2  Periodic_4_offset_2;

  // The next typedef is there for backward compatibility
  // Some users take their point type from the traits class.
  // Before this type was Point
  typedef Point_2 Point;

  typedef typename Base::Segment_2         Segment_2;
  typedef typename Base::Triangle_2        Triangle_2;

  // Delaunay specific predicates
  typedef Traits_with_offsets_adaptor<Self,
				      typename K::Side_of_oriented_circle_2>
      Side_of_oriented_circle_2;
  typedef Traits_with_offsets_adaptor<Self, typename K::Compare_distance_2>
      Compare_distance_2;
   typedef Traits_with_offsets_adaptor<Self,
 				      typename K::Side_of_bounded_disc_2>
       Side_of_bounded_disc_2;

  // Delaunay specific constructions
  typedef Traits_with_hyperbolic_offsets_adaptor<Self,
				      typename K::Construct_circumcenter_2>
      Construct_circumcenter_2;

  // Operations
  Side_of_oriented_circle_2
  side_of_oriented_circle_2_object() const {
    return Side_of_oriented_circle_2(&this->_domain);
  }
  Compare_distance_2
  compare_distance_2_object() const {
    return Compare_distance_2(&this->_domain);
  }
  Side_of_bounded_disc_2
  side_of_bounded_disc_2_object() const {
    return Side_of_bounded_disc_2(&this->_domain);
  }
  Construct_circumcenter_2
  construct_circumcenter_2_object() const {
    return Construct_circumcenter_2(&this->_domain);
  }
};

template < typename K >
class Periodic_4_hyperbolic_Delaunay_triangulation_traits_2;

} //namespace CGAL

#endif // CGAL_PERIODIC_4_HYPERBOLIC_DELAUNAY_TRIANGULATION_TRAITS_2_H
