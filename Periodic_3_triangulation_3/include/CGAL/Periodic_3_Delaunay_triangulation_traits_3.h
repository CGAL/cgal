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

#ifndef CGAL_PERIODIC_3_DELAUNAY_TRIANGULATION_TRAITS_3_H
#define CGAL_PERIODIC_3_DELAUNAY_TRIANGULATION_TRAITS_3_H

#include <CGAL/basic.h>
#include <CGAL/Periodic_3_offset_3.h>
#include <CGAL/Traits_with_offsets_adaptor.h>
#include <CGAL/Periodic_3_construct_point_3.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Periodic_3_triangulation_traits_3.h>


namespace CGAL { 

template < class Kernel, class Off = typename CGAL::Periodic_3_offset_3 >
class Periodic_3_Delaunay_triangulation_traits_base_3
  : public Periodic_3_triangulation_traits_base_3<Kernel, Off>
{
public:
  typedef Kernel                                                 K;
  typedef Off                                                    Offset;
  typedef Periodic_3_triangulation_traits_base_3<K, Offset>               Base;
  typedef Periodic_3_Delaunay_triangulation_traits_base_3< K, Offset >    Self;

  typedef typename Base::RT                   RT;
  typedef typename Base::FT                   FT;
  typedef typename Base::Point_3              Point_3;
  typedef typename Base::Vector_3             Vector_3;
  typedef typename Base::Periodic_3_offset_3  Periodic_3_offset_3;
  typedef typename Base::Iso_cuboid_3         Iso_cuboid_3;

  // The next typedef is there for backward compatibility
  // Some users take their point type from the traits class.
  // Before this type was Point
  typedef Point_3 Point;

  typedef typename Base::Segment_3         Segment_3;
  typedef typename Base::Triangle_3        Triangle_3;
  typedef typename Base::Tetrahedron_3     Tetrahedron_3;

  // Delaunay specific predicates
  typedef Traits_with_offsets_adaptor<Self,
				      typename K::Side_of_oriented_sphere_3>
      Side_of_oriented_sphere_3;
  typedef Traits_with_offsets_adaptor<Self, typename K::Compare_distance_3>
      Compare_distance_3;
   typedef Traits_with_offsets_adaptor<Self,
 				      typename K::Side_of_bounded_sphere_3>
       Side_of_bounded_sphere_3;

  // Degenerate dimension predicates
  typedef Traits_with_offsets_adaptor<Self, typename K::Coplanar_orientation_3>
      Coplanar_orientation_3;
  typedef Traits_with_offsets_adaptor<Self,
              typename K::Coplanar_side_of_bounded_circle_3>
      Coplanar_side_of_bounded_circle_3;

  // Delaunay specific constructions
  typedef Traits_with_offsets_adaptor<Self,
				      typename K::Construct_circumcenter_3>
      Construct_circumcenter_3;

  // Operations
  Side_of_oriented_sphere_3
  side_of_oriented_sphere_3_object() const {
    return Side_of_oriented_sphere_3(&this->_domain);
  }
  Compare_distance_3
  compare_distance_3_object() const {
    return Compare_distance_3(&this->_domain);
  }
  Side_of_bounded_sphere_3
  side_of_bounded_sphere_3_object() const {
    return Side_of_bounded_sphere_3(&this->_domain);
  }
  Coplanar_orientation_3
  coplanar_orientation_3_object() const {
    return Coplanar_orientation_3(&this->_domain);
  }
  Coplanar_side_of_bounded_circle_3
  coplanar_side_of_bounded_circle_3_object() const {
    return Coplanar_side_of_bounded_circle_3(&this->_domain);
  }
  Construct_circumcenter_3
  construct_circumcenter_3_object() const {
    return Construct_circumcenter_3(&this->_domain);
  }
};

template < typename K, typename Off = CGAL::Periodic_3_offset_3 >
class Periodic_3_Delaunay_triangulation_traits_3;

} //namespace CGAL

// Partial specialization for Filtered_kernel<CK>.
#include <CGAL/Filtered_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Periodic_3_Delaunay_triangulation_filtered_traits_3.h>

namespace CGAL {

// This declaration is needed to break the cyclic dependency.
template < typename K, typename Off >
class Periodic_3_Delaunay_triangulation_filtered_traits_3;

template < class K, class Off>
class Periodic_3_Delaunay_triangulation_traits_3
  : public Periodic_3_Delaunay_triangulation_traits_base_3<K, Off>
{
};

template < typename CK, typename Off >
class Periodic_3_Delaunay_triangulation_traits_3 < Filtered_kernel<CK>, Off>
  : public Periodic_3_Delaunay_triangulation_filtered_traits_3 <
  Filtered_kernel<CK>, Off >
{
public:
  typedef Filtered_kernel<CK>  Kernel;
};

template < class Off >
class Periodic_3_Delaunay_triangulation_traits_3<CGAL::Epick, Off>
  : public Periodic_3_Delaunay_triangulation_filtered_traits_3<CGAL::Epick, Off>
{
  typedef CGAL::Epick Kernel;
};

} //namespace CGAL

#endif // CGAL_PERIODIC_3_DELAUNAY_TRIANGULATION_TRAITS_3_H
