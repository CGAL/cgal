// Copyright (c) 2006-2009, 2017  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Nico Kruithof <Nico.Kruithof@sophia.inria.fr>
//                 Manuel Caroli <Manuel.Caroli@sophia.inria.fr>
//                 Mael Rouxel-Labbé

#ifndef CGAL_PERIODIC_3_DELAUNAY_TRIANGULATION_TRAITS_3_H
#define CGAL_PERIODIC_3_DELAUNAY_TRIANGULATION_TRAITS_3_H

#include <CGAL/license/Periodic_3_triangulation_3.h>

#include <CGAL/internal/Periodic_3_construct_point_3.h>
#include <CGAL/internal/Functor_with_offset_points_adaptor_3.h>
#include <CGAL/Periodic_3_offset_3.h>
#include <CGAL/Periodic_3_triangulation_traits_3.h>

#include <CGAL/basic.h>
#include <CGAL/internal/Has_boolean_tags.h>

namespace CGAL {

template < class K, class Off = typename CGAL::Periodic_3_offset_3 >
class Periodic_3_Delaunay_triangulation_traits_base_3
  : public Periodic_3_triangulation_traits_3<K, Off>
{
  typedef Periodic_3_Delaunay_triangulation_traits_base_3<K, Off>  Self;
  typedef Periodic_3_triangulation_traits_3<K, Off>                Base;

public:
  typedef K                                   Kernel;
  typedef Off                                 Offset;

  typedef typename Base::RT                   RT;
  typedef typename Base::FT                   FT;
  typedef typename Base::Point_3              Point_3;
  typedef typename Base::Periodic_3_offset_3  Periodic_3_offset_3;
  typedef typename Base::Iso_cuboid_3         Iso_cuboid_3;

  // The next typedef is there for backward compatibility
  // Some users take their point type from the traits class.
  // Before this type was Point
  typedef Point_3                             Point;

  typedef typename Base::Segment_3            Segment_3;
  typedef typename Base::Triangle_3           Triangle_3;
  typedef typename Base::Tetrahedron_3        Tetrahedron_3;

public:
  Periodic_3_Delaunay_triangulation_traits_base_3(const Iso_cuboid_3& domain,
                                                  const K& k)
    : Base(domain, k)
  { }

  // Delaunay specific predicates
  typedef Functor_with_offset_points_adaptor_3<Self, typename K::Side_of_oriented_sphere_3>
      Side_of_oriented_sphere_3;
  typedef Functor_with_offset_points_adaptor_3<Self, typename K::Compare_distance_3>
      Compare_distance_3;

  // Required for Periodic_3_Delaunay_triangulation_remove_traits
  typedef Functor_with_offset_points_adaptor_3<Self, typename K::Coplanar_orientation_3>
      Coplanar_orientation_3;
  typedef Functor_with_offset_points_adaptor_3<Self, typename K::Coplanar_side_of_bounded_circle_3>
      Coplanar_side_of_bounded_circle_3;

  // When is_Gabriel is used
   typedef Functor_with_offset_points_adaptor_3<Self, typename K::Side_of_bounded_sphere_3>
       Side_of_bounded_sphere_3;

  // Delaunay specific constructions
  typedef Functor_with_offset_points_adaptor_3<Self, typename K::Construct_circumcenter_3>
      Construct_circumcenter_3;

  using Base::construct_point_3_object;

  // Operations
  Side_of_oriented_sphere_3 side_of_oriented_sphere_3_object() const {
    return Side_of_oriented_sphere_3(this->Base::side_of_oriented_sphere_3_object(),
                                     construct_point_3_object());
  }
  Compare_distance_3 compare_distance_3_object() const {
    return Compare_distance_3(this->Base::compare_distance_3_object(),
                              construct_point_3_object());
  }
  Side_of_bounded_sphere_3 side_of_bounded_sphere_3_object() const {
    return Side_of_bounded_sphere_3(this->Base::side_of_bounded_sphere_3_object(),
                                    construct_point_3_object());
  }
  Coplanar_orientation_3 coplanar_orientation_3_object() const {
    return Coplanar_orientation_3(this->Base::coplanar_orientation_3_object(),
                                  construct_point_3_object());
  }
  Coplanar_side_of_bounded_circle_3 coplanar_side_of_bounded_circle_3_object() const {
    return Coplanar_side_of_bounded_circle_3(this->Base::coplanar_side_of_bounded_circle_3_object(),
                                             construct_point_3_object());
  }
  Construct_circumcenter_3 construct_circumcenter_3_object() const {
    return Construct_circumcenter_3(this->Base::construct_circumcenter_3_object(),
                                    construct_point_3_object());
  }
};

template < typename K,
           typename Off = CGAL::Periodic_3_offset_3,
           bool Has_filtered_predicates = internal::Has_filtered_predicates<K>::value >
class Periodic_3_Delaunay_triangulation_traits_3;

} //namespace CGAL

// Partial specialization for kernels with filtered predicates
#include <CGAL/internal/Periodic_3_Delaunay_triangulation_filtered_traits_3.h>

namespace CGAL {

template < class K, class Off>
class Periodic_3_Delaunay_triangulation_traits_3<K, Off, false>
  : public Periodic_3_Delaunay_triangulation_traits_base_3<K, Off>
{
  typedef Periodic_3_Delaunay_triangulation_traits_base_3<K, Off> Base;

public:
  typedef typename K::Iso_cuboid_3 Iso_cuboid_3;

  Periodic_3_Delaunay_triangulation_traits_3(const Iso_cuboid_3& domain = Iso_cuboid_3(0,0,0,1,1,1),
                                             const K& k = K())
    : Base(domain, k)
  { }
};

template < typename K, typename Off >
class Periodic_3_Delaunay_triangulation_traits_3<K, Off, true>
    : public Periodic_3_Delaunay_triangulation_filtered_traits_3<
    K, Off, internal::Has_static_filters<K>::value>
{
  typedef Periodic_3_Delaunay_triangulation_filtered_traits_3<
  K, Off, internal::Has_static_filters<K>::value>      Base;

public:
  typedef typename K::Iso_cuboid_3 Iso_cuboid_3;

  Periodic_3_Delaunay_triangulation_traits_3(const Iso_cuboid_3& domain = Iso_cuboid_3(0,0,0,1,1,1),
                                             const K& k = K())
    : Base(domain, k)
  { }
};

} //namespace CGAL

#endif // CGAL_PERIODIC_3_DELAUNAY_TRIANGULATION_TRAITS_3_H
