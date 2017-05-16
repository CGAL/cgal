// Copyright (c) 1999-2004,2006-2009,2013-2015,2017  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Monique Teillaud <Monique.Teillaud@inria.fr>
//                 Aymeric Pelle <Aymeric.Pelle@sophia.inria.fr>
//                 Mael Rouxel-Labbé

#ifndef CGAL_PERIODIC_3_REGULAR_TRIANGULATION_TRAITS_3_H
#define CGAL_PERIODIC_3_REGULAR_TRIANGULATION_TRAITS_3_H

#include <CGAL/license/Periodic_3_triangulation_3.h>

#include <CGAL/internal/Functor_with_offset_weighted_points_adaptor.h>
#include <CGAL/internal/Periodic_3_construct_weighted_point_3.h>
#include <CGAL/Periodic_3_offset_3.h>
#include <CGAL/Periodic_3_triangulation_traits_3.h>

#include <CGAL/basic.h>
#include <CGAL/internal/Has_boolean_tags.h>

namespace CGAL {

template <class K, class Off = typename CGAL::Periodic_3_offset_3>
class Periodic_3_regular_triangulation_traits_base_3
  : public Periodic_3_triangulation_traits_3<K, Off>
{
  typedef Periodic_3_regular_triangulation_traits_base_3<K, Off> Self;
  typedef Periodic_3_triangulation_traits_3<K, Off>              Base;

public:
  typedef K                                          Kernel;
  typedef Off                                        Offset;

  typedef typename Base::RT                          RT;
  typedef typename Base::FT                          FT;
  typedef typename Base::Point_3                     Point_3;
  typedef typename Base::Weighted_point_3            Weighted_point_3;

  typedef typename Base::Periodic_3_offset_3         Periodic_3_offset_3;
  typedef typename Base::Iso_cuboid_3                Iso_cuboid_3;

  typedef typename Base::Segment_3                   Segment_3;
  typedef typename Base::Triangle_3                  Triangle_3;
  typedef typename Base::Tetrahedron_3               Tetrahedron_3;

public:
  Periodic_3_regular_triangulation_traits_base_3(const Iso_cuboid_3& domain,
                                                 const K& k)
    : Base(domain, k)
  { }

  // Construct_weighted_point_3 with offset
  typedef Periodic_3_construct_weighted_point_3<Self, typename K::Construct_weighted_point_3>
      Construct_weighted_point_3;

  typedef Functor_with_offset_weighted_points_adaptor<Self, typename K::Power_side_of_oriented_power_sphere_3>
      Power_side_of_oriented_power_sphere_3;
  typedef Functor_with_offset_weighted_points_adaptor<Self, typename K::Compare_weighted_squared_radius_3>
      Compare_weighted_squared_radius_3;
  typedef Functor_with_offset_weighted_points_adaptor<Self, typename K::Compare_power_distance_3>
      Compare_power_distance_3;

  // Required for Periodic_3_regular_remove_traits
  typedef Functor_with_offset_weighted_points_adaptor<Self, typename K::Coplanar_orientation_3>
       Coplanar_orientation_3;

  // When dual operations are used
  typedef Functor_with_offset_weighted_points_adaptor<Self, typename K::Construct_weighted_circumcenter_3>
      Construct_weighted_circumcenter_3;

  // Operations
  using Base::construct_point_3_object;

  Construct_weighted_point_3 construct_weighted_point_3_object() const {
    return Construct_weighted_point_3(&this->_domain,
                                      this->Base::construct_weighted_point_3_object());
  }

  // construction
  Construct_weighted_circumcenter_3 construct_weighted_circumcenter_3_object() const {
    return Construct_weighted_circumcenter_3(
      this->Base::construct_weighted_circumcenter_3_object(),
      construct_point_3_object(), construct_weighted_point_3_object());
  }

  // predicates
  Power_side_of_oriented_power_sphere_3 power_side_of_oriented_power_sphere_3_object() const {
    return Power_side_of_oriented_power_sphere_3(
      this->Base::power_side_of_oriented_power_sphere_3_object(),
      construct_point_3_object(), construct_weighted_point_3_object());
  }
  Compare_power_distance_3 compare_power_distance_3_object() const {
    return Compare_power_distance_3(
      this->Base::compare_power_distance_3_object(),
      construct_point_3_object(), construct_weighted_point_3_object());
  }
  Compare_weighted_squared_radius_3 compare_weighted_squared_radius_3_object() const {
    return Compare_weighted_squared_radius_3(
      this->Base::compare_weighted_squared_radius_3_object(),
      construct_point_3_object(), construct_weighted_point_3_object());
  }
  Coplanar_orientation_3 coplanar_orientation_3_object() const {
    return Coplanar_orientation_3(
      this->Base::coplanar_orientation_3_object(),
      construct_point_3_object(), construct_weighted_point_3_object());
  }
};

template<typename K,
         typename Off = CGAL::Periodic_3_offset_3,
         bool Has_filtered_predicates = K::Has_filtered_predicates>
class Periodic_3_regular_triangulation_traits_3;

} // namespace CGAL

#include <CGAL/internal/Periodic_3_regular_triangulation_filtered_traits_3.h>

namespace CGAL
{
template<class K, class Off>
class Periodic_3_regular_triangulation_traits_3<K, Off, false /* Has_filtered_predicates */>
  : public Periodic_3_regular_triangulation_traits_base_3<K, Off>
{
  typedef Periodic_3_regular_triangulation_traits_base_3<K, Off> Base;

public:
  typedef typename K::Iso_cuboid_3 Iso_cuboid_3;

  Periodic_3_regular_triangulation_traits_3(const Iso_cuboid_3& domain = Iso_cuboid_3(0,0,0,1,1,1),
                                            const K& k = K())
    : Base(domain, k)
  { }
};

template <typename K, typename Off>
class Periodic_3_regular_triangulation_traits_3<K, Off, true /* Has_filtered_predicates */ >
    : public Periodic_3_regular_triangulation_filtered_traits_3<K, Off>
{
  typedef Periodic_3_regular_triangulation_filtered_traits_3<K, Off> Base;

public:
  typedef typename K::Iso_cuboid_3 Iso_cuboid_3;

  Periodic_3_regular_triangulation_traits_3(const Iso_cuboid_3& domain = Iso_cuboid_3(0,0,0,1,1,1),
                                            const K& k = K())
    : Base(domain, k)
  { }
};

} //namespace CGAL

#endif
