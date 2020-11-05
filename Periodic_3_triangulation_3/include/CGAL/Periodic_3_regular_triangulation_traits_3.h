// Copyright (c) 1999-2004,2006-2009,2013-2015,2017  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Monique Teillaud <Monique.Teillaud@inria.fr>
//                 Aymeric Pelle <Aymeric.Pelle@sophia.inria.fr>
//                 Mael Rouxel-Labbé

#ifndef CGAL_PERIODIC_3_REGULAR_TRIANGULATION_TRAITS_3_H
#define CGAL_PERIODIC_3_REGULAR_TRIANGULATION_TRAITS_3_H

#include <CGAL/license/Periodic_3_triangulation_3.h>

#include <CGAL/internal/Functor_with_offset_weighted_points_adaptor_3.h>
#include <CGAL/internal/Periodic_3_construct_weighted_point_3.h>
#include <CGAL/Periodic_3_offset_3.h>
#include <CGAL/Periodic_3_triangulation_traits_3.h>

#include <CGAL/basic.h>
#include <CGAL/internal/Has_boolean_tags.h>

namespace CGAL {

template <class K_,
          class Off_ = typename CGAL::Periodic_3_offset_3>
class Periodic_3_regular_triangulation_traits_base_3
  : public Periodic_3_triangulation_traits_3<K_, Off_>
{
  typedef Periodic_3_regular_triangulation_traits_base_3<K_, Off_> Self;
  typedef Periodic_3_triangulation_traits_3<K_, Off_>              Base;

public:
  typedef K_                                         Kernel;
  typedef Off_                                       Offset;

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
                                                 const Kernel& k)
    : Base(domain, k)
  { }

  // Construct_weighted_point_3 with offset
  typedef Periodic_3_construct_weighted_point_3<Self, typename Kernel::Construct_weighted_point_3>
      Construct_weighted_point_3;

  typedef Functor_with_offset_weighted_points_adaptor_3<Self, typename Kernel::Power_side_of_oriented_power_sphere_3>
      Power_side_of_oriented_power_sphere_3;
  typedef Functor_with_offset_weighted_points_adaptor_3<Self, typename Kernel::Compare_weighted_squared_radius_3>
      Compare_weighted_squared_radius_3;
  typedef Functor_with_offset_weighted_points_adaptor_3<Self, typename Kernel::Compare_power_distance_3>
      Compare_power_distance_3;

  // Undocumented, required for Is_Gabriel (Alpha shapes)
  typedef Functor_with_offset_weighted_points_adaptor_3<Self, typename Kernel::Power_side_of_bounded_power_sphere_3>
      Power_side_of_bounded_power_sphere_3;

  // Required for Periodic_3_regular_remove_traits
  typedef Functor_with_offset_weighted_points_adaptor_3<Self, typename Kernel::Coplanar_orientation_3>
       Coplanar_orientation_3;

  // When dual operations are used
  typedef Functor_with_offset_weighted_points_adaptor_3<Self, typename Kernel::Construct_weighted_circumcenter_3>
      Construct_weighted_circumcenter_3;

  // Required for Periodic_3_mesh_3
  typedef Functor_with_offset_weighted_points_adaptor_3<Self, typename Kernel::Compute_power_distance_to_power_sphere_3>
      Compute_power_distance_to_power_sphere_3;
  typedef Functor_with_offset_weighted_points_adaptor_3<Self, typename Kernel::Compute_squared_distance_3>
      Compute_squared_distance_3;

  // Operations
  Construct_weighted_point_3 construct_weighted_point_3_object() const {
    return Construct_weighted_point_3(&this->_domain,
                                      this->Base::construct_weighted_point_3_object());
  }

  // constructions
  Construct_weighted_circumcenter_3 construct_weighted_circumcenter_3_object() const {
    return Construct_weighted_circumcenter_3(
      this->Base::construct_weighted_circumcenter_3_object(),
      this->construct_point_3_object(), construct_weighted_point_3_object());
  }

  Compute_power_distance_to_power_sphere_3 compute_power_distance_to_power_sphere_3_object() const {
    return Compute_power_distance_to_power_sphere_3(
      this->Base::compute_power_distance_to_power_sphere_3_object(),
      this->construct_point_3_object(), construct_weighted_point_3_object());
  }

  Compute_squared_distance_3 compute_squared_distance_3_object() const {
    return Compute_squared_distance_3(
      this->Base::compute_squared_distance_3_object(),
      this->construct_point_3_object(), construct_weighted_point_3_object());
  }

  // predicates
  Power_side_of_bounded_power_sphere_3 power_side_of_bounded_power_sphere_3_object() const {
    return Power_side_of_bounded_power_sphere_3(
      this->Base::power_side_of_bounded_power_sphere_3_object(),
      this->construct_point_3_object(), construct_weighted_point_3_object());
  }
  Power_side_of_oriented_power_sphere_3 power_side_of_oriented_power_sphere_3_object() const {
    return Power_side_of_oriented_power_sphere_3(
      this->Base::power_side_of_oriented_power_sphere_3_object(),
      this->construct_point_3_object(), construct_weighted_point_3_object());
  }
  Compare_power_distance_3 compare_power_distance_3_object() const {
    return Compare_power_distance_3(
      this->Base::compare_power_distance_3_object(),
      this->construct_point_3_object(), construct_weighted_point_3_object());
  }
  Compare_weighted_squared_radius_3 compare_weighted_squared_radius_3_object() const {
    return Compare_weighted_squared_radius_3(
      this->Base::compare_weighted_squared_radius_3_object(),
      this->construct_point_3_object(), construct_weighted_point_3_object());
  }
  Coplanar_orientation_3 coplanar_orientation_3_object() const {
    return Coplanar_orientation_3(
      this->Base::coplanar_orientation_3_object(),
      this->construct_point_3_object(), construct_weighted_point_3_object());
  }
};

template<typename K_,
         typename Off_ = CGAL::Periodic_3_offset_3,
         bool Has_filtered_predicates_ = internal::Has_filtered_predicates<K_>::value>
class Periodic_3_regular_triangulation_traits_3;

} // namespace CGAL

#include <CGAL/internal/Periodic_3_regular_triangulation_filtered_traits_3.h>

namespace CGAL
{
template<class K_, class Off_>
class Periodic_3_regular_triangulation_traits_3<K_, Off_, false /* Has_filtered_predicates */>
  : public Periodic_3_regular_triangulation_traits_base_3<K_, Off_>
{
  typedef Periodic_3_regular_triangulation_traits_base_3<K_, Off_> Base;

public:
  typedef K_                                                       Kernel;
  typedef typename Kernel::Iso_cuboid_3                            Iso_cuboid_3;

  Periodic_3_regular_triangulation_traits_3(const Iso_cuboid_3& domain = Iso_cuboid_3(0,0,0,1,1,1),
                                            const Kernel& k = Kernel())
    : Base(domain, k)
  { }
};

template <class K_, class Off_>
class Periodic_3_regular_triangulation_traits_3<K_, Off_, true /* Has_filtered_predicates */ >
  : public Periodic_3_regular_triangulation_filtered_traits_3<K_, Off_>
{
  typedef Periodic_3_regular_triangulation_filtered_traits_3<K_, Off_> Base;

public:
  typedef K_                                                           Kernel;
  typedef typename Kernel::Iso_cuboid_3                                Iso_cuboid_3;

  Periodic_3_regular_triangulation_traits_3(const Iso_cuboid_3& domain = Iso_cuboid_3(0,0,0,1,1,1),
                                            const Kernel& k = Kernel())
    : Base(domain, k)
  { }
};

} //namespace CGAL

#endif
