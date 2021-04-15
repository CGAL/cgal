// Copyright (c) 2006-2009, 2017  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
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

template <class K_,
          class Off_ = typename CGAL::Periodic_3_offset_3 >
class Periodic_3_Delaunay_triangulation_traits_base_3
  : public Periodic_3_triangulation_traits_3<K_, Off_>
{
  typedef Periodic_3_Delaunay_triangulation_traits_base_3<K_, Off_>  Self;
  typedef Periodic_3_triangulation_traits_3<K_, Off_>                Base;

public:
  typedef K_                                  Kernel;
  typedef Off_                                Offset;

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
                                                  const Kernel& k)
    : Base(domain, k)
  { }

  // Delaunay specific predicates
  typedef Functor_with_offset_points_adaptor_3<Self, typename Kernel::Side_of_oriented_sphere_3>
      Side_of_oriented_sphere_3;
  typedef Functor_with_offset_points_adaptor_3<Self, typename Kernel::Compare_distance_3>
      Compare_distance_3;

  // Required for Periodic_3_Delaunay_triangulation_remove_traits
  typedef Functor_with_offset_points_adaptor_3<Self, typename Kernel::Coplanar_orientation_3>
      Coplanar_orientation_3;
  typedef Functor_with_offset_points_adaptor_3<Self, typename Kernel::Coplanar_side_of_bounded_circle_3>
      Coplanar_side_of_bounded_circle_3;

  // When is_Gabriel is used
   typedef Functor_with_offset_points_adaptor_3<Self, typename Kernel::Side_of_bounded_sphere_3>
       Side_of_bounded_sphere_3;

  // Delaunay specific constructions
  typedef Functor_with_offset_points_adaptor_3<Self, typename Kernel::Construct_circumcenter_3>
      Construct_circumcenter_3;

  // Operations
  Side_of_oriented_sphere_3 side_of_oriented_sphere_3_object() const {
    return Side_of_oriented_sphere_3(this->Base::side_of_oriented_sphere_3_object(),
                                     this->construct_point_3_object());
  }
  Compare_distance_3 compare_distance_3_object() const {
    return Compare_distance_3(this->Base::compare_distance_3_object(),
                              this->construct_point_3_object());
  }
  Side_of_bounded_sphere_3 side_of_bounded_sphere_3_object() const {
    return Side_of_bounded_sphere_3(this->Base::side_of_bounded_sphere_3_object(),
                                    this->construct_point_3_object());
  }
  Coplanar_orientation_3 coplanar_orientation_3_object() const {
    return Coplanar_orientation_3(this->Base::coplanar_orientation_3_object(),
                                  this->construct_point_3_object());
  }
  Coplanar_side_of_bounded_circle_3 coplanar_side_of_bounded_circle_3_object() const {
    return Coplanar_side_of_bounded_circle_3(this->Base::coplanar_side_of_bounded_circle_3_object(),
                                             this->construct_point_3_object());
  }
  Construct_circumcenter_3 construct_circumcenter_3_object() const {
    return Construct_circumcenter_3(this->Base::construct_circumcenter_3_object(),
                                    this->construct_point_3_object());
  }
};

template <typename K_,
          typename Off_ = CGAL::Periodic_3_offset_3,
          bool Has_filtered_predicates_ = internal::Has_filtered_predicates<K_>::value >
class Periodic_3_Delaunay_triangulation_traits_3;

} //namespace CGAL

// Partial specialization for kernels with filtered predicates
#include <CGAL/internal/Periodic_3_Delaunay_triangulation_filtered_traits_3.h>

namespace CGAL {

template <class K_, class Off_>
class Periodic_3_Delaunay_triangulation_traits_3<K_, Off_, false>
  : public Periodic_3_Delaunay_triangulation_traits_base_3<K_, Off_>
{
  typedef Periodic_3_Delaunay_triangulation_traits_base_3<K_, Off_> Base;

public:
  typedef K_                                                        Kernel;
  typedef typename Kernel::Iso_cuboid_3                             Iso_cuboid_3;

  Periodic_3_Delaunay_triangulation_traits_3(const Iso_cuboid_3& domain = Iso_cuboid_3(0,0,0,1,1,1),
                                             const Kernel& k = Kernel())
    : Base(domain, k)
  { }
};

template <class K_, class Off_>
class Periodic_3_Delaunay_triangulation_traits_3<K_, Off_, true>
  : public Periodic_3_Delaunay_triangulation_filtered_traits_3<
             K_, Off_, internal::Has_static_filters<K_>::value>
{
  typedef Periodic_3_Delaunay_triangulation_filtered_traits_3<
            K_, Off_, internal::Has_static_filters<K_>::value>      Base;

public:
  typedef K_                                                        Kernel;
  typedef typename Kernel::Iso_cuboid_3 Iso_cuboid_3;

  Periodic_3_Delaunay_triangulation_traits_3(const Iso_cuboid_3& domain = Iso_cuboid_3(0,0,0,1,1,1),
                                             const Kernel& k = Kernel())
    : Base(domain, k)
  { }
};

} //namespace CGAL

#endif // CGAL_PERIODIC_3_DELAUNAY_TRIANGULATION_TRAITS_3_H
