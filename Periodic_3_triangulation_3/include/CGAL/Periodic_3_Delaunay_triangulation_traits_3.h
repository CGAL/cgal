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

template <class Kernel_,
          class Offset_ = CGAL::Periodic_3_offset_3,
          class Domain_ = typename Kernel_::Iso_cuboid_3,
          class Construct_point_3_ = Default>
class Periodic_3_Delaunay_triangulation_traits_base_3
    : public Periodic_3_triangulation_traits_3<
               Kernel_, Offset_, Domain_,
               typename Default::Get<Construct_point_3_,
                                     Periodic_3_construct_point_3<Kernel_, Offset_> >::type>
{
public:
  typedef Kernel_                                                 Kernel;
  typedef Offset_                                                 Offset;
  typedef Domain_                                                 Domain;

  typedef Periodic_3_construct_point_3<Kernel, Offset>            Construct_point_3_def;
  typedef typename Default::Get<Construct_point_3_,
                                Construct_point_3_def>::type      Construct_point_3;

private:
  typedef Periodic_3_Delaunay_triangulation_traits_base_3<
            Kernel, Offset, Domain, Construct_point_3>            Self;
  typedef Periodic_3_triangulation_traits_3<
            Kernel, Offset, Domain, Construct_point_3>            Base;

public:
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
  Periodic_3_Delaunay_triangulation_traits_base_3(const Domain& domain = Domain(),
                                                  const Kernel& k = Kernel())
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

template <class K,
          class O = CGAL::Periodic_3_offset_3,
          class D = typename K::Iso_cuboid_3,
          class CP = Default,
          bool Has_filtered_predicates_ = internal::Has_filtered_predicates<K>::value >
class Periodic_3_Delaunay_triangulation_traits_3
  : public Periodic_3_Delaunay_triangulation_traits_base_3<K, O, D, CP>
{
public:
  typedef Periodic_3_Delaunay_triangulation_traits_base_3<K, O, D, CP> Base;
  Periodic_3_Delaunay_triangulation_traits_3(const D& d = D(), const K& k = K()) : Base(d, k) { }
};

#ifdef MACRO_THAT_DOESNT_EXIT_TO_MAKE_GP2T2_WORK // @tmp this is just a forward declaration normally

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

#endif // MACRO_THAT_DOESNT_EXIT_TO_MAKE_GP2T2_WORK

} //namespace CGAL

#endif // CGAL_PERIODIC_3_DELAUNAY_TRIANGULATION_TRAITS_3_H
