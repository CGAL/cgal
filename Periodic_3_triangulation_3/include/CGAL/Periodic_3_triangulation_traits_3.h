// Copyright (c) 2006-2009,2017  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Nico Kruithof <Nico.Kruithof@sophia.inria.fr>
//                 Manuel Caroli <Manuel.Caroli@sophia.inria.fr>
//                 Mael Rouxel-Labbé

#ifndef CGAL_PERIODIC_3_TRIANGULATION_TRAITS_3_H
#define CGAL_PERIODIC_3_TRIANGULATION_TRAITS_3_H

#include <CGAL/license/Periodic_3_triangulation_3.h>

#include <CGAL/internal/Periodic_3_construct_point_3.h>
#include <CGAL/internal/Functor_with_offset_points_adaptor_3.h>
#include <CGAL/Periodic_3_offset_3.h>

#include <CGAL/basic.h>
#include <CGAL/internal/Has_boolean_tags.h>
#include <CGAL/triangulation_assertions.h>

namespace CGAL {

template <class Kernel_,
          class Offset_ = CGAL::Periodic_3_offset_3,
          class Domain_ = typename Kernel_::Iso_cuboid_3,
          class Construct_point_3_ = Default>
class Periodic_3_triangulation_traits_base_3
  : public Kernel_
{
public:
  typedef Kernel_                                                 Kernel;
  typedef Offset_                                                 Offset;

  typedef Domain_                                                 Domain;

  typedef typename Kernel::RT                RT;
  typedef typename Kernel::FT                FT;
  typedef typename Kernel::Point_3           Point_3;
  typedef typename Kernel::Vector_3          Vector_3;
  typedef Offset                             Periodic_3_offset_3;
  typedef typename Kernel::Iso_cuboid_3      Iso_cuboid_3;

  // The next typedef is there for backward compatibility
  // Some users take their point type from the traits class.
  // Before this type was Point
  typedef Point_3                            Point;

  typedef typename Kernel::Segment_3         Segment_3;
  typedef typename Kernel::Triangle_3        Triangle_3;
  typedef typename Kernel::Tetrahedron_3     Tetrahedron_3;

  typedef Periodic_3_construct_point_3<Kernel, Offset>            Construct_point_3_def;
  typedef typename Default::Get<Construct_point_3_,
                                Construct_point_3_def>::type      Construct_point_3;

private:
  // Not truly "Self" since we grab the default value from CGAL::Default
  typedef Periodic_3_triangulation_traits_base_3<
            Kernel_, Offset_, Domain_, Construct_point_3>         Self;
  typedef Kernel_                                                 Base;

public:
  virtual ~Periodic_3_triangulation_traits_base_3() { }

  Periodic_3_triangulation_traits_base_3(const Domain& d = Domain(),
                                         const Kernel& k = Kernel())
    : Base(k), domain(d)
  { }

  // will be overwritten by filtered classes to create exact and approximate
  // versions of the domain
  virtual void set_domain(const Domain& d) { domain = d; }

  // Access
  const Domain& get_domain() const { return domain; }

  // Triangulation predicates
  typedef Functor_with_offset_points_adaptor_3<Self, typename Kernel::Compare_xyz_3>
      Compare_xyz_3;
  typedef Functor_with_offset_points_adaptor_3<Self, typename Kernel::Orientation_3>
      Orientation_3;

  // Triangulation constructions
  typedef Functor_with_offset_points_adaptor_3<Self, typename Kernel::Construct_segment_3>
      Construct_segment_3;
  typedef Functor_with_offset_points_adaptor_3<Self, typename Kernel::Construct_triangle_3>
      Construct_triangle_3;
  typedef Functor_with_offset_points_adaptor_3<Self, typename Kernel::Construct_tetrahedron_3>
      Construct_tetrahedron_3;

  // Operations
  Construct_point_3 construct_point_3_object() const {
    return Construct_point_3(&domain, this->Kernel::construct_point_3_object());
  }

  Compare_xyz_3 compare_xyz_3_object() const {
    return Compare_xyz_3(this->Kernel::compare_xyz_3_object(), construct_point_3_object());
  }
  Orientation_3 orientation_3_object() const {
    return Orientation_3(this->Kernel::orientation_3_object(), construct_point_3_object());
  }
  Construct_segment_3 construct_segment_3_object() const {
    return Construct_segment_3(this->Kernel::construct_segment_3_object(), construct_point_3_object());
  }
  Construct_triangle_3 construct_triangle_3_object() const {
    return Construct_triangle_3(this->Kernel::construct_triangle_3_object(), construct_point_3_object());
  }
  Construct_tetrahedron_3 construct_tetrahedron_3_object() const {
    return Construct_tetrahedron_3(this->Kernel::construct_tetrahedron_3_object(), construct_point_3_object());
  }

protected:
  Domain domain;
};

template <class K,
          class O = CGAL::Periodic_3_offset_3,
          class D = typename K::Iso_cuboid_3,
          class CP = Default,
          bool Has_filtered_predicates = internal::Has_filtered_predicates<K>::value >
class Periodic_3_triangulation_traits_3
  : public Periodic_3_triangulation_traits_base_3<K, O, D, CP> // @tmp (this is just a forward declaration normally)
{
public:
  typedef Periodic_3_triangulation_traits_base_3<K, O, D, CP> Base;
  Periodic_3_triangulation_traits_3(const D& d = D(), const K& k = K()) : Base(d, k) { }
};


#if 0 // @tmp to remove filtering

} //namespace CGAL

#include <CGAL/internal/Periodic_3_triangulation_filtered_traits_3.h>

namespace CGAL {

template <class K_, class Off_>
class Periodic_3_triangulation_traits_3<K_, Off_, false>
  : public Periodic_3_triangulation_traits_base_3<K_, Off_>
{
  typedef Periodic_3_triangulation_traits_base_3<K_, Off_> Base;

public:
  typedef K_                                               Kernel;
  typedef typename Kernel::Iso_cuboid_3                    Iso_cuboid_3;

  Periodic_3_triangulation_traits_3(const Iso_cuboid_3& domain = Iso_cuboid_3(0,0,0,1,1,1),
                                    const Kernel& k = Kernel())
    : Base(domain, k)
  { }
};

template <class K_, class Off_>
class Periodic_3_triangulation_traits_3<K_, Off_, true>
  : public Periodic_3_triangulation_filtered_traits_3<
             K_, Off_, internal::Has_static_filters<K_>::value>
{
  typedef Periodic_3_triangulation_filtered_traits_3<
            K_, Off_, internal::Has_static_filters<K_>::value>  Base;

public:
  typedef K_                                                   Kernel;
  typedef typename Kernel::Iso_cuboid_3                        Iso_cuboid_3;

  Periodic_3_triangulation_traits_3(const Iso_cuboid_3& domain = Iso_cuboid_3(0,0,0,1,1,1),
                                    const Kernel& k = Kernel())
    : Base(domain, k)
  { }
};

#endif // 

} //namespace CGAL

#endif // CGAL_PERIODIC_3_TRIANGULATION_TRAITS_3_H
