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

template <class K_,
          class Off_ = typename CGAL::Periodic_3_offset_3>
class Periodic_3_triangulation_traits_base_3
  : public K_
{
  typedef Periodic_3_triangulation_traits_base_3<K_, Off_>  Self;

public:
  typedef K_                                 Kernel;
  typedef Off_                               Offset;

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

public:
  virtual ~Periodic_3_triangulation_traits_base_3() { }

  Periodic_3_triangulation_traits_base_3(const Iso_cuboid_3& domain,
                                         const Kernel& k)
    : Kernel(k)
  {
    set_domain(domain);
  }

  // will be overwritten by filtered classes to create exact and approximate
  // versions of the domain
  virtual void set_domain(const Iso_cuboid_3& domain) {
    _domain = domain;
  }

  // Access
  const Iso_cuboid_3& get_domain() const {
    return _domain;
  }

  // Construct_point_3 with offset
  typedef Periodic_3_construct_point_3<Self, typename Kernel::Construct_point_3>
      Construct_point_3;

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
    return Construct_point_3(&_domain, this->Kernel::construct_point_3_object());
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
  Iso_cuboid_3 _domain;
};

template <class K_,
          class Off_ = CGAL::Periodic_3_offset_3,
          bool Has_filtered_predicates_ = internal::Has_filtered_predicates<K_>::value>
class Periodic_3_triangulation_traits_3;

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

} //namespace CGAL

#endif // CGAL_PERIODIC_3_TRIANGULATION_TRAITS_3_H
