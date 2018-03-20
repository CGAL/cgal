// Copyright (c) 2006-2009,2017  INRIA Sophia-Antipolis (France).
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

template < class K, class Off = typename CGAL::Periodic_3_offset_3 >
class Periodic_3_triangulation_traits_base_3
  : public K
{
  typedef Periodic_3_triangulation_traits_base_3<K, Off>  Self;
  typedef K                                               Base;

public:
  typedef K                             Kernel;
  typedef Off                           Offset;

  typedef typename K::RT                RT;
  typedef typename K::FT                FT;
  typedef typename K::Point_3           Point_3;
  typedef typename K::Vector_3          Vector_3;
  typedef Offset                        Periodic_3_offset_3;
  typedef typename K::Iso_cuboid_3      Iso_cuboid_3;

  // The next typedef is there for backward compatibility
  // Some users take their point type from the traits class.
  // Before this type was Point
  typedef Point_3                       Point;

  typedef typename K::Segment_3         Segment_3;
  typedef typename K::Triangle_3        Triangle_3;
  typedef typename K::Tetrahedron_3     Tetrahedron_3;

public:
  virtual ~Periodic_3_triangulation_traits_base_3() { }

  Periodic_3_triangulation_traits_base_3(const Iso_cuboid_3& domain,
                                         const K& k)
    : Base(k)
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
  typedef Periodic_3_construct_point_3<Self, typename K::Construct_point_3>
      Construct_point_3;

  // Triangulation predicates
  typedef Functor_with_offset_points_adaptor_3<Self, typename K::Compare_xyz_3>
      Compare_xyz_3;
  typedef Functor_with_offset_points_adaptor_3<Self, typename K::Orientation_3>
      Orientation_3;

  // Triangulation constructions
  typedef Functor_with_offset_points_adaptor_3<Self, typename K::Construct_segment_3>
      Construct_segment_3;
  typedef Functor_with_offset_points_adaptor_3<Self, typename K::Construct_triangle_3>
      Construct_triangle_3;
  typedef Functor_with_offset_points_adaptor_3<Self, typename K::Construct_tetrahedron_3>
      Construct_tetrahedron_3;

  // Operations
  Construct_point_3 construct_point_3_object() const {
    return Construct_point_3(&_domain, this->K::construct_point_3_object());
  }

  Compare_xyz_3 compare_xyz_3_object() const {
    return Compare_xyz_3(this->K::compare_xyz_3_object(), construct_point_3_object());
  }
  Orientation_3 orientation_3_object() const {
    return Orientation_3(this->K::orientation_3_object(), construct_point_3_object());
  }
  Construct_segment_3 construct_segment_3_object() const {
    return Construct_segment_3(this->K::construct_segment_3_object(), construct_point_3_object());
  }
  Construct_triangle_3 construct_triangle_3_object() const {
    return Construct_triangle_3(this->K::construct_triangle_3_object(), construct_point_3_object());
  }
  Construct_tetrahedron_3 construct_tetrahedron_3_object() const {
    return Construct_tetrahedron_3(this->K::construct_tetrahedron_3_object(), construct_point_3_object());
  }

protected:
  Iso_cuboid_3 _domain;
};

template < typename K,
           typename Off = CGAL::Periodic_3_offset_3,
           bool Has_filtered_predicates = internal::Has_filtered_predicates<K>::value >
class Periodic_3_triangulation_traits_3;

} //namespace CGAL

#include <CGAL/internal/Periodic_3_triangulation_filtered_traits_3.h>

namespace CGAL {

template < class K, class Off >
class Periodic_3_triangulation_traits_3<K, Off, false>
  : public Periodic_3_triangulation_traits_base_3<K, Off>
{
  typedef Periodic_3_triangulation_traits_base_3<K, Off> Base;

public:
  typedef typename K::Iso_cuboid_3 Iso_cuboid_3;

  Periodic_3_triangulation_traits_3(const Iso_cuboid_3& domain = Iso_cuboid_3(0,0,0,1,1,1),
                                    const K& k = K())
    : Base(domain, k)
  { }
};

template < class K, class Off >
class Periodic_3_triangulation_traits_3<K, Off, true>
  : public Periodic_3_triangulation_filtered_traits_3<
             K, Off, internal::Has_static_filters<K>::value>
{
  typedef Periodic_3_triangulation_filtered_traits_3<
            K, Off, internal::Has_static_filters<K>::value>  Base;

public:
  typedef typename K::Iso_cuboid_3 Iso_cuboid_3;

  Periodic_3_triangulation_traits_3(const Iso_cuboid_3& domain = Iso_cuboid_3(0,0,0,1,1,1),
                                    const K& k = K())
    : Base(domain, k)
  { }
};

} //namespace CGAL

#endif // CGAL_PERIODIC_3_TRIANGULATION_TRAITS_3_H
