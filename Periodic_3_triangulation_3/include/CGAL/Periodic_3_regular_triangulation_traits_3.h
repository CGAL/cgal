// Copyright (c) 1999-2004,2006-2009,2013-2015   INRIA Sophia-Antipolis (France).
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

#ifndef CGAL_PERIODIC_3_REGULAR_TRIANGULATION_TRAITS_3_H
#define CGAL_PERIODIC_3_REGULAR_TRIANGULATION_TRAITS_3_H

#include <CGAL/Periodic_3_offset_3.h>
#include <CGAL/Regular_traits_with_offsets_adaptor.h>
#include <CGAL/Periodic_3_construct_weighted_point_3.h>
#include <CGAL/Periodic_3_triangulation_traits_3.h>
#include <CGAL/triangulation_assertions.h>

#include <CGAL/basic.h>
#include <CGAL/representation_tags.h>
#include <CGAL/Kernel_traits.h>

namespace CGAL
{

template <class Kernel_, class Off = typename CGAL::Periodic_3_offset_3>
class Periodic_3_regular_triangulation_traits_base_3
  : public Periodic_3_triangulation_traits_base_3<Kernel_, Off>
{
private:
  typedef Periodic_3_regular_triangulation_traits_base_3<Kernel_, Off>     Self;
  typedef Periodic_3_triangulation_traits_base_3<Kernel_, Off>             Base;

public:
  typedef Kernel_                                    K;
  typedef Off                                        Offset;

  typedef typename Base::FT                          FT;
  typedef typename Base::Weighted_point_3            Weighted_point_3;
  typedef typename Base::Point_3                     Point_3;

  typedef typename Base::Periodic_3_offset_3         Periodic_3_offset_3;
  typedef typename Base::Iso_cuboid_3                Iso_cuboid_3;

  typedef typename Base::Segment_3                   Segment_3;
  typedef typename Base::Triangle_3                  Triangle_3;
  typedef typename Base::Tetrahedron_3               Tetrahedron_3;

  typedef Periodic_3_construct_weighted_point_3<Self, typename K::Construct_weighted_point_3>
                                          Construct_weighted_point_3;

  typedef Regular_traits_with_offsets_adaptor<Self, typename K::Power_side_of_oriented_power_sphere_3>
                                          Power_side_of_oriented_power_sphere_3;
  typedef Regular_traits_with_offsets_adaptor<Self, typename K::Compare_power_distance_3>
                                          Compare_power_distance_3;
  typedef Regular_traits_with_offsets_adaptor<Self, typename K::Construct_weighted_circumcenter_3>
                                          Construct_weighted_circumcenter_3;
  typedef Regular_traits_with_offsets_adaptor<Self, typename K::Compare_weighted_squared_radius_3>
                                          Compare_weighted_squared_radius_3;

  // Required for Periodic_3_regular_remove_traits
  typedef Regular_traits_with_offsets_adaptor<Self, typename K::Coplanar_orientation_3>
                                          Coplanar_orientation_3;

  // construction
  Construct_weighted_point_3 construct_weighted_point_3_object () const
  {
    return Construct_weighted_point_3(this->_domain);
  }

  Construct_weighted_circumcenter_3 construct_weighted_circumcenter_3_object () const
  {
    return Construct_weighted_circumcenter_3(&this->_domain);
  }

  // predicates
  Power_side_of_oriented_power_sphere_3 power_side_of_oriented_power_sphere_3_object () const
  {
    return Power_side_of_oriented_power_sphere_3(&this->_domain);
  }

  Compare_power_distance_3 compare_power_distance_3_object () const
  {
    return Compare_power_distance_3(&this->_domain);
  }

  Compare_weighted_squared_radius_3 compare_weighted_squared_radius_3_object() const
  {
    return Compare_weighted_squared_radius_3(&this->_domain);
  }

  Coplanar_orientation_3 coplanar_orientation_3_object() const {
    return Coplanar_orientation_3(&this->_domain);
  }
};

template<typename K,
         typename Off = CGAL::Periodic_3_offset_3,
         bool Has_filtered_predicates = K::Has_filtered_predicates>
class Periodic_3_regular_triangulation_traits_3;
} // namespace CGAL

// Partial specialization for Filtered_kernel<CK>.
#include <CGAL/internal/Periodic_3_regular_triangulation_filtered_traits_3.h>

namespace CGAL
{
// This declaration is needed to break the cyclic dependency.
template<typename K, typename Off>
class Periodic_3_regular_triangulation_filtered_traits_3;

template<class K, class Off, bool Has_filtered_predicates>
class Periodic_3_regular_triangulation_traits_3:
  public Periodic_3_regular_triangulation_traits_base_3<K, Off>
{ };

template < typename K, typename Off >
class Periodic_3_regular_triangulation_traits_3<K, Off, true>
: public Periodic_3_regular_triangulation_filtered_traits_3<K, Off>
{
public:
  typedef K Kernel;
};

} //namespace CGAL

#endif
