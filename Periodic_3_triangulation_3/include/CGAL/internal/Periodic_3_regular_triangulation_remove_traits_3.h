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

#ifndef CGAL_PERIODIC_3_REGULAR_TRIANGULATION_REMOVE_TRAITS_3_H
#define CGAL_PERIODIC_3_REGULAR_TRIANGULATION_REMOVE_TRAITS_3_H

#include <CGAL/basic.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Periodic_3_offset_3.h>

namespace CGAL
{

template<class Traits_, class Functor_>
class Ppoint_to_point_adaptor
{
  typedef Traits_ Traits;
  typedef Functor_ Functor;

  // `Traits::Point_3` is actually a `std::pair<Point_3, Offset>`
  // `Traits::Weighted_point_3` is actually a `std::pair<Weighted_point_3, Offset>`
  typedef typename Traits::Point_3                Point_3;
  typedef typename Traits::Weighted_point_3       Weighted_point_3;

public:
  typedef typename Functor::result_type result_type;

  Ppoint_to_point_adaptor (const Functor & functor)
      : _functor(functor)
  {
  }

  result_type operator() (const Point_3& p0, const Weighted_point_3& p1,
                          const Weighted_point_3& p2) const
  {
    return _functor(p0.first, p1.first, p2.first, p0.second, p1.second, p2.second);
  }
  result_type operator() (const Weighted_point_3& p0, const Weighted_point_3& p1) const
  {
    return _functor(p0.first, p1.first, p0.second, p1.second);
  }
  result_type operator() (const Weighted_point_3& p0, const Weighted_point_3& p1,
                          const Weighted_point_3& p2) const
  {
    return _functor(p0.first, p1.first, p2.first, p0.second, p1.second, p2.second);
  }
  result_type operator() (const Weighted_point_3& p0, const Weighted_point_3& p1,
                          const Weighted_point_3& p2, const Weighted_point_3& p3) const
  {
    return _functor(p0.first, p1.first, p2.first, p3.first, p0.second, p1.second, p2.second, p3.second);
  }
  result_type operator() (const Weighted_point_3& p0, const Weighted_point_3& p1,
                          const Weighted_point_3& p2, const Weighted_point_3& p3,
                          const Weighted_point_3& p4) const
  {
    return _functor(p0.first, p1.first, p2.first, p3.first, p4.first, p0.second, p1.second, p2.second, p3.second, p4.second);
  }

private:
  Functor _functor;
};

template<class PTTraits, class Off = typename CGAL::Periodic_3_offset_3>
class Periodic_3_regular_triangulation_remove_traits_3: public PTTraits::K
{
public:
  typedef PTTraits                      PTraits;
  typedef Off                           Offset;

private:
  typedef Periodic_3_regular_triangulation_remove_traits_3<PTraits, Offset> Self;
  typedef typename PTTraits::K                                              Base;

public:
  typedef typename PTraits::RT                                  RT;
  typedef typename PTraits::FT                                  FT;
  typedef std::pair<typename PTraits::Weighted_point_3, Offset> Weighted_point_3;
  typedef std::pair<typename PTraits::Point_3, Offset>          Point_3;
  typedef typename PTraits::Iso_cuboid_3                        Iso_cuboid_3;

  Periodic_3_regular_triangulation_remove_traits_3 (const Iso_cuboid_3& domain)
      : _traits()
  {
    _traits.set_domain(domain);
  }

  // Triangulation traits
  typedef Ppoint_to_point_adaptor<Self, typename PTraits::Compare_xyz_3> Compare_xyz_3;
  typedef Ppoint_to_point_adaptor<Self, typename PTraits::Coplanar_orientation_3> Coplanar_orientation_3;
  typedef Ppoint_to_point_adaptor<Self, typename PTraits::Orientation_3> Orientation_3;

  // Regular Triangulation traits
  typedef Ppoint_to_point_adaptor<Self, typename PTraits::Power_side_of_oriented_power_sphere_3> Power_side_of_oriented_power_sphere_3;
  typedef Ppoint_to_point_adaptor<Self, typename PTraits::Compare_power_distance_3> Compare_power_distance_3;
  typedef Ppoint_to_point_adaptor<Self, typename PTraits::Power_side_of_bounded_power_sphere_3> Power_side_of_bounded_power_sphere_3;
  typedef Ppoint_to_point_adaptor<Self, typename PTraits::Construct_weighted_circumcenter_3> Construct_weighted_circumcenter_3;
  typedef Ppoint_to_point_adaptor<Self, typename PTraits::Construct_circumcenter_3> Construct_circumcenter_3;
  typedef Ppoint_to_point_adaptor<Self, typename PTraits::Compute_squared_radius_smallest_orthogonal_sphere_3> Compute_squared_radius_smallest_orthogonal_sphere_3;
  typedef Ppoint_to_point_adaptor<Self, typename PTraits::Compute_power_product_3> Compute_power_product_3;
  typedef Ppoint_to_point_adaptor<Self, typename PTraits::Compute_power_distance_to_power_sphere_3> Compute_power_distance_to_power_sphere_3;
  typedef Ppoint_to_point_adaptor<Self, typename PTraits::Compare_weighted_squared_radius_3> Compare_weighted_squared_radius_3;

  // Operations
  Compare_xyz_3 compare_xyz_3_object () const
  {
    return Compare_xyz_3(_traits.compare_xyz_3_object());
  }
  Coplanar_orientation_3 coplanar_orientation_3_object () const
  {
    return Coplanar_orientation_3(_traits.coplanar_orientation_3_object());
  }
  Orientation_3 orientation_3_object () const
  {
    return Orientation_3(_traits.orientation_3_object());
  }
  Power_side_of_oriented_power_sphere_3 power_side_of_oriented_power_sphere_3_object () const
  {
    return Power_side_of_oriented_power_sphere_3(_traits.power_side_of_oriented_power_sphere_3_object());
  }
  Compare_power_distance_3 compare_power_distance_3_object () const
  {
    return Compare_power_distance_3(_traits.compare_power_distance_3_object());
  }
  Power_side_of_bounded_power_sphere_3 power_side_of_bounded_power_sphere_3_object () const
  {
    return Power_side_of_bounded_power_sphere_3(_traits.power_side_of_bounded_power_sphere_3_object());
  }
  Construct_weighted_circumcenter_3 construct_weighted_circumcenter_3_object () const
  {
    return Construct_weighted_circumcenter_3(_traits.construct_weighted_circumcenter_3_object());
  }
  Construct_circumcenter_3 construct_circumcenter_3_object () const
  {
    return Construct_circumcenter_3(_traits.construct_circumcenter_3_object());
  }
  Compute_squared_radius_smallest_orthogonal_sphere_3 compute_squared_radius_smallest_orthogonal_sphere_3_object () const
  {
    return Compute_squared_radius_smallest_orthogonal_sphere_3(_traits.compute_squared_radius_smallest_orthogonal_sphere_3_object());
  }
  Compute_power_product_3 compute_power_product_3_object () const
  {
    return Compute_power_product_3(_traits.compute_power_product_3_object());
  }
  Compute_power_distance_to_power_sphere_3 compute_power_distance_to_power_sphere_3_object () const
  {
    return Compute_power_distance_to_power_sphere_3(_traits.compute_power_distance_to_power_sphere_3_object());
  }
  Compare_weighted_squared_radius_3 compare_weighted_squared_radius_3_object () const
  {
    return Compare_weighted_squared_radius_3(_traits.compare_weighted_squared_radius_3_object());
  }

public:
  PTraits _traits;
};

} //namespace CGAL

#endif // CGAL_PERIODIC_3_REGULAR_TRIANGULATION_REMOVE_TRAITS_3_H
