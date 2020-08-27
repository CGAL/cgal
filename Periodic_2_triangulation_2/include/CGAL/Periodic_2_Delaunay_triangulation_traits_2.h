// Copyright (c) 1997-2013, 2017 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Nico Kruithof <Nico@nghk.nl>
//                 Mael Rouxel-Labbé

#ifndef CGAL_PERIODIC_2_DELAUNAY_TRIANGULATION_TRAITS_2_H
#define CGAL_PERIODIC_2_DELAUNAY_TRIANGULATION_TRAITS_2_H

#include <CGAL/license/Periodic_2_triangulation_2.h>

#include <CGAL/Periodic_2_triangulation_2/internal/Functor_with_offset_points_adaptor_2.h>
#include <CGAL/Periodic_2_triangulation_2/internal/traits_helpers.h>
#include <CGAL/Periodic_2_offset_2.h>
#include <CGAL/Periodic_2_triangulation_traits_2.h>

#include <CGAL/internal/Has_boolean_tags.h>

namespace CGAL {

template <class Kernel_,
          class Offset_ = CGAL::Periodic_2_offset_2,
          class Domain_ = typename Kernel_::Iso_rectangle_2>
class Periodic_2_Delaunay_triangulation_traits_base_2
  : public Periodic_2_triangulation_traits_2<Kernel_, Offset_, Domain_>
{
public:
  typedef Kernel_                                                     Kernel;
  typedef Offset_                                                     Offset;
  typedef Domain_                                                     Domain;

private:
  typedef Periodic_2_Delaunay_triangulation_traits_base_2<Kernel, Offset, Domain>        Self;
  typedef Periodic_2_triangulation_traits_2<Kernel, Offset, Domain>                      Base;

public:
  Periodic_2_Delaunay_triangulation_traits_base_2(const Domain& domain = Domain(),
                                                  const Kernel& k = Kernel())
    : Base(domain, k)
  {
    std::cout << "create 'Periodic_2_Delaunay_triangulation_traits_base_2'" << std::endl;
  }

  typedef typename P2T2::internal::Construct_point_getter<
                                      Kernel, Offset, Domain>::type   Construct_point_2;

  // Predicates
  typedef Functor_with_offset_points_adaptor_2<Self, typename Kernel::Compare_distance_2>
                                                                      Compare_distance_2;
  typedef Functor_with_offset_points_adaptor_2<Self, typename Kernel::Side_of_oriented_circle_2>
                                                                      Side_of_oriented_circle_2;
  typedef Functor_with_offset_points_adaptor_2<Self, typename Kernel::Compare_squared_radius_2>
                                                                      Compare_squared_radius_2;

  // Constructions
  typedef Functor_with_offset_points_adaptor_2<Self, typename Kernel::Construct_circumcenter_2>
                                                                      Construct_circumcenter_2;

  // Predicates
  Compare_distance_2 compare_distance_2_object() const
  {
    return Compare_distance_2(this->Base::compare_distance_2_object(),
                              this->construct_point_2_object());
  }
  Side_of_oriented_circle_2 side_of_oriented_circle_2_object() const
  {
    return Side_of_oriented_circle_2(this->Base::side_of_oriented_circle_2_object(),
                                     this->construct_point_2_object());
  }
  Compare_squared_radius_2 compare_squared_radius_2_object() const
  {
     return Compare_squared_radius_2(this->Base::compare_squared_radius_2_object(),
                                     this->construct_point_2_object());
  }

  // Constructions
  Construct_circumcenter_2 construct_circumcenter_2_object() const
  {
    return Construct_circumcenter_2(this->Base::construct_circumcenter_2_object(),
                                    this->construct_point_2_object());
  }
};

template <class K,
          class O = CGAL::Periodic_2_offset_2,
          class D = typename K::Iso_rectangle_2,
          bool Has_filtered_predicates =
            boost::mpl::and_<internal::Has_filtered_predicates<K>,
                             boost::mpl::not_<std::is_same<D, Lattice_2<K> > > >::value>
class Periodic_2_Delaunay_triangulation_traits_2;

} // namespace CGAL

// Partial specialization for kernels with filtered predicates
#include <CGAL/Periodic_2_triangulation_2/internal/Periodic_2_Delaunay_triangulation_filtered_traits_2.h>

namespace CGAL {

template <class K, class O, class D>
class Periodic_2_Delaunay_triangulation_traits_2<K, O, D, false>
  : public Periodic_2_Delaunay_triangulation_traits_base_2<K, O, D>
{
  typedef Periodic_2_Delaunay_triangulation_traits_base_2<K, O, D> Base;

public:
  Periodic_2_Delaunay_triangulation_traits_2(const D& domain  = P2T2::internal::get_default_domain<K, D>(),
                                             const K& k = K())
    : Base(domain, k)
  { }
};

// Below is filtered, but might even be statically filtered, depending on 'Has_static_filters'
template <class K, class O, class D>
class Periodic_2_Delaunay_triangulation_traits_2<K, O, D, true>
  : public Periodic_2_Delaunay_triangulation_filtered_traits_2<
             K, O, D, internal::Has_static_filters<K>::value>
{
  typedef Periodic_2_Delaunay_triangulation_filtered_traits_2<
            K, O, D, internal::Has_static_filters<K>::value>    Base;

public:
  Periodic_2_Delaunay_triangulation_traits_2(const D& domain = P2T2::internal::get_default_domain<K, D>(),
                                             const K& k = K())
    : Base(domain, k)
  { }
};

} // namespace CGAL

#endif // CGAL_PERIODIC_2_DELAUNAY_TRIANGULATION_TRAITS_2_H
