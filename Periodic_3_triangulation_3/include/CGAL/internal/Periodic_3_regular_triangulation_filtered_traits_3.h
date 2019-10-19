// Copyright (c) 1999-2004,2006-2009,2013-2015, 2017  INRIA Sophia-Antipolis (France).
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

#ifndef CGAL_PERIODIC_3_REGULAR_TRIANGULATION_FILTERED_TRAITS_3_H
#define CGAL_PERIODIC_3_REGULAR_TRIANGULATION_FILTERED_TRAITS_3_H

#include <CGAL/license/Periodic_3_triangulation_3.h>

#include <CGAL/internal/Periodic_3_triangulation_filtered_traits_3.h>
#include <CGAL/Periodic_3_regular_triangulation_traits_3.h>

#include <CGAL/basic.h>
#include <CGAL/config.h>
#include <CGAL/internal/Has_boolean_tags.h>
#include <CGAL/Interval_nt.h>
#include <CGAL/Uncertain.h>
#include <CGAL/Profile_counter.h>

namespace CGAL {

// The first template item is supposed to be a Filtered_kernel-like kernel.
template <typename K_, typename Off_>
class Periodic_3_regular_triangulation_filtered_traits_base_3
  : public Periodic_3_regular_triangulation_traits_base_3<K_, Off_>
{
  typedef Periodic_3_regular_triangulation_traits_base_3<K_, Off_>  Base;
  typedef K_                                                        Kernel;

  typedef typename Kernel::Exact_kernel                             EKernel;
  typedef typename Kernel::Approximate_kernel                       AKernel;
  typedef typename Kernel::C2E                                      C2E;
  typedef typename Kernel::C2F                                      C2F;

  // Exact traits is based on the exact kernel.
  typedef Periodic_3_regular_triangulation_traits_3<EKernel, Off_>  Exact_traits;
  // Filtering traits is based on the filtering kernel.
  typedef Periodic_3_regular_triangulation_traits_3<AKernel, Off_>  Filtering_traits;

  typedef typename Kernel::Iso_cuboid_3                             Iso_cuboid_3;

public:
  virtual ~Periodic_3_regular_triangulation_filtered_traits_base_3() { }

  Periodic_3_regular_triangulation_filtered_traits_base_3(const Iso_cuboid_3& domain,
                                                          const Kernel& k)
    :
      Base(domain, k),
      regular_traits_e(C2E()(domain)),
      regular_traits_f(C2F()(domain))
  {
    // Problem 1: above is a default initialization of the kernel in the traits.
    // Hence, if the kernel has members and we use filtered traits, then
    // the members will be default constructed here...

    // Problem 2: we have built filtered traits in P3Tfiltered_traits_3 and now
    // we also need those two...
  }

  virtual void set_domain(const Iso_cuboid_3& domain)
  {
    this->_domain = domain;
    this->set_filtrating_traits(domain);
    set_filtrating_regular_traits(domain);
  }

  void set_filtrating_regular_traits(const Iso_cuboid_3& domain)
  {
    regular_traits_e.set_domain(C2E()(domain));
    regular_traits_f.set_domain(C2F()(domain));
  }

public:
  typedef Filtered_predicate<
            typename Exact_traits::Coplanar_orientation_3,
            typename Filtering_traits::Coplanar_orientation_3,
            Offset_converter_3<C2E>,
            Offset_converter_3<C2F> >  Coplanar_orientation_3;

  typedef Filtered_predicate<
            typename Exact_traits::Power_side_of_oriented_power_sphere_3,
            typename Filtering_traits::Power_side_of_oriented_power_sphere_3,
            Offset_converter_3<C2E>,
            Offset_converter_3<C2F> >  Power_side_of_oriented_power_sphere_3;

  typedef Filtered_predicate<
            typename Exact_traits::Compare_power_distance_3,
            typename Filtering_traits::Compare_power_distance_3,
            Offset_converter_3<C2E>,
            Offset_converter_3<C2F> >  Compare_power_distance_3;

  typedef Filtered_predicate<
            typename Exact_traits::Compare_weighted_squared_radius_3,
            typename Filtering_traits::Compare_weighted_squared_radius_3,
            Offset_converter_3<C2E>,
            Offset_converter_3<C2F> >  Compare_weighted_squared_radius_3;

  Coplanar_orientation_3 coplanar_orientation_3_object() const
  {
    typename Exact_traits::Coplanar_orientation_3 pe = regular_traits_e.coplanar_orientation_3_object();
    typename Filtering_traits::Coplanar_orientation_3 pf = regular_traits_f.coplanar_orientation_3_object();

    return Coplanar_orientation_3(pe, pf);
  }

  Power_side_of_oriented_power_sphere_3 power_side_of_oriented_power_sphere_3_object() const
  {
    typename Exact_traits::Power_side_of_oriented_power_sphere_3 pe = regular_traits_e.power_side_of_oriented_power_sphere_3_object();
    typename Filtering_traits::Power_side_of_oriented_power_sphere_3 pf = regular_traits_f.power_side_of_oriented_power_sphere_3_object();

    return Power_side_of_oriented_power_sphere_3(pe, pf);
  }

  Compare_power_distance_3 compare_power_distance_3_object() const
  {
    typename Exact_traits::Compare_power_distance_3 pe = regular_traits_e.compare_power_distance_3_object();
    typename Filtering_traits::Compare_power_distance_3 pf = regular_traits_f.compare_power_distance_3_object();

    return Compare_power_distance_3(pe, pf);
  }

  Compare_weighted_squared_radius_3 compare_weighted_squared_radius_3_object() const
  {
    typename Exact_traits::Compare_weighted_squared_radius_3 pe = regular_traits_e.compare_weighted_squared_radius_3_object();
    typename Filtering_traits::Compare_weighted_squared_radius_3 pf = regular_traits_f.compare_weighted_squared_radius_3_object();

    return Compare_weighted_squared_radius_3(pe, pf);
  }

protected:
  Exact_traits regular_traits_e;
  Filtering_traits regular_traits_f;
};

} // namespace CGAL

#include <CGAL/internal/Periodic_3_regular_triangulation_statically_filtered_traits_3.h>

namespace CGAL {

template <typename K_,
          typename Off_ = typename CGAL::Periodic_3_offset_3,
          bool Has_static_filters_ = internal::Has_static_filters<K_>::value >
class Periodic_3_regular_triangulation_filtered_traits_3
  : public Periodic_3_regular_triangulation_statically_filtered_traits_3<K_, Off_>
{
  typedef Periodic_3_regular_triangulation_statically_filtered_traits_3<K_, Off_> Base;

public:
  typedef K_                                               Kernel;
  typedef typename Kernel::Iso_cuboid_3                    Iso_cuboid_3;

  Periodic_3_regular_triangulation_filtered_traits_3(const Iso_cuboid_3& domain,
                                                     const Kernel& k)
    : Base(domain, k)
  { }
};

template <typename K_, typename Off_>
class Periodic_3_regular_triangulation_filtered_traits_3<K_, Off_, false>
  : public Periodic_3_regular_triangulation_filtered_traits_base_3<K_, Off_>
{
  typedef Periodic_3_regular_triangulation_filtered_traits_base_3<K_, Off_> Base;

public:
  typedef K_                                               Kernel;
  typedef typename Kernel::Iso_cuboid_3                    Iso_cuboid_3;

  Periodic_3_regular_triangulation_filtered_traits_3(const Iso_cuboid_3& domain,
                                                     const Kernel& k)
    : Base(domain, k)
  { }
};

} //namespace CGAL

#endif // CGAL_PERIODIC_3_REGULAR_TRIANGULATION_FILTERED_TRAITS_3_H
