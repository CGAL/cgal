// Copyright (c) 2004,2006-2009,2017  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//                 Nico Kruithof <Nico.Kruithof@sophia.inria.fr>
//                 Manuel Caroli <Manuel.Caroli@sophia.inria.fr>


#ifndef CGAL_INTERNAL_PERIODIC_3_TRIANGULATION_FILTERED_TRAITS_3_H
#define CGAL_INTERNAL_PERIODIC_3_TRIANGULATION_FILTERED_TRAITS_3_H

#include <CGAL/license/Periodic_3_triangulation_3.h>

#include <CGAL/basic.h>
#include <CGAL/config.h>
#include <CGAL/Interval_nt.h>
#include <CGAL/Uncertain.h>
#include <CGAL/Profile_counter.h>
#include <CGAL/Filtered_predicate.h>
#include <CGAL/internal/Has_boolean_tags.h>

#include <CGAL/Periodic_3_triangulation_traits_3.h>

namespace CGAL
{
// The Offset_converter is parametrized by a usual kernel converter,
// and adds the conversions for Offsets.
template <typename Converter_>
struct Offset_converter_3
  : public Converter_
{
  typedef Converter_                        Base;

  typedef typename Base::Source_kernel      Source_kernel;
  typedef typename Base::Target_kernel      Target_kernel;

  typedef typename Periodic_3_triangulation_traits_base_3<Source_kernel>
                   ::Offset  Source_off;
  typedef typename Periodic_3_triangulation_traits_base_3<Source_kernel>
                   ::Point_3  Source_pt;

  typedef typename Periodic_3_triangulation_traits_base_3<Target_kernel>
                   ::Offset  Target_off;
  typedef typename Periodic_3_triangulation_traits_base_3<Target_kernel>
                   ::Point_3  Target_pt;

  using Base::operator();

  Target_off
  operator()(const Source_off &off) const
  {
    return off;
  }
};

// The first template item is supposed to be a Filtered_kernel-like kernel.
template <typename K_, typename Off_>
class Periodic_3_triangulation_filtered_traits_base_3
  : public Periodic_3_triangulation_traits_base_3<K_, Off_>
{
  typedef Periodic_3_triangulation_traits_base_3<K_, Off_> Base;
  typedef K_                                               Kernel;

  typedef typename Kernel::Exact_kernel                    EKernel;
  typedef typename Kernel::Approximate_kernel              AKernel;
  typedef typename Kernel::C2E                             C2E;
  typedef typename Kernel::C2F                             C2F;

  // Exact traits is based on the exact kernel.
  typedef Periodic_3_triangulation_traits_3<EKernel, Off_> Exact_traits;
  // Filtering traits is based on the filtering kernel.
  typedef Periodic_3_triangulation_traits_3<AKernel, Off_> Filtering_traits;

public:
  typedef typename Kernel::Iso_cuboid_3                         Iso_cuboid_3;

  virtual ~Periodic_3_triangulation_filtered_traits_base_3() { }

  Periodic_3_triangulation_filtered_traits_base_3(const Iso_cuboid_3& domain,
                                                  const Kernel& k)
    :
      Base(domain, k),
      traits_e(C2E()(domain)),
      traits_f(C2F()(domain))
  {
    // Problem: below is a default initialization of the kernel in the traits.
    // Hence, if the kernel has members and we use filtered traits, then
    // the members will be default constructed here...
  }

  virtual void set_domain(const Iso_cuboid_3& domain)
  {
    this->_domain = domain;
    set_filtrating_traits(domain);
  }

  void set_filtrating_traits(const Iso_cuboid_3& domain)
  {
    traits_e.set_domain(C2E()(domain));
    traits_f.set_domain(C2F()(domain));
  }

  typedef Filtered_predicate<
            typename Exact_traits::Compare_xyz_3,
            typename Filtering_traits::Compare_xyz_3,
            Offset_converter_3<C2E>,
            Offset_converter_3<C2F> >  Compare_xyz_3;

  typedef Filtered_predicate<
            typename Exact_traits::Orientation_3,
            typename Filtering_traits::Orientation_3,
            Offset_converter_3<C2E>,
            Offset_converter_3<C2F> >  Orientation_3;

  Compare_xyz_3 compare_xyz_3_object() const
  {
    typename Exact_traits::Compare_xyz_3 pe = traits_e.compare_xyz_3_object();
    typename Filtering_traits::Compare_xyz_3 pf = traits_f.compare_xyz_3_object();

    return Compare_xyz_3(pe, pf);
  }

  Orientation_3 orientation_3_object() const
  {
    typename Exact_traits::Orientation_3 pe = traits_e.orientation_3_object();
    typename Filtering_traits::Orientation_3 pf = traits_f.orientation_3_object();

    return Orientation_3(pe, pf);
  }

  // The following are inherited since they are constructions :
  // Construct_segment_3
  // Construct_triangle_3
  // Construct_tetrahedron_3

protected:
  Exact_traits traits_e;
  Filtering_traits traits_f;
};

template <class K_,
          class Off_ = typename CGAL::Periodic_3_offset_3,
          bool Has_static_filters_ = internal::Has_static_filters<K_>::value >
class Periodic_3_triangulation_filtered_traits_3;

} //namespace CGAL

#include <CGAL/internal/Periodic_3_triangulation_statically_filtered_traits_3.h>

namespace CGAL {

template<class K_, class Off_>
class Periodic_3_triangulation_filtered_traits_3<K_, Off_, false>
  : public Periodic_3_triangulation_filtered_traits_base_3<K_, Off_>
{
  typedef Periodic_3_triangulation_filtered_traits_base_3<K_, Off_> Base;

public:
  typedef K_                                                        Kernel;
  typedef typename Kernel::Iso_cuboid_3                             Iso_cuboid_3;

  Periodic_3_triangulation_filtered_traits_3(const Iso_cuboid_3& domain,
                                             const Kernel& k)
    : Base(domain, k)
  { }
};

template<class K_, class Off_>
class Periodic_3_triangulation_filtered_traits_3<K_, Off_, true>
  : public Periodic_3_triangulation_statically_filtered_traits_3<K_, Off_>
{
  typedef Periodic_3_triangulation_statically_filtered_traits_3<K_, Off_> Base;

public:
  typedef K_                                                              Kernel;
  typedef typename Kernel::Iso_cuboid_3                                   Iso_cuboid_3;

  Periodic_3_triangulation_filtered_traits_3(const Iso_cuboid_3& domain,
                                             const Kernel& k)
    : Base(domain, k)
  { }
};

} //namespace CGAL

#endif // CGAL_INTERNAL_PERIODIC_3_TRIANGULATION_FILTERED_TRAITS_3_H
