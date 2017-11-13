// Copyright (c) 2004,2006-2009,2017  INRIA Sophia-Antipolis (France).
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
template < typename Converter >
struct Offset_converter_3
  : public Converter
{
  typedef typename Converter::Source_kernel Source_kernel;
  typedef typename Converter::Target_kernel Target_kernel;

  typedef typename Periodic_3_triangulation_traits_base_3<Source_kernel>
                   ::Offset  Source_off;
  typedef typename Periodic_3_triangulation_traits_base_3<Source_kernel>
                   ::Point_3  Source_pt;

  typedef typename Periodic_3_triangulation_traits_base_3<Target_kernel>
                   ::Offset  Target_off;
  typedef typename Periodic_3_triangulation_traits_base_3<Target_kernel>
                   ::Point_3  Target_pt;

  using Converter::operator();

  Target_off
  operator()(const Source_off &off) const
  {
    return off;
  }
};

// The first template item is supposed to be a Filtered_kernel-like kernel.
template < typename K, typename Off >
class Periodic_3_triangulation_filtered_traits_base_3
  : public Periodic_3_triangulation_traits_base_3<K, Off>
{
  typedef Periodic_3_triangulation_traits_base_3<K, Off>   Base;

  typedef typename K::Exact_kernel                         EKernel;
  typedef typename K::Approximate_kernel                   AKernel;
  typedef typename K::C2E                                  C2E;
  typedef typename K::C2F                                  C2F;

  // Exact traits is based on the exact kernel.
  typedef Periodic_3_triangulation_traits_3<EKernel, Off>  Exact_traits;
  // Filtering traits is based on the filtering kernel.
  typedef Periodic_3_triangulation_traits_3<AKernel, Off>  Filtering_traits;

public:
  typedef typename K::Iso_cuboid_3                         Iso_cuboid_3;

  virtual ~Periodic_3_triangulation_filtered_traits_base_3() { }

  Periodic_3_triangulation_filtered_traits_base_3(const Iso_cuboid_3& domain,
                                                  const K& k)
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

template < typename K,
           typename Off = typename CGAL::Periodic_3_offset_3,
           bool Has_static_filters = internal::Has_static_filters<K>::value >
class Periodic_3_triangulation_filtered_traits_3;

} //namespace CGAL

#include <CGAL/internal/Periodic_3_triangulation_statically_filtered_traits_3.h>

namespace CGAL {

template<class K, class Off >
class Periodic_3_triangulation_filtered_traits_3<K, Off, false>
  : public Periodic_3_triangulation_filtered_traits_base_3<K, Off>
{
  typedef Periodic_3_triangulation_filtered_traits_base_3<K, Off> Base;

public:
  typedef typename K::Iso_cuboid_3 Iso_cuboid_3;

  Periodic_3_triangulation_filtered_traits_3(const Iso_cuboid_3& domain,
                                             const K& k)
    : Base(domain, k)
  { }
};

template<class K, class Off>
class Periodic_3_triangulation_filtered_traits_3<K, Off, true>
  : public Periodic_3_triangulation_statically_filtered_traits_3<K, Off>
{
  typedef Periodic_3_triangulation_statically_filtered_traits_3<K, Off> Base;

public:
  typedef typename K::Iso_cuboid_3 Iso_cuboid_3;

  Periodic_3_triangulation_filtered_traits_3(const Iso_cuboid_3& domain,
                                             const K& k)
    : Base(domain, k)
  { }
};

} //namespace CGAL

#endif // CGAL_INTERNAL_PERIODIC_3_TRIANGULATION_FILTERED_TRAITS_3_H
