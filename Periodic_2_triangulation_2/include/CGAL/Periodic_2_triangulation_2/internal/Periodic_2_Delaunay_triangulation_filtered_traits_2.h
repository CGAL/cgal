// Copyright (c) 2017 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_PERIODIC_2_DELAUNAY_TRIANGULATION_FILTERED_TRAITS_2_H
#define CGAL_PERIODIC_2_DELAUNAY_TRIANGULATION_FILTERED_TRAITS_2_H

#include <CGAL/license/Periodic_2_triangulation_2.h>

#include <CGAL/Periodic_2_triangulation_2/internal/traits_helpers.h>
#include <CGAL/Periodic_2_triangulation_2/internal/Periodic_2_triangulation_filtered_traits_2.h>
#include <CGAL/Periodic_2_Delaunay_triangulation_traits_2.h>

#include <CGAL/internal/Has_boolean_tags.h>

namespace CGAL {

// The first template item is supposed to be a Filtered_kernel-like kernel.
template <typename Kernel_,
          typename Offset_ = CGAL::Periodic_2_offset_2,
          class Domain_ = typename Kernel_::Iso_rectangle_2>
class Periodic_2_Delaunay_triangulation_filtered_traits_base_2
  : public Periodic_2_Delaunay_triangulation_traits_base_2<Kernel_, Offset_, Domain_>
{
public:
  typedef Kernel_                                                                        Kernel;
  typedef Offset_                                                                        Offset;
  typedef Domain_                                                                        Domain;

private:
  typedef Periodic_2_Delaunay_triangulation_traits_base_2<Kernel, Offset, Domain>        Base;

  typedef typename Kernel::Exact_kernel                                                  EKernel;
  typedef typename Kernel::Approximate_kernel                                            AKernel;
  typedef typename Kernel::C2E                                                           C2E;
  typedef typename Kernel::C2F                                                           C2F;

  typedef typename P2T2::internal::Exact_domain_getter<Kernel, Domain, EKernel>::type    EDomain;
  typedef typename P2T2::internal::Exact_domain_getter<Kernel, Domain, AKernel>::type    ADomain;

public:
  // Exact traits is based on the exact kernel.
  typedef Periodic_2_Delaunay_triangulation_traits_2<EKernel, Offset, EDomain>           Exact_traits;
  // Filtering traits is based on the filtering kernel.
  typedef Periodic_2_Delaunay_triangulation_traits_2<AKernel, Offset, ADomain>           Filtering_traits;

public:
  Periodic_2_Delaunay_triangulation_filtered_traits_base_2(const Domain& domain, const Kernel& k)
    :
      Base(domain, k),
      Delaunay_traits_e(P2T2::internal::convert_domain<Kernel, EKernel>(domain)),
      Delaunay_traits_f(P2T2::internal::convert_domain<Kernel, AKernel>(domain))
  {
    // @fixme
    // Problem 1: above is a default initialization of the kernel in the traits.
    // Hence, if the kernel has members and we use filtered traits, then
    // the members will be default constructed here...
    //
    // Problem-ish 2: we have built filtered traits in P2Tfiltered_traits_2 and now
    // we also need two extra exact+approx domains...

    std::cout << "create 'Periodic_2_Delaunay_triangulation_filtered_traits_base_2'" << std::endl;
  }

  void set_filtrating_Delaunay_traits(const Domain& domain)
  {
    Delaunay_traits_e.set_domain(P2T2::internal::convert_domain<Kernel, EKernel>(domain));
    Delaunay_traits_f.set_domain(P2T2::internal::convert_domain<Kernel, AKernel>(domain));
  }

  void set_domain(const Domain& domain)
  {
    Base::set_domain(domain);
    Base::set_filtrating_traits(domain);
    set_filtrating_Delaunay_traits(domain);
  }

public:
  typedef typename P2T2::internal::Construct_point_getter<
                                     Kernel, Offset, Domain>::type     Construct_point_2;

  typedef Filtered_predicate<typename Exact_traits::Compare_distance_2,
                             typename Filtering_traits::Compare_distance_2,
                             P2T2::internal::Offset_converter_2<C2E>,
                             P2T2::internal::Offset_converter_2<C2F> > Compare_distance_2;
  typedef Filtered_predicate<typename Exact_traits::Side_of_oriented_circle_2,
                             typename Filtering_traits::Side_of_oriented_circle_2,
                             P2T2::internal::Offset_converter_2<C2E>,
                             P2T2::internal::Offset_converter_2<C2F> > Side_of_oriented_circle_2;

  Compare_distance_2 compare_distance_2_object() const
  {
    typename Exact_traits::Compare_distance_2 pe = Delaunay_traits_e.compare_distance_2_object();
    typename Filtering_traits::Compare_distance_2 pf = Delaunay_traits_f.compare_distance_2_object();

    return Compare_distance_2(pe, pf);
  }

  Side_of_oriented_circle_2 side_of_oriented_circle_2_object() const
  {
    typename Exact_traits::Side_of_oriented_circle_2 pe = Delaunay_traits_e.side_of_oriented_circle_2_object();
    typename Filtering_traits::Side_of_oriented_circle_2 pf = Delaunay_traits_f.side_of_oriented_circle_2_object();

    return Side_of_oriented_circle_2(pe, pf);
  }

  // The following are inherited since they are constructions :
  // Construct_circumcenter_2

protected:
  Exact_traits Delaunay_traits_e;
  Filtering_traits Delaunay_traits_f;
};

template <class K,
          class O = CGAL::Periodic_2_offset_2,
          class D = typename K::Iso_rectangle_2,
          bool Has_static_filters_ = internal::Has_static_filters<K>::value>
class Periodic_2_Delaunay_triangulation_filtered_traits_2;

} // namespace CGAL

#include <CGAL/Periodic_2_triangulation_2/internal/Periodic_2_Delaunay_triangulation_statically_filtered_traits_2.h>

namespace CGAL {

template <class K, class O, class D>
class Periodic_2_Delaunay_triangulation_filtered_traits_2<K, O, D, false> // no static filtering
  : public Periodic_2_Delaunay_triangulation_filtered_traits_base_2<K, O, D>
{
  typedef Periodic_2_Delaunay_triangulation_filtered_traits_base_2<K, O, D>              Base;

public:
  typedef K                                                                              Kernel;
  typedef D                                                                              Domain;

  Periodic_2_Delaunay_triangulation_filtered_traits_2(const Domain& domain = P2T2::internal::get_default_domain<Kernel, Domain>(),
                                                      const Kernel& k = Kernel())
    : Base(domain, k)
  { }
};

template <class K, class O, class D>
class Periodic_2_Delaunay_triangulation_filtered_traits_2<K, O, D, true>
  : public Periodic_2_Delaunay_triangulation_statically_filtered_traits_2<K, O, D>
{
  typedef Periodic_2_Delaunay_triangulation_statically_filtered_traits_2<K, O, D>        Base;

public:
  typedef K                                                                              Kernel;
  typedef D                                                                              Domain;

  Periodic_2_Delaunay_triangulation_filtered_traits_2(const Domain& domain = P2T2::internal::get_default_domain<Kernel, Domain>(),
                                                      const Kernel& k = Kernel())
    : Base(domain, k)
  { }
};

} // namespace CGAL

#endif // CGAL_PERIODIC_2_DELAUNAY_TRIANGULATION_FILTERED_TRAITS_2_H
