// Copyright (c) 1997-2013 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Nico Kruithof <Nico@nghk.nl>,
//                 Mael Rouxel-Labbé

#ifndef CGAL_PERIODIC_2_TRIANGULATION_FILTERED_TRAITS_2_H
#define CGAL_PERIODIC_2_TRIANGULATION_FILTERED_TRAITS_2_H

#include <CGAL/license/Periodic_2_triangulation_2.h>

#include <CGAL/Periodic_2_triangulation_2/internal/traits_helpers.h>
#include <CGAL/Periodic_2_triangulation_traits_2.h>
#include <CGAL/Periodic_2_offset_2.h>

#include <CGAL/Filtered_predicate.h>

namespace CGAL {

// The argument is supposed to be a Filtered_kernel like kernel.
template <class Kernel_,
          class Offset_ = CGAL::Periodic_2_offset_2,
          class Domain_ = typename Kernel_::Iso_rectangle_2>
class Periodic_2_triangulation_filtered_traits_base_2
  : public Periodic_2_triangulation_traits_base_2<Kernel_, Offset_, Domain_>
{
public:
  typedef Kernel_                                                           Kernel;
  typedef Offset_                                                           Offset;
  typedef Domain_                                                           Domain;

private:
  typedef Periodic_2_triangulation_traits_base_2<Kernel, Offset, Domain>    Base;

  typedef typename Kernel::Exact_kernel                                     EKernel;
  typedef typename Kernel::Approximate_kernel                               AKernel;
  typedef typename Kernel::C2E                                              C2E;
  typedef typename Kernel::C2F                                              C2F;

  typedef typename P2T2::internal::Exact_domain_getter<Kernel, Domain, EKernel>::type
                                                                            EDomain;
  typedef typename P2T2::internal::Exact_domain_getter<Kernel, Domain, AKernel>::type
                                                                            ADomain;

public:
  // Exact traits is based on the exact kernel.
  typedef Periodic_2_triangulation_traits_2<EKernel, Offset, EDomain>       Exact_traits;

  // Filtering traits is based on the filtering kernel.
  typedef Periodic_2_triangulation_traits_2<AKernel, Offset, ADomain>       Filtering_traits;

public:
  Periodic_2_triangulation_filtered_traits_base_2(const Domain& domain,
                                                  const Kernel& k)
    : Base(domain, k),
      traits_e(P2T2::internal::convert_domain<Kernel, EKernel>(domain)),
      traits_f(P2T2::internal::convert_domain<Kernel, AKernel>(domain))
  {
    std::cout << "create 'Periodic_2_triangulation_filtered_traits_base_2'" << std::endl;
  }

  void set_filtrating_traits(const Domain& domain)
  {
    traits_e.set_domain(P2T2::internal::convert_domain<Kernel, EKernel>(domain));
    traits_f.set_domain(P2T2::internal::convert_domain<Kernel, AKernel>(domain));
  }

  void set_domain(const Domain& domain)
  {
    Base::set_domain(domain);
    set_filtrating_traits(domain);
  }

public:
  typedef typename P2T2::internal::Construct_point_getter<
                                     Kernel, Offset, Domain>::type          Construct_point_2;

  typedef Filtered_predicate<typename Exact_traits::Less_x_2,
                             typename Filtering_traits::Less_x_2,
                             P2T2::internal::Offset_converter_2<C2E>,
                             P2T2::internal::Offset_converter_2<C2F> >      Less_x_2;
  typedef Filtered_predicate<typename Exact_traits::Less_y_2,
                             typename Filtering_traits::Less_y_2,
                             P2T2::internal::Offset_converter_2<C2E>,
                             P2T2::internal::Offset_converter_2<C2F> >      Less_y_2;
  typedef Filtered_predicate<typename Exact_traits::Compare_x_2,
                             typename Filtering_traits::Compare_x_2,
                             P2T2::internal::Offset_converter_2<C2E>,
                             P2T2::internal::Offset_converter_2<C2F> >      Compare_x_2;
  typedef Filtered_predicate<typename Exact_traits::Compare_y_2,
                             typename Filtering_traits::Compare_y_2,
                             P2T2::internal::Offset_converter_2<C2E>,
                             P2T2::internal::Offset_converter_2<C2F> >       Compare_y_2;
  typedef Filtered_predicate<typename Exact_traits::Orientation_2,
                             typename Filtering_traits::Orientation_2,
                             P2T2::internal::Offset_converter_2<C2E>,
                             P2T2::internal::Offset_converter_2<C2F> >      Orientation_2;

  Less_x_2 less_x_2_object() const
  {
    typename Exact_traits::Less_x_2 pe = traits_e.less_x_2_object();
    typename Filtering_traits::Less_x_2 pf = traits_f.less_x_2_object();

    return Less_x_2(pe, pf);
  }
  Less_y_2 less_y_2_object() const
  {
    typename Exact_traits::Less_y_2 pe = traits_e.less_y_2_object();
    typename Filtering_traits::Less_y_2 pf = traits_f.less_y_2_object();

    return Less_y_2(pe, pf);
  }
  Compare_x_2 compare_x_2_object() const
  {
    typename Exact_traits::Compare_x_2 pe = traits_e.compare_x_2_object();
    typename Filtering_traits::Compare_x_2 pf = traits_f.compare_x_2_object();

    return Compare_x_2(pe, pf);
  }
  Compare_y_2 compare_y_2_object() const
  {
    typename Exact_traits::Compare_y_2 pe = traits_e.compare_y_2_object();
    typename Filtering_traits::Compare_y_2 pf = traits_f.compare_y_2_object();

    return Compare_y_2(pe, pf);
  }
  Orientation_2 orientation_2_object() const
  {
    typename Exact_traits::Orientation_2 pe = traits_e.orientation_2_object();
    typename Filtering_traits::Orientation_2 pf = traits_f.orientation_2_object();

    return Orientation_2(pe, pf);
  }

  // The following are inherited since they are constructions :
  // Construct_point_2
  // Construct_segment_2
  // Construct_triangle_2

protected:
  Exact_traits traits_e;
  Filtering_traits traits_f;
};

template <class K,
          class O = CGAL::Periodic_2_offset_2,
          class D = typename K::Iso_rectangle_2,
          bool Has_static_filters_ = internal::Has_static_filters<K>::value>
class Periodic_2_triangulation_filtered_traits_2;

} // namespace CGAL

#include <CGAL/Periodic_2_triangulation_2/internal/Periodic_2_triangulation_statically_filtered_traits_2.h>

namespace CGAL {

template <class K, class O, class D>
class Periodic_2_triangulation_filtered_traits_2<K, O, D, false> // no static filtering
  : public Periodic_2_triangulation_filtered_traits_base_2<K, O, D>
{
  typedef Periodic_2_triangulation_filtered_traits_base_2<K, O, D>        Base;

public:
  Periodic_2_triangulation_filtered_traits_2(const D& domain = P2T2::internal::get_default_domain<K, D>(),
                                             const K& k = K())
    : Base(domain, k)
  { }
};

template <class K, class O, class D>
class Periodic_2_triangulation_filtered_traits_2<K, O, D, true>
  : public Periodic_2_triangulation_statically_filtered_traits_2<K, O, D>
{
  typedef Periodic_2_triangulation_statically_filtered_traits_2<K, O, D>  Base;

public:
  Periodic_2_triangulation_filtered_traits_2(const D& domain = P2T2::internal::get_default_domain<K, D>(),
                                             const K& k = K())
    : Base(domain, k)
  { }
};

} // namespace CGAL

#endif // CGAL_PERIODIC_2_TRIANGULATION_FILTERED_TRAITS_2_H
