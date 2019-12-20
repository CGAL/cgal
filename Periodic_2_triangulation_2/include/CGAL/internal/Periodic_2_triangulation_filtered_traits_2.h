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

#include <CGAL/basic.h>
#include <CGAL/config.h>
#include <CGAL/Filtered_predicate.h>
#include <CGAL/Interval_nt.h>
#include <CGAL/Uncertain.h>
#include <CGAL/Profile_counter.h>

#include <CGAL/Periodic_2_triangulation_traits_2.h>

namespace CGAL
{
// The Offset_converter is parametrized by a usual kernel converter,
// and adds the conversions for Offsets.
template <typename Converter_>
struct Offset_converter_2
  : public Converter_
{
  typedef Converter_                        Base;

  typedef typename Base::Source_kernel      Source_kernel;
  typedef typename Base::Target_kernel      Target_kernel;

  typedef typename Periodic_2_triangulation_traits_base_2<Source_kernel>::Offset   Source_off;
  typedef typename Periodic_2_triangulation_traits_base_2<Source_kernel>::Point_2  Source_pt;

  typedef typename Periodic_2_triangulation_traits_base_2<Target_kernel>::Offset   Target_off;
  typedef typename Periodic_2_triangulation_traits_base_2<Target_kernel>::Point_2  Target_pt;

  using Base::operator();

  Target_off operator()(const Source_off &off) const { return off; }
};

// The argument is supposed to be a Filtered_kernel like kernel.
template <class K_, typename Off_>
class Periodic_2_triangulation_filtered_traits_base_2
  : public Periodic_2_triangulation_traits_base_2<K_, Off_>
{
  typedef Periodic_2_triangulation_traits_base_2<K_, Off_> Base;
  typedef K_                                               Kernel;

  typedef typename Kernel::Exact_kernel                    EKernel;
  typedef typename Kernel::Approximate_kernel              AKernel;
  typedef typename Kernel::C2E                             C2E;
  typedef typename Kernel::C2F                             C2F;

  // Exact traits is based on the exact kernel.
  typedef Periodic_2_triangulation_traits_2<EKernel, Off_> Exact_traits;
  // Filtering traits is based on the filtering kernel.
  typedef Periodic_2_triangulation_traits_2<AKernel, Off_> Filtering_traits;

public:
  typedef typename Kernel::Iso_rectangle_2 Iso_rectangle_2;

  virtual ~Periodic_2_triangulation_filtered_traits_base_2() { }

  Periodic_2_triangulation_filtered_traits_base_2(const Iso_rectangle_2& domain,
                                                  const Kernel& k)
    : Base(domain, k),
      traits_e(C2E()(domain)),
      traits_f(C2F()(domain))
  {
    // Problem: below is a default initialization of the kernel in the traits.
    // Hence, if the kernel has members and we use filtered traits, then
    // the members will be default constructed here...
  }

  virtual void set_domain(const Iso_rectangle_2& domain)
  {
    this->_domain = domain;
    set_filtrating_traits(domain);
  }

  void set_filtrating_traits(const Iso_rectangle_2& domain)
  {
    traits_e.set_domain(C2E()(domain));
    traits_f.set_domain(C2F()(domain));
  }

  typedef Filtered_predicate<typename Exact_traits::Less_x_2,
                             typename Filtering_traits::Less_x_2,
                             Offset_converter_2<C2E>,
                             Offset_converter_2<C2F> >  Less_x_2;
  typedef Filtered_predicate<typename Exact_traits::Less_y_2,
                             typename Filtering_traits::Less_y_2,
                             Offset_converter_2<C2E>,
                             Offset_converter_2<C2F> >  Less_y_2;
  typedef Filtered_predicate<typename Exact_traits::Compare_x_2,
                             typename Filtering_traits::Compare_x_2,
                             Offset_converter_2<C2E>,
                             Offset_converter_2<C2F> >  Compare_x_2;
  typedef Filtered_predicate<typename Exact_traits::Compare_y_2,
                             typename Filtering_traits::Compare_y_2,
                             Offset_converter_2<C2E>,
                             Offset_converter_2<C2F> >  Compare_y_2;
  typedef Filtered_predicate<typename Exact_traits::Orientation_2,
                             typename Filtering_traits::Orientation_2,
                             Offset_converter_2<C2E>,
                             Offset_converter_2<C2F> >  Orientation_2;

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

template <class K_,
          class Off_ = typename CGAL::Periodic_2_offset_2,
          bool Has_static_filters_ = internal::Has_static_filters<K_>::value >
class Periodic_2_triangulation_filtered_traits_2;

} //namespace CGAL

#include <CGAL/internal/Periodic_2_triangulation_statically_filtered_traits_2.h>

namespace CGAL {

template <class K_, class Off_>
class Periodic_2_triangulation_filtered_traits_2<K_, Off_, false>
  : public Periodic_2_triangulation_filtered_traits_base_2<K_, Off_>
{
  typedef Periodic_2_triangulation_filtered_traits_base_2<K_, Off_> Base;

public:
  typedef K_                                                        Kernel;
  typedef typename Kernel::Iso_rectangle_2                          Iso_rectangle_2;

  Periodic_2_triangulation_filtered_traits_2(const Iso_rectangle_2& domain,
                                             const Kernel& k)
    : Base(domain, k)
  { }
};

template <class K_, class Off_>
class Periodic_2_triangulation_filtered_traits_2<K_, Off_, true>
  : public Periodic_2_triangulation_statically_filtered_traits_2<K_, Off_>
{
  typedef Periodic_2_triangulation_statically_filtered_traits_2<K_, Off_> Base;

public:
  typedef K_                                                              Kernel;
  typedef typename Kernel::Iso_rectangle_2                                Iso_rectangle_2;

  Periodic_2_triangulation_filtered_traits_2(const Iso_rectangle_2& domain,
                                             const Kernel& k)
    : Base(domain, k)
  { }
};

} //namespace CGAL

#endif // CGAL_PERIODIC_2_TRIANGULATION_FILTERED_TRAITS_2_H
