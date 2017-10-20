// Copyright (c) 2017 INRIA Sophia-Antipolis (France).
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
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_PERIODIC_2_DELAUNAY_TRIANGULATION_FILTERED_TRAITS_2_H
#define CGAL_PERIODIC_2_DELAUNAY_TRIANGULATION_FILTERED_TRAITS_2_H

#include <CGAL/license/Periodic_2_triangulation_2.h>

#include <CGAL/internal/Periodic_2_triangulation_filtered_traits_2.h>
#include <CGAL/Periodic_2_Delaunay_triangulation_traits_2.h>

#include <CGAL/basic.h>
#include <CGAL/config.h>
#include <CGAL/internal/Has_boolean_tags.h>
#include <CGAL/Interval_nt.h>
#include <CGAL/Uncertain.h>
#include <CGAL/Profile_counter.h>

namespace CGAL {

// The first template item is supposed to be a Filtered_kernel-like kernel.
template < typename K, typename Off >
class Periodic_2_Delaunay_triangulation_filtered_traits_base_2
  : public Periodic_2_Delaunay_triangulation_traits_base_2<K, Off>
{
  typedef Periodic_2_Delaunay_triangulation_traits_base_2<K, Off> Base;

  // Exact traits is based on the exact kernel.
  typedef Periodic_2_Delaunay_triangulation_traits_2<typename K::Exact_kernel,
                                                     Off>      Exact_traits;
  // Filtering traits is based on the filtering kernel.
  typedef Periodic_2_Delaunay_triangulation_traits_2<typename K::Approximate_kernel,
                                                     Off>      Filtering_traits;
private:
  typedef typename K::C2E C2E;
  typedef typename K::C2F C2F;

  typedef typename K::Iso_rectangle_2 Iso_rectangle_2;

public:
  virtual ~Periodic_2_Delaunay_triangulation_filtered_traits_base_2() { }

  Periodic_2_Delaunay_triangulation_filtered_traits_base_2(const Iso_rectangle_2& domain,
                                                           const K& k)
    :
      Base(domain, k),
      Delaunay_traits_e(C2E()(domain)),
      Delaunay_traits_f(C2F()(domain))
  {
    // Problem 1: above is a default initialization of the kernel in the traits.
    // Hence, if the kernel has members and we use filtered traits, then
    // the members will be default constructed here...

    // Problem 2: we have built filtered traits in P2Tfiltered_traits_2 and now
    // we also need those two...
  }

  virtual void set_domain(const Iso_rectangle_2& domain)
  {
    this->_domain = domain;
    this->set_filtrating_traits(domain);
    set_filtrating_Delaunay_traits(domain);
  }

  void set_filtrating_Delaunay_traits(const Iso_rectangle_2& domain)
  {
    Delaunay_traits_e.set_domain(C2E()(domain));
    Delaunay_traits_f.set_domain(C2F()(domain));
  }

public:
  typedef Filtered_predicate<
            typename Exact_traits::Compare_distance_2,
            typename Filtering_traits::Compare_distance_2,
            Offset_converter_2<C2E>,
            Offset_converter_2<C2F> >  Compare_distance_2;
  typedef Filtered_predicate<
            typename Exact_traits::Side_of_oriented_circle_2,
            typename Filtering_traits::Side_of_oriented_circle_2,
            Offset_converter_2<C2E>,
            Offset_converter_2<C2F> >  Side_of_oriented_circle_2;

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

template < typename K,
           typename Off = CGAL::Periodic_2_offset_2,
           bool Has_static_filters = internal::Has_static_filters<K>::value >
class Periodic_2_Delaunay_triangulation_filtered_traits_2;

} // namespace CGAL

#include <CGAL/internal/Periodic_2_Delaunay_triangulation_statically_filtered_traits_2.h>

namespace CGAL {

template<class K, class Off>
class Periodic_2_Delaunay_triangulation_filtered_traits_2<K, Off, false>
  : public Periodic_2_Delaunay_triangulation_filtered_traits_base_2<K, Off>
{
  typedef Periodic_2_Delaunay_triangulation_filtered_traits_base_2<K, Off> Base;

public:
  typedef typename K::Iso_rectangle_2 Iso_rectangle_2;

  Periodic_2_Delaunay_triangulation_filtered_traits_2(const Iso_rectangle_2& domain,
                                                      const K& k)
    : Base(domain, k)
  { }
};

template<class K, class Off>
class Periodic_2_Delaunay_triangulation_filtered_traits_2<K, Off, true>
  : public Periodic_2_Delaunay_triangulation_statically_filtered_traits_2<K, Off>
{
  typedef Periodic_2_Delaunay_triangulation_statically_filtered_traits_2<K, Off > Base;

public:
  typedef typename K::Iso_rectangle_2 Iso_rectangle_2;

  Periodic_2_Delaunay_triangulation_filtered_traits_2(const Iso_rectangle_2& domain,
                                                      const K& k)
    : Base(domain, k)
  { }
};

} //namespace CGAL

#endif // CGAL_PERIODIC_2_DELAUNAY_TRIANGULATION_FILTERED_TRAITS_2_H
