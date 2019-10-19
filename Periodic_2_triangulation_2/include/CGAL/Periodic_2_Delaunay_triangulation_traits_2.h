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

#include <CGAL/internal/Periodic_2_construct_point_2.h>
#include <CGAL/internal/Functor_with_offset_points_adaptor_2.h>

#include <CGAL/Periodic_2_offset_2.h>
#include <CGAL/Periodic_2_triangulation_traits_2.h>

#include <CGAL/internal/Has_boolean_tags.h>

namespace CGAL {

template <class K_,
          class Off_ = typename CGAL::Periodic_2_offset_2>
class Periodic_2_Delaunay_triangulation_traits_base_2
  : public Periodic_2_triangulation_traits_2<K_, Off_>
{
  typedef Periodic_2_Delaunay_triangulation_traits_base_2<K_, Off_>  Self;
  typedef Periodic_2_triangulation_traits_2<K_, Off_>                Base;

public:
  typedef K_                                  Kernel;
  typedef Off_                                Offset;

  typedef typename Base::RT                   RT;
  typedef typename Base::FT                   FT;
  typedef typename Base::Point_2              Point_2;
  typedef typename Base::Periodic_2_offset_2  Periodic_2_offset_2;
  typedef typename Base::Iso_rectangle_2      Iso_rectangle_2;

public:
  Periodic_2_Delaunay_triangulation_traits_base_2(const Iso_rectangle_2& domain,
                                                  const Kernel& k)
    : Base(domain, k)
  { }

  typedef Functor_with_offset_points_adaptor_2<Self, typename Kernel::Compare_distance_2>
      Compare_distance_2;
  typedef Functor_with_offset_points_adaptor_2<Self, typename Kernel::Side_of_oriented_circle_2>
      Side_of_oriented_circle_2;

  typedef Functor_with_offset_points_adaptor_2<Self, typename Kernel::Construct_circumcenter_2>
      Construct_circumcenter_2;

  Compare_distance_2 compare_distance_2_object() const {
    return Compare_distance_2(this->Base::compare_distance_2_object(),
                              this->construct_point_2_object());
  }
  Side_of_oriented_circle_2 side_of_oriented_circle_2_object() const {
    return Side_of_oriented_circle_2(this->Base::side_of_oriented_circle_2_object(),
                                     this->construct_point_2_object());
  }

  Construct_circumcenter_2 construct_circumcenter_2_object() const {
    return Construct_circumcenter_2(this->Base::construct_circumcenter_2_object(),
                                    this->construct_point_2_object());
  }
};

template <class K_,
          class Off_ = CGAL::Periodic_2_offset_2,
          bool Has_filtered_predicates_ = internal::Has_filtered_predicates<K_>::value >
class Periodic_2_Delaunay_triangulation_traits_2;

} //namespace CGAL

// Partial specialization for kernels with filtered predicates
#include <CGAL/internal/Periodic_2_Delaunay_triangulation_filtered_traits_2.h>

namespace CGAL {

template <class K_, class Off_>
class Periodic_2_Delaunay_triangulation_traits_2<K_, Off_, false>
  : public Periodic_2_Delaunay_triangulation_traits_base_2<K_, Off_>
{
  typedef Periodic_2_Delaunay_triangulation_traits_base_2<K_, Off_> Base;

public:
  typedef K_                                                        Kernel;
  typedef typename Kernel::Iso_rectangle_2                          Iso_rectangle_2;

  Periodic_2_Delaunay_triangulation_traits_2(const Iso_rectangle_2& domain = Iso_rectangle_2(0,0,1,1),
                                             const Kernel& k = Kernel())
    : Base(domain, k)
  { }
};

template <typename K_, typename Off_>
class Periodic_2_Delaunay_triangulation_traits_2<K_, Off_, true>
    : public Periodic_2_Delaunay_triangulation_filtered_traits_2<
               K_, Off_, internal::Has_static_filters<K_>::value>
{
  typedef Periodic_2_Delaunay_triangulation_filtered_traits_2<
            K_, Off_, internal::Has_static_filters<K_>::value>      Base;

public:
  typedef K_                                                        Kernel;
  typedef typename Kernel::Iso_rectangle_2                          Iso_rectangle_2;

  Periodic_2_Delaunay_triangulation_traits_2(const Iso_rectangle_2& domain = Iso_rectangle_2(0,0,1,1),
                                             const Kernel& k = Kernel())
    : Base(domain, k)
  { }
};

} //namespace CGAL

#endif // CGAL_PERIODIC_2_DELAUNAY_TRIANGULATION_TRAITS_2_H
