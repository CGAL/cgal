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

#include <CGAL/internal/Functor_with_offset_points_adaptor_2.h>

#include <CGAL/Periodic_2_offset_2.h>
#include <CGAL/Periodic_2_triangulation_traits_2.h>

#include <CGAL/internal/Has_boolean_tags.h>

namespace CGAL {

template <class Kernel_,
          class Offset_ = CGAL::Periodic_2_offset_2,
          class Domain_ = typename Kernel_::Iso_rectangle_2,
          class Construct_point_2_ = Default>
class Periodic_2_Delaunay_triangulation_traits_base_2
  : public Periodic_2_triangulation_traits_2<
             Kernel_, Offset_, Domain_,
             typename Default::Get<Construct_point_2_,
                                   Periodic_2_construct_point_2<Kernel_, Offset_> >::type>
{
public:
  typedef Kernel_                                                 Kernel;
  typedef Offset_                                                 Offset;
  typedef Domain_                                                 Domain;

  typedef Periodic_2_construct_point_2<Kernel, Offset>            Construct_point_2_def;
  typedef typename Default::Get<Construct_point_2_,
                                Construct_point_2_def>::type      Construct_point_2;

private:
  typedef Periodic_2_Delaunay_triangulation_traits_base_2<
            Kernel, Offset, Domain, Construct_point_2>            Self;
  typedef Periodic_2_triangulation_traits_2<
            Kernel, Offset, Domain, Construct_point_2>            Base;

public:
  typedef typename Base::RT                                       RT;
  typedef typename Base::FT                                       FT;
  typedef typename Base::Point_2                                  Point_2;
  typedef typename Base::Periodic_2_offset_2                      Periodic_2_offset_2;
  typedef typename Base::Iso_rectangle_2                          Iso_rectangle_2;

public:
  Periodic_2_Delaunay_triangulation_traits_base_2(const Domain& domain = Domain(),
                                                  const Kernel& k = Kernel())
    : Base(domain, k)
  { }

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
  Compare_distance_2 compare_distance_2_object() const {
    return Compare_distance_2(this->Base::compare_distance_2_object(),
                              this->construct_point_2_object());
  }
  Side_of_oriented_circle_2 side_of_oriented_circle_2_object() const {
    return Side_of_oriented_circle_2(this->Base::side_of_oriented_circle_2_object(),
                                     this->construct_point_2_object());
  }
  Compare_squared_radius_2 compare_squared_radius_2_object() const {
     return Compare_squared_radius_2(this->Base::compare_squared_radius_2_object(),
                                     this->construct_point_2_object());
  }

  // Constructions
  Construct_circumcenter_2 construct_circumcenter_2_object() const {
    return Construct_circumcenter_2(this->Base::construct_circumcenter_2_object(),
                                    this->construct_point_2_object());
  }
};

template <class K,
          class O = CGAL::Periodic_2_offset_2,
          class D = typename K::Iso_rectangle_2,
          class CP = Default,
          bool Has_filtered_predicates_ = internal::Has_filtered_predicates<K>::value >
class Periodic_2_Delaunay_triangulation_traits_2
   : public Periodic_2_Delaunay_triangulation_traits_base_2<K, O, D, CP> // @tmp this is just a forward declaration normally
{
public:
  typedef Periodic_2_Delaunay_triangulation_traits_base_2<K, O, D, CP> Base;
  Periodic_2_Delaunay_triangulation_traits_2(const D& d = D(), const K& k = K()) : Base(d, k) { }
};
#if 0

} //namespace CGAL

// Partial specialization for kernels with filtered predicates
#include <CGAL/internal/Periodic_2_Delaunay_triangulation_filtered_traits_2.h>

namespace CGAL {

template <class K, class O, class D, class CP>
class Periodic_2_Delaunay_triangulation_traits_2<K, O, D, CP, false>
  : public Periodic_2_Delaunay_triangulation_traits_base_2<K, O, D, CP>
{
  typedef Periodic_2_Delaunay_triangulation_traits_base_2<K, O, D, CP> Base;

public:
  typedef K_                                                           Kernel;
  typedef typename Kernel::Domain                                      Domain;

  Periodic_2_Delaunay_triangulation_traits_2(const Domain& domain, const Kernel& k = Kernel())
    : Base(domain, k)
  { }
};

template <class K, class O, class D, class CP>
class Periodic_2_Delaunay_triangulation_traits_2<K, O, D, CP, true>
    : public Periodic_2_Delaunay_triangulation_filtered_traits_2<
               K, O, D, CP, internal::Has_static_filters<K>::value>
{
  typedef Periodic_2_Delaunay_triangulation_filtered_traits_2<
            K, O, D, CP, internal::Has_static_filters<K>::value>    Base;

public:
  typedef K                                                         Kernel;
  typedef typename Kernel::Domain                                   Domain;

  Periodic_2_Delaunay_triangulation_traits_2(const Domain& domain, const Kernel& k = Kernel())
    : Base(domain, k)
  { }
};

#endif //

} // namespace CGAL

#endif // CGAL_PERIODIC_2_DELAUNAY_TRIANGULATION_TRAITS_2_H
