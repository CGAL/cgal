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

#ifndef CGAL_PERIODIC_2_TRIANGULATION_TRAITS_2_H
#define CGAL_PERIODIC_2_TRIANGULATION_TRAITS_2_H

#include <CGAL/license/Periodic_2_triangulation_2.h>

#include <CGAL/internal/Periodic_2_construct_point_2.h>
#include <CGAL/internal/Functor_with_offset_points_adaptor_2.h>
#include <CGAL/Periodic_2_offset_2.h>

#include <CGAL/internal/Has_boolean_tags.h>
#include <CGAL/triangulation_assertions.h>

namespace CGAL {

template <class Kernel_,
          class Offset_ = CGAL::Periodic_2_offset_2,
          class Domain_ = typename Kernel_::Iso_rectangle_2,
          class Construct_point_2_ = Default>
class Periodic_2_triangulation_traits_base_2
  : public Kernel_ // @todo really shouldn't have to do that actually
{
public:
  typedef Kernel_                                                 Kernel;
  typedef Offset_                                                 Offset;
  typedef Offset                                                  Periodic_2_offset_2;

  typedef Domain_                                                 Domain;

  typedef typename Kernel::FT                                     FT;
  typedef typename Kernel::RT                                     RT;
  typedef typename Kernel::Point_2                                Point_2;
  typedef typename Kernel::Segment_2                              Segment_2;
  typedef typename Kernel::Triangle_2                             Triangle_2;
  typedef typename Kernel::Iso_rectangle_2                        Iso_rectangle_2;

  typedef Periodic_2_construct_point_2<Kernel, Offset>            Construct_point_2_def;
  typedef typename Default::Get<Construct_point_2_,
                                Construct_point_2_def>::type      Construct_point_2;

private:
  // Not truly "Self" since we grab the default value from CGAL::Default
  typedef Periodic_2_triangulation_traits_base_2<
            Kernel_, Offset_, Domain_, Construct_point_2>         Self;
  typedef Kernel_                                                 Base;

public:
  // Triangulation predicates
  typedef Functor_with_offset_points_adaptor_2<Self, typename Kernel::Less_x_2>
      Less_x_2;
  typedef Functor_with_offset_points_adaptor_2<Self, typename Kernel::Less_y_2>
      Less_y_2;
  typedef Functor_with_offset_points_adaptor_2<Self, typename Kernel::Compare_x_2>
      Compare_x_2;
  typedef Functor_with_offset_points_adaptor_2<Self, typename Kernel::Compare_y_2>
      Compare_y_2;
  typedef Functor_with_offset_points_adaptor_2<Self, typename Kernel::Orientation_2>
      Orientation_2;

  // Triangulation constructions
  typedef Functor_with_offset_points_adaptor_2<Self, typename Kernel::Construct_segment_2>
      Construct_segment_2;
  typedef Functor_with_offset_points_adaptor_2<Self, typename Kernel::Construct_triangle_2>
      Construct_triangle_2;

  // Constructor
  Periodic_2_triangulation_traits_base_2(const Domain& d = Domain(),
                                         const Kernel& k = Kernel())
    : Base(k), domain(d)
  { }

  Periodic_2_triangulation_traits_base_2(const Periodic_2_triangulation_traits_base_2& other)
    : Base(static_cast<const Base&>(other)), domain(other.get_domain())
  { }

  // Access
  void set_domain(const Domain& d) { domain = d; }
  const Domain& get_domain() const { return domain; }

  // Operations
  Construct_point_2 construct_point_2_object() const {
    return Construct_point_2(&domain, this->Kernel::construct_point_2_object());
  }

    // Predicates
  Less_x_2 less_x_2_object() const {
    return Less_x_2(this->Kernel::less_x_2_object(), construct_point_2_object());
  }
  Less_y_2 less_y_2_object() const {
    return Less_y_2(this->Kernel::less_y_2_object(), construct_point_2_object());
  }
  Compare_x_2 compare_x_2_object() const {
    return Compare_x_2(this->Kernel::compare_x_2_object(), construct_point_2_object());
  }
  Compare_y_2 compare_y_2_object() const {
    return Compare_y_2(this->Kernel::compare_y_2_object(), construct_point_2_object());
  }
  Orientation_2 orientation_2_object() const {
    return Orientation_2(this->Kernel::orientation_2_object(), construct_point_2_object());
  }

    // Constructions
  Construct_segment_2 construct_segment_2_object() const {
    return Construct_segment_2(this->Kernel::construct_segment_2_object(), construct_point_2_object());
  }
  Construct_triangle_2 construct_triangle_2_object() const {
    return Construct_triangle_2(this->Kernel::construct_triangle_2_object(), construct_point_2_object());
  }

protected:
  Domain domain;
};

// Forward declaration for the filtered traits
template <class K,
          class O = CGAL::Periodic_2_offset_2,
          class D = typename K::Iso_rectangle_2,
          class CP = Default,
          bool Has_filtered_predicates = internal::Has_filtered_predicates<K>::value >
class Periodic_2_triangulation_traits_2
  : public Periodic_2_triangulation_traits_base_2<K, O, D, CP> // @tmp (this is just a forward declaration normally)
{
public:
  typedef Periodic_2_triangulation_traits_base_2<K, O, D, CP> Base;
  Periodic_2_triangulation_traits_2(const D& d = D(), const K& k = K()) : Base(d, k) { }
};
#if 0

} // namespace CGAL

// Partial specialization for Filtered_kernel<CK>.
#include <CGAL/internal/Periodic_2_triangulation_filtered_traits_2.h>

namespace CGAL {

template <class K, class O, class D, class CP>
class Periodic_2_triangulation_traits_2<K, O, D, CP, false>
  : public Periodic_2_triangulation_traits_base_2<K, O, D, CP>
{
  typedef Periodic_2_triangulation_traits_base_2<K, O, D, CP>     Base;

public:
  typedef K                                                       Kernel;
  typedef typename Kernel::Domain                                 Domain;

  // @todo restore default construction of traits with a unit iso rectangle
  Periodic_2_triangulation_traits_2(const Domain& domain, const Kernel& k = Kernel())
    : Base(domain, k)
  { }
};

template <class K, class O, class D, class CP>
class Periodic_2_triangulation_traits_2<K, O, D, CP, true>
  : public Periodic_2_triangulation_filtered_traits_2<
             K, O, D, CP, internal::Has_static_filters<K>::value>
{
  typedef Periodic_2_triangulation_filtered_traits_2<
            K, O, D, CP, internal::Has_static_filters<K>::value>  Base;

public:
  typedef K                                                       Kernel;
  typedef typename Kernel::Domain                                 Domain;

  Periodic_2_triangulation_traits_2(const Domain& domain, const Kernel& k = Kernel())
    : Base(domain, k)
  { }
};

#endif // 

} //namespace CGAL

#endif // CGAL_PERIODIC_2_TRIANGULATION_TRAITS_2_H
