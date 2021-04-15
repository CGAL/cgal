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

template <class K_,
          class Off_ = typename CGAL::Periodic_2_offset_2 >
class Periodic_2_triangulation_traits_base_2
  : public K_
{
  typedef Periodic_2_triangulation_traits_base_2<K_, Off_>         Self;
  typedef K_                                                       Base;

public:
  typedef K_                                    Kernel;
  typedef Off_                                  Offset;

  typedef typename Kernel::FT                   FT;
  typedef typename Kernel::RT                   RT;
  typedef typename Kernel::Point_2              Point_2;
  typedef typename Kernel::Segment_2            Segment_2;
  typedef typename Kernel::Triangle_2           Triangle_2;
  typedef typename Kernel::Iso_rectangle_2      Iso_rectangle_2;

  typedef Offset                                Periodic_2_offset_2;

  typedef Periodic_2_construct_point_2<Self, typename Kernel::Construct_point_2>
      Construct_point_2;

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
  virtual ~Periodic_2_triangulation_traits_base_2() { }

  Periodic_2_triangulation_traits_base_2(const Iso_rectangle_2& domain,
                                         const Kernel& k)
    : Base(k), _domain(domain)
  { }

  // Access
  virtual void set_domain(const Iso_rectangle_2& domain) { _domain = domain; }
  Iso_rectangle_2 get_domain() const { return _domain; }

  // Operations
  Construct_point_2 construct_point_2_object() const {
    return Construct_point_2(&_domain, this->Kernel::construct_point_2_object());
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
  Iso_rectangle_2 _domain;
};


// Forward declaration for the filtered traits
template <class K_,
          class Off_ = CGAL::Periodic_2_offset_2,
          bool Has_filtered_predicates = internal::Has_filtered_predicates<K_>::value >
class Periodic_2_triangulation_traits_2;

} // namespace CGAL

// Partial specialization for Filtered_kernel<CK>.
#include <CGAL/internal/Periodic_2_triangulation_filtered_traits_2.h>

namespace CGAL
{

template <class K_, class Off_>
class Periodic_2_triangulation_traits_2<K_, Off_, false>
  : public Periodic_2_triangulation_traits_base_2<K_, Off_>
{
  typedef Periodic_2_triangulation_traits_base_2<K_, Off_>        Base;

public:
  typedef K_                                                      Kernel;
  typedef typename Kernel::Iso_rectangle_2                        Iso_rectangle_2;

  Periodic_2_triangulation_traits_2(const Iso_rectangle_2& domain = Iso_rectangle_2(0,0,1,1),
                                    const Kernel& k = Kernel())
    : Base(domain, k)
  { }
};

template <class K_, class Off_>
class Periodic_2_triangulation_traits_2<K_, Off_, true>
  : public Periodic_2_triangulation_filtered_traits_2<
             K_, Off_, internal::Has_static_filters<K_>::value>
{
  typedef Periodic_2_triangulation_filtered_traits_2<
            K_, Off_, internal::Has_static_filters<K_>::value>  Base;

public:
  typedef K_                                                   Kernel;
  typedef typename Kernel::Iso_rectangle_2                     Iso_rectangle_2;

  Periodic_2_triangulation_traits_2(const Iso_rectangle_2& domain = Iso_rectangle_2(0,0,1,1),
                                    const Kernel& k = Kernel())
    : Base(domain, k)
  { }
};

} //namespace CGAL

#endif // CGAL_PERIODIC_2_TRIANGULATION_TRAITS_2_H
