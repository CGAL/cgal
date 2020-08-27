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

#include <CGAL/Periodic_2_triangulation_2/internal/Functor_with_offset_points_adaptor_2.h>
#include <CGAL/Periodic_2_triangulation_2/internal/traits_helpers.h>
#include <CGAL/Periodic_2_offset_2.h>

#include <CGAL/internal/Has_boolean_tags.h>

#include <boost/mpl/and.hpp>
#include <boost/mpl/not.hpp>

#include <type_traits>

namespace CGAL {

// default template parameterss are defined in traits_helper.h's forward declaration
template <class Kernel_, class Offset_, class Domain_>
class Periodic_2_triangulation_traits_base_2
  : public Kernel_ // @todo should just be a member of the class
{
public:
  typedef Kernel_                                                           Kernel;
  typedef Offset_                                                           Offset;
  typedef Offset                                                            Periodic_2_offset_2;

  typedef Domain_                                                           Domain;

  typedef typename Kernel::FT                                               FT;
  typedef typename Kernel::RT                                               RT;
  typedef typename Kernel::Point_2                                          Point_2;
  typedef typename Kernel::Segment_2                                        Segment_2;
  typedef typename Kernel::Triangle_2                                       Triangle_2;
  typedef typename Kernel::Iso_rectangle_2                                  Iso_rectangle_2;

private:
  typedef Periodic_2_triangulation_traits_base_2<Kernel, Offset, Domain>    Self;
  typedef Kernel_                                                           Base;

public:
  typedef typename P2T2::internal::Construct_point_getter<
                                     Kernel, Offset, Domain>::type          Construct_point_2;

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
  { std::cout << "create 'Periodic_2_triangulation_traits_base_2'" << std::endl; }

  Periodic_2_triangulation_traits_base_2(const Periodic_2_triangulation_traits_base_2& other)
    : Base(static_cast<const Base&>(other)), domain(other.get_domain())
  { }

  // Access
  void set_domain(const Domain& d) { domain = d; }
  const Domain& get_domain() const { return domain; }

  // Operations
  Construct_point_2 construct_point_2_object() const
  {
    return Construct_point_2(&domain, this->Kernel::construct_point_2_object());
  }

    // Predicates
  Less_x_2 less_x_2_object() const
  {
    return Less_x_2(this->Kernel::less_x_2_object(), construct_point_2_object());
  }
  Less_y_2 less_y_2_object() const
  {
    return Less_y_2(this->Kernel::less_y_2_object(), construct_point_2_object());
  }
  Compare_x_2 compare_x_2_object() const
  {
    return Compare_x_2(this->Kernel::compare_x_2_object(), construct_point_2_object());
  }
  Compare_y_2 compare_y_2_object() const
  {
    return Compare_y_2(this->Kernel::compare_y_2_object(), construct_point_2_object());
  }
  Orientation_2 orientation_2_object() const
  {
    return Orientation_2(this->Kernel::orientation_2_object(), construct_point_2_object());
  }

  // Constructions
  Construct_segment_2 construct_segment_2_object() const
  {
    return Construct_segment_2(this->Kernel::construct_segment_2_object(), construct_point_2_object());
  }
  Construct_triangle_2 construct_triangle_2_object() const
  {
    return Construct_triangle_2(this->Kernel::construct_triangle_2_object(), construct_point_2_object());
  }

protected:
  Domain domain;
};

// Forward declaration for the filtered traits

// @exact disable filtering for Lattice-based Periodic
template <class K,
          class O = CGAL::Periodic_2_offset_2,
          class D = typename K::Iso_rectangle_2,
          bool Has_filtered_predicates =
            boost::mpl::and_<internal::Has_filtered_predicates<K>,
                             boost::mpl::not_<std::is_same<D, Lattice_2<K> > > >::value>
class Periodic_2_triangulation_traits_2;

} // namespace CGAL

// Partial specialization for kernels with filtered predicates
#include <CGAL/Periodic_2_triangulation_2/internal/Periodic_2_triangulation_filtered_traits_2.h>

namespace CGAL {

template <class K, class O, class D>
class Periodic_2_triangulation_traits_2<K, O, D, false> // no filtering
  : public Periodic_2_triangulation_traits_base_2<K, O, D>
{
  typedef Periodic_2_triangulation_traits_base_2<K, O, D>     Base;

public:
  Periodic_2_triangulation_traits_2(const D& domain = P2T2::internal::get_default_domain<K, D>(),
                                    const K& k = K())
    : Base(domain, k)
  { }
};

// Below is filtered, but might even be statically filtered, depending on 'Has_static_filters'
template <class K, class O, class D>
class Periodic_2_triangulation_traits_2<K, O, D, true>
  : public Periodic_2_triangulation_filtered_traits_2<K, O, D, internal::Has_static_filters<K>::value>
{
  typedef Periodic_2_triangulation_filtered_traits_2<
            K, O, D, internal::Has_static_filters<K>::value>  Base;

public:
  Periodic_2_triangulation_traits_2(const D& domain = P2T2::internal::get_default_domain<K, D>(),
                                    const K& k = K())
    : Base(domain, k)
  { }
};

} // namespace CGAL

#endif // CGAL_PERIODIC_2_TRIANGULATION_TRAITS_2_H
