// Copyright (c) 1997-2013 INRIA Sophia-Antipolis (France).
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
// $URL$
// $Id$
//
// Author(s)     : Nico Kruithof <Nico@nghk.nl>

#include <CGAL/Periodic_2_triangulation_traits_2.h>

#ifndef CGAL_PERIODIC_2_TRIANGULATION_FILTERED_TRAITS_2_H
#define CGAL_PERIODIC_2_TRIANGULATION_FILTERED_TRAITS_2_H

#include <string>
#include <CGAL/basic.h>
#include <CGAL/config.h>
#include <CGAL/Filtered_predicate.h>
#include <CGAL/Interval_nt.h>
#include <CGAL/Uncertain.h>
#include <CGAL/Profile_counter.h>
#include <CGAL/Periodic_2_triangulation_traits_2.h>

namespace CGAL
{
// This template class is a wrapper that implements the filtering for any
// predicate (dynamic filters with IA).

// TODO :
// - each predicate in the default kernel should define a tag that says if it
//   wants to be filtered or not (=> all homogeneous predicate define this
//   tag).  We could even test-suite that automatically.  It makes a strong
//   new requirement on the kernel though...
//   Could be done with a traits mechanism ?
//   A default template could use the current IA, but other tags or whatever
//   could specify no filtering at all, or static filtering...
// - same thing for constructions => virtual operator() ?
// - similarly, constructions should have a tag saying if they can throw or
//   not, or we let all this up to the compiler optimizer to figure out ?
// - Some caching could be done at the Point_2 level.


template <class EP, class AP, class C2E, class C2A, bool Protection = true>
class Filtered_periodic_predicate_2
  : public Filtered_predicate <EP, AP, C2E, C2A, Protection>
{
  typedef Filtered_predicate<EP, AP, C2E, C2A, Protection> Base;
public:
  Filtered_periodic_predicate_2() : Base() {}

  // These constructors are used for constructive predicates.
  // You should try to avoid constructive predicates, as they will construct
  // the exact values systematically (in the ctor), rather than lazily.
  template <class OE, class OA>
  Filtered_periodic_predicate_2(const OE * oe, const OA * oa) : Base(EP(oe), AP(oa)) {}


};

// The Offset_converter is parametrized by a usual kernel converter,
// and adds the conversions for Offsets.
template < typename Converter >
struct Offset_converter_2
  : public Converter
{
  typedef typename Converter::Source_kernel Source_kernel;
  typedef typename Converter::Target_kernel Target_kernel;

  typedef typename Periodic_2_triangulation_traits_base_2<Source_kernel>
  ::Offset_2 Source_off;
  /*   typedef typename Periodic_2_triangulation_traits_base_2<Source_kernel> */
  /*                    ::Point_2  Source_pt; */

  typedef typename Periodic_2_triangulation_traits_base_2<Target_kernel>
  ::Offset_2 Target_off;
  /*   typedef typename Periodic_2_triangulation_traits_base_2<Target_kernel> */
  /*                    ::Point_2  Target_pt; */


  using Converter::operator();

  Target_off
  operator()(const Source_off &off) const
  {
    return off;
  }
};

// The argument is supposed to be a Filtered_kernel like kernel.
template < typename K, typename Off >
class Periodic_2_triangulation_filtered_traits_base_2
  : public Periodic_2_triangulation_traits_base_2<K, Off>
{
  typedef Periodic_2_triangulation_traits_base_2<K, Off> Base;

  // Exact traits is based on the exact kernel.
  typedef Periodic_2_triangulation_traits_2<typename K::Exact_kernel, Off>
  Exact_traits;
  // Filtering traits is based on the filtering kernel.
  typedef Periodic_2_triangulation_traits_2<typename K::Approximate_kernel, Off>
  Filtering_traits;
private:
  typedef typename K::C2E C2E;
  typedef typename K::C2F C2F;

  typedef typename C2E::Target_kernel::Iso_rectangle_2 Exact_iso_rectangle_2;
  typedef typename C2F::Target_kernel::Iso_rectangle_2 Approximate_iso_rectangle_2;

public:
  typedef typename K::Iso_rectangle_2 Iso_rectangle_2;

  Periodic_2_triangulation_filtered_traits_base_2() {}


  void set_domain(const Iso_rectangle_2& domain)
  {
    C2E c2e;
    C2F c2f;
    this->_domain = domain;
    this->_domain_e = c2e(this->_domain);
    this->_domain_f = c2f(this->_domain);
  }

  typedef Filtered_periodic_predicate_2 <
  typename Exact_traits::Less_x_2,
           typename Filtering_traits::Less_x_2,
           Offset_converter_2<C2E>,
           Offset_converter_2<C2F> >  Less_x_2;
  typedef Filtered_periodic_predicate_2 <
  typename Exact_traits::Less_y_2,
           typename Filtering_traits::Less_y_2,
           Offset_converter_2<C2E>,
           Offset_converter_2<C2F> >  Less_y_2;
  typedef Filtered_periodic_predicate_2 <
  typename Exact_traits::Orientation_2,
           typename Filtering_traits::Orientation_2,
           Offset_converter_2<C2E>,
           Offset_converter_2<C2F> >  Orientation_2;
  typedef Filtered_periodic_predicate_2 <
  typename Exact_traits::Side_of_oriented_circle_2,
           typename Filtering_traits::Side_of_oriented_circle_2,
           Offset_converter_2<C2E>,
           Offset_converter_2<C2F> >  Side_of_oriented_circle_2;
  typedef Filtered_periodic_predicate_2 <
  typename Exact_traits::Compare_distance_2,
           typename Filtering_traits::Compare_distance_2,
           Offset_converter_2<C2E>,
           Offset_converter_2<C2F> >  Compare_distance_2;

  Less_x_2 less_x_2_object() const
  {
    return Less_x_2(&_domain_e, &_domain_f);
  }
  Less_y_2 less_y_2_object() const
  {
    return Less_y_2(&_domain_e, &_domain_f);
  }
  Orientation_2 orientation_2_object() const
  {
    return Orientation_2(&_domain_e, &_domain_f);
  }
  Side_of_oriented_circle_2 side_of_oriented_circle_2_object() const
  {
    return Side_of_oriented_circle_2(&_domain_e, &_domain_f);
  }
  Compare_distance_2 compare_distance_2_object() const
  {
    return Compare_distance_2(&_domain_e, &_domain_f);
  }

  // The following are inherited since they are constructions :
  // Compute_squared_radius_2 !!!
  // Construct_center_2
  // Construct_circumcenter_2
  // Construct_bisector_2
  // Construct_point_2
  // Construct_segment_2
  // Construct_triangle_2
  // Construct_direction_2
  // Construct_ray_2

protected:
  Exact_iso_rectangle_2 _domain_e;
  Approximate_iso_rectangle_2 _domain_f;
};

} //namespace CGAL

#include <CGAL/Periodic_2_triangulation_statically_filtered_traits_2.h>

namespace CGAL
{

template < typename K, typename Off = typename CGAL::Periodic_2_offset_2, bool Has_static_filters = K::Has_static_filters >
class Periodic_2_triangulation_filtered_traits_2;

template < typename K, typename Off, bool Has_static_filters >
class Periodic_2_triangulation_filtered_traits_2
  : public Periodic_2_triangulation_filtered_traits_base_2<K, Off> {};

template < typename K, typename Off >
class Periodic_2_triangulation_filtered_traits_2<K, Off, true>
  : public Periodic_2_triangulation_statically_filtered_traits_2 <
  Periodic_2_triangulation_filtered_traits_base_2<K, Off> > {};

} //namespace CGAL

#endif // CGAL_PERIODIC_2_TRIANGULATION_FILTERED_TRAITS_2_H
