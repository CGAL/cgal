// Copyright (c) 2006,2007,2008,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s): Ron Wein          <wein@post.tau.ac.il>
//            Baruch Zukerman   <baruchzu@post.tau.ac.il>
//            Efi Fogel         <efif@post.tau.ac.il>
//
#ifndef CGAL_ARRANGEMENT_ON_SURFACE_2_GLOBAL_H
#define CGAL_ARRANGEMENT_ON_SURFACE_2_GLOBAL_H

#include <CGAL/license/Arrangement_on_surface_2.h>

/*! \file
 * Global insertion functions for the Arrangement_2 class.
 */

#include <boost/type_traits.hpp>
#include <boost/mpl/if.hpp>
#include <boost/type_traits.hpp>
#include <list>

#include <CGAL/Arrangement_on_surface_2.h>
#include <CGAL/Arr_accessor.h>
#include <CGAL/Arrangement_zone_2.h>
#include <CGAL/Arrangement_2/Arr_compute_zone_visitor.h>
#include <CGAL/Arrangement_2/Arr_do_intersect_zone_visitor.h>
#include <CGAL/Arrangement_2/Arr_traits_adaptor_2.h>
#include <CGAL/Arr_point_location_result.h>
#include <CGAL/No_intersection_surface_sweep_2.h>
#include <CGAL/Surface_sweep_2/Arr_insertion_ss_visitor.h>
#include <CGAL/Surface_sweep_2/Arr_no_intersection_insertion_ss_visitor.h>
#include <CGAL/Surface_sweep_2.h>
#include <CGAL/Surface_sweep_2/Surface_sweep_2_utils.h>
#include <CGAL/Surface_sweep_2/Do_interior_intersect_visitor.h>
#include <CGAL/Surface_sweep_2/Arr_construction_ss_visitor.h>
#include <CGAL/Surface_sweep_2/Arr_insertion_ss_visitor.h>
#include <CGAL/Surface_sweep_2/Arr_construction_event.h>
#include <CGAL/Surface_sweep_2/Arr_construction_subcurve.h>
#include <CGAL/Surface_sweep_2/Arr_insertion_traits_2.h>
#include <CGAL/Surface_sweep_2/No_overlap_event_base.h>
#include <CGAL/Surface_sweep_2/No_overlap_subcurve.h>
#include <CGAL/assertions.h>
#include <CGAL/use.h>
#include <CGAL/IO/Arr_iostream.h>

namespace CGAL {

namespace Ss2 = Surface_sweep_2;

//-----------------------------------------------------------------------------
// Insert a curve into the arrangement (incremental insertion).
// The inserted curve may intersect the existing arrangement.
//
// The last parameter is used to resolve ambiguity between this function and
// do_intersect of X_monotone_curve_2 in case that X_monotone_curve_2 and
// Curve_2 are the same class.
// The last parameter should be boost::false_type but we used a
// workaround since it didn't compile in FC3_g++-3.4.4 with the error of:
//
// error: no matching function for call to `do_intersect(Arrangement_2<>&,
// const Arr_segment_2&, const Arr_walk_along_line_point_location<>&,
// mpl_::bool_< true>)'
//
template <typename GeometryTraits_2, typename TopologyTraits,
          typename PointLocation, typename ZoneVisitor>
void insert(Arrangement_on_surface_2<GeometryTraits_2, TopologyTraits>& arr,
            const typename GeometryTraits_2::Curve_2& c,
            const PointLocation& pl, ZoneVisitor &visitor,
            boost::is_same<int, double>::type)
{
  typedef GeometryTraits_2                              Gt2;
  typedef TopologyTraits                                Tt;

  typedef Arrangement_on_surface_2<Gt2, Tt>             Arr;
  typedef ZoneVisitor                                   Zone_visitor;

  typedef typename Gt2::Point_2                         Point_2;
  typedef typename Gt2::X_monotone_curve_2              X_monotone_curve_2;
  typedef boost::variant<Point_2, X_monotone_curve_2>   Make_x_monotone_result;

  // Obtain an arrangement accessor.
  Arr_accessor<Arr> arr_access(arr);

  // Initialize a zone-computation object an a visitor that performs the
  // incremental insertion.
  Arrangement_zone_2<Arr, Zone_visitor> arr_zone(arr, &visitor);

  // Break the input curve into x-monotone subcurves and isolated points.
  std::list<Make_x_monotone_result> x_objects;
  const auto* traits = arr.geometry_traits();
  traits->make_x_monotone_2_object()(c, std::back_inserter(x_objects));

  // Insert each x-monotone curve into the arrangement.
  for (const auto& x_obj : x_objects) {
    // Act according to the type of the current object.
    const auto* x_curve = boost::get<X_monotone_curve_2>(&x_obj);
    if (x_curve != nullptr) {
      // Inserting an x-monotone curve:
      // Initialize the zone-computation object with the given curve.
      arr_zone.init(*x_curve, pl);

      // Notify the arrangement observers that a global operation is about to
      // take place.
      arr_access.notify_before_global_change();

      // Insert the current x-monotone curve into the arrangement.
      arr_zone.compute_zone();

      // Notify the arrangement observers that the global operation has been
      // completed.
      arr_access.notify_after_global_change();
      continue;
    }
    const auto* iso_p = boost::get<Point_2>(&x_obj);
    CGAL_assertion(iso_p != nullptr);

    // Inserting a point into the arrangement:
    insert_point(arr, *iso_p, pl);
  }
}

//-----------------------------------------------------------------------------
// Insert an x-monotone curve into the arrangement (incremental insertion).
// The inserted x-monotone curve may intersect the existing arrangement.
//
// The last parameter is used to resolve ambiguity between this function and
// do_intersect of Curve_2 in case that X_monotone_curve_2 and Curve_2 are the
// same class. The last parameter should be boost::true_type but we used a
// workaround since it didn't compile in FC3_g++-3.4.4 with the error of:
//
// error: no matching function for call to `do_intersect(Arrangement_2<>&,
// const Arr_segment_2&, const Arr_walk_along_line_point_location<>&,
// mpl_::bool_< true>)'
//
//
template <typename GeometryTraits_2, typename TopologyTraits,
          typename PointLocation, typename ZoneVisitor>
void insert(Arrangement_on_surface_2<GeometryTraits_2, TopologyTraits>& arr,
            const typename GeometryTraits_2::X_monotone_curve_2& c,
            const PointLocation& pl, ZoneVisitor &visitor,
            boost::is_same<int, int>::type)
{
  typedef GeometryTraits_2                              Gt2;
  typedef TopologyTraits                                Tt;

  typedef Arrangement_on_surface_2<Gt2, Tt>             Arr;
  typedef ZoneVisitor                                   Zone_visitor;

  // Obtain an arrangement accessor.
  Arr_accessor<Arr> arr_access(arr);

  // Initialize a zone-computation object an a visitor that performs the
  // incremental insertion.
  Arrangement_zone_2<Arr, Zone_visitor> arr_zone(arr, &visitor);

  // Initialize the zone-computation object with the given curve.
  arr_zone.init(c, pl);

  // Notify the arrangement observers that a global operation is about to
  // take place.
  arr_access.notify_before_global_change();

  // Insert the x-monotone curve into the arrangement.
  arr_zone.compute_zone();

  // Notify the arrangement observers that the global operation has been
  // completed.
  arr_access.notify_after_global_change();
}

//-----------------------------------------------------------------------------
// Common interface for the insert of the Curve_2 and X_monotone_curve_2
template <typename GeometryTraits_2, typename TopologyTraits, typename Curve,
          typename PointLocation, typename ZoneVisitor>
void insert(Arrangement_on_surface_2<GeometryTraits_2, TopologyTraits>& arr,
            const Curve& c, const PointLocation& pl, ZoneVisitor &visitor)
{
  typedef GeometryTraits_2                              Gt2;
  typedef typename Gt2::X_monotone_curve_2              X_monotone_curve_2;
  typedef typename boost::is_same<Curve, X_monotone_curve_2>::type
    Is_x_monotone;

  insert(arr, c, pl, visitor, Is_x_monotone());
}

// In some compilers there is a template deduction disambiguity between this
// function and the function receiving two InputIterator.
// For now the solution is to add a dummy variable at the end (referring
// to point-location). Maybe the proper solution is to use boost::enable_if
// together with appropriate tag.
template <typename GeometryTraits_2, typename TopologyTraits, typename Curve,
          typename PointLocation>
void insert(Arrangement_on_surface_2<GeometryTraits_2, TopologyTraits>& arr,
            const Curve& c, const PointLocation& pl,
            typename PointLocation::Point_2*)
{
  typedef TopologyTraits                                Tt;

  typedef typename Tt::Zone_insertion_visitor           Zone_visitor;

  Zone_visitor visitor;
  insert(arr, c, pl, visitor);
}

//-----------------------------------------------------------------------------
// Insert a curve/x-monotone curve into the arrangement (incremental
// insertion).
// The inserted x-monotone curve may intersect the existing arrangement.
// Overloaded version with no point location object - using the default point
// location.
//
template <typename GeometryTraits_2, typename TopologyTraits, typename Curve>
void insert(Arrangement_on_surface_2<GeometryTraits_2, TopologyTraits>& arr,
            const Curve& c)
{
  typedef TopologyTraits                                Tt;

  // Create a default point-location object and use it to insert the curve.
  typename Tt::Default_point_location_strategy def_pl(arr);

  insert(arr, c, def_pl);
}

/*! Insert a range of x-monotone curves into an empty arrangement
 * \param arr the resulting arrangement
 * \param begin the beginning of the curve range
 * \param end past-the-end curve range
 */
template <typename GeometryTraits_2, typename TopologyTraits,
          typename InputIterator>
void
insert_empty(Arrangement_on_surface_2<GeometryTraits_2, TopologyTraits>& arr,
             InputIterator begin_xcurves, InputIterator end_xcurves)
{
  typedef GeometryTraits_2                              Gt2;
  typedef TopologyTraits                                Tt;

  // Arrangement types
  typedef Arrangement_on_surface_2<Gt2, Tt>             Arr;
  typedef typename Arr::Allocator                       Allocator;

  // Surface sweep types
  typedef Arr_construction_event<Gt2, Arr, Allocator>   C_event;
  typedef Arr_construction_subcurve<Gt2, C_event, Allocator>
                                                        C_curve;
  typedef typename Tt::template Construction_helper<C_event, C_curve>
                                                        C_helper;
  typedef Arr_construction_ss_visitor<C_helper>         C_visitor;

  typedef typename C_visitor::Geometry_traits_2         Cgt2;

  const Gt2* geom_traits = arr.geometry_traits();
  C_visitor visitor(&arr);

  /* We would like to avoid copy construction of the geometry traits class.
   * Copy construction is undesired, because it may results with data
   * duplication or even data loss.
   *
   * If the type C_visitor::Geometry_traits_2 is the same as the type
   * GeometryTraits_2, use a reference to GeometryTraits_2 to avoid constructing
   * a new one.  Otherwise, instantiate a local variable of the former and
   * provide the later as a single parameter to the constructor.
   *
   * Use the form 'A a(*b);' and not ''A a = b;' to handle the case where A has
   * only an implicit constructor, (which takes *b as a parameter).
   */
  typename boost::mpl::if_<boost::is_same<Gt2, Cgt2>, const Cgt2&, Cgt2>::type
    traits(*geom_traits);

  // Define a surface-sweep instance and perform the sweep:
  Ss2::Surface_sweep_2<C_visitor> surface_sweep(&traits, &visitor);
  surface_sweep.sweep(begin_xcurves, end_xcurves);
}

/*! Insert a range of x-monotone curves and a range of isolated points into
 * an empty arrangement
 * \param arr the resulting arrangement
 * \param begin_xcurves the beginning of the curve range
 * \param end_xcurves past-the-end curve range
 * \param begin_points the beginning of the point range
 * \param end_points past-the-end point range
 */
template <typename GeometryTraits_2, typename TopologyTraits,
          typename XcInputIterator, typename PInputIterator>
void insert_empty(Arrangement_on_surface_2<GeometryTraits_2, TopologyTraits>&
                  arr,
                  XcInputIterator begin_xcurves, XcInputIterator end_xcurves,
                  PInputIterator begin_points, PInputIterator end_points)
{
  typedef GeometryTraits_2                              Gt2;
  typedef TopologyTraits                                Tt;

  // Arrangement types
  typedef Arrangement_on_surface_2<Gt2, Tt>             Arr;
  typedef typename Arr::Allocator                       Allocator;

  // Surface sweep types
  typedef Arr_construction_event<Gt2, Arr, Allocator>   C_event;
  typedef Arr_construction_subcurve<Gt2, C_event, Allocator>
                                                        C_curve;
  typedef typename Tt::template Construction_helper<C_event, C_curve>
                                                        C_helper;
  typedef Arr_construction_ss_visitor<C_helper>         C_visitor;

  typedef typename C_visitor::Geometry_traits_2         Cgt2;

  const Gt2* geom_traits = arr.geometry_traits();
  C_visitor visitor(&arr);

  /* We would like to avoid copy construction of the geometry traits class.
   * Copy construction is undesired, because it may results with data
   * duplication or even data loss.
   *
   * If the type C_visitor::Geometry_traits_2 is the same as the type
   * GeometryTraits_2, use a reference to GeometryTraits_2 to avoid constructing
   * a new one.  Otherwise, instantiate a local variable of the former and
   * provide the later as a single parameter to the constructor.
   *
   * Use the form 'A a(*b);' and not ''A a = b;' to handle the case where A has
   * only an implicit constructor, (which takes *b as a parameter).
   */
  typename boost::mpl::if_<boost::is_same<Gt2, Cgt2>, const Cgt2&, Cgt2>::type
    traits(*geom_traits);

  // Define a surface-sweep instance and perform the sweep.
  Ss2::Surface_sweep_2<C_visitor> surface_sweep(&traits, &visitor);
  surface_sweep.sweep(begin_xcurves, end_xcurves, begin_points, end_points);
}

/*! Insert a range of x-monotone curves into a non-empty arrangement
 * \param arr the resulting arrangement
 * \param begin the beginning of the curve range
 * \param end past-the-end curve range
 */
template <typename GeometryTraits_2, typename TopologyTraits,
          typename XcInputIterator, typename PInputIterator>
void insert_non_empty(Arrangement_on_surface_2<GeometryTraits_2,
                      TopologyTraits>& arr,
                      XcInputIterator begin_xcurves, XcInputIterator end_xcurves,
                      PInputIterator begin_points, PInputIterator end_points)
{
  typedef GeometryTraits_2                              Gt2;
  typedef TopologyTraits                                Tt;

  // Arrangement types
  typedef Arrangement_on_surface_2<Gt2, Tt>             Arr;
  typedef typename Arr::Allocator                       Allocator;

  // Surface sweep types
  typedef Arr_insertion_traits_2<Gt2, Arr>              Igt2;
  typedef Arr_construction_event<Igt2, Arr, Allocator>  I_event;
  typedef Arr_construction_subcurve<Igt2, I_event, Allocator>
                                                        I_curve;
  typedef typename Tt::template Insertion_helper<I_event, I_curve>
                                                        I_helper;
  typedef Arr_insertion_ss_visitor<I_helper>            I_visitor;
  typedef typename Igt2::X_monotone_curve_2             Ex_x_monotone_curve_2;
  typedef typename Igt2::Point_2                        Ex_point_2;

  const Gt2* geom_traits = arr.geometry_traits();
  I_visitor visitor(&arr);

  /* We would like to avoid copy construction of the geometry traits class.
   * Copy construction is undesired, because it may results with data
   * duplication or even data loss.
   *
   * If the type Igt2 is the same as the type
   * GeometryTraits_2, use a reference to GeometryTraits_2 to avoid constructing
   * a new one.  Otherwise, instantiate a local variable of the former and
   * provide the later as a single parameter to the constructor.
   *
   * Use the form 'A a(*b);' and not ''A a = b;' to handle the case where A has
   * only an implicit constructor, (which takes *b as a parameter).
   */
  typename boost::mpl::if_<boost::is_same<Gt2, Igt2>, const Igt2&, Igt2>::type
    traits(*geom_traits);

  // Create a set of existing as well as new curves and points.
  std::list<Ex_x_monotone_curve_2> ex_cvs;
  std::list<Ex_point_2> ex_pts;

  Ss2::prepare_for_sweep(arr,
                         begin_xcurves, end_xcurves,   // the x-monotone curves
                         begin_points, end_points,     // the points (if any)
                         std::back_inserter(ex_cvs),
                         std::back_inserter(ex_pts),
                         &traits);

  // Define a basic surface-sweep instance and perform the sweep.
  Ss2::Surface_sweep_2<I_visitor> surface_sweep(&traits, &visitor);
  surface_sweep.sweep(ex_cvs.begin(), ex_cvs.end(),ex_pts.begin(), ex_pts.end());
}

//-----------------------------------------------------------------------------
// Insert a range of curves into the arrangement (aggregated insertion).
// The inserted curves may intersect one another and may also intersect the
// existing arrangement.
//
// The last parameter is used to resolve ambiguity between this function and
// do_intersect of X_monotone_curve_2 in case that X_monotone_curve_2 and
// Curve_2 are the same class.
// The last parameter should be boost::false_type but we used a
// workaround since it didn't compile in FC3_g++-3.4.4 with the error of:
//
// error: no matching function for call to `do_intersect(Arrangement_2<>&,
// const Arr_segment_2&, const Arr_walk_along_line_point_location<>&,
// mpl_::bool_< true>)'
//
template <typename GeometryTraits_2, typename TopologyTraits,
          typename InputIterator>
void insert(Arrangement_on_surface_2<GeometryTraits_2, TopologyTraits>& arr,
            InputIterator begin, InputIterator end,
            boost::is_same<int, double>::type)
{
  typedef GeometryTraits_2                              Gt2;
  typedef TopologyTraits                                Tt;

  typedef Arrangement_on_surface_2<Gt2, Tt>             Arr;
  typedef typename Gt2::Point_2                         Point_2;
  typedef typename Gt2::X_monotone_curve_2              X_monotone_curve_2;

  // Obtain an arrangement accessor.
  Arr_accessor<Arr> arr_access(arr);

  // Notify the arrangement observers that a global operation is about to
  // take place.
  arr_access.notify_before_global_change();

  // Subdivide the input curves into x-monotone subcurves and isolated points.
  const Gt2* geom_traits = arr.geometry_traits();
  std::list<X_monotone_curve_2> xcurves;
  std::list<Point_2> iso_points;

  Ss2::make_x_monotone(begin, end,
                       std::back_inserter(xcurves),
                       std::back_inserter(iso_points),
                       geom_traits);

  if (arr.is_empty()) insert_empty(arr, xcurves.begin(), xcurves.end(),
                                   iso_points.begin(), iso_points.end());
  else insert_non_empty(arr, xcurves.begin(), xcurves.end(),
                        iso_points.begin(), iso_points.end());

  // Notify the arrangement observers that the global operation has been
  // completed.
  arr_access.notify_after_global_change();
}

//-----------------------------------------------------------------------------
// Insert a range of x-monotone curves into the arrangement (aggregated
// insertion). The inserted x-monotone curves may intersect one another and
// may also intersect the existing arrangement.
//
// The last parameter is used to resolve ambiguity between this function and
// insert of Curve_2 in case that X_monotone_curve_2 and Curve_2 are the
// same class. The last parameter should be boost::true_type but we used a
// workaround since it didn't compile in FC3_g++-3.4.4 with the error of:
//
// error: no matching function for call to `do_intersect(Arrangement_2<>&,
// const Arr_segment_2&, const Arr_walk_along_line_point_location<>&,
// mpl_::bool_< true>)'
//
template <typename GeometryTraits_2, typename TopologyTraits,
          typename InputIterator>
void insert(Arrangement_on_surface_2<GeometryTraits_2, TopologyTraits>& arr,
            InputIterator begin, InputIterator end,
            boost::is_same<int, int>::type)
{
  typedef GeometryTraits_2                              Gt2;
  typedef TopologyTraits                                Tt;

  typedef Arrangement_on_surface_2<Gt2, Tt>             Arr;

  // Obtain an arrangement accessor.
  Arr_accessor<Arr> arr_access(arr);

  // Notify the arrangement observers that a global operation is about to
  // take place.
  arr_access.notify_before_global_change();

  // Choose the operation depending on whether the input arrangement is
  // empty (then we construct it from scratch), or not (where we just insert
  // the new curves).
  if (arr.is_empty()) insert_empty(arr, begin, end);
  else {
    // The arrangement is not empty: use the insertion visitor.
    std::list<typename Gt2::Point_2> empty;
    insert_non_empty(arr, begin, end, empty.begin(), empty.end());
  }

  // Notify the arrangement observers that the global operation has been
  // completed.
  arr_access.notify_after_global_change();
}

//-----------------------------------------------------------------------------
// Common interface for the inserts of the Curve_2 and X_monotone_curve_2
template <typename GeometryTraits_2, typename TopologyTraits,
          typename InputIterator>
void insert(Arrangement_on_surface_2<GeometryTraits_2,TopologyTraits>& arr,
            InputIterator begin, InputIterator end)
{
  typedef GeometryTraits_2                              Gt2;

  typedef typename Gt2::X_monotone_curve_2              X_monotone_curve_2;
  typedef typename std::iterator_traits<InputIterator>::value_type
                                                        Iterator_value_type;

  typedef typename boost::is_same<Iterator_value_type,X_monotone_curve_2>::type
                                                        Is_x_monotone;

  return insert(arr, begin, end, Is_x_monotone());
}

//-----------------------------------------------------------------------------
// Insert an x-monotone curve into the arrangement (incremental insertion)
// when the location of the left endpoint of the curve is known and is
// given as an isertion hint.
// The inserted x-monotone curve may intersect the existing arrangement.
//
template <typename GeometryTraits_2, typename TopologyTraits>
void insert(Arrangement_on_surface_2<GeometryTraits_2, TopologyTraits>& arr,
            const typename GeometryTraits_2::X_monotone_curve_2& c,
            typename Arr_point_location_result<
              Arrangement_on_surface_2<GeometryTraits_2,
                                       TopologyTraits> >::type obj)
{
  typedef GeometryTraits_2                              Gt2;
  typedef TopologyTraits                                Tt;

  typedef Arrangement_on_surface_2<Gt2, Tt>             Arr;
  typedef typename Tt::Zone_insertion_visitor           Zone_visitor;

  // Obtain an arrangement accessor.
  Arr_accessor<Arr> arr_access(arr);

  // Define a zone-computation object an a visitor that performs the
  // incremental insertion.
  Zone_visitor visitor;
  Arrangement_zone_2<Arr, Zone_visitor> arr_zone(arr, &visitor);

  // Initialize the zone-computation object with the given curve.
  arr_zone.init_with_hint(c, obj);

  // Notify the arrangement observers that a global operation is about to
  // take place.
  arr_access.notify_before_global_change();

  // Insert the x-monotone curve into the arrangement.
  arr_zone.compute_zone();

  // Notify the arrangement observers that the global operation has been
  // completed.
  arr_access.notify_after_global_change();
}

// ----------------------------------------------------------------------------
// backward compatibility functions.
/* DEPRECATED use insert() instead */
template <typename GeometryTraits_2, typename TopologyTraits,
          typename PointLocation>
CGAL_DEPRECATED void insert_x_monotone_curve
(Arrangement_on_surface_2<GeometryTraits_2, TopologyTraits>& arr,
 const typename GeometryTraits_2::X_monotone_curve_2& c,
 const PointLocation& pl)
{
  insert(arr, c, pl);
}

/* DEPRECATED use insert() instead */
template <typename GeometryTraits_2, typename TopologyTraits>
CGAL_DEPRECATED void insert_x_monotone_curve
(Arrangement_on_surface_2<GeometryTraits_2, TopologyTraits>& arr,
 const typename GeometryTraits_2::X_monotone_curve_2& c)
{
  insert(arr, c);
}

/* DEPRECATED use insert() instead */
template <typename GeometryTraits_2, typename TopologyTraits,
          typename InputIterator>
CGAL_DEPRECATED void insert_x_monotone_curves
(Arrangement_on_surface_2<GeometryTraits_2, TopologyTraits>& arr,
 InputIterator begin, InputIterator end)
{
  insert(arr, begin, end);
}

/* DEPRECATED use insert() instead */
template <typename GeometryTraits_2, typename TopologyTraits>
CGAL_DEPRECATED void insert_x_monotone_curve
(Arrangement_on_surface_2<GeometryTraits_2, TopologyTraits>& arr,
 const typename GeometryTraits_2::X_monotone_curve_2& c,
 typename Arr_point_location_result<
   Arrangement_on_surface_2<GeometryTraits_2, TopologyTraits> >::type obj)
{
  insert(arr, c, obj);
}

/* DEPRECATED use insert() instead */
template <typename GeometryTraits_2, typename TopologyTraits,
          typename PointLocation>
CGAL_DEPRECATED
void insert_curve(Arrangement_on_surface_2<GeometryTraits_2, TopologyTraits>&
                  arr,
                  const typename GeometryTraits_2::Curve_2& c,
                  const PointLocation& pl)
{
  insert(arr, c, pl);
}

/* DEPRECATED use insert() instead */
template <typename GeometryTraits_2, typename TopologyTraits>
CGAL_DEPRECATED
void insert_curve(Arrangement_on_surface_2<GeometryTraits_2, TopologyTraits>&
                  arr,
                  const typename GeometryTraits_2::Curve_2& c)
{
  insert(arr, c);
}

/* DEPRECATED use insert() instead */
template <typename GeometryTraits_2, typename TopologyTraits,
          typename InputIterator>
CGAL_DEPRECATED
void insert_curves(Arrangement_on_surface_2<GeometryTraits_2, TopologyTraits>&
                   arr,
                   InputIterator begin, InputIterator end)
{
  insert(arr, begin, end);
}

//-----------------------------------------------------------------------------
// Insert an x-monotone curve into the arrangement, such that the curve
// interior does not intersect with any existing edge or vertex in the
// arragement (incremental insertion).
//
template <typename GeometryTraits_2, typename TopologyTraits,
          typename PointLocation>
typename Arrangement_on_surface_2<GeometryTraits_2, TopologyTraits>::
Halfedge_handle
insert_non_intersecting_curve
(Arrangement_on_surface_2<GeometryTraits_2, TopologyTraits>& arr,
 const typename GeometryTraits_2::X_monotone_curve_2& c,
 const PointLocation& pl)
{
  typedef GeometryTraits_2                              Gt2;
  typedef TopologyTraits                                Tt;

  typedef Arrangement_on_surface_2<Gt2, Tt>             Arr;

  typedef Arr_traits_basic_adaptor_2<typename Arr::Geometry_traits_2>
                                                        Traits_adaptor_2;
  typedef typename Arr::Vertex_const_handle             Vertex_const_handle;
  typedef typename Arr::Halfedge_const_handle           Halfedge_const_handle;
  typedef typename Arr::Face_const_handle               Face_const_handle;
  CGAL_USE_TYPE(Halfedge_const_handle);

  const Traits_adaptor_2* geom_traits =
    static_cast<const Traits_adaptor_2*>(arr.geometry_traits());
  Arr_accessor<Arr> arr_access(arr);

  // Check whether the left end has boundary conditions, and locate it in the
  // arrangement accordingly.
  auto bx1 = geom_traits->parameter_space_in_x_2_object()(c, ARR_MIN_END);
  auto by1 = geom_traits->parameter_space_in_y_2_object()(c, ARR_MIN_END);
  const Vertex_const_handle* vh1 = nullptr;

  typedef Arr_point_location_result<Arr>        Pl_result;

  typename Pl_result::type obj1;
  if ((bx1 == ARR_INTERIOR) && (by1 == ARR_INTERIOR)) {
    // We have a normal left endpoint with no boundary conditions:
    // use a point-location query.
    obj1 = pl.locate(geom_traits->construct_min_vertex_2_object()(c));

    // The endpoint must not lie on an existing edge, but may coincide with
    // and existing vertex vh1.
    CGAL_precondition_msg(boost::get<Halfedge_const_handle>(&obj1) == nullptr,
                          "The curve must not intersect an existing edge.");

  }
  else {
    // We have a left end with boundary conditions. Use the accessor to locate
    // the feature that contains it.
    obj1 = arr_access.locate_curve_end(c, ARR_MIN_END, bx1, by1);
    CGAL_precondition_msg(boost::get<Halfedge_const_handle>(&obj1) == nullptr,
                          "The curve must not overlap an existing edge.");
  }
  vh1 = Pl_result::template assign<Vertex_const_handle>(&obj1);

  // Check whether the right end has boundary conditions, and locate it in the
  // arrangement accordingly.
  auto bx2 = geom_traits->parameter_space_in_x_2_object()(c, ARR_MAX_END);
  auto by2 = geom_traits->parameter_space_in_y_2_object()(c, ARR_MAX_END);
  const Vertex_const_handle* vh2 = nullptr;

  typename Pl_result::type obj2;
  if ((bx2 == ARR_INTERIOR) && (by2 == ARR_INTERIOR)) {
    // We have a normal right endpoint with no boundary conditions:
    // use a point-location query.
    obj2 = pl.locate(geom_traits->construct_max_vertex_2_object()(c));

    // The endpoint must not lie on an existing edge, but may coincide with
    // and existing vertex vh2.
    CGAL_precondition_msg(boost::get<Halfedge_const_handle>(&obj2) == nullptr,
                          "The curve must not intersect an existing edge.");
  }
  else {
    // We have a right end with boundary conditions. Use the accessor to locate
    // the feature that contains it.
    // std::cout << "before locate_curve_end()"
    //           << ", bx2: " << bx2
    //           << ", by2: " << by2
    //           << std::endl;
    obj2 = arr_access.locate_curve_end(c, ARR_MAX_END, bx2, by2);
    CGAL_precondition_msg(boost::get<Halfedge_const_handle>(&obj2) == nullptr,
                          "The curve must not overlap an existing edge.");
  }
  vh2 = Pl_result::template assign<Vertex_const_handle>(&obj2);

  // Notify the arrangement observers that a global operation is about to
  // take place.
  arr_access.notify_before_global_change();

  // Check whether the located features containing the curve endpoints
  // are vertices or faces, and use the proper specialized insertion function
  // accordingly.
  typename Arr::Halfedge_handle new_he;

  if (vh1 != nullptr) {
    if (vh2 != nullptr) {
      // Both endpoints are associated with a existing vertices.
      // In this case insert_at_vertices() already returns a halfedge
      // directed from left to right.
      new_he = arr.insert_at_vertices(c,
                                      arr.non_const_handle(*vh1),
                                      arr.non_const_handle(*vh2));
    }
    else {
      // Only the left endpoint is associated with an existing vertex.
      // In this case insert_from_left_vertex() returns a halfedge directed
      // to the new vertex it creates, so it is already directed from left to
      // right.
      new_he = arr.insert_from_left_vertex(c, arr.non_const_handle(*vh1));
    }
  }
  else {
    if (vh2 != nullptr) {
      // Only the right endpoint is associated with an existing vertex.
      // In this case insert_from_left_vertex() returns a halfedge directed
      // to the new vertex it creates, so it is directed from right to left
      // and we take its twin halfedge instead.
      new_he = arr.insert_from_right_vertex(c, arr.non_const_handle(*vh2));
      new_he = new_he->twin();
    }
    else {
      // Both endpoints are not associated with existing vertices, so
      // we must insert the curve in the interior of a face.
      // In this case insert_in_face_interior() already returns a halfedge
      // directed from left to right.
      const Face_const_handle* fh1 = boost::get<Face_const_handle>(&obj1);
      const Face_const_handle* fh2 = boost::get<Face_const_handle>(&obj2);

      // std::cout << arr << std::endl;
      // std::cout << "(*fh1)->number_of_outer_ccbs(): "
      //           << (*fh1)->number_of_outer_ccbs() << std::endl;
      // std::cout << "(*fh2)->number_of_outer_ccbs(): "
      //           << (*fh2)->number_of_outer_ccbs() << std::endl;

      CGAL_assertion_msg
        ((fh1 != nullptr) && (fh2 != nullptr) && ((*fh1) == (*fh2)),
         "The curve intersects the interior of existing edges.");

      if ((fh1 != nullptr) && (fh2 != nullptr) && (*fh1 == *fh2)) {
        new_he = arr.insert_in_face_interior(c, arr.non_const_handle (*fh1));
      }
    }
  }

  // Notify the arrangement observers that the global operation has been
  // completed.
  arr_access.notify_after_global_change();

  // Return the resulting halfedge from the insertion operation.
  return new_he;
}

//-----------------------------------------------------------------------------
// Insert an x-monotone curve into the arrangement, such that the curve
// interior does not intersect with any existing edge or vertex in the
// arragement (incremental insertion).
// Overloaded version with no point location object.
//
template <typename GeometryTraits_2, typename TopologyTraits>
typename Arrangement_on_surface_2<GeometryTraits_2, TopologyTraits>::
Halfedge_handle
insert_non_intersecting_curve
(Arrangement_on_surface_2<GeometryTraits_2, TopologyTraits>& arr,
 const typename GeometryTraits_2::X_monotone_curve_2& c)
{
  typedef TopologyTraits                                Tt;

  // Create a default point-location object and use it to insert the curve.
  typename Tt::Default_point_location_strategy def_pl(arr);

  return (insert_non_intersecting_curve(arr, c, def_pl));
}

/*! Insert a range of x-monotone curves into an empty arrangement
 * \param arr the resulting arrangement
 * \param begin the beginning of the curve range
 * \param end past-the-end curve range
 */
template <typename GeometryTraits_2, typename TopologyTraits,
          typename InputIterator>
void non_intersecting_insert_empty(Arrangement_on_surface_2<GeometryTraits_2,
                                                            TopologyTraits>& arr,
                                   InputIterator begin_xcurves,
                                   InputIterator end_xcurves)
{
  typedef GeometryTraits_2                              Gt2;
  typedef TopologyTraits                                Tt;

  // Arrangement types
  typedef Arrangement_on_surface_2<Gt2, Tt>             Arr;
  typedef typename Arr::Allocator                       Allocator;

  // Surface sweep types

  // The forth parameter of Arr_construction_subcurve is the base class of
  // Arr_construction_subcurve. By default Arr_construction_subcurve derives
  // from Default_subcurve, which is suitable for general curves
  // including overlapping curves. Here we bypass Default_subcurve and
  // force Arr_construction_subcurve to derive directly from
  // No_overlap_subcurve (which is the base class of Default_subcurve).
  typedef Arr_construction_event<Gt2, Arr, Allocator,
                                 Ss2::No_overlap_event_base,
                                 Ss2::No_overlap_subcurve>
                                                        Nxc_event;
  typedef Arr_construction_subcurve<Gt2, Nxc_event, Allocator,
                                    Ss2::No_overlap_subcurve>
                                                        Nxc_curve;
  typedef typename Tt::template No_intersection_construction_helper<Nxc_event,
                                                                    Nxc_curve>
                                                        Nxc_helper;
  typedef Arr_construction_ss_visitor<Nxc_helper>       Nxc_visitor;

  const Gt2* traits = arr.geometry_traits();
  Nxc_visitor visitor(&arr);

  // Define a basic surface-sweep instance (which is not supposed to handle
  // insersections) and perform the sweep.
  Ss2::No_intersection_surface_sweep_2<Nxc_visitor>
    surface_sweep(traits, &visitor);
  surface_sweep.sweep(begin_xcurves, end_xcurves);
}

/*! Insert a range of x-monotone curves into an empty arrangement
 * \param arr the resulting arrangement
 * \param begin the beginning of the curve range
 * \param end past-the-end curve range
 */
template <typename GeometryTraits_2, typename TopologyTraits,
          typename XcInputIterator, typename PInputIterator>
void non_intersecting_insert_empty(Arrangement_on_surface_2<GeometryTraits_2,
                                                            TopologyTraits>& arr,
                                   XcInputIterator begin_xcurves,
                                   XcInputIterator end_xcurves,
                                   PInputIterator begin_points,
                                   PInputIterator end_points)
{
  typedef GeometryTraits_2                              Gt2;
  typedef TopologyTraits                                Tt;

  // Arrangement types
  typedef Arrangement_on_surface_2<Gt2, Tt>             Arr;
  typedef typename Arr::Allocator                       Allocator;

  // Surface sweep types
  // Type definition for the no-intersection construction surface-sweep visitor.
  // The third parameter of Arr_construction_subcurve is the base class of
  // Arr_construction_subcurve. By default Arr_construction_subcurve derives
  // from Default_subcurve, which is suitable for general curves
  // including overlapping curves. Here we bypass Default_subcurve and
  // force Arr_construction_subcurve to derive directly from
  // No_overlap_subcurve (which is the base class of Default_subcurve).
  typedef Arr_construction_event<Gt2, Arr, Allocator,
                                 Ss2::No_overlap_event_base,
                                 Ss2::No_overlap_subcurve>
                                                        Nxc_event;
  typedef Arr_construction_subcurve<Gt2, Nxc_event, Allocator,
                                    Ss2::No_overlap_subcurve>
                                                        Nxc_curve;
  typedef typename Tt::template No_intersection_construction_helper<Nxc_event,
                                                                    Nxc_curve>
                                                        Nxc_helper;
  typedef Arr_construction_ss_visitor<Nxc_helper>       Nxc_visitor;

  const Gt2* traits = arr.geometry_traits();
  Nxc_visitor visitor(&arr);

  // Define a basic surface-sweep instance (which is not supposed to handle
  // insersections) and perform the sweep.
  Ss2::No_intersection_surface_sweep_2<Nxc_visitor>
    surface_sweep(traits, &visitor);
  surface_sweep.sweep(begin_xcurves, end_xcurves, begin_points, end_points);
}

/*! Insert a range of x-monotone curves into a non-empty arrangement
 * \param arr the resulting arrangement
 * \param begin the beginning of the curve range
 * \param end past-the-end curve range
 */
template <typename GeometryTraits_2, typename TopologyTraits,
          typename XcInputIterator, typename PInputIterator>
void
non_intersecting_insert_non_empty(Arrangement_on_surface_2<GeometryTraits_2,
                                                           TopologyTraits>& arr,
                                  XcInputIterator begin_xcurves,
                                  XcInputIterator end_xcurves,
                                  PInputIterator begin_points,
                                  PInputIterator end_points)
{
  typedef GeometryTraits_2                              Gt2;
  typedef TopologyTraits                                Tt;

  // Arrangement types
  typedef Arrangement_on_surface_2<Gt2, Tt>             Arr;
  typedef typename Arr::Allocator                       Allocator;

  // Surface sweep types
  typedef Arr_basic_insertion_traits_2<Gt2, Arr>        Igt2;

  // The third parameter of Arr_construction_subcurve is the base class of
  // Arr_construction_subcurve. By default Arr_construction_subcurve derives
  // from Default_subcurve, which is suitable for general curves
  // including overlapping curves. Here we bypass Default_subcurve and
  // force Arr_construction_subcurve to derive directly from
  // No_overlap_subcurve (which is the base class of Default_subcurve).
  typedef Arr_construction_event<Igt2, Arr, Allocator,
                                 Ss2::No_overlap_event_base,
                                 Ss2::No_overlap_subcurve>
                                                        Nxi_event;
  typedef Arr_construction_subcurve<Igt2, Nxi_event, Allocator,
                                    Ss2::No_overlap_subcurve>
                                                        Nxi_curve;
  // typedef typename Tt::template No_intersection_insertion_event<Allocator>
  //                                                       Nxi_event;
  // typedef typename Tt::template No_intersection_insertion_curve<Nxi_event>
  //                                                       Nxi_curve;
  typedef typename Tt::template No_intersection_insertion_helper<Nxi_event,
                                                                 Nxi_curve>
                                                        Nxi_Helper;
  typedef Arr_no_intersection_insertion_ss_visitor<Nxi_Helper>
                                                        Nxi_visitor;
  typedef typename Igt2::X_monotone_curve_2             Ex_x_monotone_curve_2;
  typedef typename Igt2::Point_2                        Ex_point_2;

  const Gt2* geom_traits = arr.geometry_traits();
  Nxi_visitor visitor(&arr);

  /* We would like to avoid copy construction of the geometry traits class.
   * Copy construction is undesired, because it may results with data
   * duplication or even data loss.
   *
   * If the type Nxi_visitor::Geometry_traits_2 is the same as the type
   * GeometryTraits_2, use a reference to GeometryTraits_2 to avoid constructing
   * a new one.  Otherwise, instantiate a local variable of the former and
   * provide the later as a single parameter to the constructor.
   *
   * Use the form 'A a(*b);' and not ''A a = b;' to handle the case where A has
   * only an implicit constructor, (which takes *b as a parameter).
   */
  typename boost::mpl::if_<boost::is_same<Gt2, Igt2>, const Igt2&, Igt2>::type
    traits(*geom_traits);

  // Create a set of existing as well as new curves and points.
  std::list<Ex_x_monotone_curve_2> ex_cvs;
  std::list<Ex_point_2> ex_pts;

  Ss2::prepare_for_sweep(arr,
                         begin_xcurves, end_xcurves,   // the x-monotone curves
                         begin_points, end_points,     // the points (if any)
                         std::back_inserter(ex_cvs),
                         std::back_inserter(ex_pts),
                         &traits);

  // Define a basic surface-sweep instance and perform the sweep.
  Ss2::No_intersection_surface_sweep_2<Nxi_visitor>
    surface_sweep(&traits, &visitor);
  surface_sweep.sweep(ex_cvs.begin(), ex_cvs.end(),
                      ex_pts.begin(), ex_pts.end());
}

//-----------------------------------------------------------------------------
// Insert a range of pairwise interior-disjoint x-monotone curves into
// the arrangement, such that the curve interiors do not intersect with
// any existing edge or vertex in the arragement (aggregated insertion).
//
template <typename GeometryTraits_2, typename TopologyTraits,
          typename InputIterator>
void insert_non_intersecting_curves
(Arrangement_on_surface_2<GeometryTraits_2, TopologyTraits>& arr,
 InputIterator begin, InputIterator end)
{
  typedef GeometryTraits_2                              Gt2;
  typedef TopologyTraits                                Tt;

  typedef Arrangement_on_surface_2<Gt2, Tt>             Arr;

  // Obtain an arrangement accessor.
  Arr_accessor<Arr> arr_access(arr);

  // Notify the arrangement observers that a global operation is about to
  // take place.
  arr_access.notify_before_global_change();

  // Choose the operation depending on whether the input arrangement is
  // empty (then we construct it from scratch), or not (where we just insert
  // the new curves).
  if (arr.is_empty()) non_intersecting_insert_empty(arr, begin, end);
  else {
    std::list<typename Gt2::Point_2> empty;
    non_intersecting_insert_non_empty(arr, begin, end,
                                      empty.begin(), empty.end());
  }

  // Notify the arrangement observers that the global operation has been
  // completed.
  arr_access.notify_after_global_change();
}

//-----------------------------------------------------------------------------
// Remove an edge from the arrangement. In case it is possible to merge
// the edges incident to the end-vertices of the removed edge after its
// deletion, the function performs these merges as well.
//
template <typename GeometryTraits_2, typename TopologyTraits>
typename Arrangement_on_surface_2<GeometryTraits_2, TopologyTraits>::Face_handle
remove_edge(Arrangement_on_surface_2<GeometryTraits_2, TopologyTraits>& arr,
            typename Arrangement_on_surface_2<GeometryTraits_2,
            TopologyTraits>::Halfedge_handle e)
{
  typedef GeometryTraits_2                              Gt2;
  typedef TopologyTraits                                Tt;

  typedef Arrangement_on_surface_2<Gt2, Tt>             Arr;
  typedef Arr_traits_adaptor_2<Gt2>                     Traits_adaptor_2;

  // Notify the arrangement observers that a global operation is about to
  // take place.
  Arr_accessor<Arr> arr_access(arr);

  arr_access.notify_before_global_change();

  // Keep track of the end-vertices of the edge we are about to remove.
  typename Arr::Vertex_handle v_ends[2];
  bool is_removed[2];

  v_ends[0] = e->source();
  is_removed[0] =
    (v_ends[0]->is_at_open_boundary() || (v_ends[0]->degree() == 1));
  v_ends[1] = e->target();
  is_removed[1] =
    (v_ends[1]->is_at_open_boundary() || (v_ends[1]->degree() == 1));

  // Remove the edge from the arrangement.
  typename Arr::Face_handle    face = arr.remove_edge (e);

  // Examine the end-vertices: If a vertex has now two incident edges, and the
  // curves associated with these edges can be merged, merge the two edges and
  // remove the vertex.

  const Traits_adaptor_2* traits =
    static_cast<const Traits_adaptor_2*>(arr.geometry_traits());

  typename Arr::Halfedge_around_vertex_circulator circ;
  typename Arr::Halfedge_handle e1, e2;
  for (size_t i = 0; i < 2; i++) {
    if (! is_removed[i] && v_ends[i]->degree() == 2) {
      // Get the two edges incident to the end-vertex.
      circ = v_ends[i]->incident_halfedges();
      e1 = circ;
      ++circ;
      e2 = circ;

      // Check if it is possible to merge the two edges.
      if (traits->are_mergeable_2_object() (e1->curve(), e2->curve())) {
        // Merge the two curves.
        typename Gt2::X_monotone_curve_2   cv;
        traits->merge_2_object()(e1->curve(), e2->curve(), cv);

        // Merge the two edges in the arrangement.
        arr.merge_edge(e1, e2, cv);
      }
    }
  }

  // Notify the arrangement observers that the global operation has been
  // completed.
  arr_access.notify_after_global_change();

  // Return the face remaining after the removal of the edge.
  return face;
}

//-----------------------------------------------------------------------------
// Insert a vertex that corresponds to a given point into the arrangement.
// The inserted point may lie on any existing arrangement feature.
//
template <typename GeometryTraits_2, typename TopologyTraits,
          typename PointLocation>
typename Arrangement_on_surface_2<GeometryTraits_2, TopologyTraits>::Vertex_handle
insert_point(Arrangement_on_surface_2<GeometryTraits_2, TopologyTraits>& arr,
             const typename GeometryTraits_2::Point_2& p,
             const PointLocation& pl)
{
  typedef GeometryTraits_2                              Gt2;
  typedef TopologyTraits                                Tt;

  typedef Arrangement_on_surface_2<Gt2, Tt>             Arr;

  typedef typename Arr::Vertex_const_handle             Vertex_const_handle;
  typedef typename Arr::Halfedge_const_handle           Halfedge_const_handle;
  typedef typename Arr::Face_const_handle               Face_const_handle;

  // Act according to the type of arrangement feature that contains the point.
  typename Arr::Vertex_handle vh_for_p;

  // Locate the given point in the arrangement.
  auto obj = pl.locate(p);

  // Notify the arrangement observers that a global operation is about to
  // take place.
  Arr_accessor<Arr> arr_access(arr);

  arr_access.notify_before_global_change();

  const Face_const_handle* fh = boost::get<Face_const_handle>(&obj);
  if (fh != nullptr) {
    // p lies inside a face: Insert it as an isolated vertex it the interior of
    // this face.
    vh_for_p = arr.insert_in_face_interior(p, arr.non_const_handle (*fh));
  }
  else {
    const Halfedge_const_handle* hh = boost::get<Halfedge_const_handle>(&obj);
    if (hh != nullptr) {
      // p lies in the interior of an edge: Split this edge to create a new
      // vertex associated with p.
      typename Gt2::X_monotone_curve_2 sub_cv1, sub_cv2;
      typename Arr::Halfedge_handle split_he;

      const auto* gt = arr.geometry_traits();
      gt->split_2_object()((*hh)->curve(), p, sub_cv1, sub_cv2);
      split_he = arr.split_edge(arr.non_const_handle(*hh), sub_cv1, sub_cv2);

      // The new vertex is the target of the returned halfedge.
      vh_for_p = split_he->target();
    }
    else {
      // p lies on an existing vertex, so we just update this vertex.
      const Vertex_const_handle* vh = boost::get<Vertex_const_handle>(&obj);
      CGAL_assertion(vh != nullptr);
      vh_for_p = arr.modify_vertex (arr.non_const_handle (*vh), p);
    }
  }

  // Notify the arrangement observers that the global operation has been
  // completed.
  arr_access.notify_after_global_change();

  // Return a handle for the vertex associated with p.
  return vh_for_p;
}

//-----------------------------------------------------------------------------
// Insert a vertex that corresponds to a given point into the arrangement.
// The inserted point may lie on any existing arrangement feature.
//
template <typename GeometryTraits_2, typename TopologyTraits>
typename Arrangement_on_surface_2<GeometryTraits_2, TopologyTraits>::
Vertex_handle
insert_point(Arrangement_on_surface_2<GeometryTraits_2, TopologyTraits>& arr,
             const typename GeometryTraits_2::Point_2& p)
{
  typedef TopologyTraits                                Tt;

  // Create a default point-location object and use it to insert the point.
  typename Tt::Default_point_location_strategy def_pl(arr);

  return insert_point(arr, p, def_pl);
}

//-----------------------------------------------------------------------------
// Remove a vertex from the arrangement.
//
template <typename GeometryTraits_2, typename TopologyTraits>
bool remove_vertex(Arrangement_on_surface_2<GeometryTraits_2, TopologyTraits>&
                   arr,
                   typename Arrangement_on_surface_2<
                     GeometryTraits_2, TopologyTraits>::Vertex_handle v)
{
  typedef GeometryTraits_2                              Gt2;
  typedef TopologyTraits                                Tt;

  typedef Arrangement_on_surface_2<Gt2, Tt>             Arr;
  typedef Arr_traits_adaptor_2<Gt2>                     Traits_adaptor_2;

  // Notify the arrangement observers that a global operation is about to
  // take place.
  Arr_accessor<Arr> arr_access(arr);

  arr_access.notify_before_global_change();

  // Act according to the number of edges incident to v.
  bool removed = false;

  if (v->is_isolated()) {
    // In case v is an isolated vertex, simply remove it.
    arr.remove_isolated_vertex (v);
    removed = true;
  }
  else if (v->degree() == 2) {
    // If the vertex has now two incident edges, and the curves associated
    // with these edges can be merged, merge the two edges and remove the
    // vertex.
    const Traits_adaptor_2* traits =
      static_cast<const Traits_adaptor_2*>(arr.geometry_traits());
    typename Arr::Halfedge_around_vertex_circulator circ;
    typename Arr::Halfedge_handle e1, e2;

    circ = v->incident_halfedges();
    e1 = circ;
    ++circ;
    e2 = circ;

    if (traits->are_mergeable_2_object() (e1->curve(), e2->curve())) {
      // Merge the two curves.
      typename Gt2::X_monotone_curve_2   cv;
      traits->merge_2_object()(e1->curve(), e2->curve(), cv);

      // Merge the two edges in the arrangement.
      arr.merge_edge(e1, e2, cv);
      removed = true;
    }
  }

  // Notify the arrangement observers that the global operation has been
  // completed.
  arr_access.notify_after_global_change();

  // Return an indication whether the vertex has been removed or not.
  return removed;
}

//-----------------------------------------------------------------------------
// Check the validity of the arrangement. In particular, check that the
// edegs are disjoint-interior, and the holes are located in their proper
// position.
//
template <typename GeometryTraits_2, typename TopologyTraits>
bool
is_valid(const Arrangement_on_surface_2<GeometryTraits_2, TopologyTraits>& arr)
{
  typedef GeometryTraits_2                              Gt2;
  typedef TopologyTraits                                Tt;

  // Arrangement types (iterator and circulator types).
  typedef Arrangement_on_surface_2<Gt2, Tt>             Arr;
  typedef typename Arr::Allocator                       Allocator;
  typedef typename Arr::Edge_const_iterator             Edge_const_iterator;
  typedef typename Arr::Halfedge_const_handle           Halfedge_const_handle;
  typedef typename Arr::Inner_ccb_const_iterator        Inner_ccb_const_iterator;
  typedef typename Arr::Face_const_iterator             Face_const_iterator;
  typedef typename Arr::Face_const_handle               Face_const_handle;
  typedef typename Arr::Vertex_const_handle             Vertex_const_handle;
  typedef typename Arr::Isolated_vertex_const_iterator
    Isolated_vertex_const_iterator;
  typedef typename Arr::Halfedge_around_vertex_const_circulator
    Halfedge_around_vertex_const_circulator;

  typedef typename Gt2::X_monotone_curve_2              X_monotone_curve_2;

  // The surface-sweep types:
  typedef Ss2::Do_interior_intersect_visitor<Gt2, Allocator>
                                                        Visitor;
  typedef Ss2::Surface_sweep_2<Visitor>                 Surface_sweep_2;

  // First use the internal validity check.
  if (!arr.is_valid()) return (false);

  // Perform a sweep over all subcurves associated with arrangement edges.
  std::vector<X_monotone_curve_2> curves_vec(arr.number_of_edges());
  unsigned int i = 0;

  Edge_const_iterator eit;
  for (eit = arr.edges_begin(); eit != arr.edges_end(); ++eit, i++)
    curves_vec[i] = eit->curve();

  Visitor visitor;
  const Gt2* traits = arr.geometry_traits();
  Surface_sweep_2 surface_sweep(traits, &visitor);
  visitor.sweep_xcurves(curves_vec.begin(), curves_vec.end());
  bool are_edges_disjoint = (! visitor.found_intersection());

  if (!are_edges_disjoint) {
    CGAL_warning_msg(are_edges_disjoint,
                     "Arrangement edges are not disjoint in their interior.");
    return false;
  }

  // Check that the holes and isolated vertices are located where they should.
  // At the same time, we prepare a vector that consists of all isolated
  // vertices and all leftmost vertices from every hole.
  std::list<std::pair<Vertex_const_handle, Face_const_handle> > vf_list;

  typename Gt2::Compare_xy_2 compare_xy = traits->compare_xy_2_object();
  Face_const_iterator fit;
  Face_const_handle fh;
  Inner_ccb_const_iterator ic_it;
  Halfedge_const_handle ccb;
  Isolated_vertex_const_iterator iv_it;
  Vertex_const_handle left_v;
  bool is_first;

  for (fit = arr.faces_begin(); fit != arr.faces_end(); ++fit) {
    // Check all holes in the current face.
    fh = fit;
    for (ic_it = fh->inner_ccbs_begin(); ic_it != fh->inner_ccbs_end(); ++ic_it)
    {
      ccb = *ic_it;
      is_first = true;

      do {
        if (ccb->face() != fit)
          return (false);

        if (is_first ||
            compare_xy (ccb->target()->point(), left_v->point()) == SMALLER)
        {
          left_v = ccb->target();
          is_first = false;
        }

        ccb = ccb->next();

      } while (ccb != *ic_it);

      vf_list.push_back (std::make_pair (left_v, fh));
    }

    // Check all isolated vertices in the current face.
    for (iv_it = fh->isolated_vertices_begin();
         iv_it != fh->isolated_vertices_end(); ++iv_it)
    {
      if (iv_it->face() != fit) return (false);

      vf_list.push_back (std::make_pair (Vertex_const_handle(iv_it), fh));
    }
  }

  // Shoot a vertical ray from each vertex we have collected downward, and
  // check that this vertex is really contained in the proper face.
  auto comp_y_at_x_right = traits->compare_y_at_x_right_2_object();
  auto comp_y_at_x_left = traits->compare_y_at_x_left_2_object();

  typename Tt::Default_point_location_strategy def_pl(arr);
  const Halfedge_const_handle invalid_he;

  Face_const_handle in_face;
  for (auto vf_iter = vf_list.begin(); vf_iter != vf_list.end(); ++vf_iter) {
    // Perform ray-shooting from the current vertex.
    Vertex_const_handle curr_v = vf_iter->first;
    auto obj = def_pl.ray_shoot_down(curr_v->point());

    // if (CGAL::assign(he_below, obj)) {
    if (auto* he_below_p = boost::get<Halfedge_const_handle>(&obj)) {
      // Hit an edge; take the incident face of the halfedge directed to the
      // right.
      auto he_below = *he_below_p;
      in_face = (he_below->direction() == ARR_RIGHT_TO_LEFT) ?
        he_below->twin()->face() : he_below->face();
    }
    else if (auto* v_below_p = boost::get<Vertex_const_handle>(&obj)) {
      auto v_below = *v_below_p;
      // Hit a vertex.
      if (v_below->is_isolated()) in_face = v_below->face();
      else {
        // Get the first halfedge around v_below that is directed from left to
        // right and the first halfedge that is directed from right to left.
        Halfedge_around_vertex_const_circulator circ =
          v_below->incident_halfedges();
        Halfedge_around_vertex_const_circulator first = circ;
        Halfedge_const_handle he_left;  // A halfedge to the left of v_below.
        Halfedge_const_handle he_right; // A halfedge to the right of v_below.
        do {
          if (circ->direction() == ARR_LEFT_TO_RIGHT) he_left = circ;
          else {
            he_right = circ;
            if ((he_left != invalid_he) && (he_right != invalid_he)) break;
          }
        } while(++circ != first);

        CGAL_assertion((he_left != invalid_he) || (he_right != invalid_he));

        if (he_left != invalid_he && he_right != invalid_he) {
          while (he_left->direction() == ARR_LEFT_TO_RIGHT)
            he_left = he_left->next()->twin();

          he_left = he_left->twin()->prev();
          CGAL_assertion (he_left->direction() == ARR_LEFT_TO_RIGHT);
          in_face = he_left->face();
        }
        else if (he_left != invalid_he) {
          Comparison_result res;
          Halfedge_const_handle he_curr = he_left;

          // as long as we have next he_left halfedge which is above
          do {
            he_left = he_curr;
            he_curr = he_left->next()->twin();
            res = comp_y_at_x_left(he_curr->curve(), he_left->curve(),
                                   v_below->point());
          } while(res == LARGER);
          in_face = he_left->face();
        }
        else {
          Comparison_result res;
          Halfedge_const_handle he_curr = he_right;
          do {
            // as long as we have he_right halfedge which is below
            he_right = he_curr;
            he_curr = he_right->next()->twin();
            res = comp_y_at_x_right(he_curr->curve(),
                                    he_right->curve(),
                                    v_below->point());
          } while(res == SMALLER);
          in_face = he_right->face();
        }
      }
    }
    else {
      auto* in_face_p = boost::get<Face_const_handle>(&obj);
      CGAL_assertion(in_face_p);
      in_face = *in_face_p;
      // Hit nothing (an unbounded face is returned).
      CGAL_assertion(in_face->is_unbounded());
    }

    if (vf_iter->second != in_face) {
      CGAL_warning_msg (false,
                        "An inner component is located in the wrong face.");
      return false;
    }
  }

  // If we reached here, the arrangement is valid:
  return true;
}

//-----------------------------------------------------------------------------
// Compute the zone of the given x-monotone curve in the existing arrangement.
// Meaning, it output the arrangment's vertices, edges and faces that the
// x-monotone curve intersects.
template <typename GeometryTraits_2, typename TopologyTraits,
  typename OutputIterator, typename PointLocation>
OutputIterator
zone(Arrangement_on_surface_2<GeometryTraits_2, TopologyTraits>& arr,
     const typename GeometryTraits_2::X_monotone_curve_2& c,
     OutputIterator oi,
     const PointLocation& pl)
{
  typedef GeometryTraits_2                              Gt2;
  typedef TopologyTraits                                Tt;

  // Obtain an arrangement accessor.
  typedef Arrangement_on_surface_2<Gt2, Tt>             Arr;

  // Define a zone-computation object an a visitor that performs the
  // intersection check.
  typedef Arr_compute_zone_visitor<Arr, OutputIterator> Zone_visitor;

  Zone_visitor visitor(oi);
  Arrangement_zone_2<Arr, Zone_visitor> arr_zone(arr, &visitor);

  arr_zone.init(c, pl);
  arr_zone.compute_zone();

  return oi;
}

//-----------------------------------------------------------------------------
// Compute the zone of the given x-monotone curve in the existing arrangement.b
// Overloaded version with no point location object - the walk point-location
// strategy is used as default.
//
template <typename GeometryTraits_2, typename TopologyTraits,
          typename OutputIterator>
OutputIterator
zone(Arrangement_on_surface_2<GeometryTraits_2, TopologyTraits>& arr,
     const typename GeometryTraits_2::X_monotone_curve_2& c,
     OutputIterator oi)
{
  typedef TopologyTraits                                Tt;

  // Create a default point-location object and use it to insert the curve.
  typename Tt::Default_point_location_strategy def_pl(arr);

  //insert the curve using the walk point location
  zone(arr, c, oi, def_pl);
  return oi;
}


//-----------------------------------------------------------------------------
// Checks whether the given x-monotone curve intersects the existing arrangement.
// The last parameter is used to resolve ambiguity between this function and
// do_intersect of Curve_2 in case that X_monotone_curve_2 and Curve_2 are the
// same class. The last parameter should be boost::true_type but we used a
// workaround since it didn't compile in FC3_g++-3.4.4 with the error of:
//
// error: no matching function for call to `do_intersect(Arrangement_on_surface_2<>&,
// const Arr_segment_2&, const Arr_walk_along_line_point_location<>&, mpl_::bool_< true>)'
//
template <typename GeometryTraits_2, typename TopologyTraits,
          typename PointLocation>
bool
do_intersect(Arrangement_on_surface_2<GeometryTraits_2, TopologyTraits>& arr,
             const typename GeometryTraits_2::X_monotone_curve_2& c,
             const PointLocation& pl, boost::is_same<int, int>::type)
{
  typedef GeometryTraits_2                              Gt2;
  typedef TopologyTraits                                Tt;

  // Obtain an arrangement accessor.
  typedef Arrangement_on_surface_2<Gt2, Tt>             Arr;

  // Define a zone-computation object an a visitor that performs the
  // intersection check.
  typedef Arr_do_intersect_zone_visitor<Arr>            Zone_visitor;

  Zone_visitor visitor;
  Arrangement_zone_2<Arr, Zone_visitor> arr_zone(arr, &visitor);

  arr_zone.init(c, pl);
  arr_zone.compute_zone();

  return visitor.do_intersect();
}

//-----------------------------------------------------------------------------
// Checks if the given curve intersects the existing arrangement.
// The last parameter is used to resolve ambiguity between this function and
// do_intersect of X_monotone_curve_2 in case that X_monotone_curve_2 and
// Curve_2 are the same class.
// The last parameter should be boost::false_type but we used a
// workaround since it didn't compile in FC3_g++-3.4.4 with the error of:
//
// error: no matching function for call to
// `do_intersect(Arrangement_on_surface_2<>&,
// const Arr_segment_2&, const Arr_walk_along_line_point_location<>&, mpl_::bool_< true>)'
//
template <typename GeometryTraits_2, typename TopologyTraits,
          typename PointLocation>
bool
do_intersect(Arrangement_on_surface_2<GeometryTraits_2, TopologyTraits>& arr,
             const typename GeometryTraits_2::X_monotone_curve_2& c,
             const PointLocation& pl, boost::is_same<int, double>::type)
{
  typedef GeometryTraits_2                              Gt2;
  typedef TopologyTraits                                Tt;

  // Obtain an arrangement accessor.
  typedef Arrangement_on_surface_2<Gt2, Tt>             Arr;

  // Break the input curve into x-monotone subcurves and isolated points.
  typedef Arr_traits_adaptor_2<Gt2>                     Traits_adaptor_2;

  typedef typename Gt2::Point_2                         Point_2;
  typedef typename Gt2::X_monotone_curve_2              X_monotone_curve_2;
  typedef boost::variant<Point_2, X_monotone_curve_2>   Make_x_monotone_result;
  typedef typename Arr::Face_const_handle               Face_const_handle;

  const Traits_adaptor_2* traits =
    static_cast<const Traits_adaptor_2*>(arr.geometry_traits());

  std::list<Make_x_monotone_result> x_objects;
  traits->make_x_monotone_2_object()(c, std::back_inserter(x_objects));

  // Insert each x-monotone curve into the arrangement.
  for (const auto& x_obj : x_objects) {
    // Act according to the type of the current object.
    const X_monotone_curve_2* x_curve = boost::get<X_monotone_curve_2>(&x_obj);
    if (x_curve != nullptr) {
      // Check if the x-monotone subcurve intersects the arrangement.
      if (do_intersect(arr, *x_curve, pl) == true) return true;
      continue;
    }

    const Point_2* iso_p = boost::get<Point_2>(&x_obj);
    CGAL_assertion(iso_p != nullptr);

    // Check whether the isolated point lies inside a face (otherwise,
    // it conincides with a vertex or an edge).
    auto obj = pl.locate(*iso_p);
    if (boost::get<Face_const_handle>(&x_obj) != nullptr) return true;
  }

  // If we reached here, the curve does not intersect the arrangement.
  return false;
}

//-----------------------------------------------------------------------------
// Common interface for the do_intersect of the Curve_2 and X_monotone_curve_2
template <typename GeometryTraits_2, typename TopologyTraits, typename Curve,
          typename PointLocation>
bool
do_intersect(Arrangement_on_surface_2<GeometryTraits_2, TopologyTraits>& arr,
             const Curve& c, const PointLocation& pl)
{
  typedef GeometryTraits_2                              Gt2;

  typedef typename Gt2::X_monotone_curve_2              X_monotone_curve_2;

  typedef typename boost::is_same<Curve, X_monotone_curve_2>::type
                                                        Is_x_monotone;

  return do_intersect(arr, c, pl, Is_x_monotone());
}

//-----------------------------------------------------------------------------
// Checks whether the given curve intersects the existing arrangement.
// Overloaded version with no point location object - the walk point-location
// strategy is used as default.
template <typename GeometryTraits_2, typename TopologyTraits, typename Curve>
bool
do_intersect(Arrangement_on_surface_2<GeometryTraits_2, TopologyTraits>& arr,
             const Curve& c)
{
  typedef TopologyTraits                                Tt;

  // Create a default point-location object and use it to insert the curve.
  typename Tt::Default_point_location_strategy def_pl(arr);

  // check whether the curve intersects the arrangement using the walk point
  // location.
  return do_intersect(arr, c, def_pl);
}

} // namespace CGAL

#endif
