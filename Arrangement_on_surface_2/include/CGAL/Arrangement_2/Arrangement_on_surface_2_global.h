// Copyright (c) 2006,2007,2008,2009,2010,2011 Tel-Aviv University (Israel).
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
//
// Author(s)     : Ron Wein          <wein@post.tau.ac.il>
//                 Baruch Zukerman   <baruchzu@post.tau.ac.il>
//                 Efi Fogel         <efif@post.tau.ac.il>
//
#ifndef CGAL_ARRANGEMENT_ON_SURFACE_2_GLOBAL_H
#define CGAL_ARRANGEMENT_ON_SURFACE_2_GLOBAL_H

#include <CGAL/license/Arrangement_on_surface_2.h>


/*! \file
 * Global insertion functions for the Arrangement_2 class.
 */

#include <boost/type_traits.hpp>
#include <boost/mpl/if.hpp>

#include <CGAL/Arrangement_on_surface_2.h>
#include <CGAL/Arr_accessor.h>
#include <CGAL/Basic_sweep_line_2.h>
#include <CGAL/Sweep_line_2.h>
#include <CGAL/Arrangement_zone_2.h>
#include <CGAL/Arrangement_2/Arr_compute_zone_visitor.h>
#include <CGAL/Arrangement_2/Arr_do_intersect_zone_visitor.h>
#include <CGAL/Arrangement_2/Arr_traits_adaptor_2.h>
#include <CGAL/Sweep_line_2/Sweep_line_2_utils.h>
#include <CGAL/Sweep_line_2/Sweep_line_2_visitors.h>
#include <CGAL/assertions.h>
#include <CGAL/use.h>

#include <CGAL/IO/Arr_iostream.h>

#include <boost/type_traits.hpp>

#include <list>

namespace CGAL {

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
template <class GeomTraits, class TopTraits, class PointLocation,
  class ZoneVisitor>
void insert (Arrangement_on_surface_2<GeomTraits, TopTraits>& arr,
             const typename GeomTraits::Curve_2& c,
             const PointLocation& pl, ZoneVisitor &visitor,
             boost::is_same<int, double>::type)
{
  typedef Arrangement_on_surface_2<GeomTraits, TopTraits>  Arr;
  typedef ZoneVisitor                                      Zone_visitor;

  // Obtain an arrangement accessor.
  Arr_accessor<Arr>                      arr_access (arr);

  // Initialize a zone-computation object an a visitor that performs the
  // incremental insertion.
  Arrangement_zone_2<Arr, Zone_visitor>  arr_zone (arr, &visitor);

  // Break the input curve into x-monotone subcurves and isolated points.
  std::list<CGAL::Object>                         x_objects;
  std::list<CGAL::Object>::const_iterator         obj_iter;
  const typename GeomTraits::X_monotone_curve_2  *x_curve;
  const typename GeomTraits::Point_2             *iso_p;

  arr.geometry_traits()->make_x_monotone_2_object()
    (c, std::back_inserter (x_objects));

  // Insert each x-monotone curve into the arrangement.
  for (obj_iter = x_objects.begin(); obj_iter != x_objects.end(); ++obj_iter)
  {
    // Act according to the type of the current object.
    x_curve =
      object_cast<typename GeomTraits::X_monotone_curve_2> (&(*obj_iter));

    if (x_curve != NULL)
    {
      // Inserting an x-monotone curve:
      // Initialize the zone-computation object with the given curve.
      arr_zone.init (*x_curve, pl);

      // Notify the arrangement observers that a global operation is about to
      // take place.
      arr_access.notify_before_global_change();

      // Insert the current x-monotone curve into the arrangement.
      arr_zone.compute_zone();

      // Notify the arrangement observers that the global operation has been
      // completed.
      arr_access.notify_after_global_change();
    }
    else
    {
      iso_p = object_cast<typename GeomTraits::Point_2> (&(*obj_iter));
      CGAL_assertion (iso_p != NULL);

      // Inserting a point into the arrangement:
      insert_point (arr, *iso_p, pl);
    }
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
template <class GeomTraits, class TopTraits, class PointLocation,
  class ZoneVisitor>
void insert(Arrangement_on_surface_2<GeomTraits, TopTraits>& arr,
            const typename GeomTraits::X_monotone_curve_2& c,
            const PointLocation& pl, ZoneVisitor &visitor,
            boost::is_same<int, int>::type)
{
  typedef Arrangement_on_surface_2<GeomTraits, TopTraits>  Arr;
  typedef ZoneVisitor                                      Zone_visitor;

  // Obtain an arrangement accessor.
  Arr_accessor<Arr>                      arr_access (arr);

  // Initialize a zone-computation object an a visitor that performs the
  // incremental insertion.
  Arrangement_zone_2<Arr, Zone_visitor>  arr_zone (arr, &visitor);

  // Initialize the zone-computation object with the given curve.

/*  Needs to be deleted!!!
    {
    std::cout << std::endl;
    std::cout << "xxxxxxxxxxxx" << std::endl;
    std::cout << "c: " << c << std::endl;
    std::cout << std::endl;
    typename Arr::Vertex_const_iterator  vit;
    std::cout << arr.number_of_vertices() << " vertices:" << std::endl;
    for (vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit)
    {
      std::cout << "(" << vit->point() << ")";
      if (vit->is_isolated())
        std::cout << " - Isolated." << std::endl;
      else
        std::cout << " - degree " << vit->degree() << std::endl;
    }
    std::cout << std::endl;
    std::cout << "xxxxxxxxxxxx" << std::endl;
  }
*/
  arr_zone.init (c, pl);

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
template <class GeomTraits, class TopTraits, class Curve, class PointLocation,
  class ZoneVisitor>
void insert (Arrangement_on_surface_2<GeomTraits, TopTraits>& arr,
             const Curve& c, const PointLocation& pl, ZoneVisitor &visitor)
{
  typedef typename GeomTraits::X_monotone_curve_2       X_monotone_curve_2;

  typedef typename boost::is_same<Curve, X_monotone_curve_2>::type
    Is_x_monotone;

  insert(arr, c, pl, visitor, Is_x_monotone());
}

// In some compilers there is a template deduction disambiguity between this
// function and the function receiving two InputIterator.
// For now the solution is to add a dummy variable at the end (referring
// to point-location). Maybe the proper solution is to use boost::enable_if
// together with appropriate tag.
template <class GeomTraits, class TopTraits, class Curve, class PointLocation>
void insert (Arrangement_on_surface_2<GeomTraits, TopTraits>& arr,
             const Curve& c, const PointLocation& pl,
	     typename PointLocation::Point_2*)
{
  typedef typename TopTraits::Zone_insertion_visitor       Zone_visitor;

  Zone_visitor visitor;
  insert (arr, c, pl, visitor);
}

//-----------------------------------------------------------------------------
// Insert a curve/x-monotone curve into the arrangement (incremental
// insertion).
// The inserted x-monotone curve may intersect the existing arrangement.
// Overloaded version with no point location object - using the default point
// location.
//
template <class GeomTraits, class TopTraits, class Curve>
void insert(Arrangement_on_surface_2<GeomTraits, TopTraits>& arr,
            const Curve& c)
{
  // Create a default point-location object and use it to insert the curve.
  typename TopTraits::Default_point_location_strategy    def_pl (arr);

  insert (arr, c, def_pl);
}

/*! Insert a range of x-monotone curves into an empty arrangement
 * \param arr the resulting arrangement
 * \param begin the begining of the curve range
 * \param end past-the-end curve range
 */
template <typename GeomTraits, typename TopTraits, typename InputIterator>
void insert_empty(Arrangement_on_surface_2<GeomTraits, TopTraits>& arr,
                  InputIterator begin_xcurves, InputIterator end_xcurves)
{
  typedef typename TopTraits::Sweep_line_construction_visitor
                                                        Construct_visitor;
  typedef typename Construct_visitor::Traits_2          Construct_traits;

  const GeomTraits * geom_traits = arr.geometry_traits();
  Construct_visitor visitor(&arr);

  /* We would like to avoid copy construction of the geometry traits class.
   * Copy construction is undesired, because it may results with data
   * duplication or even data loss.
   *
   * If the type Construct_visitor::Traits_2 is the same as the type
   * GeomTraits, use a reference to GeomTraits to avoid constructing a new one.
   * Otherwise, instantiate a local variable of the former and provide
   * the later as a single parameter to the constructor.
   *
   * Use the form 'A a(*b);' and not ''A a = b;' to handle the case where A has
   * only an implicit constructor, (which takes *b as a parameter).
   */
  typename boost::mpl::if_<boost::is_same<GeomTraits, Construct_traits>,
                           const Construct_traits&, Construct_traits>::type
    traits(*geom_traits);

  // Define a sweep-line instance and perform the sweep:
  Sweep_line_2<typename Construct_visitor::Traits_2, Construct_visitor,
               typename Construct_visitor::Subcurve,
               typename Construct_visitor::Event>
    sweep_line(&traits, &visitor);
  sweep_line.sweep(begin_xcurves, end_xcurves);
}

/*! Insert a range of x-monotone curves and a range of isolated points into
 * an empty arrangement
 * \param arr the resulting arrangement
 * \param begin_xcurves the begining of the curve range
 * \param end_xcurves past-the-end curve range
 * \param begin_points the begining of the point range
 * \param end_points past-the-end point range
 */
template <typename GeomTraits, typename TopTraits,
          typename XcInputIterator, typename PInputIterator>
void insert_empty(Arrangement_on_surface_2<GeomTraits, TopTraits>& arr,
                  XcInputIterator begin_xcurves, XcInputIterator end_xcurves,
                  PInputIterator begin_points, PInputIterator end_points)
{
  typedef typename TopTraits::Sweep_line_construction_visitor
                                                        Construct_visitor;
  typedef typename Construct_visitor::Traits_2          Construct_traits;

  const GeomTraits * geom_traits = arr.geometry_traits();
  Construct_visitor visitor(&arr);

  /* We would like to avoid copy construction of the geometry traits class.
   * Copy construction is undesired, because it may results with data
   * duplication or even data loss.
   *
   * If the type Construct_visitor::Traits_2 is the same as the type
   * GeomTraits, use a reference to GeomTraits to avoid constructing a new one.
   * Otherwise, instantiate a local variable of the former and provide
   * the later as a single parameter to the constructor.
   *
   * Use the form 'A a(*b);' and not ''A a = b;' to handle the case where A has
   * only an implicit constructor, (which takes *b as a parameter).
   */
  typename boost::mpl::if_<boost::is_same<GeomTraits, Construct_traits>,
                           const Construct_traits&, Construct_traits>::type
    traits(*geom_traits);

  // Define a sweep-line instance and perform the sweep.
  Sweep_line_2<typename Construct_visitor::Traits_2, Construct_visitor,
               typename Construct_visitor::Subcurve,
               typename Construct_visitor::Event>
    sweep_line(&traits, &visitor);
  sweep_line.sweep(begin_xcurves, end_xcurves, begin_points, end_points);
}

/*! Insert a range of x-monotone curves into a non-empty arrangement
 * \param arr the resulting arrangement
 * \param begin the begining of the curve range
 * \param end past-the-end curve range
 */
template <typename GeomTraits, typename TopTraits,
          typename XcInputIterator, typename PInputIterator>
void insert_non_empty(Arrangement_on_surface_2<GeomTraits, TopTraits>& arr,
                      XcInputIterator begin_xcurves,
                      XcInputIterator end_xcurves,
                      PInputIterator begin_points, PInputIterator end_points)
{
  typedef typename TopTraits::Sweep_line_insertion_visitor
                                                        Insert_visitor;
  typedef typename Insert_visitor::Traits_2             Insert_traits;
  typedef typename Insert_visitor::Traits_2::X_monotone_curve_2
                                                        Ex_x_monotone_curve_2;
  typedef typename Insert_visitor::Traits_2::Point_2    Ex_point_2;

  const GeomTraits * geom_traits = arr.geometry_traits();
  Insert_visitor visitor(&arr);

  /* We would like to avoid copy construction of the geometry traits class.
   * Copy construction is undesired, because it may results with data
   * duplication or even data loss.
   *
   * If the type Construct_visitor::Traits_2 is the same as the type
   * GeomTraits, use a reference to GeomTraits to avoid constructing a new one.
   * Otherwise, instantiate a local variable of the former and provide
   * the later as a single parameter to the constructor.
   *
   * Use the form 'A a(*b);' and not ''A a = b;' to handle the case where A has
   * only an implicit constructor, (which takes *b as a parameter).
   */
  typename boost::mpl::if_<boost::is_same<GeomTraits, Insert_traits>,
                           const Insert_traits&, Insert_traits>::type
    traits(*geom_traits);

  // Create a set of existing as well as new curves and points.
  std::list<Ex_x_monotone_curve_2> ex_cvs;
  std::list<Ex_point_2> ex_pts;

  prepare_for_sweep(arr,
                    begin_xcurves, end_xcurves,   // the x-monotone curves
                    begin_points, end_points,     // the points (if any)
                    std::back_inserter(ex_cvs),
                    std::back_inserter(ex_pts),
                    &traits);

  // Define a basic sweep-line instance and perform the sweep.
  Sweep_line_2<typename Insert_visitor::Traits_2, Insert_visitor,
               typename Insert_visitor::Subcurve,
               typename Insert_visitor::Event>
    sweep_line(&traits, &visitor);
  sweep_line.sweep(ex_cvs.begin(), ex_cvs.end(),ex_pts.begin(), ex_pts.end());
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
template <class GeomTraits, class TopTraits, class InputIterator>
void insert (Arrangement_on_surface_2<GeomTraits, TopTraits>& arr,
             InputIterator begin, InputIterator end,
             boost::is_same<int, double>::type)
{
  typedef Arrangement_on_surface_2<GeomTraits, TopTraits>   Arr;
  typedef typename GeomTraits::Point_2                      Point_2;
  typedef typename GeomTraits::X_monotone_curve_2           X_monotone_curve_2;

  // Obtain an arrangement accessor.
  Arr_accessor<Arr>                  arr_access (arr);

  // Notify the arrangement observers that a global operation is about to
  // take place.
  arr_access.notify_before_global_change();

  // Subdivide the input curves into x-monotone subcurves and isolated points.
  const GeomTraits * geom_traits = arr.geometry_traits();
  std::list<X_monotone_curve_2>      xcurves;
  std::list<Point_2>                 iso_points;

  make_x_monotone(begin, end,
                  std::back_inserter(xcurves), std::back_inserter(iso_points),
                  geom_traits);

  if (arr.is_empty())
    insert_empty(arr, xcurves.begin(), xcurves.end(),
                 iso_points.begin(), iso_points.end());
  else
    insert_non_empty(arr, xcurves.begin(), xcurves.end(),
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
template <class GeomTraits, class TopTraits, class InputIterator>
void insert(Arrangement_on_surface_2<GeomTraits, TopTraits>& arr,
            InputIterator begin, InputIterator end,
            boost::is_same<int, int>::type)
{
  typedef Arrangement_on_surface_2<GeomTraits, TopTraits>   Arr;

  // Obtain an arrangement accessor.
  Arr_accessor<Arr>                  arr_access (arr);

  // Notify the arrangement observers that a global operation is about to
  // take place.
  arr_access.notify_before_global_change();

  // Choose the operation depending on whether the input arrangement is
  // empty (then we construct it from scratch), or not (where we just insert
  // the new curves).
  if (arr.is_empty())
    insert_empty(arr, begin, end);
  else {
    // The arrangement is not empty: use the insertion visitor.
    std::list<typename GeomTraits::Point_2> empty;
    insert_non_empty(arr, begin, end, empty.begin(), empty.end());
  }

  // Notify the arrangement observers that the global operation has been
  // completed.
  arr_access.notify_after_global_change();
}

//-----------------------------------------------------------------------------
// Common interface for the inserts of the Curve_2 and X_monotone_curve_2
template <class GeomTraits, class TopTraits, class InputIterator>
void insert (Arrangement_on_surface_2<GeomTraits,TopTraits>& arr,
             InputIterator begin, InputIterator end)
{
  typedef typename GeomTraits::X_monotone_curve_2       X_monotone_curve_2;
  typedef typename std::iterator_traits<InputIterator>::value_type
    Iterator_value_type;

  typedef typename boost::is_same<Iterator_value_type,X_monotone_curve_2>::type
    Is_x_monotone;

  return insert (arr, begin, end, Is_x_monotone());
}

//-----------------------------------------------------------------------------
// Insert an x-monotone curve into the arrangement (incremental insertion)
// when the location of the left endpoint of the curve is known and is
// given as an isertion hint.
// The inserted x-monotone curve may intersect the existing arrangement.
//
template <class GeomTraits, class TopTraits>
void insert (Arrangement_on_surface_2<GeomTraits, TopTraits>& arr,
             const typename GeomTraits::X_monotone_curve_2& c,
             const Object& obj)
{
  typedef Arrangement_on_surface_2<GeomTraits, TopTraits>  Arr;
  typedef typename TopTraits::Zone_insertion_visitor       Zone_visitor;

  // Obtain an arrangement accessor.
  Arr_accessor<Arr>                      arr_access (arr);

  // Define a zone-computation object an a visitor that performs the
  // incremental insertion.
  Zone_visitor                           visitor;
  Arrangement_zone_2<Arr, Zone_visitor>  arr_zone (arr, &visitor);

  // Initialize the zone-computation object with the given curve.
  arr_zone.init_with_hint (c, obj);

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
template <class GeomTraits, class TopTraits, class PointLocation>
CGAL_DEPRECATED void insert_x_monotone_curve
(Arrangement_on_surface_2<GeomTraits, TopTraits>& arr,
 const typename GeomTraits::X_monotone_curve_2& c,
 const PointLocation& pl)
{
  insert(arr, c, pl);
}

/* DEPRECATED use insert() instead */
template <class GeomTraits, class TopTraits>
CGAL_DEPRECATED void insert_x_monotone_curve
(Arrangement_on_surface_2<GeomTraits, TopTraits>& arr,
 const typename GeomTraits::X_monotone_curve_2& c)
{
  insert(arr, c);
}

/* DEPRECATED use insert() instead */
template <class GeomTraits, class TopTraits, class InputIterator>
CGAL_DEPRECATED void insert_x_monotone_curves
(Arrangement_on_surface_2<GeomTraits, TopTraits>& arr,
 InputIterator begin, InputIterator end)
{
  insert(arr, begin, end);
}

/* DEPRECATED use insert() instead */
template <class GeomTraits, class TopTraits>
CGAL_DEPRECATED void insert_x_monotone_curve
(Arrangement_on_surface_2<GeomTraits, TopTraits>& arr,
 const typename GeomTraits::X_monotone_curve_2& c,
 const Object& obj)
{
  insert(arr, c, obj);
}

/* DEPRECATED use insert() instead */
template <class GeomTraits, class TopTraits, class PointLocation>
CGAL_DEPRECATED
void insert_curve(Arrangement_on_surface_2<GeomTraits, TopTraits>& arr,
                  const typename GeomTraits::Curve_2& c,
                  const PointLocation& pl)
{
  insert(arr, c, pl);
}

/* DEPRECATED use insert() instead */
template <class GeomTraits, class TopTraits>
CGAL_DEPRECATED
void insert_curve(Arrangement_on_surface_2<GeomTraits, TopTraits>& arr,
                  const typename GeomTraits::Curve_2& c)
{
  insert(arr, c);
}

/* DEPRECATED use insert() instead */
template <class GeomTraits, class TopTraits, class InputIterator>
CGAL_DEPRECATED
void insert_curves (Arrangement_on_surface_2<GeomTraits, TopTraits>& arr,
                    InputIterator begin, InputIterator end)
{
  insert(arr, begin, end);
}

//-----------------------------------------------------------------------------
// Insert an x-monotone curve into the arrangement, such that the curve
// interior does not intersect with any existing edge or vertex in the
// arragement (incremental insertion).
//
template <class GeomTraits, class TopTraits, class PointLocation>
typename Arrangement_on_surface_2<GeomTraits, TopTraits>::Halfedge_handle
insert_non_intersecting_curve
(Arrangement_on_surface_2<GeomTraits, TopTraits>& arr,
 const typename GeomTraits::X_monotone_curve_2& c,
 const PointLocation& pl)
{
  typedef Arrangement_on_surface_2<GeomTraits, TopTraits>  Arr;

  typedef Arr_traits_basic_adaptor_2<typename Arr::Geometry_traits_2>
                                                         Traits_adaptor_2;
  typedef typename Arr::Vertex_const_handle              Vertex_const_handle;
  typedef typename Arr::Halfedge_const_handle            Halfedge_const_handle;
  CGAL_USE_TYPE(Halfedge_const_handle);

  const Traits_adaptor_2* geom_traits =
    static_cast<const Traits_adaptor_2*> (arr.geometry_traits());
  Arr_accessor<Arr> arr_access(arr);

  // Check whether the left end has boundary conditions, and locate it in the
  // arrangement accordingly.
  const Arr_parameter_space  bx1 =
    geom_traits->parameter_space_in_x_2_object()(c, ARR_MIN_END);
  const Arr_parameter_space  by1 =
    geom_traits->parameter_space_in_y_2_object()(c, ARR_MIN_END);
  CGAL::Object obj1;
  const Vertex_const_handle* vh1 = NULL;

  if ((bx1 == ARR_INTERIOR) && (by1 == ARR_INTERIOR)) {
    // We have a normal left endpoint with no boundary conditions:
    // use a point-location query.
    obj1 = pl.locate(geom_traits->construct_min_vertex_2_object()(c));

    // The endpoint must not lie on an existing edge, but may coincide with
    // and existing vertex vh1.
    CGAL_precondition_msg
      (object_cast<Halfedge_const_handle>(&obj1) == NULL,
       "The curve must not intersect an existing edge.");

    vh1 = object_cast<Vertex_const_handle>(&obj1);
  }
  else {
    // We have a left end with boundary conditions. Use the accessor to locate
    // the feature that contains it.
    obj1 = arr_access.locate_curve_end(c, ARR_MIN_END, bx1, by1);

    CGAL_precondition_msg
      (object_cast<Halfedge_const_handle>(&obj1) == NULL,
       "The curve must not overlap an existing edge.");

    vh1 = object_cast<Vertex_const_handle>(&obj1);
  }

  // Check whether the right end has boundary conditions, and locate it in the
  // arrangement accordingly.
  const Arr_parameter_space  bx2 =
    geom_traits->parameter_space_in_x_2_object()(c, ARR_MAX_END);
  const Arr_parameter_space  by2 =
    geom_traits->parameter_space_in_y_2_object()(c, ARR_MAX_END);
  CGAL::Object obj2;
  const Vertex_const_handle* vh2 = NULL;

  if ((bx2 == ARR_INTERIOR) && (by2 == ARR_INTERIOR)) {
    // We have a normal right endpoint with no boundary conditions:
    // use a point-location query.
    obj2 = pl.locate(geom_traits->construct_max_vertex_2_object()(c));

    // The endpoint must not lie on an existing edge, but may coincide with
    // and existing vertex vh2.
    CGAL_precondition_msg
      (object_cast<Halfedge_const_handle>(&obj2) == NULL,
       "The curve must not intersect an existing edge.");

    vh2 = object_cast<Vertex_const_handle>(&obj2);
  }
  else {
    // We have a right end with boundary conditions. Use the accessor to locate
    // the feature that contains it.
    // std::cout << "before locate_curve_end()"
    //           << ", bx2: " << bx2
    //           << ", by2: " << by2
    //           << std::endl;
    obj2 = arr_access.locate_curve_end(c, ARR_MAX_END, bx2, by2);

    CGAL_precondition_msg
      (object_cast<Halfedge_const_handle>(&obj2) == NULL,
       "The curve must not overlap an existing edge.");

    vh2 = object_cast<Vertex_const_handle>(&obj2);
  }

  // Notify the arrangement observers that a global operation is about to
  // take place.
  arr_access.notify_before_global_change();

  // Check whether the located features containing the curve endpoints
  // are vertices or faces, and use the proper specialized insertion function
  // accordingly.
  typename Arr::Halfedge_handle new_he;

  if (vh1 != NULL) {
    if (vh2 != NULL) {
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
    if (vh2 != NULL) {
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
      const typename Arr::Face_const_handle* fh1 =
        object_cast<typename Arr::Face_const_handle>(&obj1);
      const typename Arr::Face_const_handle* fh2 =
        object_cast<typename Arr::Face_const_handle>(&obj2);

      // std::cout << arr << std::endl;
      // std::cout << "(*fh1)->number_of_outer_ccbs(): "
      //           << (*fh1)->number_of_outer_ccbs() << std::endl;
      // std::cout << "(*fh2)->number_of_outer_ccbs(): "
      //           << (*fh2)->number_of_outer_ccbs() << std::endl;

      CGAL_assertion_msg
        ((fh1 != NULL) && (fh2 != NULL) && ((*fh1) == (*fh2)),
         "The curve intersects the interior of existing edges.");

      if ((fh1 != NULL) && (fh2 != NULL) && (*fh1 == *fh2)) {
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
template <class GeomTraits, class TopTraits>
typename Arrangement_on_surface_2<GeomTraits, TopTraits>::Halfedge_handle
insert_non_intersecting_curve
    (Arrangement_on_surface_2<GeomTraits, TopTraits>& arr,
     const typename GeomTraits::X_monotone_curve_2& c)
{
  // Create a default point-location object and use it to insert the curve.
  typename TopTraits::Default_point_location_strategy  def_pl (arr);

  return (insert_non_intersecting_curve (arr, c, def_pl));
}

/*! Insert a range of x-monotone curves into an empty arrangement
 * \param arr the resulting arrangement
 * \param begin the begining of the curve range
 * \param end past-the-end curve range
 */
template <typename GeomTraits, typename TopTraits, typename InputIterator>
void non_intersecting_insert_empty(Arrangement_on_surface_2<GeomTraits,
                                                            TopTraits>& arr,
                                   InputIterator begin_xcurves,
                                   InputIterator end_xcurves)
{
  const GeomTraits * geom_traits = arr.geometry_traits();
  typedef typename TopTraits::Sweep_line_non_intersecting_construction_visitor
                                                            Construct_visitor;
  Construct_visitor visitor(&arr);
  typename Construct_visitor::Traits_2 traits(*geom_traits);

  // Define a basic sweep-line instance (which is not supposed to handle
  // insersections) and perform the sweep.
  Basic_sweep_line_2<typename Construct_visitor::Traits_2, Construct_visitor,
                     typename Construct_visitor::Subcurve,
                     typename Construct_visitor::Event>
    sweep_line(&traits, &visitor);
  sweep_line.sweep(begin_xcurves, end_xcurves);
}

/*! Insert a range of x-monotone curves into an empty arrangement
 * \param arr the resulting arrangement
 * \param begin the begining of the curve range
 * \param end past-the-end curve range
 */
template <typename GeomTraits, typename TopTraits,
          typename XcInputIterator, typename PInputIterator>
void non_intersecting_insert_empty(Arrangement_on_surface_2<GeomTraits,
                                                            TopTraits>& arr,
                                   XcInputIterator begin_xcurves,
                                   XcInputIterator end_xcurves,
                                   PInputIterator begin_points,
                                   PInputIterator end_points)
{
  typedef typename TopTraits::Sweep_line_non_intersecting_construction_visitor
                                                        Construct_visitor;
  typedef typename Construct_visitor::Traits_2          Construct_traits;

  const GeomTraits * geom_traits = arr.geometry_traits();
  Construct_visitor visitor(&arr);

  /* We would like to avoid copy construction of the geometry traits class.
   * Copy construction is undesired, because it may results with data
   * duplication or even data loss.
   *
   * If the type Construct_visitor::Traits_2 is the same as the type
   * GeomTraits, use a reference to GeomTraits to avoid constructing a new one.
   * Otherwise, instantiate a local variable of the former and provide
   * the later as a single parameter to the constructor.
   *
   * Use the form 'A a(*b);' and not ''A a = b;' to handle the case where A has
   * only an implicit constructor, (which takes *b as a parameter).
   */
  typename boost::mpl::if_<boost::is_same<GeomTraits, Construct_traits>,
                           const Construct_traits&, Construct_traits>::type
    traits(*geom_traits);

  // Define a basic sweep-line instance (which is not supposed to handle
  // insersections) and perform the sweep.
  Basic_sweep_line_2<typename Construct_visitor::Traits_2, Construct_visitor,
                     typename Construct_visitor::Subcurve,
                     typename Construct_visitor::Event>
    sweep_line(&traits, &visitor);
  sweep_line.sweep(begin_xcurves, end_xcurves, begin_points, end_points);
}

/*! Insert a range of x-monotone curves into a non-empty arrangement
 * \param arr the resulting arrangement
 * \param begin the begining of the curve range
 * \param end past-the-end curve range
 */
template <typename GeomTraits, typename TopTraits,
          typename XcInputIterator, typename PInputIterator>
void non_intersecting_insert_non_empty(Arrangement_on_surface_2<GeomTraits,
                                                                TopTraits>& arr,
                                       XcInputIterator begin_xcurves,
                                       XcInputIterator end_xcurves,
                                       PInputIterator begin_points,
                                       PInputIterator end_points)
{
  typedef typename TopTraits::Sweep_line_non_intersecting_insertion_visitor
                                                        Insert_visitor;
  typedef typename Insert_visitor::Traits_2             Insert_traits;
  typedef typename Insert_visitor::Traits_2::X_monotone_curve_2
                                                        Ex_x_monotone_curve_2;
  typedef typename Insert_visitor::Traits_2::Point_2    Ex_point_2;

  const GeomTraits * geom_traits = arr.geometry_traits();
  Insert_visitor visitor(&arr);

  /* We would like to avoid copy construction of the geometry traits class.
   * Copy construction is undesired, because it may results with data
   * duplication or even data loss.
   *
   * If the type Construct_visitor::Traits_2 is the same as the type
   * GeomTraits, use a reference to GeomTraits to avoid constructing a new one.
   * Otherwise, instantiate a local variable of the former and provide
   * the later as a single parameter to the constructor.
   *
   * Use the form 'A a(*b);' and not ''A a = b;' to handle the case where A has
   * only an implicit constructor, (which takes *b as a parameter).
   */
  typename boost::mpl::if_<boost::is_same<GeomTraits, Insert_traits>,
                           const Insert_traits&, Insert_traits>::type
    traits(*geom_traits);

  // Create a set of existing as well as new curves and points.
  std::list<Ex_x_monotone_curve_2> ex_cvs;
  std::list<Ex_point_2> ex_pts;

  prepare_for_sweep(arr,
                    begin_xcurves, end_xcurves,   // the x-monotone curves
                    begin_points, end_points,     // the points (if any)
                    std::back_inserter(ex_cvs),
                    std::back_inserter(ex_pts),
                    &traits);

  // Define a basic sweep-line instance and perform the sweep.
  Basic_sweep_line_2<typename Insert_visitor::Traits_2, Insert_visitor,
                     typename Insert_visitor::Subcurve,
                       typename Insert_visitor::Event>
    sweep_line(&traits, &visitor);
  sweep_line.sweep(ex_cvs.begin(), ex_cvs.end(), ex_pts.begin(), ex_pts.end());
}

//-----------------------------------------------------------------------------
// Insert a range of pairwise interior-disjoint x-monotone curves into
// the arrangement, such that the curve interiors do not intersect with
// any existing edge or vertex in the arragement (aggregated insertion).
//
template <class GeomTraits, class TopTraits, class InputIterator>
void insert_non_intersecting_curves
  (Arrangement_on_surface_2<GeomTraits, TopTraits>& arr,
   InputIterator begin, InputIterator end)
{
  typedef Arrangement_on_surface_2<GeomTraits, TopTraits>   Arr;

  // Obtain an arrangement accessor.
  Arr_accessor<Arr>                  arr_access (arr);

  // Notify the arrangement observers that a global operation is about to
  // take place.
  arr_access.notify_before_global_change();

  // Choose the operation depending on whether the input arrangement is
  // empty (then we construct it from scratch), or not (where we just insert
  // the new curves).
  if (arr.is_empty())
    non_intersecting_insert_empty(arr, begin, end);
  else {
    std::list<typename GeomTraits::Point_2> empty;
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
template <class GeomTraits, class TopTraits>
typename Arrangement_on_surface_2<GeomTraits, TopTraits>::Face_handle
remove_edge
(Arrangement_on_surface_2<GeomTraits, TopTraits>& arr,
 typename Arrangement_on_surface_2<GeomTraits, TopTraits>::Halfedge_handle e)
{
  typedef Arrangement_on_surface_2<GeomTraits, TopTraits>   Arr;
  typedef Arr_traits_adaptor_2<GeomTraits>                  Traits_adaptor_2;

  // Notify the arrangement observers that a global operation is about to
  // take place.
  Arr_accessor<Arr>    arr_access (arr);

  arr_access.notify_before_global_change();

  // Keep track of the end-vertices of the edge we are about to remove.
  typename Arr::Vertex_handle  v_ends[2];
  bool                         is_removed[2];

  v_ends[0] = e->source();
  is_removed[0] =
      (v_ends[0]->is_at_open_boundary() || v_ends[0]->degree() == 1);
  v_ends[1] = e->target();
  is_removed[1] =
      (v_ends[1]->is_at_open_boundary() || v_ends[1]->degree() == 1);

  // Remove the edge from the arrangement.
  typename Arr::Face_handle    face = arr.remove_edge (e);

  // Examine the end-vertices: If a vertex has now two incident edges, and the
  // curves associated with these edges can be merged, merge the two edges and
  // remove the vertex.

  const Traits_adaptor_2 * traits =
    static_cast<const Traits_adaptor_2*> (arr.geometry_traits());

  typename Arr::Halfedge_around_vertex_circulator  circ;
  typename Arr::Halfedge_handle                    e1, e2;
  int                                              i;

  for (i = 0; i < 2; i++)
  {
    if (! is_removed[i] && v_ends[i]->degree() == 2)
    {
      // Get the two edges incident to the end-vertex.
      circ = v_ends[i]->incident_halfedges();
      e1 = circ;
      ++circ;
      e2 = circ;

      // Check if it is possible to merge the two edges.
      if (traits->are_mergeable_2_object() (e1->curve(), e2->curve()))
      {
        // Merge the two curves.
        typename GeomTraits::X_monotone_curve_2   cv;
        traits->merge_2_object() (e1->curve(), e2->curve(),
                                  cv);

        // Merge the two edges in the arrangement.
         arr.merge_edge (e1, e2, cv);
      }
    }
  }

  // Notify the arrangement observers that the global operation has been
  // completed.
  arr_access.notify_after_global_change();

  // Return the face remaining after the removal of the edge.
  return (face);
}

//-----------------------------------------------------------------------------
// Insert a vertex that corresponds to a given point into the arrangement.
// The inserted point may lie on any existing arrangement feature.
//
template <class GeomTraits, class TopTraits, class PointLocation>
typename Arrangement_on_surface_2<GeomTraits, TopTraits>::Vertex_handle
insert_point (Arrangement_on_surface_2<GeomTraits, TopTraits>& arr,
              const typename GeomTraits::Point_2& p,
              const PointLocation& pl)
{
  typedef Arrangement_on_surface_2<GeomTraits, TopTraits>   Arr;

  // Act according to the type of arrangement feature that contains the point.
  const typename Arr::Face_const_handle      *fh;
  const typename Arr::Halfedge_const_handle  *hh;
  const typename Arr::Vertex_const_handle    *vh;
  typename Arr::Vertex_handle                 vh_for_p;

  // Locate the given point in the arrangement.
  CGAL::Object        obj = pl.locate (p);

  // Notify the arrangement observers that a global operation is about to
  // take place.
  Arr_accessor<Arr>    arr_access (arr);

  arr_access.notify_before_global_change();

  if ((fh = object_cast<typename Arr::Face_const_handle>(&obj)) != NULL)
  {
    // p lies inside a face: Insert it as an isolated vertex it the interior of
    // this face.
    vh_for_p = arr.insert_in_face_interior (p,
                                            arr.non_const_handle (*fh));
  }
  else if ((hh =
            object_cast<typename Arr::Halfedge_const_handle>(&obj)) != NULL)
  {
    // p lies in the interior of an edge: Split this edge to create a new
    // vertex associated with p.
    typename GeomTraits::X_monotone_curve_2   sub_cv1, sub_cv2;
    typename Arr::Halfedge_handle             split_he;

    arr.geometry_traits()->split_2_object() ((*hh)->curve(), p,
                                             sub_cv1, sub_cv2);

    split_he = arr.split_edge (arr.non_const_handle (*hh),
                               sub_cv1, sub_cv2);

    // The new vertex is the target of the returned halfedge.
    vh_for_p = split_he->target();
  }
  else
  {
    // In this case p lies on an existing vertex, so we just update this
    // vertex.
    vh = object_cast<typename Arr::Vertex_const_handle>(&obj);
    CGAL_assertion (vh != NULL);

    vh_for_p = arr.modify_vertex (arr.non_const_handle (*vh), p);
  }

  // Notify the arrangement observers that the global operation has been
  // completed.
  arr_access.notify_after_global_change();

  // Return a handle for the vertex associated with p.
  return (vh_for_p);
}

//-----------------------------------------------------------------------------
// Insert a vertex that corresponds to a given point into the arrangement.
// The inserted point may lie on any existing arrangement feature.
//
template <class GeomTraits, class TopTraits>
typename Arrangement_on_surface_2<GeomTraits, TopTraits>::Vertex_handle
insert_point (Arrangement_on_surface_2<GeomTraits, TopTraits>& arr,
              const typename GeomTraits::Point_2& p)
{
  // Create a default point-location object and use it to insert the point.
  typename TopTraits::Default_point_location_strategy    def_pl (arr);

  return (insert_point (arr, p, def_pl));
}

//-----------------------------------------------------------------------------
// Remove a vertex from the arrangement.
//
template <class GeomTraits, class TopTraits>
bool remove_vertex
    (Arrangement_on_surface_2<GeomTraits, TopTraits>& arr,
     typename Arrangement_on_surface_2<GeomTraits, TopTraits>::Vertex_handle v)
{
  typedef Arrangement_on_surface_2<GeomTraits, TopTraits>   Arr;
  typedef Arr_traits_adaptor_2<GeomTraits>                  Traits_adaptor_2;

  // Notify the arrangement observers that a global operation is about to
  // take place.
  Arr_accessor<Arr>    arr_access (arr);

  arr_access.notify_before_global_change();

  // Act according to the number of edges incident to v.
  bool      removed = false;

  if (v->is_isolated())
  {
    // In case v is an isolated vertex, simply remove it.
    arr.remove_isolated_vertex (v);
    removed = true;
  }
  else if (v->degree() == 2)
  {
    // If the vertex has now two incident edges, and the curves associated
    // with these edges can be merged, merge the two edges and remove the
    // vertex.
    const Traits_adaptor_2 * traits =
      static_cast<const Traits_adaptor_2*>(arr.geometry_traits());
    typename Arr::Halfedge_around_vertex_circulator  circ;
    typename Arr::Halfedge_handle                    e1, e2;

    circ = v->incident_halfedges();
    e1 = circ;
    ++circ;
    e2 = circ;

    if (traits->are_mergeable_2_object() (e1->curve(), e2->curve()))
    {
      // Merge the two curves.
      typename GeomTraits::X_monotone_curve_2   cv;
      traits->merge_2_object() (e1->curve(), e2->curve(),
                                cv);

      // Merge the two edges in the arrangement.
      arr.merge_edge (e1, e2, cv);
      removed = true;
    }
  }

  // Notify the arrangement observers that the global operation has been
  // completed.
  arr_access.notify_after_global_change();

  // Return an indication whether the vertex has been removed or not.
  return (removed);
}

//-----------------------------------------------------------------------------
// Check the validity of the arrangement. In particular, check that the
// edegs are disjoint-interior, and the holes are located in their proper
// position.
//
template <class GeomTraits, class TopTraits>
bool is_valid (const Arrangement_on_surface_2<GeomTraits, TopTraits>& arr)
{
  // First use the internal validity check.
  if(!arr.is_valid())
    return (false);

  typedef Arrangement_on_surface_2<GeomTraits, TopTraits>   Arr;
  typedef GeomTraits                                        Geometry_traits_2;
  typedef typename Geometry_traits_2::X_monotone_curve_2    X_monotone_curve_2;

  // Define the sweep-line types:
  typedef Sweep_line_do_curves_x_visitor<Geometry_traits_2> Visitor;
  typedef Sweep_line_2<Geometry_traits_2, Visitor>          Sweep_line_2;

  // Define the arrangement iterator and circulator types:
  typedef typename Arr::Edge_const_iterator           Edge_const_iterator;
  typedef typename Arr::Halfedge_const_handle         Halfedge_const_handle;
  typedef typename Arr::Inner_ccb_const_iterator      Inner_ccb_const_iterator;
  typedef typename Arr::Face_const_iterator           Face_const_iterator;
  typedef typename Arr::Face_const_handle             Face_const_handle;
  typedef typename Arr::Vertex_const_handle           Vertex_const_handle;
  typedef typename Arr::Isolated_vertex_const_iterator
                                       Isolated_vertex_const_iterator;
  typedef typename Arr::Halfedge_around_vertex_const_circulator
                                       Halfedge_around_vertex_const_circulator;

  // Perform a sweep over all subcurves associated with arrangement edges.
  std::vector<X_monotone_curve_2>   curves_vec (arr.number_of_edges());
  Edge_const_iterator               eit;
  unsigned int                      i = 0;

  for (eit = arr.edges_begin(); eit != arr.edges_end(); ++eit, i++)
    curves_vec[i] = eit->curve();

  Visitor            visitor;
  const Geometry_traits_2 *traits = arr.geometry_traits();
  Sweep_line_2       sweep_line (traits, &visitor);

  visitor.sweep_xcurves (curves_vec.begin(), curves_vec.end());

  bool               are_edges_disjoint = (! visitor.found_intersection());

  if (!are_edges_disjoint)
  {
    CGAL_warning_msg (are_edges_disjoint,
                      "Arrangement edges are not disjoint in their interior.");
    return (false);
  }

  // Check that the holes and isolated vertices are located where they should.
  // At the same time, we prepare a vector that consists of all isolated
  // vertices and all leftmost vertices from every hole.
  std::list<std::pair<Vertex_const_handle, Face_const_handle> >  vf_list;

  typename Geometry_traits_2::Compare_xy_2  compare_xy =
    traits->compare_xy_2_object();
  Face_const_iterator               fit;
  Face_const_handle                 fh;
  Inner_ccb_const_iterator          ic_it;
  Halfedge_const_handle             ccb;
  Isolated_vertex_const_iterator    iv_it;
  Vertex_const_handle               left_v;
  bool                              is_first;

  for (fit = arr.faces_begin(); fit != arr.faces_end(); ++fit)
  {
    // Check all holes in the current face.
    fh = fit;
    for (ic_it = fh->inner_ccbs_begin();
         ic_it != fh->inner_ccbs_end(); ++ic_it)
    {
      ccb = *ic_it;
      is_first = true;

      do
      {
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
      if (iv_it->face() != fit)
        return (false);

      vf_list.push_back (std::make_pair (Vertex_const_handle(iv_it), fh));
    }
  }

  // Shoot a vertical ray from each vertex we have collected downward, and
  // check that this vertex is really contained in the proper face.
  typename Geometry_traits_2::Compare_y_at_x_right_2  comp_y_at_x_right =
                                      traits->compare_y_at_x_right_2_object();
  typename Geometry_traits_2::Compare_y_at_x_left_2   comp_y_at_x_left =
                                      traits->compare_y_at_x_left_2_object();

  typename std::list<std::pair<Vertex_const_handle,
                               Face_const_handle> >::iterator    vf_iter;
  typename TopTraits::Default_point_location_strategy            def_pl (arr);
  Vertex_const_handle                                            curr_v;
  Object                                                         obj;
  Halfedge_const_handle                                          he_below;
  Vertex_const_handle                                            v_below;
  Face_const_handle                                              in_face;
  Halfedge_around_vertex_const_circulator                        first, circ;
  bool                                                           assign_ok;
  const Halfedge_const_handle                                    invalid_he;

  for (vf_iter = vf_list.begin(); vf_iter != vf_list.end(); ++vf_iter)
  {
    // Perform ray-shooting from the current vertex.
    curr_v = vf_iter->first;
    obj = def_pl.ray_shoot_down(curr_v->point());

    if (CGAL::assign(he_below, obj))
    {
      // Hit an edge - take the incident face of the halfedge directed to the
      // right.
      if (he_below->direction() == ARR_RIGHT_TO_LEFT)
        he_below = he_below->twin();

      in_face = he_below->face();
    }
    else if (CGAL::assign(v_below, obj))
    {
      // Hit a vertex.
      if(v_below->is_isolated())
      {
        in_face = v_below->face();
      }
      else
      {
        // Get the first halfedge around v_below that is directed from left to
        // right and the first halfedge that is directed from right to left.
        first = circ = v_below->incident_halfedges();
        Halfedge_const_handle he_left;  // A halfedge to the left of v_below.
        Halfedge_const_handle he_right; // A halfedge to the right of v_below.

        do
        {
          if (circ->direction() == ARR_LEFT_TO_RIGHT)
          {
            he_left = circ;
          }
          else
          {
            he_right = circ;
            if (he_left != invalid_he && he_right != invalid_he)
              break;
          }
          ++circ;

        } while(circ != first);

        CGAL_assertion (he_left != invalid_he || he_right != invalid_he);

        if (he_left != invalid_he && he_right != invalid_he)
        {
          while (he_left->direction() == ARR_LEFT_TO_RIGHT)
            he_left = he_left->next()->twin();

          he_left = he_left->twin()->prev();
          CGAL_assertion (he_left->direction() == ARR_LEFT_TO_RIGHT);
          in_face = he_left->face();
        }
        else if (he_left != invalid_he)
        {
          Comparison_result     res;
          Halfedge_const_handle he_curr = he_left;

          do // as long as we have next he_left halfedge which is above
          {
            he_left = he_curr;
            he_curr = he_left->next()->twin();
            res = comp_y_at_x_left (he_curr->curve(),
                                    he_left->curve(),
                                    v_below->point());
          } while(res == LARGER);
          in_face = he_left->face();
        }
        else
        {
          Comparison_result     res;
          Halfedge_const_handle he_curr = he_right;

          do // as long as we have he_right halfedge which is below
          {
            he_right = he_curr;
            he_curr = he_right->next()->twin();
            res = comp_y_at_x_right (he_curr->curve(),
                                     he_right->curve(),
                                     v_below->point());
          } while(res == SMALLER);
          in_face = he_right->face();
        }
      }
    }
    else
    {
      // Hit nothing (an unbounded face is returned).
      assign_ok = CGAL::assign(in_face, obj);

      CGAL_assertion (assign_ok && in_face->is_unbounded());

      if (! assign_ok)
        return (false);
    }

    if (vf_iter->second != in_face)
    {
      CGAL_warning_msg (false,
                        "An inner component is located in the wrong face.");
      return (false);
    }
  }

  // If we reached here, the arrangement is valid:
  return (true);
}

//-----------------------------------------------------------------------------
// Compute the zone of the given x-monotone curve in the existing arrangement.
// Meaning, it output the arrangment's vertices, edges and faces that the
// x-monotone curve intersects.
template <class GeomTraits, class TopTraits,
  class OutputIterator, class PointLocation>
OutputIterator zone (Arrangement_on_surface_2<GeomTraits, TopTraits>& arr,
                     const typename GeomTraits::X_monotone_curve_2& c,
                     OutputIterator oi,
                     const PointLocation& pl)
{
  // Obtain an arrangement accessor.
  typedef Arrangement_on_surface_2<GeomTraits,TopTraits>
    Arrangement_on_surface_2;

  // Define a zone-computation object an a visitor that performs the
  // intersection check.
  typedef Arr_compute_zone_visitor<Arrangement_on_surface_2, OutputIterator>
    Zone_visitor;

  Zone_visitor                                     visitor (oi);
  Arrangement_zone_2<Arrangement_on_surface_2, Zone_visitor>
    arr_zone (arr, &visitor);

  arr_zone.init (c, pl);
  arr_zone.compute_zone();

  return (oi);
}

//-----------------------------------------------------------------------------
// Compute the zone of the given x-monotone curve in the existing arrangement.b
// Overloaded version with no point location object - the walk point-location
// strategy is used as default.
//
template <class GeomTraits, class TopTraits, class OutputIterator>
OutputIterator zone (Arrangement_on_surface_2<GeomTraits, TopTraits>& arr,
                     const typename GeomTraits::X_monotone_curve_2& c,
                     OutputIterator oi)
{
  // Create a default point-location object and use it to insert the curve.
  typename TopTraits::Default_point_location_strategy    def_pl (arr);

  //insert the curve using the walk point location
  zone (arr, c, oi, def_pl);
  return oi;
}


//-----------------------------------------------------------------------------
// Checks if the given x-monotone curve intersects the existing arrangement.
// The last parameter is used to resolve ambiguity between this function and
// do_intersect of Curve_2 in case that X_monotone_curve_2 and Curve_2 are the
// same class. The last parameter should be boost::true_type but we used a
// workaround since it didn't compile in FC3_g++-3.4.4 with the error of:
//
// error: no matching function for call to `do_intersect(Arrangement_on_surface_2<>&,
// const Arr_segment_2&, const Arr_walk_along_line_point_location<>&, mpl_::bool_< true>)'
//
template <class GeomTraits, class TopTraits, class PointLocation>
bool do_intersect (Arrangement_on_surface_2<GeomTraits, TopTraits>& arr,
                   const typename GeomTraits::X_monotone_curve_2& c,
                   const PointLocation& pl, boost::is_same<int, int>::type)
{
  // Obtain an arrangement accessor.
  typedef Arrangement_on_surface_2<GeomTraits,TopTraits>
    Arrangement_on_surface_2;

  // Define a zone-computation object an a visitor that performs the
  // intersection check.
  typedef Arr_do_intersect_zone_visitor<Arrangement_on_surface_2>  Zone_visitor;

  Zone_visitor                                     visitor;
  Arrangement_zone_2<Arrangement_on_surface_2, Zone_visitor>
    arr_zone (arr, &visitor);

  arr_zone.init (c, pl);
  arr_zone.compute_zone();

  return (visitor.do_intersect());
}

//-----------------------------------------------------------------------------
// Checks if the given curve intersects the existing arrangement.
// The last parameter is used to resolve ambiguity between this function and
// do_intersect of X_monotone_curve_2 in case that X_monotone_curve_2 and
// Curve_2 are the same class.
// The last parameter should be boost::false_type but we used a
// workaround since it didn't compile in FC3_g++-3.4.4 with the error of:
//
// error: no matching function for call to `do_intersect(Arrangement_on_surface_2<>&,
// const Arr_segment_2&, const Arr_walk_along_line_point_location<>&, mpl_::bool_< true>)'
//
template <class GeomTraits, class TopTraits, class PointLocation>
bool do_intersect (Arrangement_on_surface_2<GeomTraits, TopTraits>& arr,
                   const typename GeomTraits::X_monotone_curve_2& c,
                   const PointLocation& pl, boost::is_same<int, double>::type)
{
  // Obtain an arrangement accessor.
  typedef Arrangement_on_surface_2<GeomTraits,TopTraits>
    Arrangement_on_surface_2;

  // Break the input curve into x-monotone subcurves and isolated points.
  typedef Arr_traits_adaptor_2<GeomTraits>                   Traits_adaptor_2;

  const Traits_adaptor_2 * traits =
    static_cast<const Traits_adaptor_2*> (arr.geometry_traits());

  std::list<CGAL::Object>                        x_objects;
  std::list<CGAL::Object>::const_iterator        obj_iter;
  const typename GeomTraits::X_monotone_curve_2  *x_curve;
  const typename GeomTraits::Point_2             *iso_p;

  traits->make_x_monotone_2_object() (c,
                                      std::back_inserter (x_objects));

  // Insert each x-monotone curve into the arrangement.
  for (obj_iter = x_objects.begin(); obj_iter != x_objects.end(); ++obj_iter)
  {
    // Act according to the type of the current object.
    x_curve = object_cast<typename GeomTraits::X_monotone_curve_2>
      (&(*obj_iter));
    if (x_curve != NULL)
    {
      // Check if the x-monotone subcurve intersects the arrangement.
      if (do_intersect(arr, *x_curve, pl) == true)
        return true;
    }
    else
    {
      iso_p = object_cast<typename GeomTraits::Point_2> (&(*obj_iter));
      CGAL_assertion (iso_p != NULL);

      // Check whether the isolated point lies inside a face (otherwise,
      // it conincides with a vertex or an edge).
      CGAL::Object  obj = pl.locate (*iso_p);

      return (object_cast<typename
              Arrangement_on_surface_2::Face_const_handle>(&obj) != NULL);
    }
  }

  // If we reached here, the curve does not intersect the arrangement.
  return (false);
}

//-----------------------------------------------------------------------------
// Common interface for the do_intersect of the Curve_2 and X_monotone_curve_2
template <class GeomTraits, class TopTraits, class Curve, class PointLocation>
bool do_intersect (Arrangement_on_surface_2<GeomTraits, TopTraits>& arr,
                   const Curve& c, const PointLocation& pl)
{
  typedef typename GeomTraits::X_monotone_curve_2       X_monotone_curve_2;

  typedef typename boost::is_same<Curve, X_monotone_curve_2>::type
    Is_x_monotone;

  return do_intersect(arr, c, pl, Is_x_monotone());
}

//-----------------------------------------------------------------------------
// Checks if the given curve intersects the existing arrangement.
// Overloaded version with no point location object - the walk point-location
// strategy is used as default.
template <class GeomTraits, class TopTraits, class Curve>
bool do_intersect (Arrangement_on_surface_2<GeomTraits, TopTraits>& arr,
                   const Curve& c)
{
  // Create a default point-location object and use it to insert the curve.
  typename TopTraits::Default_point_location_strategy    def_pl (arr);

  // check if the curve intersects the arrangement using the walk point
  // location.
  return do_intersect (arr, c, def_pl);
}


} //namespace CGAL

#endif
