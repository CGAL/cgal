// Copyright (c) 2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
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
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>
//                 Ron Wein <wein@post.tau.ac.il>
//                 Efi Fogel <efif@post.tau.ac.il>

#ifndef CGAL_ARR_BATCHED_POINT_LOCATION_H
#define CGAL_ARR_BATCHED_POINT_LOCATION_H

#include <CGAL/license/Arrangement_on_surface_2.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Arrangement_on_surface_2.h>
#include <CGAL/No_intersection_surface_sweep_2.h>
#include <CGAL/Arr_point_location/Arr_batched_point_location_traits_2.h>
#include <CGAL/Surface_sweep_2/No_overlap_event.h>
#include <CGAL/Surface_sweep_2/No_overlap_subcurve.h>
#include <CGAL/Surface_sweep_2/Arr_batched_pl_ss_visitor.h>

#include <vector>
#include <boost/mpl/if.hpp>
#include <boost/type_traits.hpp>

namespace CGAL {

namespace Ss2 = Surface_sweep_2;

/*! Issue a batched point-location query on an arrangement given an input
 * range of points.
 * \param arr The arrangement.
 * \param points_begin An iterator for the range of query points.
 * \param points_end A past-the-end iterator for the range of query points.
 * \param oi Output: An output iterator for the query results.
 * \pre The value-type of PointsIterator is Arrangement::Point_2,
 *      and the value-type of OutputIterator is is pair<Point_2, Result>,
 *      where Result is either
 *       (i) Object or
 *      (ii) boost::optional<boost::variant<Vertex_const_handle,
 *                                          Halfedge_const_handle,
 *                                          Face_const_handle> >.
 *      It represents the arrangement feature containing the point.
 */
template <typename GeometryTraits_2, typename TopologyTraits,
          typename PointsIterator, typename OutputIterator>
OutputIterator
locate(const Arrangement_on_surface_2<GeometryTraits_2, TopologyTraits>& arr,
       PointsIterator points_begin, PointsIterator points_end,
       OutputIterator oi)
{
  typedef GeometryTraits_2                              Gt2;
  typedef TopologyTraits                                Tt;
  typedef OutputIterator                                Output_iterator;

  // Arrangement types:
  typedef Arrangement_on_surface_2<Gt2, Tt>             Arr;
  typedef typename Arr::Halfedge_const_handle           Halfedge_const_handle;
  typedef typename Arr::Vertex_const_iterator           Vertex_const_iterator;
  typedef typename Arr::Edge_const_iterator             Edge_const_iterator;
  typedef typename Arr::Vertex_const_handle             Vertex_const_handle;
  typedef typename Arr::Halfedge_const_handle           Halfedge_const_handle;
  typedef typename Arr::Allocator                       Allocator;

  // Surface sweep types
  typedef Arr_batched_point_location_traits_2<Arr>      Bgt2;
  typedef Ss2::No_overlap_event<Bgt2, Allocator>        Bpl_event;
  typedef Ss2::No_overlap_subcurve<Bgt2, Bpl_event, Allocator>
                                                        Bpl_curve;
  typedef typename Tt::template Batched_point_location_helper<Bpl_event,
                                                              Bpl_curve>
                                                        Bpl_helper;
  typedef Arr_batched_pl_ss_visitor<Bpl_helper, Output_iterator>
                                                        Bpl_visitor;
  typedef typename Bgt2::X_monotone_curve_2             Bpl_x_monotone_curve_2;
  typedef typename Bgt2::Point_2                        Bpl_point_2;

  // Go over all arrangement edges and collect their associated x-monotone
  // curves. To each curve we attach a halfedge handle going from right to
  // left.
  std::vector<Bpl_x_monotone_curve_2> xcurves_vec(arr.number_of_edges());
  Edge_const_iterator eit;
  std::size_t i(0);
  for (eit = arr.edges_begin(); eit != arr.edges_end(); ++eit) {
    // Associate each x-monotone curve with the halfedge that represent it
    // that is directed from right to left.
    Halfedge_const_handle he =
      (eit->direction() == ARR_RIGHT_TO_LEFT) ? eit : eit->twin();
    xcurves_vec[i++] = Bpl_x_monotone_curve_2(eit->curve(), he);
  }

  // Go over all isolated vertices and collect their points. To each point
  // we attach its vertex handle.
  std::vector<Bpl_point_2> iso_pts_vec(arr.number_of_isolated_vertices());
  Vertex_const_iterator vit;
  i = 0;
  for (vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit) {
    if (vit->is_isolated()) {
      Vertex_const_handle iso_v = vit;
      iso_pts_vec[i++] = Bpl_point_2(vit->point(), iso_v);
    }
  }

  // Obtain a extended traits-class object.
  const Gt2* geom_traits = arr.geometry_traits();

  /* We would like to avoid copy construction of the geometry traits class.
   * Copy construction is undesired, because it may results with data
   * duplication or even data loss.
   *
   * If the type Bgt2 is the same as the type Gt2, use a reference to Gt2
   * to avoid constructing a new one.  Otherwise, instantiate a local variable
   * of the former and provide the later as a single parameter to the
   * constructor.
   *
   * Use the form 'A a(*b);' and not ''A a = b;' to handle the case where A has
   * only an implicit constructor, (which takes *b as a parameter).
   */
  typename boost::mpl::if_<boost::is_same<Gt2, Bgt2>, const Bgt2&, Bgt2>::type
    ex_traits(*geom_traits);

  // Define the sweep-line visitor and perform the sweep.
  Bpl_visitor visitor(&arr, oi);
  Ss2::No_intersection_surface_sweep_2<Bpl_visitor>
    surface_sweep(&ex_traits, &visitor);
  surface_sweep.sweep(xcurves_vec.begin(), xcurves_vec.end(), // Curves.
                      iso_pts_vec.begin(), iso_pts_vec.end(), // Action points.
                      points_begin, points_end);              // Query points.

  return oi;
}

} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif
