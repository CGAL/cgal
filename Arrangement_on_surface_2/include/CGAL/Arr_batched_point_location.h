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
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>
//                 Ron Wein <wein@post.tau.ac.il>
//                 Efi Fogel <efif@post.tau.ac.il>

#ifndef CGAL_ARR_BATCHED_POINT_LOCATION_H
#define CGAL_ARR_BATCHED_POINT_LOCATION_H

#include <CGAL/Arrangement_on_surface_2.h>
#include <CGAL/Basic_sweep_line_2.h>

#include <vector>
#include <boost/mpl/if.hpp>
#include <boost/type_traits.hpp>

namespace CGAL {

/*!
 * Issue a batched point-location query on an arrangement given an input
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
template<typename GeomTraits, typename TopTraits,
         typename PointsIterator, typename OutputIterator> 
OutputIterator
locate(const Arrangement_on_surface_2<GeomTraits, TopTraits>& arr,
       PointsIterator points_begin, PointsIterator points_end,
       OutputIterator oi)
{
  // Arrangement types:
  typedef Arrangement_on_surface_2<GeomTraits, TopTraits>  Arr;
  typedef typename TopTraits::template
          Sweep_line_batched_point_location_visitor<OutputIterator>
                                                           Bpl_visitor;

  typedef typename Arr::Halfedge_const_handle          Halfedge_const_handle;
  typedef typename Arr::Vertex_const_iterator          Vertex_const_iterator;
  typedef typename Arr::Edge_const_iterator            Edge_const_iterator;
  typedef typename Arr::Vertex_const_handle            Vertex_const_handle;
  typedef typename Arr::Halfedge_const_handle          Halfedge_const_handle;

  typedef typename Bpl_visitor::Traits_2               Bpl_traits_2;
  typedef typename Bpl_traits_2::X_monotone_curve_2    Bpl_x_monotone_curve_2;
  typedef typename Bpl_traits_2::Point_2               Bpl_point_2;

  // Go over all arrangement edges and collect their associated x-monotone
  // curves. To each curve we attach a halfedge handle going from right to
  // left.
  std::vector<Bpl_x_monotone_curve_2>  xcurves_vec(arr.number_of_edges());
  Edge_const_iterator                  eit;
  unsigned int                         i = 0;
  for (eit = arr.edges_begin(); eit != arr.edges_end(); ++eit) {
    // Associate each x-monotone curve with the halfedge that represent it
    // that is directed from right to left.
    Halfedge_const_handle he =
      (eit->direction() == ARR_RIGHT_TO_LEFT) ? eit : eit->twin();
    xcurves_vec[i++] = Bpl_x_monotone_curve_2(eit->curve(), he);
  }

  // Go over all isolated vertices and collect their points. To each point
  // we attach its vertex handle.
  std::vector<Bpl_point_2>    iso_pts_vec(arr.number_of_isolated_vertices());
  Vertex_const_iterator       vit;
  i = 0;
  for (vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit) {
    if (vit->is_isolated()) {
      Vertex_const_handle iso_v = vit;
      iso_pts_vec[i++] = Bpl_point_2(vit->point(), iso_v);
    }
  }
    
  // Obtain a extended traits-class object.
  GeomTraits* geom_traits = const_cast<GeomTraits*>(arr.geometry_traits());

  /* We would like to avoid copy construction of the geometry traits class.
   * Copy construction is undesired, because it may results with data
   * duplication or even data loss.
   *
   * If the type Bpl_traits_2 is the same as the type
   * GeomTraits, use a reference to GeomTraits to avoid constructing a new one.
   * Otherwise, instantiate a local variable of the former and provide
   * the later as a single parameter to the constructor.
   * 
   * Use the form 'A a(*b);' and not ''A a = b;' to handle the case where A has
   * only an implicit constructor, (which takes *b as a parameter).
   */
  typename boost::mpl::if_<boost::is_same<GeomTraits, Bpl_traits_2>,
                           const Bpl_traits_2&, Bpl_traits_2>::type
    ex_traits(*geom_traits);

  // Define the sweep-line visitor and perform the sweep.
  Bpl_visitor   visitor(&arr, oi);
  Basic_sweep_line_2<typename Bpl_visitor::Traits_2, Bpl_visitor,
                     typename Bpl_visitor::Subcurve,
                     typename Bpl_visitor::Event>
    sweep_line(&ex_traits, &visitor);
  sweep_line.sweep(xcurves_vec.begin(), xcurves_vec.end(), // Curves.
                   iso_pts_vec.begin(), iso_pts_vec.end(), // Action points.
                   points_begin, points_end);              // Query points.

  return oi;
}

} //namespace CGAL

#endif
