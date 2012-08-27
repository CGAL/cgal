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
// Author(s)     : Ron Wein <wein@post.tau.ac.il>

#ifndef CGAL_ARR_VERTICAL_DECOMPOSITION_2_H
#define CGAL_ARR_VERTICAL_DECOMPOSITION_2_H

#include <CGAL/Arrangement_on_surface_2.h>
#include <CGAL/Basic_sweep_line_2.h>

#include <vector>
#include <boost/mpl/if.hpp>
#include <boost/type_traits.hpp>

namespace CGAL {

/*!
 * Perform a vertical decomposition of an arrangement, by performing a
 * "batched vertical ray-shooting" query from all arrangement vertices.
 * \param arr The arrangement.
 * \param oi Output: An output iterator of the vertices, each paired with
 *                   a pair of arrangement features that lie below and above
 *                   it, respectively.
 *                   The vertices are sorted by increasing xy-order.
 * \return A past-the-end iterator for the ordered arrangement vertices.
 * \pre The value-type of OutputIterator is
 *      pair<Vertex_const_handle, pair<Object, Object> >, where
 *      the Object represents a handle to an arrangement feature.
 */
template<typename GeomTraits, typename TopTraits,
         typename OutputIterator>
OutputIterator
decompose(const Arrangement_on_surface_2<GeomTraits, TopTraits>& arr,
          OutputIterator oi)
{
  // Arrangement types:
  typedef Arrangement_on_surface_2<GeomTraits, TopTraits> Arrangement_2;
  typedef typename TopTraits::template
    Sweep_line_vertical_decomposition_visitor<OutputIterator>
                                                          Vd_visitor;

  typedef typename Arrangement_2::Vertex_const_iterator   Vertex_const_iterator;
  typedef typename Arrangement_2::Edge_const_iterator     Edge_const_iterator;
  typedef typename Arrangement_2::Vertex_const_handle     Vertex_const_handle;
  typedef typename Arrangement_2::Halfedge_const_handle   Halfedge_const_handle;
  typedef typename Arrangement_2::X_monotone_curve_2      X_monotone_curve_2;
  typedef typename Arrangement_2::Point_2                 Point_2;
  typedef typename Vd_visitor::Traits_2                   Vd_traits_2;
  typedef typename Vd_traits_2::X_monotone_curve_2        Vd_x_monotone_curve_2;
  typedef typename Vd_traits_2::Point_2                   Vd_point_2;

  // Go over all arrangement edges and collect their associated x-monotone
  // curves. To each curve we attach a halfedge handle going from right to
  // left.
  std::vector<Vd_x_monotone_curve_2>  xcurves_vec (arr.number_of_edges());
  Edge_const_iterator                 eit;
  Halfedge_const_handle               he;
  unsigned int                        i = 0;

  for (eit = arr.edges_begin(); eit != arr.edges_end(); ++eit, ++i) 
  {
    // Associate each x-monotone curve with the halfedge that represents it
    // and is directed from right to left.
    if (eit->direction() == ARR_RIGHT_TO_LEFT)
      he = eit;
    else
      he = eit->twin();
    //attempt to solve compile problem in one of the tests. created the
    // tmp_curve instead of passing eit->curve() as a parmeter to the function
    X_monotone_curve_2 tmp_curve = eit->curve();
    xcurves_vec[i] = Vd_x_monotone_curve_2 (tmp_curve, he);
  }

  // Go over all isolated vertices and collect their points. To each point
  // we attach its vertex handle.
  std::vector<Vd_point_2>     iso_pts_vec (arr.number_of_isolated_vertices());
  Vertex_const_iterator       vit;
  Vertex_const_handle         iso_v;

  i = 0;
  for (vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit)
  {
    // Associate isolated point with the vertex that represents it.
    if (vit->is_isolated())
    {
      iso_v = vit;
      //attempt to solve compile problem in one of the tests. created the
      // tmp_curve instead of passing eit->curve() as a parmeter to the
      // function  
      Point_2 tmp_point = vit->point();
      iso_pts_vec[i] = Vd_point_2 (tmp_point, iso_v);
      ++i;
    }
  }

  // Obtain a extended traits-class object.
  const GeomTraits * geom_traits = arr.geometry_traits();

  /* We would like to avoid copy construction of the geometry traits class.
   * Copy construction is undesired, because it may results with data
   * duplication or even data loss.
   *
   * If the type Vd_traits_2 is the same as the type
   * GeomTraits, use a reference to GeomTraits to avoid constructing a new one.
   * Otherwise, instantiate a local variable of the former and provide
   * the later as a single parameter to the constructor.
   * 
   * Use the form 'A a(*b);' and not ''A a = b;' to handle the case where A has
   * only an implicit constructor, (which takes *b as a parameter).
   */
  typename boost::mpl::if_<boost::is_same<GeomTraits, Vd_traits_2>,
                           const Vd_traits_2&, Vd_traits_2>::type
    ex_traits(*geom_traits);

  // Define the sweep-line visitor and perform the sweep.
  Vd_visitor    visitor (&arr, &oi);
  Basic_sweep_line_2<typename Vd_visitor::Traits_2,
                     Vd_visitor,
                     typename Vd_visitor::Subcurve,
                     typename Vd_visitor::Event>
    sweep_line (&ex_traits, &visitor);

  sweep_line.sweep (xcurves_vec.begin(), xcurves_vec.end(),  // Curves.
                    iso_pts_vec.begin(), iso_pts_vec.end()); // Action points.

  // Return a past-the-end iterator.
  return (oi);
}


} //namespace CGAL

#endif
