// Copyright (c) 2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Ron Wein <wein@post.tau.ac.il>

#ifndef CGAL_ARR_VERTICAL_DECOMPOSITION_2_H
#define CGAL_ARR_VERTICAL_DECOMPOSITION_2_H

#include <CGAL/license/Arrangement_on_surface_2.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Arrangement_on_surface_2.h>
#include <CGAL/No_intersection_surface_sweep_2.h>
#include <CGAL/Arr_point_location/Arr_batched_point_location_traits_2.h>
#include <CGAL/Surface_sweep_2/No_overlap_event.h>
#include <CGAL/Surface_sweep_2/No_overlap_subcurve.h>
#include <CGAL/Surface_sweep_2/Arr_vert_decomp_ss_visitor.h>

#include <vector>
#include <boost/mpl/if.hpp>
#include <boost/type_traits.hpp>

namespace CGAL {

namespace Ss2 = Surface_sweep_2;

/*! Perform a vertical decomposition of an arrangement, by performing a
 * "batched vertical ray-shooting" query from all arrangement vertices.
 * \param arr The arrangement.
 * \param oi An output iterator of the vertices, each paired with a pair of
 *           arrangement features that lie below and above it, respectively.
 *           The vertices are sorted by increasing xy-order.
 *           The OutputIterator dereferences the type \c
 *           pair<Vertex_const_handle, pair<Vert_type, Vert_type> >, where
 *           \c Vert_type is an optional handle to an arrangement feature.
 * \return A past-the-end iterator for the ordered arrangement vertices.
 */
template <typename GeometryTraits_2, typename TopologyTraits,
          typename OutputIterator>
OutputIterator
decompose(const Arrangement_on_surface_2<GeometryTraits_2, TopologyTraits>& arr,
          OutputIterator oi)
{
  typedef GeometryTraits_2                              Gt2;
  typedef TopologyTraits                                Tt;
  typedef OutputIterator                                Output_iterator;

  // Arrangement types:
  typedef Arrangement_on_surface_2<Gt2, Tt>             Arr;
  typedef typename Arr::Vertex_const_iterator           Vertex_const_iterator;
  typedef typename Arr::Edge_const_iterator             Edge_const_iterator;
  typedef typename Arr::Vertex_const_handle             Vertex_const_handle;
  typedef typename Arr::Halfedge_const_handle           Halfedge_const_handle;
  typedef typename Arr::X_monotone_curve_2              X_monotone_curve_2;
  typedef typename Arr::Point_2                         Point_2;
  typedef typename Arr::Allocator                       Allocator;

  // Surface sweep types:
  typedef Arr_batched_point_location_traits_2<Arr>      Vgt2;
  typedef Ss2::No_overlap_event<Vgt2, Allocator>        Vd_event;
  typedef Ss2::No_overlap_subcurve<Vgt2, Vd_event, Allocator>
                                                        Vd_curve;
  typedef typename Tt::template
    Vertical_decomposition_helper<Vd_event, Vd_curve>   Vd_helper;
  typedef Arr_vert_decomp_ss_visitor<Vd_helper, Output_iterator>
                                                        Vd_visitor;
  typedef typename Vgt2::X_monotone_curve_2             Vd_x_monotone_curve_2;
  typedef typename Vgt2::Point_2                        Vd_point_2;

  // Go over all arrangement edges and collect their associated x-monotone
  // curves. To each curve we attach a halfedge handle going from right to
  // left.
  std::vector<Vd_x_monotone_curve_2> xcurves_vec(arr.number_of_edges());

  size_t i(0);
  Edge_const_iterator eit;
  for (eit = arr.edges_begin(); eit != arr.edges_end(); ++eit) {
    // Associate each x-monotone curve with the halfedge that represents it
    // and is directed from right to left.
    Halfedge_const_handle he = (eit->direction() == ARR_RIGHT_TO_LEFT) ?
      eit : eit->twin();
    //attempt to solve compile problem in one of the tests. created the
    // tmp_curve instead of passing eit->curve() as a parmeter to the function
    X_monotone_curve_2 tmp_curve = eit->curve();
    xcurves_vec[i++] = Vd_x_monotone_curve_2(tmp_curve, he);
  }

  // Go over all isolated vertices and collect their points. To each point
  // we attach its vertex handle.
  std::vector<Vd_point_2> iso_pts_vec(arr.number_of_isolated_vertices());
  i = 0;
  Vertex_const_iterator vit;
  for (vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit) {
    // Associate isolated point with the vertex that represents it.
    if (vit->is_isolated()) {
      Vertex_const_handle iso_v = vit;
      //attempt to solve compile problem in one of the tests. created the
      // tmp_curve instead of passing eit->curve() as a parmeter to the
      // function
      Point_2 tmp_point = vit->point();
      iso_pts_vec[i++] = Vd_point_2(tmp_point, iso_v);
    }
  }

  // Obtain a extended traits-class object.
  const Gt2* geom_traits = arr.geometry_traits();

  /* We would like to avoid copy construction of the geometry traits class.
   * Copy construction is undesired, because it may results with data
   * duplication or even data loss.
   *
   * If the type Vgt2 is the same as the type Gt2, use a
   * reference to Gt2 to avoid constructing a new one.  Otherwise,
   * instantiate a local variable of the former and provide the later as a
   * single parameter to the constructor.
   *
   * Use the form 'A a(*b);' and not ''A a = b;' to handle the case where A has
   * only an implicit constructor, (which takes *b as a parameter).
   */
  typename boost::mpl::if_<boost::is_same<Gt2, Vgt2>, const Vgt2&, Vgt2>::type
    ex_traits(*geom_traits);

  // Define the sweep-line visitor and perform the sweep.
  Vd_visitor visitor(&arr, &oi);
  Ss2::No_intersection_surface_sweep_2<Vd_visitor>
    surface_sweep(&ex_traits, &visitor);
  surface_sweep.sweep(xcurves_vec.begin(), xcurves_vec.end(),  // Curves.
                      iso_pts_vec.begin(), iso_pts_vec.end()); // Action points.

  // Return a past-the-end iterator.
  return oi;
}

} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif
