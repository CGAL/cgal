// Copyright (c) 2006  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s)     : Ron Wein <wein@post.tau.ac.il>

#ifndef CGAL_ARR_VERTICAL_DECOMPOSITION_2_H
#define CGAL_ARR_VERTICAL_DECOMPOSITION_2_H

#include <CGAL/Basic_sweep_line_2.h>
#include <CGAL/Sweep_line_2/Sweep_line_subcurve.h>
#include <CGAL/Sweep_line_2/Sweep_line_event.h>
#include <CGAL/Arr_point_location/Arr_vertical_decomposition_visitor.h>
#include <CGAL/Arr_point_location/Arr_batched_point_location_traits_2.h>
#include <vector>

CGAL_BEGIN_NAMESPACE

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
template<typename Traits, typename Dcel, typename OutputIterator>
OutputIterator decompose (const Arrangement_2<Traits,Dcel>& arr, 
                          OutputIterator oi)
{
  // Arrangement types:
  typedef Arrangement_2<Traits,Dcel>                    Arrangement_2;
  typedef typename Arrangement_2::Traits_2              Traits_2;
  typedef typename Traits_2::X_monotone_curve_2         Base_x_monotone_curve_2;

  typedef typename Arrangement_2::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Arrangement_2::Vertex_const_iterator Vertex_const_iterator;
  typedef typename Arrangement_2::Edge_const_iterator   Edge_const_iterator;
  typedef typename Arrangement_2::Size                  Size;

  // Define the meta-traits class:
  typedef Arr_batched_point_location_traits_2<Arrangement_2>
                                                        Meta_traits_2;

  typedef typename Meta_traits_2::X_monotone_curve_2    X_monotone_curve_2;
  typedef typename Meta_traits_2::Point_2               Point_2;

  // Define the sweep-line visitor:
  typedef Arr_vertical_decomposition_visitor<Meta_traits_2,
                                             Arrangement_2,
                                             OutputIterator>  Visitor;
  
  typedef Basic_sweep_line_2<Meta_traits_2, Visitor>    Sweep_line;

  // Go over all arrangement edges.
  std::vector<X_monotone_curve_2>  xcurves_vec (arr.number_of_edges());
  Edge_const_iterator              eit;
  Size                             i = 0;

  for (eit = arr.edges_begin(); eit != arr.edges_end(); ++eit, ++i) 
  {
    // Associate each x-monotone curve with the halfedge that represents it
    // and is directed from right to left.
    if (eit->direction() == LEFT_TO_RIGHT)
      xcurves_vec[i] = X_monotone_curve_2(eit->curve(),eit);
    else
      xcurves_vec[i] = X_monotone_curve_2(eit->curve(),eit->twin());
  }

  // Go over all arrangement vertices.
  Vertex_const_iterator   vit;
  std::vector<Point_2>    iso_points;

  for (vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit)
  {
    // Associate isolated point with the vertex that represents it.
    if (vit->is_isolated())
      iso_points.push_back (Point_2 (vit->point(), vit));
  }
  
  // Perform the sweep and fill the vertical mapping.
  Visitor          visitor (arr, oi);
  Meta_traits_2    meta_tr (*(arr.get_traits()));
  Sweep_line       sweep_line (&meta_tr ,&visitor);
  
  sweep_line.sweep (xcurves_vec.begin(),
                    xcurves_vec.end(),
                    iso_points.begin(),  
                    iso_points.end()); 

  // Return a past-the-end iterator.
  return (visitor.output_iterator());
}


CGAL_END_NAMESPACE

#endif
