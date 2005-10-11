// Copyright (c) 2005  Tel-Aviv University (Israel).
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef CGAL_ARR_BATCHED_POINT_LOCATION_H
#define CGAL_ARR_BATCHED_POINT_LOCATION_H

#include <CGAL/Basic_sweep_line_2.h>
#include <CGAL/Sweep_line_2/Sweep_line_subcurve.h>
#include <CGAL/Sweep_line_2/Sweep_line_event.h>
#include <CGAL/Arr_point_location/Arr_batched_point_location_visitor.h>
#include <CGAL/Arr_point_location/Arr_batched_point_location_meta_traits.h>
#include <vector>
#include <iostream>

CGAL_BEGIN_NAMESPACE

/*!
 * Perform a batched point-location operation on an arrangement.
 * \param arr The arrangement.
 * \param points_begin An iterator for the range of query points.
 * \param points_end A past-the-end iterator for the range of query points.
 * \param oi Output: An output iterator for the query results.
 * \pre The value-type of PointsIterator is Arrangement::Point_2,
 *      and the value-type of OutputIterator is pair<Point_2,Object>, where
 *      the Object represents the arrangement feature containing the points.
 */
template<class Arrangement, class PointsIterator, class OutputIterator> 
OutputIterator locate(const Arrangement& arr,
                      PointsIterator points_begin,
                      PointsIterator points_end,
                      OutputIterator oi)
{
  // Arrangement types:
  typedef typename Arrangement::Traits_2               Traits_2;
  typedef typename Traits_2::X_monotone_curve_2        Base_X_monotone_curve_2;
  typedef typename Arrangement::Halfedge_const_handle  Halfedge_const_handle;
  typedef typename Arrangement::Vertex_const_iterator  Vertex_const_iterator;
  typedef typename Arrangement::Edge_const_iterator    Edge_const_iterator;
  typedef typename Arrangement::Size                   Size;


  // Define meta-traits class for the batched point location:
  typedef Arr_batched_point_location_meta_traits<Traits_2, Arrangement>
                                                       Meta_traits_2;

  typedef typename Meta_traits_2::X_monotone_curve_2   X_monotone_curve_2;
  typedef typename Meta_traits_2::Point_2              Point_2;

  // Define the sweep-line visitor:
  typedef Arr_batched_point_location_visitor<Meta_traits_2,
                                             OutputIterator,
                                             Arrangement>  Visitor;
  
  typedef Basic_sweep_line_2<Meta_traits_2, Visitor>       Sweep_line;

  // Go over all arrangement edges.
  std::vector<X_monotone_curve_2>  xcurves_vec (arr.number_of_edges());
  Edge_const_iterator              eit;
  Size                             i = 0;

  for (eit = arr.edges_begin(); eit != arr.edges_end(); ++eit, ++i) 
  {
    // Associate each x-monotone curve with the halfedge that represent it
    // that is directed from right to left.
    if (eit->direction() == LARGER)
      xcurves_vec[i] = X_monotone_curve_2(eit->curve(),eit);
    else
      xcurves_vec[i] = X_monotone_curve_2(eit->curve(),eit->twin());
  }

  std::vector<Point_2>    iso_points;
  for(Vertex_const_iterator v_itr = arr.vertices_begin();
      v_itr != arr.vertices_end();
      ++v_itr)
  {
    if(v_itr->is_isolated())
    iso_points.push_back(Point_2(v_itr->point(), v_itr));
  }
  
  // Perform the sweep, while initializing it with all query points as event
  // points.
  Visitor     visitor (oi, arr);
  Sweep_line  sweep_line (&visitor);
  
  sweep_line.sweep(xcurves_vec.begin(),
		               xcurves_vec.end(),
                   iso_points.begin(),  
                   iso_points.end(), 
		               points_begin, 
		               points_end);
  
  return visitor.get_output_iterator();  // return a past_end iterator
}


CGAL_END_NAMESPACE

#endif
