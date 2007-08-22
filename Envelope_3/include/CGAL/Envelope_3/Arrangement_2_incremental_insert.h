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
// $URL$
// $Id$
// 
//
// Author(s)     : Michal Meyerovitch     <gorgymic@post.tau.ac.il>
//                 Baruch Zukerman        <baruchzu@post.tau.ac.il>

#ifndef CGAL_ARRANGEMENT_2_INCREMENTAL_INSERT_H
#define CGAL_ARRANGEMENT_2_INCREMENTAL_INSERT_H

/*! \file
 * Global incremental insertion function for the Arrangement_2 class
 * with Zone Visitor parameter.
 */

#include <CGAL/Arrangement_2.h>
#include <CGAL/Arrangement_2/Arr_traits_adaptor_2.h>
#include <CGAL/Arr_accessor.h>
#include <CGAL/Arrangement_zone_2.h>
#include <CGAL/Arr_topology_traits/Arr_planar_inc_insertion_zone_visitor.h>

#include <list>
#include <map>

CGAL_BEGIN_NAMESPACE


//-----------------------------------------------------------------------------
// Insert a curve into the arrangement (incremental insertion).
// The inserted x-monotone curve may intersect the existing arrangement.
// Use the given zone visitor.
//
template <class Traits, class Dcel, class PointLocation, class ZoneVisitor>
void insert_x_monotone_curve (Arrangement_2<Traits,Dcel>& arr, 
                   const typename Traits::X_monotone_curve_2& c,
                   const PointLocation& pl,
                   ZoneVisitor& visitor)
{


  // Obtain an arrangement accessor.
  typedef Arrangement_2<Traits,Dcel>                     Arrangement_2;

  Arr_accessor<Arrangement_2>                      arr_access (arr);
  Arrangement_zone_2<Arrangement_2, ZoneVisitor>   arr_zone (arr, &visitor);

  arr_zone.init (c, pl);

  // Notify the arrangement observers that a global operation is about to
  // take place.
  arr_access.notify_before_global_change();

  // Insert the current x-monotone curve into the arrangement.
  arr_zone.compute_zone();

  // Notify the arrangement observers that the global operation has been
  // completed.
  arr_access.notify_after_global_change();
                                                   
  //// Break the input curve into x-monotone subcurves and isolated points.
  //typedef Arr_traits_adaptor_2<Traits>                   Traits_wrapper_2;

  //Traits_wrapper_2   *traits =
  //                      static_cast<Traits_wrapper_2*>(arr.get_traits());

  //std::list<Object>                    x_objects;
  //std::list<Object>::const_iterator    obj_iter;
  //typename Traits::X_monotone_curve_2  x_curve;
  //typename Traits::Point_2             iso_p;
  //Object                               obj;
  //bool                                 assign_success;

 /* traits->make_x_monotone_2_object() (c,
                                      std::back_inserter (x_objects));*/

  //// Insert each x-monotone curve into the arrangement.
  //for (obj_iter = x_objects.begin(); obj_iter != x_objects.end(); ++obj_iter)
  //{
  //  // Act according to the type of the current object.
  //  if (assign (x_curve, *obj_iter))
  //  {
  //    // Inserting an x-monotone curve:
  //    // Initialize the zone-computation object with the given curve.
  //    arr_zone.init (x_curve, pl);

  //    // Notify the arrangement observers that a global operation is about to
  //    // take place.
  //    arr_access.notify_before_global_change();

  //    // Insert the current x-monotone curve into the arrangement.
  //    arr_zone.compute_zone();

  //    // Notify the arrangement observers that the global operation has been
  //    // completed.
  //    arr_access.notify_after_global_change();
  //  }
  //  else
  //  {
  //    assign_success = assign (iso_p, *obj_iter);

  //    CGAL_assertion (assign_success);
  //    if (! assign_success)
  //      continue;

  //    // Inserting a point into the arrangement:
  //    //insert_vertex (arr, iso_p, pl);
  //    // we use the version with the visitor
  //    insert_point(arr, iso_p, pl, visitor);
  //    
  //  }
  //}

  return;
}

//-----------------------------------------------------------------------------
// Insert a vertex that corresponds to a given point into the arrangement.
// The inserted point may lie on any existing arrangement feature.
// Use the given visitor to actually change the arrangement.
// The visitor's methods should return the Vertex_handle of the new vertex,
// if one was created, or an invalid handle, if a vertex wasn't created.
//
template <class Traits, class Dcel, class PointLocation, class Visitor>
typename Arrangement_2<Traits,Dcel>::Vertex_handle
insert_point (Arrangement_2<Traits,Dcel>& arr,
               const typename Traits::Point_2& p,
               const PointLocation& pl,
               Visitor& visitor)
{
  // Obtain an arrangement accessor.
  typedef Arrangement_2<Traits,Dcel>                     Arrangement_2;
  // Act according to the type of arrangement feature that contains the point.
  const typename Arrangement_2::Face_const_handle      *fh;
  const typename Arrangement_2::Halfedge_const_handle  *hh;
  const typename Arrangement_2::Vertex_const_handle    *vh;
  typename Arrangement_2::Vertex_handle                 vh_for_p;

  Arr_accessor<Arrangement_2>                           arr_access (arr);

  // Locate the given point in the arrangement.
  CGAL::Object        obj = pl.locate (p);

  // Notify the arrangement observers that a global operation is about to
  // take place.
  arr_access.notify_before_global_change();

  visitor.init(&arr);
  
  if ((fh = object_cast<typename Arrangement_2::Face_const_handle>(&obj))
      != NULL)
  {
    vh_for_p = visitor.found_point_in_face(p, arr.non_const_handle (*fh));
  }
  else if ((hh =
	    object_cast<typename Arrangement_2::Halfedge_const_handle>(&obj))
	   != NULL)
  {
    vh_for_p = visitor.found_point_on_edge(p , arr.non_const_handle (*hh));
  }
  else
  {
    // In this case p lies on an existing vertex, so we just update this
    // vertex.
    vh = object_cast<typename Arrangement_2::Vertex_const_handle>(&obj);
    CGAL_assertion (vh != NULL);
    vh_for_p = visitor.found_point_on_vertex(p, arr.non_const_handle (*vh));
  }

  // Notify the arrangement observers that the global operation has been
  // completed.
  arr_access.notify_after_global_change();


  // Return a handle for the vertex associated with p.
  return (vh_for_p);

}

CGAL_END_NAMESPACE

#endif
