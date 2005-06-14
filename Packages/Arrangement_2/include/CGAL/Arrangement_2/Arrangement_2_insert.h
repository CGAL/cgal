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
// Author(s)     : Ron Wein          <wein@post.tau.ac.il>
//                 Baruch Zukerman   <baruchzu@post.tau.ac.il>
//                 (based on old version by Eyal Flato)
#ifndef CGAL_ARRANGEMENT_2_INSERT_H
#define CGAL_ARRANGEMENT_2_INSERT_H

/*! \file
 * Global insertion functions for the Arrangement_2 class.
 */

#include <CGAL/Arrangement_2/Arr_traits_wrapper_2.h>
#include <CGAL/Arrangement_2/Arr_accessor.h>
#include <CGAL/Arrangement_zone_2.h>
#include <CGAL/Arrangement_2/Arr_inc_insertion_zone_visitor.h>
#include <CGAL/Sweep_line_2/Arr_aggregate_insert.h>
#include <CGAL/Sweep_line_2/Arr_non_x_aggregate_insert.h>
#include <list>
#include <map>

CGAL_BEGIN_NAMESPACE

//-----------------------------------------------------------------------------
// Insert a curve into the arrangement (incremental insertion).
// The inserted x-monotone curve may intersect the existing arrangement.
//
template <class Traits, class Dcel, class PointLocation>
void insert (Arrangement_2<Traits,Dcel>& arr, 
	     const PointLocation& pl,
	     const typename Traits::Curve_2& c)
{
  // Obtain an arrangement accessor.
  typedef Arrangement_2<Traits,Dcel>                     Arrangement_2;

  Arr_accessor<Arrangement_2>                      arr_access (arr);

  // Define a zone-computation object an a visitor that performs the
  // incremental insertion.
  typedef Arr_inc_insertion_zone_visitor<Arrangement_2>  Zone_visitor;

  Zone_visitor                                     visitor;
  Arrangement_zone_2<Arrangement_2, Zone_visitor>  arr_zone (arr, &visitor);

  // Break the input curve into x-monotone subcurves.
  typedef Arr_traits_wrapper_2<Traits>                   Traits_wrapper_2;

  const Traits_wrapper_2   *traits =
                        static_cast<const Traits_wrapper_2*>(arr.get_traits());

  typedef std::list<typename Traits::X_monotone_curve_2> Curves_list;
  Curves_list                             x_curves;
  typename Curves_list::const_iterator    x_iter; 
  Object                                  obj;

  traits->make_x_monotone_2_object() (c,
                                      std::back_inserter (x_curves));

  // Insert each x-monotone curve into the arrangement.
  for (x_iter = x_curves.begin(); x_iter != x_curves.end(); ++x_iter)
  {
    // Initialize the zone-computation object with the given curve.
    arr_zone.init (*x_iter, pl);

    // Notify the arrangement observers that a global operation is about to 
    // take place.
    arr_access.notify_before_global_change();

    // Insert the current x-monotone curve into the arrangement.
    arr_zone.compute_zone();

    // Notify the arrangement observers that the global operation has been
    // completed.
    arr_access.notify_after_global_change();
  }

  return;
}

//-----------------------------------------------------------------------------
// Insert a range of curves into the arrangement (aggregated insertion). 
// The inserted curves may intersect one another and may also intersect the 
// existing arrangement.
//
template <class Traits, class Dcel, class InputIterator>
void insert (Arrangement_2<Traits,Dcel>& arr,
	     InputIterator begin, InputIterator end)
{
  // Notify the arrangement observers that a global operation is about to 
  // take place.
  typedef Arrangement_2<Traits,Dcel>                     Arrangement_2;

  Arr_accessor<Arrangement_2>                      arr_access (arr);

  arr_access.notify_before_global_change();

  // Perform the aggregated insertion.
  Arr_aggregate_insert<Arrangement_2>  agg_insert_obj (arr.get_traits(), &arr);
  agg_insert_obj.insert_curves (begin, end);

  // Notify the arrangement observers that the global operation has been
  // completed.
  arr_access.notify_after_global_change();

  return;
}

//-----------------------------------------------------------------------------
// Insert an x-monotone curve into the arrangement (incremental insertion).
// The inserted x-monotone curve may intersect the existing arrangement.
//
template <class Traits, class Dcel, class PointLocation>
void insert_x_monotone (Arrangement_2<Traits,Dcel>& arr,
			const PointLocation& pl,
			const typename Traits::X_monotone_curve_2& c)
{
  // Obtain an arrangement accessor.
  typedef Arrangement_2<Traits,Dcel>                     Arrangement_2;

  Arr_accessor<Arrangement_2>                      arr_access (arr);

  // Define a zone-computation object an a visitor that performs the
  // incremental insertion.
  typedef Arr_inc_insertion_zone_visitor<Arrangement_2>  Zone_visitor;
  Zone_visitor                                     visitor;
  Arrangement_zone_2<Arrangement_2, Zone_visitor>  arr_zone (arr, &visitor);

  // Initialize the zone-computation object with the given curve.
  arr_zone.init (*x_iter, pl);

  // Notify the arrangement observers that a global operation is about to 
  // take place.
  arr_access.notify_before_global_change();

  // Insert the x-monotone curve into the arrangement.
  arr_zone.compute_zone();

  // Notify the arrangement observers that the global operation has been
  // completed.
  arr_access.notify_after_global_change();

  return;
}

//-----------------------------------------------------------------------------
// Insert a range of x-monotone curves into the arrangement (aggregated
// insertion). The inserted x-monotone curves may intersect one another and
// may also intersect the existing arrangement.
//
template <class Traits, class Dcel, class InputIterator>
void insert_x_monotone (Arrangement_2<Traits,Dcel>& arr,
			InputIterator begin, InputIterator end)
{
  // Notify the arrangement observers that a global operation is about to 
  // take place.
  typedef Arrangement_2<Traits,Dcel>                     Arrangement_2;

  Arr_accessor<Arrangement_2>                      arr_access (arr);

  arr_access.notify_before_global_change();

  // Perform the aggregated insertion.
  Arr_aggregate_insert<Arrangement_2>  agg_insert_obj (arr.get_traits(), &arr);
  agg_insert_obj.insert_x_curves (begin, end);

  // Notify the arrangement observers that the global operation has been
  // completed.
  arr_access.notify_after_global_change();

  return;
}

//-----------------------------------------------------------------------------
// Insert an x-monotone curve into the arrangement, such that the curve
// interior does not intersect with any existing edge or vertex in the
// arragement (incremental insertion).
//
template <class Traits, class Dcel, class PointLocation>
typename Arrangement_2<Traits,Dcel>::Halfedge_handle
insert_non_intersecting (Arrangement_2<Traits,Dcel>& arr,
			 const PointLocation& pl,
			 const typename Traits::X_monotone_curve_2& c)
{
  // Locate the curve endpoints.
  Object   obj1 =
            pl.locate (arr.get_traits()->construct_min_vertex_2_object() (c));
  Object   obj2 =
            pl.locate (arr.get_traits()->construct_max_vertex_2_object() (c));

  // Notify the arrangement observers that a global operation is about to 
  // take place.
  typedef Arrangement_2<Traits,Dcel>                     Arrangement_2;

  Arr_accessor<Arrangement_2>                      arr_access (arr);

  arr_access.notify_before_global_change();

  // In principal, the result of the point-location queries should not be
  // halfedges, because this function does not allow the inserted curve to
  // intersect the interior of any existing halfedge.
  CGAL_precondition_code (
    typename Arrangement_2::Halfedge_const_handle  hh;
  );
  
  CGAL_precondition_msg (!(assign (hh, obj1)) && !(assign (hh, obj2)),
                      "The curve intersects the interior of existing edges.");
  
  // Check whether the located features containing the curve endpoints
  // are vertices or faces, and use the proper specialized insertion function
  // accordingly.
  typename Arrangement_2::Vertex_const_handle  vh1;
  typename Arrangement_2::Vertex_const_handle  vh2;
  typename Arrangement_2::Halfedge_handle      res;

  if (assign (vh1, obj1))
  {
    if (assign (vh2, obj2))
    {
      // Both endpoints are associated with a existing vertices.
      res = arr.insert_at_vertices (c,
                                    arr.non_const_handle (vh1),
                                    arr.non_const_handle (vh2));
    }
    else
    {
      // Only the left endpoint is associated with an existing vertex.
      res = arr.insert_from_left_vertex (c,
					 arr.non_const_handle (vh1));
    }
  }
  else
  {
    if (assign (vh2, obj2))
    {
      // Only the right endpoint is associated with an existing vertex.
      res = arr.insert_from_right_vertex (c,
					  arr.non_const_handle (vh2));
    }
    else
    {
      // Both endpoints are not associated with existing vertices, so
      // we must insert the curve in the interior of a face.
      typename Arrangement_2::Face_const_handle  fh1;
      typename Arrangement_2::Face_const_handle  fh2;

      bool    succ1 = assign (fh1, obj1); 
      bool    succ2 = assign (fh2, obj2);

      CGAL_precondition_msg (succ1 && succ2 && fh1 == fh2,
                      "The curve intersects the interior of existing edges.");

      if (succ1 && succ2 && fh1 == fh2)
      {
        res = arr.insert_in_face_interior (c,
                                           arr.non_const_handle (fh1));
      }
    }
  }

  // Notify the arrangement observers that the global operation has been
  // completed.
  arr_access.notify_after_global_change();

  // Return the resulting halfedge from the insertion operation.
  return (res);
}

//-----------------------------------------------------------------------------
// Insert a range of pairwise interior-disjoint x-monotone curves into
// the arrangement, such that the curve interiors do not intersect with
// any existing edge or vertex in the arragement (aggregated insertion).
//
template <class Traits, class Dcel, class InputIterator>
void insert_non_intersecting (Arrangement_2<Traits,Dcel>& arr,
			      InputIterator begin, InputIterator end)
{
  // Notify the arrangement observers that a global operation is about to 
  // take place.
  typedef Arrangement_2<Traits,Dcel>                     Arrangement_2;

  Arr_accessor<Arrangement_2>                      arr_access (arr);

  arr_access.notify_before_global_change();

  // Perform the aggregated insertion.
  Arr_non_x_aggregate_insert<Arrangement_2>  agg_insert_obj (arr.get_traits(), 
							     &arr);
  agg_insert_obj.insert_curves(begin, end);

  // Notify the arrangement observers that the global operation has been
  // completed.
  arr_access.notify_after_global_change();

  return;
}

//-----------------------------------------------------------------------------
// Remove an edge from the arrangement. In case it is possible to merge
// the edges incident to the end-vertices of the removed edge after its
// deletion, the function performs these merges as well.
//
template <class Traits, class Dcel>
typename Arrangement_2<Traits,Dcel>::Face_handle
remove_edge (Arrangement_2<Traits,Dcel>& arr,
	     typename Arrangement_2<Traits,Dcel>::Halfedge_handle e)
{
  // Notify the arrangement observers that a global operation is about to 
  // take place.
  typedef Arrangement_2<Traits,Dcel>                     Arrangement_2;

  Arr_accessor<Arrangement_2>                      arr_access (arr);

  arr_access.notify_before_global_change();

  // Keep track of the end-vertices of the edge we are about to remove.
  typename Arrangement_2::Vertex_handle  v_ends[2];

  v_ends[0] = e.source();
  v_ends[1] = e.target();

  // Remove the edge from the arrangement.
  typename Arrangement_2::Face_handle    face = arr.remove_edge (e);

  // Examine the end-vertices: If a vertex has now two incident edges, and the
  // curves associated with these edges can be merged, merge the two edges and
  // remove the vertex.
  typedef Arr_traits_wrapper_2<Traits>                   Traits_wrapper_2;

  const Traits_wrapper_2                *traits =
                        static_cast<const Traits_wrapper_2*>(arr.get_traits());

  typename Arrangement_2::Halfedge_around_vertex_circulator  circ;
  typename Arrangement_2::Halfedge_handle                    e1, e2;
  int                                                        i;

  for (i = 0; i < 2; i++)
  {
    if (v_ends[i].degree() == 2)
    {
      // Get the two edges incident to the end-vertex.
      circ = v_ends[i].incident_halfedges();
      e1 = *circ;
      ++circ;
      e2 = *circ;

      // Check if it is possible to merge the two edges.
      if (traits->are_mergeable_2_object() (e1.curve(), e2.curve()))
      {
        // Merge the two curves.
        typename Traits::X_monotone_curve_2   cv;
        traits->merge_2_object() (e1.curve(), e2.curve(),
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

CGAL_END_NAMESPACE

#endif
