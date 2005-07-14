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

#include <CGAL/Arrangement_2.h>
#include <CGAL/Arrangement_2/Arr_traits_wrapper_2.h>
#include <CGAL/Arrangement_2/Arr_accessor.h>
#include <CGAL/Arrangement_zone_2.h>
#include <CGAL/Arrangement_2/Arr_inc_insertion_zone_visitor.h>
#include <CGAL/Arr_walk_along_line_point_location.h>
#include <CGAL/Sweep_line_2/Arr_construction.h>
#include <CGAL/Sweep_line_2/Arr_addition.h>
#include <CGAL/Sweep_line_2/Arr_non_x_construction.h>
#include <CGAL/Sweep_line_2/Arr_non_x_addition.h>
#include <list>
#include <map>

CGAL_BEGIN_NAMESPACE

//-----------------------------------------------------------------------------
// Insert a curve into the arrangement (incremental insertion).
// The inserted x-monotone curve may intersect the existing arrangement.
//
template <class Traits, class Dcel, class PointLocation>
void insert (Arrangement_2<Traits,Dcel>& arr, 
             const typename Traits::Curve_2& c,
             const PointLocation& pl)
{
  // Obtain an arrangement accessor.
  typedef Arrangement_2<Traits,Dcel>                     Arrangement_2;

  Arr_accessor<Arrangement_2>                      arr_access (arr);

  // Define a zone-computation object an a visitor that performs the
  // incremental insertion.
  typedef Arr_inc_insertion_zone_visitor<Arrangement_2>  Zone_visitor;

  Zone_visitor                                     visitor;
  Arrangement_zone_2<Arrangement_2, Zone_visitor>  arr_zone (arr, &visitor);

  // Break the input curve into x-monotone subcurves and isolated points.
  typedef Arr_traits_wrapper_2<Traits>                   Traits_wrapper_2;

  Traits_wrapper_2   *traits =
                        static_cast<Traits_wrapper_2*> (arr.get_traits());

  std::list<CGAL::Object>                     x_objects;
  std::list<CGAL::Object>::const_iterator     obj_iter;
  const typename Traits::X_monotone_curve_2  *x_curve;
  const typename Traits::Point_2             *iso_p;

  traits->make_x_monotone_2_object() (c,
                                      std::back_inserter (x_objects));

  // Insert each x-monotone curve into the arrangement.
  for (obj_iter = x_objects.begin(); obj_iter != x_objects.end(); ++obj_iter)
  {
    // Act according to the type of the current object.
    x_curve = object_cast<typename Traits::X_monotone_curve_2> (&(*obj_iter));
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
      iso_p = object_cast<typename Traits::Point_2> (&(*obj_iter));
      CGAL_assertion (iso_p != NULL);

      // Inserting a point into the arrangement:
      insert_vertex (arr, *iso_p, pl);
    }
  }

  return;
}

//-----------------------------------------------------------------------------
// Insert a curve into the arrangement (incremental insertion).
// The inserted x-monotone curve may intersect the existing arrangement.
// Overloaded version with no point location object - the walk point-location
// strategy is used as default.
//
template <class Traits, class Dcel>
void insert (Arrangement_2<Traits,Dcel>& arr,
             const typename Traits::Curve_2& c)
{
  typedef Arrangement_2<Traits, Dcel>                          Arrangement_2;
  typedef Arr_walk_along_line_point_location<Arrangement_2>    Walk_pl;
  
  // create walk point location object
  Walk_pl    walk_pl(arr);

  //insert the curve using the walk point location
  insert(arr, c, walk_pl);
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

  if(arr.is_empty())
  {
    // Perform the aggregated insertion.
    Arr_construction<Arrangement_2>              arr_construct (arr);
    arr_construct.insert_curves (begin, end);
  }
  else
  {
    Arr_addition<Arrangement_2>                  arr_adder(arr);
    arr_adder.insert_curves (begin, end);
  }

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
                        const typename Traits::X_monotone_curve_2& c,
                        const PointLocation& pl)
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
  arr_zone.init (c, pl);

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
// Insert an x-monotone curve into the arrangement (incremental insertion).
// The inserted x-monotone curve may intersect the existing arrangement.
// Overloaded version with no point location object - the walk point-location
// strategy is used as default.
//
template <class Traits, class Dcel, class PointLocation>
void insert_x_monotone (Arrangement_2<Traits,Dcel>& arr,
                        const typename Traits::X_monotone_curve_2& c)
{
  typedef Arrangement_2<Traits, Dcel>                          Arrangement_2;
  typedef Arr_walk_along_line_point_location<Arrangement_2>    Walk_pl;
  
  // create walk point location object
  Walk_pl    walk_pl(arr);

  insert_x_monotone(arr, c, walk_pl);
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
  if(arr.is_empty())
  {
    // Perform the aggregated insertion.
    Arr_construction<Arrangement_2>              arr_construct (arr);
    arr_construct.insert_x_curves (begin, end);
  }
  else
  {
    Arr_addition<Arrangement_2>                  arr_adder(arr);
    arr_adder.insert_x_curves (begin, end);
  }
 
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
                         const typename Traits::X_monotone_curve_2& c,
                         const PointLocation& pl)
{
  // Locate the curve endpoints.
  CGAL::Object   obj1 =
    pl.locate (arr.get_traits()->construct_min_vertex_2_object() (c));
  CGAL::Object   obj2 =
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
  
  CGAL_precondition_msg
    (object_cast<typename Arrangement_2::Halfedge_const_handle>(&obj1) == NULL
     &&
     object_cast<typename Arrangement_2::Halfedge_const_handle>(&obj2) == NULL,
     "The curve intersects the interior of existing edges.");
  
  // Check whether the located features containing the curve endpoints
  // are vertices or faces, and use the proper specialized insertion function
  // accordingly.
  const typename Arrangement_2::Vertex_const_handle  *vh1;
  const typename Arrangement_2::Vertex_const_handle  *vh2;
  typename Arrangement_2::Halfedge_handle             res;

  vh1 = object_cast<typename Arrangement_2::Vertex_const_handle> (&obj1);
  vh2 = object_cast<typename Arrangement_2::Vertex_const_handle> (&obj2);

  if (vh1 != NULL)
  {
    if (vh2 != NULL)
    {
      // Both endpoints are associated with a existing vertices.
      // In this case insert_at_vertices() already returns a halfedge directed
      // from left to right.
      res = arr.insert_at_vertices (c,
                                    arr.non_const_handle (*vh1),
                                    arr.non_const_handle (*vh2));
    }
    else
    {
      // Only the left endpoint is associated with an existing vertex.
      // In this case insert_from_left_vertex() returns a halfedge directed to 
      // the new vertex it creates, so it is already directed from left to
      // right.
      res = arr.insert_from_left_vertex (c,
                                         arr.non_const_handle (*vh1));
    }
  }
  else
  {
    if (vh2 != NULL)
    {
      // Only the right endpoint is associated with an existing vertex.
      // In this case insert_from_left_vertex() returns a halfedge directed to
      // the new vertex it creates, so it is directed from right to left and
      // we take its twin halfedge instead.
      res = arr.insert_from_right_vertex (c,
                                          arr.non_const_handle (*vh2));
      res = res->twin();
    }
    else
    {
      // Both endpoints are not associated with existing vertices, so
      // we must insert the curve in the interior of a face.
      // In this case insert_in_face_interior() already returns a halfedge
      // directed from left to right.
      const typename Arrangement_2::Face_const_handle  *fh1;
      const typename Arrangement_2::Face_const_handle  *fh2;

      fh1 = object_cast<typename Arrangement_2::Face_const_handle> (&obj1);
      fh2 = object_cast<typename Arrangement_2::Face_const_handle> (&obj2);

      CGAL_assertion_msg (fh1 != NULL && fh2 != NULL && *fh1 == *fh2,
                      "The curve intersects the interior of existing edges.");

      if (fh1 != NULL && fh2 != NULL && *fh1 == *fh2)
      {
	res = arr.insert_in_face_interior (c,
					   arr.non_const_handle (*fh1));
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
// Insert an x-monotone curve into the arrangement, such that the curve
// interior does not intersect with any existing edge or vertex in the
// arragement (incremental insertion).
// Overloaded version with no point location object - the walk point-location
// strategy is used as default.
//
template <class Traits, class Dcel>
typename Arrangement_2<Traits,Dcel>::Halfedge_handle
insert_non_intersecting (Arrangement_2<Traits,Dcel>& arr,
                         const typename Traits::X_monotone_curve_2& c)
{
  typedef Arrangement_2<Traits, Dcel>                          Arrangement_2;
  typedef Arr_walk_along_line_point_location<Arrangement_2>    Walk_pl;
  
  // create walk point location object
  Walk_pl    walk_pl(arr);

  return (insert_non_intersecting(arr, c, walk_pl));
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
  if(arr.is_empty())
  {
    // Perform the aggregated insertion.
    Arr_non_x_construction<Arrangement_2>  non_x_construct (arr);
    non_x_construct.insert_curves(begin, end);
  }
  else
  {
    Arr_non_x_addition<Arrangement_2>      non_x_adder(arr);
    non_x_adder.insert_curves (begin, end);
  }
  
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

  v_ends[0] = e->source();
  v_ends[1] = e->target();

  // Remove the edge from the arrangement.
  typename Arrangement_2::Face_handle    face = arr.remove_edge (e);

  // Examine the end-vertices: If a vertex has now two incident edges, and the
  // curves associated with these edges can be merged, merge the two edges and
  // remove the vertex.
  typedef Arr_traits_wrapper_2<Traits>                   Traits_wrapper_2;

  Traits_wrapper_2                *traits =
                         static_cast<Traits_wrapper_2*>(arr.get_traits());

  typename Arrangement_2::Halfedge_around_vertex_circulator  circ;
  typename Arrangement_2::Halfedge_handle                    e1, e2;
  int                                                        i;

  for (i = 0; i < 2; i++)
  {
    if (v_ends[i]->degree() == 2)
    {
      // Get the two edges incident to the end-vertex.
      circ = v_ends[i]->incident_halfedges();
      e1 = circ->handle();
      ++circ;
      e2 = circ->handle();

      // Check if it is possible to merge the two edges.
      if (traits->are_mergeable_2_object() (e1->curve(), e2->curve()))
      {
        // Merge the two curves.
        typename Traits::X_monotone_curve_2   cv;
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
template <class Traits, class Dcel, class PointLocation>
typename Arrangement_2<Traits,Dcel>::Vertex_handle
insert_vertex (Arrangement_2<Traits,Dcel>& arr,
               const typename Traits::Point_2& p,
               const PointLocation& pl)
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

  if ((fh = object_cast<typename Arrangement_2::Face_const_handle>(&obj)) 
      != NULL) 
  {
    // p lies inside a face: Insert it as an isolated vertex it the interior of
    // this face.
    vh_for_p = arr.insert_isolated_vertex (p,
                                           arr.non_const_handle (*fh));
  }
  else if ((hh = 
	    object_cast<typename Arrangement_2::Halfedge_const_handle>(&obj)) 
	   != NULL) 
  {
    // p lies in the interior of an edge: Split this edge to create a new
    // vertex associated with p.
    typename Traits::X_monotone_curve_2                   sub_cv1, sub_cv2;
    typename Arrangement_2::Halfedge_handle  split_he;
   
    arr.get_traits()->split_2_object() ((*hh)->curve(), p,
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
    vh = object_cast<typename Arrangement_2::Vertex_const_handle>(&obj);
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
template <class Traits, class Dcel>
typename Arrangement_2<Traits,Dcel>::Vertex_handle
insert_vertex (Arrangement_2<Traits,Dcel>& arr,
               const typename Traits::Point_2& p)
{
  typedef Arrangement_2<Traits, Dcel>                          Arrangement_2;
  typedef Arr_walk_along_line_point_location<Arrangement_2>    Walk_pl;
  
  // create walk point location object
  Walk_pl    walk_pl(arr);

  return (insert_vertex(arr, p, walk_pl));
}

//-----------------------------------------------------------------------------
// Remove a vertex from the arrangement.
//
template <class Traits, class Dcel>
bool remove_vertex (Arrangement_2<Traits,Dcel>& arr,
                    typename Arrangement_2<Traits,Dcel>::Vertex_handle v)
{
  // Notify the arrangement observers that a global operation is about to 
  // take place.
  typedef Arrangement_2<Traits,Dcel>                     Arrangement_2;
  typedef Arr_traits_wrapper_2<Traits>                   Traits_wrapper_2;

  Arr_accessor<Arrangement_2>                      arr_access (arr);

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
    Traits_wrapper_2                *traits =
                        static_cast<Traits_wrapper_2*>(arr.get_traits());

    typename Arrangement_2::Halfedge_around_vertex_circulator  circ;
    typename Arrangement_2::Halfedge_handle                    e1, e2;

    circ = v->incident_halfedges();
    e1 = circ->handle();
    ++circ;
    e2 = circ->handle();

    if (traits->are_mergeable_2_object() (e1->curve(), e2->curve()))
    {
      // Merge the two curves.
      typename Traits::X_monotone_curve_2   cv;
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

CGAL_END_NAMESPACE

#endif
