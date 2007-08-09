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
// Author(s)     : Ron Wein          <wein@post.tau.ac.il>
//                 Baruch Zukerman   <baruchzu@post.tau.ac.il>
//                 Efi Fogel         <efif@post.tau.ac.il>
//                 Ophir Setter      <ophirset@post.tau.ac.il>
//
#ifndef CGAL_ARRANGEMENT_2_INSERT_H
#define CGAL_ARRANGEMENT_2_INSERT_H

/*! \file
 * Global insertion functions for the Arrangement_2 class.
 */

#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_accessor.h>
#include <CGAL/Arrangement_zone_2.h>
#include <CGAL/Arr_walk_along_line_point_location.h>
#include <CGAL/Arrangement_2/Arr_traits_adaptor_2.h>
#include <CGAL/Arrangement_2/Arr_inc_insertion_zone_visitor.h>
#include <CGAL/Arrangement_2/Arr_compute_zone_visitor.h>
#include <CGAL/Arrangement_2/Arr_do_intersect_zone_visitor.h>
#include <CGAL/Sweep_line_2/Arr_construction.h>
#include <CGAL/Sweep_line_2/Arr_addition.h>
#include <CGAL/Sweep_line_2/Arr_non_x_construction.h>
#include <CGAL/Sweep_line_2/Arr_non_x_addition.h>
#include <CGAL/Sweep_line_2/Sweep_line_2_visitors.h>
#include <CGAL/Sweep_line_2.h>
#include <CGAL/Arr_naive_point_location.h>

#include <boost/type_traits.hpp>

#include <set>
#include <list>
#include <map>

CGAL_BEGIN_NAMESPACE

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
template <class Traits, class Dcel, class PointLocation>
void insert_curve (Arrangement_2<Traits,Dcel>& arr, 
                   const typename Traits::Curve_2& c,
                   const PointLocation& pl, boost::is_same<int, double>::type)
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
  typedef Arr_traits_adaptor_2<Traits>                   Traits_adaptor_2;

  Traits_adaptor_2   *traits =
                        static_cast<Traits_adaptor_2*> (arr.get_traits());

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
      insert_point (arr, *iso_p, pl);
    }
  }

  return;
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
template <class Traits, class Dcel, class PointLocation>
void insert_curve (Arrangement_2<Traits,Dcel>& arr,
                   const typename Traits::X_monotone_curve_2& c,
                   const PointLocation& pl, 
                   boost::is_same<int, int>::type)
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
// Common interface for the insert_curve of the Curve_2 and X_monotone_curve_2
template <class Traits, class Dcel, class Curve, class PointLocation>
void insert_curve (Arrangement_2<Traits,Dcel>& arr,
                   const Curve& c, const PointLocation& pl)
{
  typedef typename Traits::X_monotone_curve_2       X_monotone_curve_2;
  
  typedef typename boost::is_same<Curve, X_monotone_curve_2>::type
    Is_x_monotone;
  
  return insert_curve(arr, c, pl, Is_x_monotone());
}

//-----------------------------------------------------------------------------
// Insert a curve/x-monotone curve into the arrangement (incremental 
// insertion).
// The inserted x-monotone curve may intersect the existing arrangement.
// Overloaded version with no point location object - the walk point-location
// strategy is used as default.
//
template <class Traits, class Dcel, class Curve>
void insert_curve (Arrangement_2<Traits,Dcel>& arr,
                   const Curve& c)
{
  typedef Arrangement_2<Traits, Dcel>                          Arrangement_2;
  typedef Arr_walk_along_line_point_location<Arrangement_2>    Walk_pl;
  
  // create walk point location object
  Walk_pl    walk_pl(arr);

  //insert the curve using the walk point location
  insert_curve (arr, c, walk_pl);
  return;
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
template <class Traits, class Dcel, class InputIterator>
void insert_curves (Arrangement_2<Traits,Dcel>& arr,
                    InputIterator begin, InputIterator end,
                    boost::is_same<int, double>::type)
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
// Insert a range of x-monotone curves into the arrangement (aggregated
// insertion). The inserted x-monotone curves may intersect one another and
// may also intersect the existing arrangement.
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
template <class Traits, class Dcel, class InputIterator>
void insert_curves (Arrangement_2<Traits,Dcel>& arr,
                    InputIterator begin, InputIterator end,
                    boost::is_same<int, int>::type)
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
// Common interface for the insert_curves of the Curve_2 and X_monotone_curve_2
template <class Traits, class Dcel, class InputIterator>
void insert_curves (Arrangement_2<Traits,Dcel>& arr,
                    InputIterator begin, InputIterator end)
{
  typedef typename Traits::X_monotone_curve_2       X_monotone_curve_2;
  typedef typename std::iterator_traits<InputIterator>::value_type 
    Iterator_value_type;

  typedef typename boost::is_same<Iterator_value_type, 
    X_monotone_curve_2>::type    Is_x_monotone;
  
  return insert_curves (arr, begin, end, Is_x_monotone());
}

//-----------------------------------------------------------------------------
// Insert an x-monotone curve into the arrangement (incremental insertion)
// when the location of the left endpoint of the curve is known and is
// given as an isertion hint.
// The inserted x-monotone curve may intersect the existing arrangement.
//
template <class Traits, class Dcel>
void insert_curve (Arrangement_2<Traits,Dcel>& arr,
                   const typename Traits::X_monotone_curve_2& c,
                   const Object& obj)
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
  arr_zone.init_with_hint (c, obj);

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

// ----------------------------------------------------------------------------
// backward compatibility functions.
template <class Traits, class Dcel, class PointLocation>
void insert_x_monotone_curve (Arrangement_2<Traits,Dcel>& arr,
                              const typename Traits::X_monotone_curve_2& c,
                              const PointLocation& pl)
{
  insert_curve(arr, c, pl);
}
template <class Traits, class Dcel>
void insert_x_monotone_curve (Arrangement_2<Traits,Dcel>& arr,
                              const typename Traits::X_monotone_curve_2& c)
{
  insert_curve(arr, c);
}
template <class Traits, class Dcel, class InputIterator>
void insert_x_monotone_curves (Arrangement_2<Traits,Dcel>& arr,
                               InputIterator begin, InputIterator end)
{
  insert_curves(arr, begin, end);
}
template <class Traits, class Dcel>
void insert_x_monotone_curve (Arrangement_2<Traits,Dcel>& arr,
                              const typename Traits::X_monotone_curve_2& c,
                              const Object& obj)
{
  insert_curve(arr, c, obj);
}

//-----------------------------------------------------------------------------
// Insert an x-monotone curve into the arrangement, such that the curve
// interior does not intersect with any existing edge or vertex in the
// arragement (incremental insertion).
//
template <class Traits, class Dcel, class PointLocation>
typename Arrangement_2<Traits,Dcel>::Halfedge_handle
insert_non_intersecting_curve (Arrangement_2<Traits,Dcel>& arr,
                               const typename Traits::X_monotone_curve_2& c,
                               const PointLocation& pl)
{
  typedef Arrangement_2<Traits,Dcel>                     Arrangement_2;
  typedef Arr_traits_basic_adaptor_2<typename Arrangement_2::Traits_2>
                                                         Traits_adaptor_2;
  typedef typename Arrangement_2::Vertex_const_handle    Vertex_const_handle;
  typedef typename Arrangement_2::Halfedge_const_handle  Halfedge_const_handle;

  Arr_accessor<Arrangement_2>       arr_access (arr);
  const Traits_adaptor_2           *traits =
    static_cast<const Traits_adaptor_2*> (arr.get_traits());

  // Check whether the left end of c lies at infinity, or whether it is a
  // normal endpoint, and locate it in the arrangement accordingly.
  const Boundary_type  inf_x1 = traits->boundary_in_x_2_object() (c, MIN_END);
  const Boundary_type  inf_y1 = traits->boundary_in_y_2_object() (c, MIN_END);
  CGAL::Object         obj1;
  const Halfedge_const_handle *fict_hh1 = NULL;
  const Vertex_const_handle   *vh1 = NULL;

  if (inf_x1 == NO_BOUNDARY && inf_y1 == NO_BOUNDARY)
  {
    // We have a normal left endpoint.
    obj1 = pl.locate (arr.get_traits()->construct_min_vertex_2_object() (c));

    // The endpoint must not lie on an existing edge, but may coincide with
    // and existing vertex vh1.
    CGAL_precondition_msg
      (object_cast<Halfedge_const_handle> (&obj1) == NULL,
       "The curve must not intersect an existing edge.");

    vh1 = object_cast<Vertex_const_handle> (&obj1);
  }
  else
  {
    // We have an unbounded left end.
    obj1 = arr_access.locate_unbounded_end (c, MIN_END);
    
    // The unbounded end should lie on a fictitious edge.
    CGAL_precondition_msg
        (object_cast<Vertex_const_handle> (&obj1) == NULL,
         "The curve must not overlap an existing edge.");
  
    fict_hh1 = object_cast<Halfedge_const_handle> (&obj1);
    CGAL_assertion (fict_hh1 != NULL);
  }

  // Check whether the right end of c lies at infinity, or whether it is a
  // normal endpoint, and locate it in the arrangement accordingly.
  const Boundary_type  inf_x2 = traits->boundary_in_x_2_object() (c, MAX_END);
  const Boundary_type  inf_y2 = traits->boundary_in_y_2_object() (c, MAX_END);
  CGAL::Object         obj2;
  const Halfedge_const_handle *fict_hh2 = NULL;
  const Vertex_const_handle   *vh2 = NULL;

  if (inf_x2 == NO_BOUNDARY && inf_y2 == NO_BOUNDARY)
  {
    // We have a normal right endpoint.
    obj2 = pl.locate (arr.get_traits()->construct_max_vertex_2_object() (c));

    // The endpoint must not lie on an existing edge, but may coincide with
    // and existing vertex vh2.
    CGAL_precondition_msg
      (object_cast<Halfedge_const_handle> (&obj2) == NULL,
       "The curve must not intersect an existing edge.");

    vh2 = object_cast<Vertex_const_handle> (&obj2);
  }
  else
  {
    // We have an unbounded right end.
    obj2 = arr_access.locate_unbounded_end (c, MAX_END);

    // The unbounded end should lie on a fictitious edge.
    CGAL_precondition_msg
        (object_cast<Vertex_const_handle> (&obj2) == NULL,
         "The curve must not overlap an existing edge.");
  
    fict_hh2 = object_cast<Halfedge_const_handle> (&obj2);
    CGAL_assertion (fict_hh2 != NULL);
  }

  // Notify the arrangement observers that a global operation is about to 
  // take place.
  arr_access.notify_before_global_change();

  // Check whether the located features containing the curve endpoints
  // are vertices or faces, and use the proper specialized insertion function
  // accordingly.
  typename Arrangement_2::Halfedge_handle             res;

  if (fict_hh1 == NULL && fict_hh2 == NULL)
  {
    // Both endpoints are finite.
    if (vh1 != NULL)
    {
      if (vh2 != NULL)
      {
        // Both endpoints are associated with a existing vertices.
        // In this case insert_at_vertices() already returns a halfedge
        // directed from left to right.
        res = arr.insert_at_vertices (c,
                                      arr.non_const_handle (*vh1),
                                      arr.non_const_handle (*vh2));
      }
      else
      {
        // Only the left endpoint is associated with an existing vertex.
        // In this case insert_from_left_vertex() returns a halfedge directed
        // to the new vertex it creates, so it is already directed from left to
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
        // In this case insert_from_left_vertex() returns a halfedge directed
        // to the new vertex it creates, so it is directed from right to left
        // and we take its twin halfedge instead.
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

        CGAL_assertion_msg 
          (fh1 != NULL && fh2 != NULL && *fh1 == *fh2,
           "The curve intersects the interior of existing edges.");

        if (fh1 != NULL && fh2 != NULL && *fh1 == *fh2)
        {
          res = arr.insert_in_face_interior (c,
                                             arr.non_const_handle (*fh1));
        }
      }
    }
  }
  else if (fict_hh1 != NULL && fict_hh2 == NULL)
  {
    // The left end is unbounded and the right endpoint is bounded.
    if (vh2 != NULL)
    {
      res = arr.insert_from_right_vertex (c,
                                          arr.non_const_handle (*vh2));
      res = res->twin();
    }
    else
    {
      res = arr.insert_in_face_interior (c,
                                         arr.non_const_handle (*fict_hh1));
    }
  }
  else if (fict_hh1 == NULL && fict_hh2 != NULL)
  {
    // The left endpoint is bounded and the right end is unbounded.
    if (vh1 != NULL)
    {
      res = arr.insert_from_left_vertex (c,
                                         arr.non_const_handle (*vh1));
    }
    else
    {
      res = arr.insert_in_face_interior (c,
                                         arr.non_const_handle (*fict_hh2));
    }
  }
  else
  {
    // Both curve ends are unbounded. In this case the two fictitious
    // halfedges must belong to the same unbounded face.
    typename Arrangement_2::Face_const_handle  fh1 = (*fict_hh1)->face();
    typename Arrangement_2::Face_const_handle  fh2 = (*fict_hh2)->face();

    CGAL_assertion_msg 
      (fh1 == fh2,
       "The curve intersects the interior of existing edges.");
    
    if (fh1 == fh2)
    {
      res = arr.insert_in_face_interior (c,
                                         arr.non_const_handle (fh1));
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
insert_non_intersecting_curve (Arrangement_2<Traits,Dcel>& arr,
                               const typename Traits::X_monotone_curve_2& c)
{
  typedef Arrangement_2<Traits, Dcel>                          Arrangement_2;
  typedef Arr_walk_along_line_point_location<Arrangement_2>    Walk_pl;
  
  // create walk point location object
  Walk_pl    walk_pl(arr);

  return (insert_non_intersecting_curve (arr, c, walk_pl));
}

//-----------------------------------------------------------------------------
// Insert a range of pairwise interior-disjoint x-monotone curves into
// the arrangement, such that the curve interiors do not intersect with
// any existing edge or vertex in the arragement (aggregated insertion).
//
template <class Traits, class Dcel, class InputIterator>
void insert_non_intersecting_curves (Arrangement_2<Traits,Dcel>& arr,
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
  bool                                   is_removed[2];

  v_ends[0] = e->source();
  is_removed[0] = (v_ends[0]->is_at_infinity() || v_ends[0]->degree() == 1);
  v_ends[1] = e->target();
  is_removed[1] = (v_ends[1]->is_at_infinity() || v_ends[1]->degree() == 1);

  // Remove the edge from the arrangement.
  typename Arrangement_2::Face_handle    face = arr.remove_edge (e);

  // Examine the end-vertices: If a vertex has now two incident edges, and the
  // curves associated with these edges can be merged, merge the two edges and
  // remove the vertex.
  typedef Arr_traits_adaptor_2<Traits>                   Traits_adaptor_2;

  Traits_adaptor_2                *traits =
                         static_cast<Traits_adaptor_2*>(arr.get_traits());

  typename Arrangement_2::Halfedge_around_vertex_circulator  circ;
  typename Arrangement_2::Halfedge_handle                    e1, e2;
  int                                                        i;

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
insert_point (Arrangement_2<Traits,Dcel>& arr,
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
    vh_for_p = arr.insert_in_face_interior (p,
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
insert_point (Arrangement_2<Traits,Dcel>& arr,
              const typename Traits::Point_2& p)
{
  typedef Arrangement_2<Traits, Dcel>                          Arrangement_2;
  typedef Arr_walk_along_line_point_location<Arrangement_2>    Walk_pl;
  
  // create walk point location object
  Walk_pl    walk_pl(arr);

  return (insert_point (arr, p, walk_pl));
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
  typedef Arr_traits_adaptor_2<Traits>                   Traits_adaptor_2;

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
    Traits_adaptor_2                *traits =
                        static_cast<Traits_adaptor_2*>(arr.get_traits());

    typename Arrangement_2::Halfedge_around_vertex_circulator  circ;
    typename Arrangement_2::Halfedge_handle                    e1, e2;

    circ = v->incident_halfedges();
    e1 = circ;
    ++circ;
    e2 = circ;

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

//-----------------------------------------------------------------------------
// Check the validity of the arrangement. In particular, check that the
// edegs are disjoint-interior, and the holes are located in their proper
// position.
//
template <class Traits_, class Dcel_>
bool is_valid (const Arrangement_2<Traits_,Dcel_>& arr)
{
  // First use the internal validity check.
  if(!arr.is_valid())
    return (false);

  typedef Traits_                                       Traits_2;
  typedef Arrangement_2<Traits_,Dcel_>                  Arrangement_2;
  typedef typename Arrangement_2::X_monotone_curve_2    X_monotone_curve_2;
  
  // Define the sweep-line types:
  typedef Sweep_line_do_curves_x_visitor<Traits_2>      Visitor;
  typedef Sweep_line_2<Traits_2, Visitor>               Sweep_line_2;

  // Define the arrangement iterator and circulator types:
  typedef typename Arrangement_2::Edge_const_iterator   Edge_const_iterator;
  typedef typename Arrangement_2::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Arrangement_2::Inner_ccb_const_iterator  Hole_const_iterator;
  typedef typename Arrangement_2::Face_const_iterator   Face_const_iterator;
  typedef typename Arrangement_2::Face_const_handle     Face_const_handle;
  typedef typename Arrangement_2::Vertex_const_handle   Vertex_const_handle;
  typedef typename Arrangement_2::Isolated_vertex_const_iterator
                                              Isolated_vertex_const_iterator;
  typedef typename Arrangement_2::Ccb_halfedge_const_circulator 
                                                 Ccb_halfedge_const_circulator;
  typedef typename Arrangement_2::Halfedge_around_vertex_const_circulator 
                                       Halfedge_around_vertex_const_circulator;
  
  // Perform a sweep over all subcurves associated with arrangement edges.
  std::vector<X_monotone_curve_2>   curves_vec(arr.number_of_edges());
  Edge_const_iterator               eit;
  unsigned int                      i = 0;
  
  for (eit = arr.edges_begin(); eit != arr.edges_end(); ++eit, i++)
    curves_vec[i] = eit->curve();

  Visitor       visitor;
  Traits_2     *traits = const_cast<Traits_2 *> (arr.get_traits());
  Sweep_line_2  sweep_line (traits, &visitor);

  visitor.sweep_xcurves (curves_vec.begin(), curves_vec.end());
  
  bool          are_edges_disjoint = !visitor.found_x();

  if (!are_edges_disjoint)
  {
    CGAL_warning_msg (are_edges_disjoint,
                      "Edges are not disjoint in their interior.");
    return (false);
  }

  // Check that the holes and isolated vertices are located where they should.
  // At the same time, we prepare a vector that consists of all isolated
  // vertices and all leftmost vertices from every hole. 
  std::list<std::pair<Vertex_const_handle, Face_const_handle> >  vf_list;

  typename Traits_2::Compare_xy_2   compare_xy = traits->compare_xy_2_object();
  Face_const_iterator               fit;
  Face_const_handle                 fh;
  Hole_const_iterator              hoit;
  Halfedge_const_handle             ccb;
  Isolated_vertex_const_iterator  ivit;
  Vertex_const_handle               left_v;
  bool                              is_first;

  for (fit = arr.faces_begin(); fit != arr.faces_end(); ++fit)
  {
    // Check all holes in the current face.
    fh = fit;
    for(hoit = fh->holes_begin(); hoit != fh->holes_end(); ++hoit)
    {
      ccb = *hoit;
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

      } while (ccb != *hoit);

      vf_list.push_back(std::make_pair (left_v, fh));
    }

    // Check all isolated vertices in the current face.
    for(ivit = fh->isolated_vertices_begin();
        ivit != fh->isolated_vertices_end(); ++ivit)
    {
      if (ivit->face() != fit)
        return (false);

      vf_list.push_back (std::make_pair (Vertex_const_handle(ivit), fh));
    }
  }

  // Shoot a vertical ray from each vertex we have collected downward, and
  // check that this vertex is really contained in the proper face.
  typename Traits_2::Compare_y_at_x_right_2  comp_y_at_x_right =
                                      traits->compare_y_at_x_right_2_object();
  typename Traits_2::Compare_y_at_x_left_2   comp_y_at_x_left =
                                      traits->compare_y_at_x_left_2_object();

  typename std::list<std::pair<Vertex_const_handle,
                               Face_const_handle> >::iterator    vf_iter;
  Arr_naive_point_location<Arrangement_2>                        pl(arr);
  Vertex_const_handle                                            curr_v;
  Object                                                         obj;
  Halfedge_const_handle                                          he_below;
  Vertex_const_handle                                            v_below;
  Face_const_handle                                              in_face;
  Halfedge_around_vertex_const_circulator                        first, circ;
  bool                                                           assign_ok;

  for (vf_iter = vf_list.begin(); vf_iter != vf_list.end(); ++vf_iter)
  {
    // Perform ray-shooting from the current vertex.
    curr_v = vf_iter->first;
    obj = pl.ray_shoot_down(curr_v->point());

    if (CGAL::assign(he_below, obj))
    {
      // Hit an edge - take the incident face of the halfedge directed to the
      // right.
      if (he_below->direction() == LARGER)
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
        // Get the first halfedge aroung v_below that is directed from left to
        // right and the first halfedge that is directed from right to left.
        first = circ = v_below->incident_halfedges();
        Halfedge_const_handle he_left;  // A halfedge to the left of v_below.
        Halfedge_const_handle he_right; // A halfedge to the right of v_below.

        do
        {
          if (circ->direction() == SMALLER)
          {
            he_left = circ;
          }
          else
          {
            he_right = circ;
            if((he_left != Halfedge_const_handle()) &&
               (he_right != Halfedge_const_handle()))
              break;
          }
          ++circ;
        
        } while(circ != first);

        CGAL_assertion (he_left != Halfedge_const_handle() || 
                        he_right != Halfedge_const_handle()); 

        if ((he_left != Halfedge_const_handle()) &&
            (he_right != Halfedge_const_handle()))
        {
          while (he_left -> direction() == SMALLER)
          {
            he_left = he_left->next()->twin();
          }
          he_left = he_left->twin()->prev();
          CGAL_assertion(he_left->direction() == SMALLER);
          in_face = he_left->face();
        }
        else
        {
          if (he_left != Halfedge_const_handle())
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
                        "Found a hole that is located in the wrong face.");
      return (false);
    }
  }

  // If we reached here, the arrangement is valid:
  return true;
}



//-----------------------------------------------------------------------------
// Compute the zone of the given x-monotone curve in the existing arrangement.
// Meaning, it output the arrangment's vertices, edges and faces that the 
// x-monotone curve intersects.
template <class Traits, class Dcel, class OutputIterator, class PointLocation>
OutputIterator zone (Arrangement_2<Traits,Dcel>& arr, 
                     const typename Traits::X_monotone_curve_2& c,
                     OutputIterator oi,
                     const PointLocation& pl)
{
  // Obtain an arrangement accessor.
  typedef Arrangement_2<Traits,Dcel>                     Arrangement_2;

  // Define a zone-computation object an a visitor that performs the
  // intersection check.
  typedef Arr_compute_zone_visitor<Arrangement_2, OutputIterator>  
    Zone_visitor;
  
  Zone_visitor                                     visitor (oi);
  Arrangement_zone_2<Arrangement_2, Zone_visitor>  arr_zone (arr, &visitor);

  arr_zone.init (c, pl);
  arr_zone.compute_zone();

  return (oi);
}

//-----------------------------------------------------------------------------
// Compute the zone of the given x-monotone curve in the existing arrangement.b
// Overloaded version with no point location object - the walk point-location
// strategy is used as default.
//
template <class Traits, class Dcel, class OutputIterator>
OutputIterator zone (Arrangement_2<Traits,Dcel>& arr, 
                     const typename Traits::X_monotone_curve_2& c,
                     OutputIterator oi)
{
  typedef Arrangement_2<Traits, Dcel>                          Arrangement_2;
  typedef Arr_walk_along_line_point_location<Arrangement_2>    Walk_pl;
  
  // create walk point location object
  Walk_pl    walk_pl(arr);

  //insert the curve using the walk point location
  zone (arr, c, oi, walk_pl);
  return oi;
}


//-----------------------------------------------------------------------------
// Checks if the given x-monotone curve intersects the existing arrangement.
// The last parameter is used to resolve ambiguity between this function and 
// do_intersect of Curve_2 in case that X_monotone_curve_2 and Curve_2 are the 
// same class. The last parameter should be boost::true_type but we used a 
// workaround since it didn't compile in FC3_g++-3.4.4 with the error of:
//
// error: no matching function for call to `do_intersect(Arrangement_2<>&, 
// const Arr_segment_2&, const Arr_walk_along_line_point_location<>&, mpl_::bool_< true>)'
//
template <class Traits, class Dcel, class PointLocation>
bool do_intersect (Arrangement_2<Traits,Dcel>& arr, 
                   const typename Traits::X_monotone_curve_2& c,
                   const PointLocation& pl, boost::is_same<int, int>::type)
{
  // Obtain an arrangement accessor.
  typedef Arrangement_2<Traits,Dcel>                     Arrangement_2;

  // Define a zone-computation object an a visitor that performs the
  // intersection check.
  typedef Arr_do_intersect_zone_visitor<Arrangement_2>  Zone_visitor;
  
  Zone_visitor                                     visitor;
  Arrangement_zone_2<Arrangement_2, Zone_visitor>  arr_zone (arr, &visitor);

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
// error: no matching function for call to `do_intersect(Arrangement_2<>&, 
// const Arr_segment_2&, const Arr_walk_along_line_point_location<>&, mpl_::bool_< true>)'
//
template <class Traits, class Dcel, class PointLocation>
bool do_intersect (Arrangement_2<Traits,Dcel>& arr, 
                   const typename Traits::Curve_2& c,
                   const PointLocation& pl, boost::is_same<int, double>::type)
{
  // Obtain an arrangement accessor.
  typedef Arrangement_2<Traits,Dcel>                     Arrangement_2;

  // Break the input curve into x-monotone subcurves and isolated points.
  typedef Arr_traits_adaptor_2<Traits>                   Traits_adaptor_2;

  Traits_adaptor_2   *traits =
                        static_cast<Traits_adaptor_2*> (arr.get_traits());

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
      // Check if the x-monotone subcurve intersects the arrangement.
      if (do_intersect(arr, *x_curve, pl) == true)
        return true;
    }
    else
    {
      iso_p = object_cast<typename Traits::Point_2> (&(*obj_iter));
      CGAL_assertion (iso_p != NULL);

      // Check whether the isolated point lies inside a face (otherwise,
      // it conincides with a vertex or an edge).
      CGAL::Object  obj = pl.locate (*iso_p);

      return (object_cast<typename Arrangement_2::Face_const_handle>(&obj) !=
              NULL);
    }
  }

  // If we reached here, the curve does not intersect the arrangement.
  return (false);
}

//-----------------------------------------------------------------------------
// Common interface for the do_intersect of the Curve_2 and X_monotone_curve_2
template <class Traits, class Dcel, class Curve, class PointLocation>
bool do_intersect (Arrangement_2<Traits,Dcel>& arr, const Curve& c, 
                   const PointLocation& pl)
{
  typedef typename Traits::X_monotone_curve_2       X_monotone_curve_2;
  
  typedef typename boost::is_same<Curve, X_monotone_curve_2>::type
    Is_x_monotone;
  
  return do_intersect(arr, c, pl, Is_x_monotone());
}

//-----------------------------------------------------------------------------
// Checks if the given curve intersects the existing arrangement.
// Overloaded version with no point location object - the walk point-location
// strategy is used as default.
template <class Traits, class Dcel, class Curve>
bool do_intersect (Arrangement_2<Traits, Dcel>& arr, 
                   const Curve& c)
{
  typedef Arrangement_2<Traits, Dcel>                          Arrangement_2;
  typedef Arr_walk_along_line_point_location<Arrangement_2>    Walk_pl;
  
  // create walk point location object
  Walk_pl    walk_pl(arr);
  
  // check if the curve intersects the arrangement using the walk point 
  // location.
  return do_intersect (arr, c, walk_pl);
}


CGAL_END_NAMESPACE

#endif
