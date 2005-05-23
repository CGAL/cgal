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
#include <CGAL/Sweep_line_2/Arr_aggregate_insert.h>
#include <CGAL/Sweep_line_2/Arr_non_x_aggregate_insert.h>
#include <list>
#include <map>

CGAL_BEGIN_NAMESPACE

/*!
 * \class Auxiliary class for inserting an intersecting x-monotone curve into
 * an arrangement.
 */
template <class Arrangement_>
class _Arr_x_monotone_curve_inserter
{
public:

  typedef Arrangement_                            Arrangement_2;
  typedef typename Arrangement_2::Traits_2        Traits_2;

  typedef typename Arrangement_2::Vertex_handle   Vertex_handle;
  typedef typename Arrangement_2::Halfedge_handle Halfedge_handle;
  typedef typename Arrangement_2::Face_handle     Face_handle;

  typedef typename Traits_2::Point_2              Point_2;
  typedef typename Traits_2::X_monotone_curve_2   X_monotone_curve_2;
                   
protected:

  typedef Arr_traits_wrapper_2<Traits_2>          Traits_wrapper_2;

  typedef std::pair<Point_2, unsigned int>        Intersect_point_2;
  typedef std::list<Object>                       Intersect_list;
  typedef std::map<const X_monotone_curve_2*,
                   Intersect_list>                Intersect_map;
  typedef typename Intersect_map::value_type      Intersect_map_entry;
  typedef typename Intersect_map::iterator        Intersect_map_iterator;

  // Data members:
  Arrangement_2&          arr;      // The associated arrangement.
  const Traits_wrapper_2 *traits;   // Its associated traits object.

  Intersect_map           inter_map;   // Stores all computed intersections.

  X_monotone_curve_2  cv;              // The current portion of the
                                       // inserted curve.
  Point_2             left_pt;         // Its current left endpoint.
  Point_2             right_pt;        // Its right endpoint.
  Vertex_handle       curr_v;          // The current incident vertex of cv.
  Halfedge_handle     prev_he;         // The predecessor of cv around this
                                       // vertex. Its incident face is the
                                       // next face that cv will penetrate.
  bool                prev_he_valid;   // Indicates whether prev_he is valid.
  
  Point_2             intersect_p;     // The next intersection point.
  unsigned int        ip_mult;         // Its multiplicity 
                                       // (0 in case of an overlap).
  bool                found_intersect; // Have we found an intersection
                                       // (or an overlap).
  X_monotone_curve_2  overlap_cv;      // The currently discovered overlap.
  bool                found_overlap;   // Have we found an overlap.
  Halfedge_handle     intersect_he;    // The halfedge that intersects cv
                                       // (or overlaps it).
  
  Halfedge_handle     inserted_he;     // The last inserted halfedge.

public:

  /*!
   * Constructor.
   * \param _arr The arrangement to insert curves into.
   */
  _Arr_x_monotone_curve_inserter (Arrangement_2& _arr) :
    arr (_arr)
  {
    traits = static_cast<const Traits_wrapper_2*> (arr.get_traits());
  }

  /*!
   * Insert an intersecting x-monotone curve into the arrangement.
   * \param _cv The x-monotone curve to be inserted.
   * \param obj An object containing the left endpoint of the curve.
   * \return A handle to the rightmost edge in the zone of the
   *         inserted curve.
   */
  Halfedge_handle insert (const X_monotone_curve_2& _cv,
			  Object obj)
  {
    // Set the curve and its endpoints.
    cv = _cv;
    left_pt = traits->construct_min_vertex_2_object() (cv);
    right_pt = traits->construct_max_vertex_2_object() (cv);

    // Initialize flags.
    bool    done = false;

    prev_he_valid = false;
    found_overlap = false;

    // Locate the arrangement feature containing the left endpoint of the
    // curve.
    Vertex_handle          source_v;
    X_monotone_curve_2     sub_cv1, sub_cv2;
    Halfedge_handle        split_he;
    X_monotone_curve_2     icv;
  
    typename Arrangement_2::Vertex_const_handle    vh;
    typename Arrangement_2::Halfedge_const_handle  hh;
    typename Arrangement_2::Face_const_handle      fh;
                       
    if (assign (vh, obj))
    {
      // The left endpoint coincides with an existing vertex. We simply set 
      // the current vertex appropriately.
      curr_v = arr.non_const_handle (vh);
    }
    else if (assign (hh, obj))
    {
      // The left endpoint lies in the interior of an existing halfedge.
      // In this case we split the halfedge and set the source vertex as the
      // split point. We first split the curve associated with the halfedge.
      traits->split_2_object() (hh.curve(),
                                left_pt,
                                sub_cv1, sub_cv2);

      // Split the halfedge containing left_pt.
      split_he = arr.split_edge (arr.non_const_handle (hh),
                                 sub_cv1, sub_cv2);

      // The vertex associated with left_pt is now the target of the returned
      // halfedge.
      curr_v = arr.non_const_handle (split_he.target());
    }
    else
    {
      // The left endpoint lies inside a face.
      assign (fh, obj);

      // Insert the curve at the interior of the face.
      // This function also sets the curr_v vertex.
      done = _insert_in_face_interior (arr.non_const_handle(fh));

      // In case we have just discovered an overlap, insert the overlapping
      // subcurve as well.
      if (! done && found_overlap)
        done = _insert_overlap();
    }

    // Insert the curve (or what is remaining of it) into the arrangement,
    // starting from the source vertex we have computed.
    while (! done)
    {
      // Check if we know the face the curve is going to penetrate now.
      if (prev_he_valid)
      {
        // prev_he is valid, and its target vertex that represents the left
        // endpoint of cv.
	source_v = prev_he.target();
      }
      else
      {
        // The vertex that representing the left endpoint of cv is curr_v.
        source_v = curr_v;

        // Locate the curve around the source vertex - that is, find a halfedge
        // prev_he such that cv should be placed between prev_he and its
        // current successor around the vertex, going in a clockwise order.
        found_overlap = _find_prev_around_vertex ();

        if (found_overlap)
        {
          // In this case cv overlaps the curve associated with intersect_he.
          // Compute the overlapping subcurve to the right of curr_v.
          obj = _compute_next_intersection (intersect_he);

          bool  assign_success = assign (overlap_cv, obj);

          CGAL_assertion (assign_success);
          if (assign_success)
          {
            // Remove the overlap from the map.
            _remove_next_intersection (intersect_he);

            // Insert the overlap and continue to the end of the loop.
            done = _insert_overlap ();
            continue;
          }
        }
        else
        {
          // If we have no overlap, we now know the identity of the next face:
          // It is the incident face of prev_he.
          prev_he_valid = true;
        }
      }

      // At this point we can isert the curve starting from the prev_he
      // inside its incident face.
      done = _insert_from_face_boundary (source_v);

      // In case we have just discovered an overlap, insert the overlapping
      // subcurve as well.
      if (! done && found_overlap)
      {
        // Remove the overlap from the map.
        _remove_next_intersection (intersect_he);

        done = _insert_overlap();
      }
    }

    // Clear the intersections map and the list of halfedges that originated
    // from cv.
    inter_map.clear();

    return (inserted_he);
  }

private:

  /*!
   * Find a face containing the insertec curve around the curr_v vertex.
   * Sets prev_he to be the predecessor of cv around the vertex.
   * In case an overlap occurs, sets intersect_he to be the overlapping edge.
   * \return (true) if cv overlaps with the curve associated with cv;
   *         (false) if there is no overlap.
   */
  bool _find_prev_around_vertex ()
  {
    // Go over the incident halfedges of v, going in a clockwise order.
    typename Arrangement_2::Halfedge_around_vertex_circulator he_first;
    typename Arrangement_2::Halfedge_around_vertex_circulator he_curr;
    bool                                                      cv_equals_curr;
    typename Arrangement_2::Halfedge_around_vertex_circulator he_next;
    bool                                                      cv_equals_next;
    bool                                                      is_between;

    he_first = curr_v.incident_halfedges();
    he_curr = he_first;
    he_next = he_curr;
    ++he_next;

    if (he_curr == he_next)
    {
      // In case there is just a single incident halfedge around curr_v,
      // we should insert cv right after this halfedge.
      prev_he = *he_curr;

      // Note that cv extends to the right of curr_v. In case the single
      // halfedge also extends to the right of curr_v (its source is to
      // the right), check if an overlap occurs.
      if ((traits->compare_xy_2_object()
                                      (curr_v.point(),
                                       prev_he.source().point()) == SMALLER) &&
          (traits->compare_y_at_x_right_2_object() (prev_he.curve(), cv,
                                                    curr_v.point()) == EQUAL))
      {
        // Mark that an overlap occurs:
        intersect_he = *he_curr;
        return (true);
      }

      // We have no overlap - mark that prev_he is valid.
      prev_he_valid = true;
      return (false);
    }

    // Find the face containing cv around the vertex.
    typename Traits_wrapper_2::Is_between_cw_2  is_between_cw =
                                             traits->is_between_cw_2_object();

    do
    {
      // Check if it is possible to insert cv in between the current curve
      // and the next curve, going in a clockwise direction around v.
      is_between = is_between_cw (cv,
                                  (*he_curr).curve(), (*he_next).curve(),
                                  curr_v.point(),
                                  cv_equals_curr, cv_equals_next);

      // Check the case of overlaps:
      if (cv_equals_curr)
      {
        // cv overlaps with the curve of he_curr:
        intersect_he = *he_curr;
        return (true);
      }
      else if (cv_equals_next)
      {
        // cv overlaps with the curve of he_next:
        intersect_he = *he_next;
        return (true);
      }

      if (is_between)
      {
        // We can conclude that cv should be placed between he_curr and
        // he_next (in a clockwise order), and no overlap occurs.
        prev_he = *he_curr;
        prev_he_valid = true;
        return (false);
      }

      // Proceed to the next halfedges around the vertex.
      ++he_curr;
      ++he_next;

    } while (he_curr != he_first);

    // We should never reach here:
    CGAL_assertion (false);
    prev_he_valid = false;
    return (false);
  }

  /*!
   * Get the next intersection of cv with the given halfedge.
   * \param he A handle to the halfedge.
   * \return An object representing the next intersection: Intersect_point_2
   *         in case of a simple intersection point, X_monotone_curve_2 in
   *         case of an overlap, and an empty object if there is no
   *         intersection.
   */
  Object _compute_next_intersection (Halfedge_handle he)
  {
    // Get a pointer to the curve associated with the halfedge.
    const X_monotone_curve_2  *p_curve = &(he.curve());

    // Try to locate the intersections with this curve in the intersections
    // map.
    Intersect_map_iterator     iter = inter_map.find (p_curve);

    if (iter != inter_map.end())
    {
      // The intersections with the curve have already been computed.
      // Retrieve the intersections list from the map.
      Intersect_list&          inter_list = iter->second;

      if (inter_list.empty())
        return Object();

      // Locate the first intersection that lies to the right of left_pt.
      Intersect_point_2        ip;
      bool                     valid_intersection;

      while (! inter_list.empty())
      {
        // Compare that current object with left_pt.
        if (assign (ip, inter_list.front()))
        {
          valid_intersection = 
            (traits->compare_xy_2_object() (ip.first, left_pt) == LARGER);
        }
        else
        {
          X_monotone_curve_2   icv;

          assign (icv, inter_list.front());
          valid_intersection = (traits->compare_xy_2_object() 
                                (traits->construct_min_vertex_2_object()(icv), 
                                 left_pt) != SMALLER);
        }

        if (valid_intersection)
          // Found an intersection to left_pt's right.
          return (inter_list.front());

        // Discard the current intersection, which lies to left_pt's left.
        inter_list.pop_front();
      }

      // If we reached here, the list of intersections is empty:
      return Object();
    }

    // The intersections with the curve have not been computed yet, so we
    // have to compute them now.
    Intersect_list           inter_list;

    traits->intersect_2_object() (cv, he.curve(),
                                  std::back_inserter(inter_list));


    // Discard all intersection lying to the left of left_pt.
    Intersect_point_2        ip;
    bool                     valid_intersection;

    while (! inter_list.empty())
    {
      // Compare that current object with left_pt.
      if (assign (ip, inter_list.front()))
      {
        valid_intersection = 
          (traits->compare_xy_2_object() (ip.first, left_pt) == LARGER);
      }
      else
      {
        X_monotone_curve_2   icv;
        
        assign (icv, inter_list.front());
        valid_intersection = (traits->compare_xy_2_object() 
                              (traits->construct_min_vertex_2_object() (icv), 
                               left_pt) != SMALLER);
      }

      if (valid_intersection)
        // Found an intersection to left_pt's right.
        break;

      // Discard the current intersection, which lies to left_pt's left.
      inter_list.pop_front();
    }

    // Insert the list of valid intersections into the map.
    Intersect_map_entry        entry (p_curve, inter_list);

    inter_map.insert (entry);

    // Return the first intersection object computed (may be empty).
    if (inter_list.empty())
      return Object();
    else
      return (inter_list.front());
  }

  /*!
   * Remove the next intersection of cv with the given halfedge from the map.
   * \param he A handle to the halfedge.
   * \pre The list of intersections with the curve of he has already been
   *      computed, and it is not empty.
   */
  void _remove_next_intersection (Halfedge_handle he)
  {
    // Get a pointer to the curve associated with the halfedge.
    const X_monotone_curve_2  *p_curve = &(he.curve());

    // Locate the intersections with this curve in the intersections map.
    Intersect_map_iterator     iter = inter_map.find (p_curve);

    CGAL_assertion (iter != inter_map.end());
    CGAL_assertion (! iter->second.empty());

    // Remove the first object in the list of intersections.
    iter->second.pop_front();
    return;
  }

  /*!
   * Compute the (lexicographically) leftmost intersection of the inserted
   * curve with the boundary of a given face in the arrangement.
   * The function computes sets intersect_p, intersect_he (or alternatively
   * overlap_cv and intersect_he) and set the flags found_intersect and
   * found_overlap accordingly.
   * \param face A handle to the face.
   */
  void _leftmost_intersection_with_face_boundary (Face_handle face)
  {
    // Mark that we have not found any intersection (or overlap) yet.
    found_intersect = false;
    found_overlap = false;

    // Go over the outer boundary of the face (if one exists), and try to
    // locate intersections of cv with the edges along the boundary.
    typename Traits_wrapper_2::Compare_xy_2            compare_xy =
                                      traits->compare_xy_2_object();
    typename Traits_wrapper_2::Is_in_x_range_2         is_in_x_range =
                                      traits->is_in_x_range_2_object();
    typename Traits_wrapper_2::Construct_min_vertex_2  min_vertex =
                                      traits->construct_min_vertex_2_object();

    typename Arrangement_2::Ccb_halfedge_circulator  he_first;
    typename Arrangement_2::Ccb_halfedge_circulator  he_curr;

    Object                  obj;
    X_monotone_curve_2      icv;
    Intersect_point_2       int_p;
    Point_2                 ip;

    if (! face.is_unbounded())
    {
      // Get circulators for the outer boundary of the face.
      he_first = face.outer_ccb();
      he_curr = he_first;

      do
      {
        // If we already have an intersection point, compare it to the
        // endpoints of the curve associated with the current halfedge,
        //  in order to filter unnecessary intersection computations.
        if (found_intersect &&
            compare_xy ((*he_curr).source().point(), intersect_p) == LARGER &&
            compare_xy ((*he_curr).target().point(), intersect_p) == LARGER)
        {
          // The current x-monotone curve lies entirely to the right of
          // ip_left, so its intersection with cv (if any) cannot lie to
          // the left of this point. We therefore do not need to compute
          // this intersection.
          ++he_curr;
          continue;
        }

        // Check whether the two curves overlap in their x-range (in order
        // to avoid unnecessary intersection computations).
        if (! is_in_x_range (cv, (*he_curr).curve()))
        {
          // In case there is no overlap, the two x-monotone curve obviously
          // do not intersect.
          ++he_curr;
          continue;
        }

        // Compute the next intersection of cv and the current halfedge.
        obj = _compute_next_intersection ((*he_curr));

        if (! obj.is_empty())
        {
          // We have found an intersection (either a simple point or an
          // overlapping x-monotone curve).
          if (assign (int_p, obj))
          {
            ip = int_p.first;

            // Found a simple intersection point. Check if it is the leftmost
            // intersection point so far.
            if (! found_intersect ||
                compare_xy (ip, intersect_p) == SMALLER)
            {
              // Store the leftmost intersection point and the halfedge handle.
              intersect_p = ip;
	      ip_mult = int_p.second;
              intersect_he = *he_curr;
              found_overlap = false;
            }
          }
          else
          {
            // We have located an overlapping curve. Assign ip as its left
            // endpoint.
            assign (icv, obj);
            ip = min_vertex (icv);

            // Check if this endpoint it is the leftmost intersection point so
            // far.
            if (! found_intersect ||
                compare_xy (ip, intersect_p) == SMALLER)
            {
              // Store the leftmost intersection point and the halfedge handle.
              intersect_p = ip;
	      ip_mult = 0;
              overlap_cv = icv;
              intersect_he = *he_curr;
              found_overlap = true;
            }
          }

          // Mark that we found an intersection.
          found_intersect = true;
        }

        // Move to the next edge along the outer boundary,
        ++he_curr;
      
      } while (he_curr != he_first);

    } // End: if (! face.is_unbounded())

    // Go over the boundary of the holes inside the face (if there exist any),
    // and try to locate intersections of cv with the edges along the boundary
    // of each hole.
    typename Arrangement_2::Holes_iterator   hole_it;
  
    for (hole_it = face.holes_begin(); hole_it != face.holes_end(); ++hole_it)
    {
      // Get circulators for the boundary of the current hole.
      he_first = *hole_it;
      he_curr = he_first;

      do
      {
        // If we already have an intersection point, compare it to the
        // endpoints of the curve associated with the current halfedge,
        //  in order to filter unnecessary intersection computations.
        if (found_intersect &&
            compare_xy ((*he_curr).source().point(), intersect_p) == LARGER &&
            compare_xy ((*he_curr).target().point(), intersect_p) == LARGER)
        {
          // The current x-monotone curve lies entirely to the right of
          // ip_left, so its intersection with cv (if any) cannot lie to
          // the left of this point. We therefore do not need to compute
          // this intersection.
          ++he_curr;
          continue;
        }

        // Check whether the two curves overlap in their x-range (in order
        // to avoid unnecessary intersection computations).
        if (! is_in_x_range (cv, (*he_curr).curve()))
        {
          // In case there is no overlap, the two x-monotone curve obviously
          // do not intersect.
          ++he_curr;
          continue;        
        }

        // Compute the next intersection of cv and the current halfedge.
        obj = _compute_next_intersection ((*he_curr));

        if (! obj.is_empty())
        {
          // We have found an intersection (either a simple point or an
          // overlapping x-monotone curve).
          if (assign (int_p, obj))
          {
            ip = int_p.first;

            // Found a simple intersection point. Check if it is the leftmost
            // intersection point so far.
            if (! found_intersect ||
                compare_xy (ip, intersect_p) == SMALLER)
            {
              // Store the leftmost intersection point and the halfedge
              // handle.
              intersect_p = ip;
	      ip_mult = int_p.second;
              intersect_he = *he_curr;
              found_overlap = false;
            }
          }
          else
          {
            // We have located an overlapping curve. Assign ip as its left
            // endpoint.
            assign (icv, obj);
            ip = min_vertex (icv);

            // Check if this endpoint it is the leftmost intersection point
            // so far.
            if (! found_intersect ||
                compare_xy (ip, intersect_p) == SMALLER)
            {
              // Store the leftmost intersection point and the halfedge
              // handle.
              intersect_p = ip;
	      ip_mult = 0;
              overlap_cv = icv;
              intersect_he = *he_curr;
              found_overlap = true;
            }
          }

          // Mark that we found an intersection.
          found_intersect = true;
        }

        // Move to the next edge along the outer boundary,
        ++he_curr;

      } while (he_curr != he_first);

    } // End: traversal of the holes inside the face.

    // Remove the next intersection associated with intersect_he, as we have
    // now reported it and do not want to encounter it again.
    if (found_intersect)
    {
      _remove_next_intersection (intersect_he);
    }

    return;
  }

  /*!
   * Get an arrangement vertex representing a given point on an existing edge.
   * This function sets curr_v to be this vertex.
   * \param pt The given point.
   * \param he A halfedge containing this point.
   * \return (true) if we had to split he in order to create a new vertex;
   *         (false) if the vertex had already existed in the arrangement.
   */
  bool _get_vertex_on_halfedge (const Point_2& pt,
                                Halfedge_handle he)
  {
    // Check if pt corresponds to one of he's end vertices.
    if (traits->equal_2_object() (pt, he.source().point()))
    {
      curr_v = he.source();
      return (false);
    }
    else if (traits->equal_2_object() (pt, he.target().point()))
    {
      curr_v = he.target();
      return (false);
    }

    // In this case, pt is contained in the interior of the curve associated
    // with the given halfedge. Split this curve at pt.
    X_monotone_curve_2   sub_cv1, sub_cv2;
    Halfedge_handle      split_he;

    traits->split_2_object() (he.curve(),
                              pt,
                              sub_cv1, sub_cv2);

    // Split the halfedge.
    split_he = arr.split_edge (he,
                               sub_cv1, sub_cv2);

    // The vertex associated with pt is now the target of the returned
    // halfedge.
    curr_v = split_he.target();
    return (true);
  }

  /*!
   * Locate the previous halfedge for the insertion of the given subcurve
   * around the vertex curr_v.
   * \param cv_ins The curve to be inserted, whose right endpoint coincides 
   *               with v.
   * \param prev Output: The previous halfedge around the vertex (for the
   *                     insertion of cv at v).
   * \param next Output: The next halfedge (for the insertion of the
   *                     continuation of the curve around v).
   * \pre There are exactly two halfedges incident to v.
   */
  void _find_prev_around_new_vertex (const X_monotone_curve_2& cv_ins,
				     Halfedge_handle& prev,
				     Halfedge_handle& next)
  {
    CGAL_assertion (traits->equal_2_object()
		    (traits->construct_max_vertex_2_object() (cv_ins),
		     curr_v.point()));

    // Get the two halfedges incident to the vertex.
    typename Arrangement_2::Halfedge_around_vertex_circulator  circ;
    Halfedge_handle                                            he1, he2;

    circ = curr_v.incident_halfedges();
    he1 = *circ;
    ++circ;
    he2 = *circ;

    CGAL_assertion (++circ == curr_v.incident_halfedges());

    // Check which halfedge is defined to the left of the vertex, and check
    // whether cv_ins lies above or below it.
    Comparison_result       he_res =
      traits->compare_xy_2_object() (he1.source().point(), curr_v.point());
    const Halfedge_handle&  he_left = (he_res == SMALLER) ? he1 : he2;
    const Halfedge_handle&  he_right = (he_res == SMALLER) ? he2 : he1;
    Comparison_result       pos_res = 
      traits->compare_y_position_2_object() (cv_ins, he_left.curve());

    if (pos_res == SMALLER)
    {
      // If cv_ins lies below the existing left halfegde, we should take the
      // right halfedge to be its predecessor and the left one to be its 
      // successor.
      prev = he_right;
      next = he_left;
    }
    else
    {
      CGAL_assertion (pos_res != EQUAL);

      // If cv_ins lies above the existing left halfegde, we should take the
      // left halfedge to be its predecessor and the right one to be its 
      // successor.
      prev = he_left;
      next = he_right;
    }
    
    return;
  }
  
  /*!
   * Insert an x-monotone curve whose left endpoint lies inside a given face
   * of the arrangement.
   * This function updates cv and its left endpoint and also sets the curr_v,
   * prev_he and prev_he_valid.
   * In case of overlaps, it sets also overlap_cv and intersect_he.
   * \param face The given face.
   * \return (true) if we are done with the inserion process;
   *         (false) if we still have a remaining portion of cv to continue
   *         with.
   */
  bool _insert_in_face_interior (Face_handle face)
  {
    // Find the first intersection of the curve with the face boundary.
    _leftmost_intersection_with_face_boundary (face);

    if (! found_intersect)
    {
      // In case we found no intersection with the face boundary, we can
      // insert the entire curve as a new hole inside the given face (as we
      // know its left endpoint is located in this face).
      inserted_he = arr.insert_in_face_interior (cv, face);

      // Mark that we are done with the insertion process.
      prev_he_valid = false;
      return (true);
    }

    // In this case found_intersect is true and intersect_he is the edge that
    // cv next intersects (or overlaps). If found_overlap is also true,
    // then overlap_cv is set and intersect_p is the left endpoint of the
    // overlapping subcurve. Otherwise, intersect_p is a simple intersection
    // point.
    bool                  done = false;
    X_monotone_curve_2    sub_cv1, sub_cv2;
    bool                  new_vertex_created;

    if (traits->equal_2_object() (intersect_p, right_pt))
    {
      // If the intersection point is cv's right endpoint, the interior of cv
      // does not intersect any existing halfedge. In this case, we only have
      // to insert cv to the arrangement and we are done.
      sub_cv1 = cv;
      done = true;
    }
    else
    {
      // Split cv at the intersection point.
      traits->split_2_object() (cv,
                                intersect_p,
                                sub_cv1, sub_cv2);

      // Set cv to be the remaining portion.
      left_pt = intersect_p;
      cv = sub_cv2;
    }

    // Compute an arrangement vertex representing the intersection point (and
    // set curr_v to be this vertex).
    new_vertex_created = _get_vertex_on_halfedge (intersect_p,
                                                  intersect_he);

    if (new_vertex_created)
    {
      // In case a new vertex has been created, we can obtain the next face
      // that cv penetrates by comparing it to the two existing halfedges
      // around curr_v (which we have just split).
      Halfedge_handle    prev_he_around_v, next_he_around_v;

      _find_prev_around_new_vertex (sub_cv1, 
				    prev_he_around_v, next_he_around_v);

      inserted_he = arr.insert_from_vertex (sub_cv1, prev_he_around_v);

      // In case cv crosses the halfedge (the multiplicity of the intersection
      // is odd), assign a valid prev_he handle.
      if ((ip_mult % 2) == 1)
      {
	prev_he = next_he_around_v;
	prev_he_valid = true;
      }
      else
      {
	prev_he_valid = false;
      }
    }
    else
    {
      // We intersect an existing vertex and will have to locate cv around it.
      inserted_he = arr.insert_from_vertex (sub_cv1, curr_v);

      prev_he_valid = false;
    }

    return (done);
  }

  /*!
   * Insert an x-monotone curve whose left endpoint lies on the boundary of a
   * given face of the arrangement - the incident face of prev_he.
   * This function updates cv and its left endpoint and also sets the curr_v,
   * prev_he and prev_he_valid.
   * In case of overlaps, it sets also overlap_cv and intersect_he.
   * \param source_v The vertex representing the current left endpoint of cv.
   * \return (true) if we are done with the inserion process;
   *         (false) if we still have a remaining portion of cv to continue
   *         with.
   */
  bool _insert_from_face_boundary (Vertex_handle source_v)
  {
    // The halfedge prev_he is directed from the vertex representing the
    // left endpoint of cv, and its incident face is the one cv now
    // penetrates.
    Face_handle      face = prev_he.face();

    // Find the first intersection of the curve with the face boundary.
    _leftmost_intersection_with_face_boundary (face);

    if (! found_intersect)
    {
      // In case we found no intersection with the face boundary, we can
      // insert the curve from the source vertex, right after prev_he
      inserted_he = arr.insert_from_vertex (cv, prev_he);

      // Mark that we are done with the insertion process.
      prev_he_valid = false;
      return (true);
    }

    // In this case found_intersect is true and intersect_he is the edge that
    // cv next intersects (or overlaps). If found_overlap is also true,
    // then overlap_cv is set and intersect_p is the left endpoint of the
    // overlapping subcurve. Otherwise, intersect_p is a simple intersection
    // point.
    bool                  done = false;
    X_monotone_curve_2    sub_cv1, sub_cv2;
    bool                  new_vertex_created;

    if (traits->equal_2_object() (intersect_p, right_pt))
    {
      // If the intersection point is cv's right endpoint, the interior of cv
      // does not intersect any existing halfedge. In this case, we only have
      // to insert cv to the arrangement and we are done.
      sub_cv1 = cv;
      done = true;
    }
    else
    {
      // Split cv at the intersection point.
      traits->split_2_object() (cv,
                                intersect_p,
                                sub_cv1, sub_cv2);

      // Set cv to be the remaining portion.
      left_pt = intersect_p;
      cv = sub_cv2;
    }

    // Compute an arrangement vertex representing the intersection point (and
    // set curr_v to be this vertex).
    new_vertex_created = _get_vertex_on_halfedge (intersect_p,
                                                  intersect_he);

    // In case a new vertex has been created, we can obtain the next face that
    // cv penetrates by taking the twin halfedge of the halfedge it intersects.
    if (new_vertex_created)
    {
      // Check if we have just split the halfedge that prev_he refers to.
      // If so, prev_he's target is now the new vertex, and we have to proceed
      // to the next halfedge (whose target is source_v).
      if (intersect_he == prev_he)
	prev_he = prev_he.next();

      // In case a new vertex has been created, we can obtain the next face
      // that cv penetrates by comparing it to the two existing halfedges
      // around curr_v (which we have just split).
      Halfedge_handle    prev_he_around_v, next_he_around_v;

      _find_prev_around_new_vertex (sub_cv1,
				    prev_he_around_v, next_he_around_v);

      // The predecessor halfedge of the left end-vertex of sub_cv1 is prev_he,
      // while the right predecessor is prev_he_around_cv:
      inserted_he = arr.insert_at_vertices (sub_cv1, 
					    prev_he, prev_he_around_v);

      // In case cv crosses the halfedge (the multiplicity of the intersection
      // is odd), assign a valid prev_he handle.
      if ((ip_mult % 2) == 1)
      {
	prev_he = next_he_around_v;
	prev_he_valid = true;
      }
      else
      {
	prev_he_valid = false;
      }
    }
    else
    {
      // Insert the left sub_cv1 between the source vertex (reprsenting
      // the left endpoint) and the vertex we have just obtained. We will
      // later have to locate the face cv penetrates next around this existing
      // vertex.
      inserted_he = arr.insert_at_vertices (sub_cv1, source_v, curr_v);
      prev_he_valid = false;
    }

    return (done);
  }

  /*!
   * Insert an overlapping subcurve overlap_cv of cv and the curve currently
   * associated with intersect_he.
   * This function updates cv and its left endpoint and also sets the curr_v,
   * while invalidating prev_he.
   * \return (true) if we are done with the inserion process;
   *         (false) if we still have a remaining portion of cv to continue
   *         with.
   */
  bool _insert_overlap ()
  {
    // Get the right endpoints of overlap_cv and the right end vertex of
    // the overlapping halfedge intersect_he.
    Point_2  p_right_cv = traits->construct_max_vertex_2_object() (overlap_cv);
    Vertex_handle v_right;

    if (traits->compare_xy_2_object()
        (intersect_he.source().point(),
         intersect_he.target().point()) == SMALLER)
    {
      v_right = intersect_he.target();
    }
    else
    {
      v_right = intersect_he.source();
    }

    // Compare the two right endpoints. Note that overlap_cv cannot extend to
    // the right longer than the halfedge it overlaps.
    const Point_2&         p_right_he = v_right.point();
    X_monotone_curve_2     sub_cv1, sub_cv2;
    Halfedge_handle        split_he;

    Comparison_result      res = traits->compare_xy_2_object() (p_right_cv,
                                                                p_right_he);

    CGAL_assertion (res != LARGER);

    if (res == EQUAL)
    {
      // The overlap is with the entire halfedge. In this case we replace the
      // curve currently assoicated with intersect_he by icv.
      arr.modify_edge (intersect_he, overlap_cv);

      // Move the current incident vertex of cv to the right endpoint of the
      // overlapping subcurve.
      curr_v = v_right;
    }
    else
    {
      // In this case intersect_he overlaps just a portion of prev_he. We split
      // the curve associated with this halfedge at the right endpoint of
      // overlap_cv.
      traits->split_2_object() (intersect_he.curve(),
                                p_right_cv,
                                sub_cv1, sub_cv2);

      // Split the halfedge. The left part is associated with the overlapping
      // subcurve, and its right part with the portion split from the original
      // curve.
      split_he = arr.split_edge (intersect_he,
                                 overlap_cv, sub_cv2);

      // The vertex associated with left_pt is now the target of the returned
      // halfedge.
      curr_v = arr.non_const_handle (split_he.target());
    }

    // Mark that we have dealt with the overlap.
    found_overlap = false;

    // Mark that prev_he is no longer valid, and that we have to locate the
    // next incident face of cv around curr_v.
    prev_he_valid = false;

    // Compare the right endpoint of cv to the right endpoint of cv.
    if (traits->equal_2_object() (p_right_cv, right_pt))
    {
      // We reached cv's right endpoint, so we are done with the insertion
      // process.
      return (true);
    }

    // Split cv at right endpoint of the overlapping curve.
    traits->split_2_object() (cv,
                              p_right_cv,
                              sub_cv1, sub_cv2);

    // Set cv to be the remaining portion.
    left_pt = p_right_cv;
    cv = sub_cv2;

    // We are still not done with the insertion process:
    return (false);
  }
};

//-----------------------------------------------------------------------------
// Global functions for the various insertion operations:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Insert a curve into the arrangement (incremental insertion).
// The inserted x-monotone curve may intersect the existing arrangement.
//
template <class Arrangement, class PointLocation>
void arr_insert (Arrangement& arr, const PointLocation& pl,
                 const typename Arrangement::Traits_2::Curve_2& c)
{
  // Obtain an arrangement accessor.
  Arr_accessor<Arrangement>                           arr_access (arr);

  // Create an auxiliary inserter object.
  _Arr_x_monotone_curve_inserter<Arrangement>         arr_inserter (arr);

  // Break the input curve into x-monotone subcurves.
  typedef Arr_traits_wrapper_2<typename Arrangement::Traits_2>  
                                                              Traits_wrapper_2;

  const Traits_wrapper_2   *traits =
                        static_cast<const Traits_wrapper_2*>(arr.get_traits());

  typedef std::list<typename Arrangement::X_monotone_curve_2> Curves_list;
  Curves_list                             x_curves;
  typename Curves_list::const_iterator    x_iter; 
  Object                                  obj;

  traits->make_x_monotone_2_object() (c,
                                      std::back_inserter (x_curves));

  // Insert each x-monotone curve into the arrangement.
  for (x_iter = x_curves.begin(); x_iter != x_curves.end(); ++x_iter)
  {
    // Locate the left endpoint of the current x-monotone curve.
    obj = pl.locate (traits->construct_min_vertex_2_object() (*x_iter));

    // Notify the arrangement observers that a global operation is about to 
    // take place.
    arr_access.notify_before_global_change();

    // Insert the current x-monotone curve into the arrangement.
    arr_inserter.insert (*x_iter, obj);

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
template <class Arrangement, class InputIterator>
void arr_insert (Arrangement& arr,
                 InputIterator begin, InputIterator end)
{
  // Notify the arrangement observers that a global operation is about to 
  // take place.
  Arr_accessor<Arrangement>                             arr_access (arr);

  arr_access.notify_before_global_change();

  // Perform the aggregated insertion.
  Arr_aggregate_insert<Arrangement>  agg_insert_obj (arr.get_traits(), &arr);
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
template <class Arrangement, class PointLocation>
void arr_insert_x_monotone (Arrangement& arr, const PointLocation& pl,
			    const typename Arrangement::X_monotone_curve_2& c)
{
  // Obtain an arrangement accessor.
  Arr_accessor<Arrangement>                           arr_access (arr);

  // Create an auxiliary inserter object and insert the x-monotone curve.
  _Arr_x_monotone_curve_inserter<Arrangement>         arr_inserter (arr);

  // Locate the left endpoint of the x-monotone curve.
  Object   obj = pl.locate (traits->construct_min_vertex_2_object() (*x_iter));

  // Notify the arrangement observers that a global operation is about to 
  // take place.
  arr_access.notify_before_global_change();

  // Insert the x-monotone curve into the arrangement.
  arr_inserter.insert (c, obj);

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
template <class Arrangement, class InputIterator>
void arr_insert_x_monotone (Arrangement& arr,
			    InputIterator  begin , InputIterator  end )
{
  // Notify the arrangement observers that a global operation is about to 
  // take place.
  Arr_accessor<Arrangement>                             arr_access (arr);

  arr_access.notify_before_global_change();

  // Perform the aggregated insertion.
  Arr_aggregate_insert<Arrangement>  agg_insert_obj (arr.get_traits(), &arr);
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
template <class Arrangement, class PointLocation>
typename Arrangement::Halfedge_handle
arr_insert_non_intersecting
                (Arrangement& arr, const PointLocation& pl,
                 const typename Arrangement::X_monotone_curve_2& c)
{
  // Locate the curve endpoints.
  Object   obj1 =
            pl.locate (arr.get_traits()->construct_min_vertex_2_object() (c));
  Object   obj2 =
            pl.locate (arr.get_traits()->construct_max_vertex_2_object() (c));

  // Notify the arrangement observers that a global operation is about to 
  // take place.
  Arr_accessor<Arrangement>                             arr_access (arr);

  arr_access.notify_before_global_change();

  // In principal, the result of the point-location queries should not be
  // halfedges, because this function does not allow the inserted curve to
  // intersect the interior of any existing halfedge.
  CGAL_precondition_code (
    typename Arrangement::Halfedge_const_handle  hh;
  );
  
  CGAL_precondition_msg (!(assign (hh, obj1)) && !(assign (hh, obj2)),
                      "The curve intersects the interior of existing edges.");
  
  // Check whether the located features containing the curve endpoints
  // are vertices or faces, and use the proper specialized insertion function
  // accordingly.
  typename Arrangement::Vertex_const_handle    vh1;
  typename Arrangement::Vertex_const_handle    vh2;
  typename Arrangement::Halfedge_handle        res;

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
      // Only the first endpoint is associated with an existing vertex.
      res = arr.insert_from_vertex (c,
                                    arr.non_const_handle (vh1));
    }
  }
  else
  {
    if (assign (vh2, obj2))
    {
      // Only the second endpoint is associated with an existing vertex.
      res = arr.insert_from_vertex (c,
                                    arr.non_const_handle (vh2));
    }
    else
    {
      // Both endpoints are not associated with existing vertices, so
      // we must insert the curve in the interior of a face.
      typename Arrangement::Face_const_handle      fh1;
      typename Arrangement::Face_const_handle      fh2;

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
template <class Arrangement, class InputIterator>
void arr_insert_non_intersecting
                (Arrangement& arr,
                 InputIterator begin, InputIterator end)
{
  // Notify the arrangement observers that a global operation is about to 
  // take place.
  Arr_accessor<Arrangement>                             arr_access (arr);

  arr_access.notify_before_global_change();

  // Perform the aggregated insertion.
  Arr_non_x_aggregate_insert<Arrangement>  agg_insert_obj (arr.get_traits(), 
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
template <class Arrangement>
typename Arrangement::Face_handle
arr_remove_edge (Arrangement& arr,
                 typename Arrangement::Halfedge_handle e)
{
  // Notify the arrangement observers that a global operation is about to 
  // take place.
  Arr_accessor<Arrangement>                             arr_access (arr);

  arr_access.notify_before_global_change();

  // Keep track of the end-vertices of the edge we are about to remove.
  typename Arrangement::Vertex_handle  v_ends[2];

  v_ends[0] = e.source();
  v_ends[1] = e.target();

  // Remove the edge from the arrangement.
  typename Arrangement::Face_handle    face = arr.remove_edge (e);

  // Examine the end-vertices: If a vertex has now two incident edges, and the
  // curves associated with these edges can be merged, merge the two edges and
  // remove the vertex.
  typedef Arr_traits_wrapper_2<typename Arrangement::Traits_2>  
                                                              Traits_wrapper_2;

  const Traits_wrapper_2                *traits =
                        static_cast<const Traits_wrapper_2*>(arr.get_traits());

  typename Arrangement::Halfedge_around_vertex_circulator  circ;
  typename Arrangement::Halfedge_handle                    e1, e2;
  int                                                      i;

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
        typename Arrangement::X_monotone_curve_2   cv;
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
