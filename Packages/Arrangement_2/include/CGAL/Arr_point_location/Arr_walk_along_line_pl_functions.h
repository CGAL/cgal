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
// Author(s)     : Ron Wein   <wein@post.tau.ac.il>
//                 (based on old version by Oren Nechushtan
//                                      and Iddo Hanniel)
#ifndef CGAL_ARR_WALK_ALONG_LINE_POINT_LOCATION_FUNCTIONS_H
#define CGAL_ARR_WALK_ALONG_LINE_POINT_LOCATION_FUNCTIONS_H

/*! \file
 * Member-function definitions for the
 * Arr_walk_along_line_point_location<Arrangement> class.
 */

CGAL_BEGIN_NAMESPACE

//-----------------------------------------------------------------------------
// Locate the arrangement feature containing the given point.
//
template <class Arrangement>
Object Arr_walk_along_line_point_location<Arrangement>::locate
    (const Point_2& p) const
{
  // Start from the unbounded face, and an invalid halfedge representing
  // the closest edge to p from above it so far.
  typename Traits_wrapper_2::Equal_2  equal = traits->equal_2_object();
  Holes_const_iterator   holes_it;
  Face_const_handle      face = p_arr->unbounded_face();
  Halfedge_const_handle  closest_he;
  bool                   is_in_face;
  bool                   is_on_edge;
  bool                   found_containing_hole;

  do
  {
    // Go over the holes in the current face.
    found_containing_hole = false;
    for (holes_it = face.holes_begin();
	 holes_it != face.holes_end() && !found_containing_hole;
	 ++holes_it)
    {
      // Check if the point is inside the current hole.
      is_in_face = _is_in_connected_component (p, *holes_it,
					       true,        // Shoot up.
					       true,        // Including p.
					       closest_he, is_on_edge);

      // Check if the query point is located on the returned edge.
      if (is_on_edge)
      {
	// Check if the point is located on one of the edge endpoints.
	if (equal (p, closest_he.source().point()))
	{
	  // The query point is located on the source vertex:
	  return (make_object (closest_he.source()));
	}
	else if (equal (p, closest_he.target().point()))
	{
	  // The query point is located on the target vertex:
	  return (make_object (closest_he.target()));
	}

	// The query point is located in the edge interior:
	return (make_object (closest_he));
      }

      // Check if the point is contained in the interior of the current hole.
      if (is_in_face)
      {
	// Move inside the faces that constitute the hole, the first one
	// being incident face of the twin of closest halfedge found so far.
	CGAL_assertion (face != closest_he.twin().face());

	face = closest_he.twin().face();

	// Perform a vertical walk along the faces of the hole until locating
	// a face that contains the qeury point.
	do
	{
	  CGAL_assertion_code (
	    Halfedge_const_handle  old_closest_he = closest_he;
	  );
 
	  is_in_face = _is_in_connected_component (p, face.outer_ccb(),
						   true,    // Shoot up.
						   true,    // Including p.
						   closest_he, is_on_edge);

	  // Check if the query point was located on an edge.
	  if (is_on_edge)
	  {
	    // Check if the point is located on one of the edge endpoints.
	    if (equal (p, closest_he.source().point()))
	    {
	      // The query point is located on the source vertex:
	      return (make_object (closest_he.source()));
	    }
	    else if (equal (p, closest_he.target().point()))
	    {
	      // The query point is located on the target vertex:
	      return (make_object (closest_he.target()));
	    }

	    // The query point is located in the edge iterior:
	    return (make_object (closest_he));
	  }
	  
	  // If the point is not contained in the face, move to the neighboring
	  // face from below, using the closest halfedge located so far.
	  if (! is_in_face)
	  {
	    CGAL_assertion (old_closest_he != closest_he);
	    face = closest_he.twin().face();
	  }

	} while (! is_in_face);

	// We have located a face in the hole that contains the query point.
	// This will break the internal loop on holes, and we shall proceed
	// for another iteration of the external loop, trying to locate p in
	// one of the hole of this new face.
	found_containing_hole = true;
      }
    } // End loop on the current face's holes.

  } while (found_containing_hole);

  // If we reached here, we did not locate the query point in any of the holes
  // inside the current face, so we conclude it is contained in this face.
  // However, we first have to check whether the query point coincides with
  // any of the isolated vertices contained inside this face.
  Isolated_vertices_const_iterator   iso_verts_it;

  for (iso_verts_it = face.isolated_vertices_begin();
       iso_verts_it != face.isolated_vertices_end(); ++iso_verts_it)
  {
    if (equal (p, (*iso_verts_it).point()))
      return (make_object (*iso_verts_it));
  }

  // The query point is contained in the face interior:
  return (make_object (face));
}

//-----------------------------------------------------------------------------
// Locate the arrangement feature which a vertical ray emanating from the
// given point hits.
//
template <class Arrangement>
Object Arr_walk_along_line_point_location<Arrangement>::
_vertical_ray_shoot (const Point_2& p,
		     bool shoot_up) const
{
  // Start from the unbounded face, and an invalid halfedge representing
  // the closest edge to p from above it so far.
  typename Traits_wrapper_2::Is_vertical_2        is_vertical =
                                            traits->is_vertical_2_object();

  Holes_const_iterator   holes_it;
  Face_const_handle      face = p_arr->unbounded_face();
  Halfedge_const_handle  closest_he;
  bool                   is_in_face;
  bool                   is_on_edge;
  bool                   found_containing_hole;

  do
  {
    // Go over the holes in the current face.
    found_containing_hole = false;
    for (holes_it = face.holes_begin();
	 holes_it != face.holes_end() && !found_containing_hole;
	 ++holes_it)
    {
      // Check if the point is inside the current hole.
      is_in_face = _is_in_connected_component (p, *holes_it,
					       shoot_up,
					       false,     // Not including p.
					       closest_he, is_on_edge);

      // Check if the query point is located on the returned edge.
      // This can happen only if the edge is vertical.
      if (is_on_edge)
      {
	CGAL_assertion (is_vertical (closest_he.curve()));
	return (make_object (closest_he));
      }

      // Check if the point is contained in the interior of the current hole.
      if (is_in_face)
      {
	// Move inside the faces that constitute the hole, the first one
	// being incident face of the twin of closest halfedge found so far.
	CGAL_assertion (face != closest_he.twin().face());

	face = closest_he.twin().face();

	// Perform a vertical walk along the faces of the hole until locating
	// a face that contains the qeury point.
	do
	{
	  CGAL_assertion_code (
	    Halfedge_const_handle  old_closest_he = closest_he;
	  );
 
	  is_in_face = _is_in_connected_component (p, face.outer_ccb(),
						   shoot_up,
						   false,   // Not including p.
						   closest_he, is_on_edge);

	  // Check if the query point was located on an edge.
	  // This can happen only if the edge is vertical.
	  if (is_on_edge)
	  {
	    CGAL_assertion (is_vertical (closest_he.curve()));
	    return (make_object (closest_he));
	  }
	  
	  // If the point is not contained in the face, move to the neighboring
	  // face from below (or above, if we shoot downward), using the
	  // closest halfedge located so far.
	  if (! is_in_face)
	  {
	    CGAL_assertion (old_closest_he != closest_he);
	    face = closest_he.twin().face();
	  }

	} while (! is_in_face);

	// We have located a face in the hole that contains the query point.
	// This will break the internal loop on holes, and we shall proceed
	// for another iteration of the external loop, trying to locate p in
	// one of the holes of this new face.
	found_containing_hole = true;
      }
    } // End loop on the current face's holes.

  } while (found_containing_hole);

  // Check whether one of the isolated vertices in the face containing p lies
  // above (or below) it, closer than the closest halfdge we have located.
  typename Traits_wrapper_2::Compare_x_2          compare_x =
                                             traits->compare_x_2_object();
  typename Traits_wrapper_2::Compare_xy_2         compare_xy =
                                             traits->compare_xy_2_object();
  typename Traits_wrapper_2::Compare_y_at_x_2     compare_y_at_x =
                                             traits->compare_y_at_x_2_object();

  const Comparison_result point_above_under = (shoot_up ? SMALLER : LARGER);

  Isolated_vertices_const_iterator   iso_verts_it;
  Vertex_const_handle                closest_iso_v;
  const Vertex_const_handle          invalid_v;
  const Halfedge_const_handle        invalid_he;

  for (iso_verts_it = face.isolated_vertices_begin();
       iso_verts_it != face.isolated_vertices_end(); ++iso_verts_it)
  {
    // The current isolated vertex should have the same x-coordinate as the
    // query point in order to be below or above it.
    if (compare_x (p, (*iso_verts_it).point()) != EQUAL)
      continue;

    // Make sure the isolated vertex is above the query point (if we shoot up)
    // or below it (if we shoot down).
    if (compare_xy (p, (*iso_verts_it).point()) != point_above_under)
      continue;

    // Check if the current isolated vertex lies closer to the query point than
    // the closest feature so far.
    if (closest_iso_v == invalid_v)
    {
      // Compare the current isolated vertex with the closest halfedge.
      if (closest_he == invalid_he ||
          compare_y_at_x ((*iso_verts_it).point(),
                          closest_he.curve()) == point_above_under)
      {
        closest_iso_v = *iso_verts_it;
      }
    }
    else if (compare_xy ((*iso_verts_it).point(),
                         closest_iso_v.point()) == point_above_under)
    {
      closest_iso_v = *iso_verts_it;
    }
  }

  if (closest_iso_v != invalid_v)
  {
    // The first object we encounter when we shoot a vertical ray from p is
    // an isolated vertex:
    return (make_object (closest_iso_v));
  }

  // If we reached here, closest_he is the closest edge from above (below)
  // the query point.
  if (closest_he == invalid_he)
  {
    // We did not encounter any edge above (below) the query point:
    return Object();
  }

  // Check if one of closest_he's end vertices lies directly above (below) the
  // query point, and if so, return this vertex.
  if (! is_vertical (closest_he.curve()))
  {
    if (traits->compare_x_2_object() (closest_he.source().point(), p) == EQUAL)
      return (make_object (closest_he.source()));

    if (traits->compare_x_2_object() (closest_he.target().point(), p) == EQUAL)
      return (make_object (closest_he.target()));
  }
  else
  {
    // The entire vertical segment is above (below) the query point: Return the
    // endpoint closest to it.
    const bool    is_directed_up =
      (traits->compare_xy_2_object() (closest_he.source().point(),
				      closest_he.target().point()) == SMALLER);

    if ((shoot_up && is_directed_up) ||
	(! shoot_up && ! is_directed_up))
    {
      return (make_object (closest_he.source()));
    }
    else
    {
      return (make_object (closest_he.target()));
    }
  }

  // The interior of the edge is closest to the query point:
  return (make_object (closest_he));
}

//-----------------------------------------------------------------------------
// Find the closest feature to p (and lying above or below it) along the
// boundary of the given connected component.
//
template <class Arrangement>
bool Arr_walk_along_line_point_location<Arrangement>::
_is_in_connected_component (const Point_2& p,
			    Ccb_halfedge_const_circulator circ,
			    bool shoot_up,
			    bool inclusive,
			    Halfedge_const_handle& closest_he,
			    bool& is_on_edge) const
{
  // Set the results for comparison acording to the ray direction.
  const Comparison_result point_above_under = (shoot_up ? SMALLER : LARGER);
  const Comparison_result curve_above_under = (shoot_up ? LARGER : SMALLER);

  // Keep a counter of the number of halfedges of the connected component's
  // boundary that intersect an upward (or downward, if shoot_up is false)
  // vertical ray emanating from the query point p (except for some degenerate
  // cases that are explained below).
  unsigned int              n_ray_intersections = 0;

  typename Traits_wrapper_2::Equal_2              equal = 
                                            traits->equal_2_object();
  typename Traits_wrapper_2::Is_vertical_2        is_vertical =
                                            traits->is_vertical_2_object();
  typename Traits_wrapper_2::Compare_x_2          compare_x =
                                            traits->compare_x_2_object();
  typename Traits_wrapper_2::Compare_y_at_x_2     compare_y_at_x =
                                            traits->compare_y_at_x_2_object();
  typename Traits_wrapper_2::Compare_y_position_2 compare_y_position =
                                        traits->compare_y_position_2_object();

  // Start from the first non-vertical segment in the connected component.
  Ccb_halfedge_const_circulator  first = circ;
  bool                           found_non_vertical = false;

  do
  {
    // Stop if we found a non-vertical curve.
    if (! is_vertical ((*first).curve()))
    {
      found_non_vertical = true;
      break;
    }

    if (inclusive)
    {
      // Check if the current vertical curve contains the query point.
      if (compare_x ((*first).source().point(), p) == EQUAL &&
	  compare_y_at_x (p, (*first).curve()) == EQUAL)
      {
	closest_he = *first;
	is_on_edge = true;
	return (true);
      }
    }
    else
    {
      // Check if the current vertical curve contains the query point in its
      // interior.
      if (compare_x ((*first).source().point(), p) == EQUAL &&
	  compare_y_at_x (p, (*first).curve()) == EQUAL &&
	  ! equal ((*first).source().point(), p) &&
	  ! equal ((*first).target().point(), p))
      {
	closest_he = *first;
	is_on_edge = true;
	return (true);
      }
    }

    // Move to the next curve.
    ++first;

  } while (first != circ);

  if (! found_non_vertical)
  {
    // In this case the entire component is comprised of vertical segments,
    // so it has any empty interior and p cannot lie inside it.
    return (false);
  }

  // Go over all curves of the boundary, starting from the non-vertical curve
  // we have located, and count those which are above p.
  const Halfedge_const_handle    invalid_he;
  Ccb_halfedge_const_circulator  curr = first;
  Comparison_result              source_res, target_res;
  Comparison_result              res;

  do
  {
    // Ignore the current edge if p is not in the x-range of its curve.
    source_res = compare_x ((*curr).source().point(), p);
    target_res = compare_x ((*curr).target().point(), p);

    if (source_res == target_res && source_res != EQUAL)
    {
      ++curr;
      continue;
    }

    // Check whether p lies above or below the curve.
    res = compare_y_at_x (p, (*curr).curve());

    if (res == EQUAL)
    {
      // The current edge contains the query point. If the seach is inclusive
      // we return the edge. Otherwise, we return it only if it is vertical,
      // and contains p in its interior.
      if (inclusive)
      {
	closest_he = *curr;
	is_on_edge = true;
	return (true);
      }
      else
      {
	if (is_vertical ((*curr).curve()) &&
	    ! equal ((*curr).source().point(), p) &&
	    ! equal ((*curr).target().point(), p))
	{
	  closest_he = *curr;
	  is_on_edge = true;
	  return (true);
	}
      }
    }

    // If the point is above the current edge (or below it, if we shoot down),
    // move to the next edge.
    if (res != point_above_under)
    {
      ++curr;
      continue;
    }

    // Note that we do not have count intersections (actually these are 
    // overlaps) of the vertical ray we shoot with vertical segments along
    // the boundary.
    if ( ! is_vertical ((*curr).curve())) 
    {
      // The current curve is not vertical. Check the query point is in the
      // semi-open x-range (source, target] of this curve and lies below it.
      if (source_res != EQUAL)
      {
	// If the current curve is closer to p than the closest curve found
	// so far, assign its halfedge as the closest one.
	if (closest_he == invalid_he ||
	    compare_y_position (closest_he.curve(), 
				(*curr).curve()) == curve_above_under)
	{
	  closest_he = *curr;
	}

        // In the degenerate case where p lies below the target vertex of
        // the current halfedge, we have to be a bit careful:
        if (target_res == EQUAL)
        {
          // Locate the next halfedge along the boundary that does not
          // contain a vertical segment.
          Halfedge_const_handle  next_non_vert = *curr;

          do
          {
            next_non_vert = next_non_vert.next();

            CGAL_assertion (next_non_vert != *curr);

          } while (is_vertical (next_non_vert.curve()));

          // In case the source of the current curve and the target of
          // the next non-vertical curve lie on opposite sides of the
          // ray we shoot from p (case (a)), we have to count an
          // intersection. Otherwise, we have a "tangency" with the ray
          // (case (b)) and it is not necessary to count it.
          //
          //            +--------+              +                 .
          //            |   next                 \ next           .
          //            |                         \               .
          //            +                          +              .
          //           /                          /               .
          //     curr /                     curr /                .
          //         /                          /                 .
          //        +  (.)p                    +  (.)p            .
          //
          //          (a)                        (b)
          //
	  target_res = compare_x (next_non_vert.target().point(), p);

          CGAL_assertion (source_res != EQUAL && target_res != EQUAL);

          if (source_res != target_res)
            n_ray_intersections++;
        }
        else 
        {
          // In this case p lies under the interior of the current x-montone
          // curve, so the vertical ray we shoot intersects it exactly once.
          n_ray_intersections++;
        }
      }
    }
    else
    {
      // Assign an edge associated with a vertical curve as the closest edge
      // only if the vertical curve has an endpoint closer to p than the
      // closest edge so far.
      if (closest_he == invalid_he ||
	  (compare_y_at_x ((*curr).source().point(),
			   closest_he.curve()) == point_above_under) ||
	  (compare_y_at_x ((*curr).target().point(),
			   closest_he.curve()) == point_above_under))
      {
	closest_he = *curr;
      }
    }

    // Proceed to the next halfedge along the component boundary.
    ++curr;

  } while (curr != first);

  // The query point lies inside the connected components if and only if the
  // ray we shoot from it intersects the boundary an odd number of time.
  is_on_edge = false;
  return ((n_ray_intersections % 2) != 0);
}

CGAL_END_NAMESPACE

#endif
