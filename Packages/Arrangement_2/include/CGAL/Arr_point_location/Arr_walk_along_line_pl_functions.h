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
  Face_const_handle     face = p_arr->unbounded_face();
  Face_const_handle     last_face;    
  Halfedge_const_handle closest_he;
  Holes_const_iterator  holes_it;
  Locate_type           locate_type = WPL_UFACE;

  while (face != last_face)
  {
    // Mark the last face visited.
    last_face = face;
    
    // Go over the holes in the current face.
    for (holes_it = face.holes_begin(); 
	 holes_it != face.holes_end() && last_face == face; 
	 ++holes_it)
    {
      if (_find_closest_feature (p, *holes_it,
				 true,           // Shoot up.
				 true,           // Include the query point.
				 closest_he, locate_type))
      {
        switch (locate_type)
	{
	case WPL_VERTEX:
	  // Return the target vertex of the closest edge.
	  return (make_object (closest_he.target()));

	case WPL_EDGE:
	  // Return the closest edge (which contains p in its interior).
	  return (make_object (closest_he));

	case WPL_FACE:
	  {
	    // Perform the walk from the face we have located.
	    _walk_along_line (p,
			      true,           // Shoot up.
			      true,           // Include the query point.
			      closest_he, locate_type);
	    
	    switch (locate_type)
	    {
	    case WPL_VERTEX:
	      // Return the target vertex of the closest edge.
	      return (make_object (closest_he.target()));
	      
	    case WPL_EDGE:
	      // Return the closest edge (which contains p in its interior).
	      return (make_object (closest_he));
	      
	    case WPL_FACE:
	      // Proceed to the next face.
	      face = closest_he.face();
	      break;
	      
	    default:
	      CGAL_assertion (false);
	    }
	  }
          break;

	case WPL_UFACE:
	  // Stop the walk:
          break;
        }
      }
    } // End loop on the current face's holes.
  }

  // Check if we should return the unbounded face:
  if (locate_type == WPL_UFACE)
    return (make_object (p_arr->unbounded_face()));
    
  // Return a handle to the face that contains p in its interior.
  return (make_object (face));
}

//-----------------------------------------------------------------------------
// Locate the arrangement feature which a vertical ray emanating from the
// given point hits.
//
template <class Arrangement>
Object Arr_walk_along_line_point_location<Arrangement>::_vertical_ray_shoot
    (const Point_2& p,
     bool shoot_up) const
{
  // Start from the unbounded face, and an invalid halfedge representing
  // the closest edge to p from above it so far. 
  Face_const_handle     face = p_arr->unbounded_face();
  Face_const_handle     last_face;    
  Halfedge_const_handle closest_he;
  Holes_const_iterator  holes_it;
  Locate_type           locate_type = WPL_UFACE;

  while (f != last_face)
  {
    // Mark the last face visited.
    last_face = face;
    
    // Go over the holes in the current face.
    for (holes_it = face.holes_begin(); 
	 holes_it != face.holes_end() && last_face == face; 
	 ++holes_it)
    {
      if (_find_closest_feature (p, *holes_it,
				 shoot_up,   // Shoot up or down.
				 false,      // Do not include the query point.
				 closest_he, locate_type))
      {
        switch (locate_type)
	{
	case WPL_VERTEX:
	  // Return the target vertex of the closest edge.
	  return (make_object (closest_he.target()));

	case WPL_EDGE:
	  {
	    // Perform the walk from face incident to closest_he.
	    _walk_along_line (p,
			      shoot_up,   // Shoot up or down.
			      false,      // Do not include the query point.
			      closest_he, locate_type);
	    
	    switch (locate_type)
	    {
	    case WPL_VERTEX:
	      // Proceed to the next face.
	      face = closest_he.twin().face();
	      break;
	      
	    case WPL_EDGE:
	      // Proceed to the next face.
	      face = closest_he.face();
	      break;
	      
	    case WPL_UFACE:
	      break;
	      
	    default:
	      CGAL_assertion (false);
	    }
	  }
          break;

	default:
	  CGAL_assertion (false);
        }
      }
    } // End loop on the current face's holes.
  }

  // If we reached here, we have encountered no vertex or that above p
  // (or below it, if we shoot down), so we return the unbounded face.
  return (make_object (p_arr->unbounded_face()));
}

//-----------------------------------------------------------------------------
// Walk along a the vertical line enamating from p from the given halfedge.
//
template <class Arrangement>
void Arr_walk_along_line_point_location<Arrangement>::_walk_along_line
    (const Point_2& p,
     bool shoot_up, bool inclusive,
     Halfedge_const_handle& closest_he,
     Locate_type& locate_type) const     
{
  Face_const_handle   face = 
    (inclusive || locate_type != WPL_VERTEX) ? closest_he.face() : 
                                               closest_he.twin().face();
  Face_const_handle   last_face;

  do
  {
    last_face = face;

    if (! face.is_unbounded()) 
    {
      _find_closest_feature (p, face.outer_ccb(),
			     shoot_up, inclusive,
			     closest_he, locate_type);
    }
        
    face = (inclusive || locate_type != WPL_VERTEX) ? closest_he.face() : 
                                                      closest_he.twin().face();
  
  } while((inclusive == (locate_type == WPL_UFACE)) && last_face != face);

  return;
}

//-----------------------------------------------------------------------------
// Find the closest feature to p (and lying above or below it) along the
// boundary of the given connected component.
//
template <class Arrangement>
bool Arr_walk_along_line_point_location<Arrangement>::_find_closest_feature
    (const Point_2& p,
     Ccb_halfedge_const_circulator circ,
     bool shoot_up, bool inclusive,
     Halfedge_const_handle& closest_he,
     Locate_type& locate_type) const
{
  // Set the results for comparison acording to the ray direction.
  const Comparison_result point_above_under = (shoot_up ? SMALLER : LARGER);
  const Comparison_result curve_above_under = (shoot_up ? LARGER : SMALLER);

  // The inclusive flag indicated whether the vertical ray includes its source:
  // If not (in case of vertical ray shooting) we should find the edge or
  // vertex right above (or below) the query point p.
  // If it does (in case of point location) and p lies on a vertex of on an
  // edge, we return this feature.

  // Check whether we have already encountered an intersection of the ray with
  // a previous halfedge (the parameter he is a handle to a valid edge).
  // We use he to store the halfedge closest to p so far.
  Halfedge_const_handle  invalid_halfedge;
  bool                   closest_he_valid = (closest_he != invalid_halfedge);
 
  // This flag indicates whether p is inside the connected component, as
  // determined by the number of times the vertical ray intersects the boundary
  // egdes (in case of an odd number p is inside, in case of an even number
  // p is outside). At the moment we have no intersections yet, so:
  bool    inside_cc = false;

  // Go over all edges in the connected components.
  typename Traits_wrapper_2::Compare_xy_2           compare_xy = 
                                      traits->compare_xy_2_object();
  typename Traits_wrapper_2::Is_in_x_range_2        is_in_x_range = 
                                      traits->is_in_x_range_2_object();
  typename Traits_wrapper_2::Compare_y_at_x_2       compare_y_at_x = 
                                      traits->compare_y_at_x_2_object();
  typename Traits_wrapper_2::Is_vertical_2          is_vertical = 
                                      traits->is_vertical_2_object();
  typename Traits_wrapper_2::Compare_y_position_2   compare_y_position =
                                      traits->compare_y_position_2_object();
  Ccb_halfedge_const_circulator curr = circ;
  Comparison_result             res1, res2;
  Comparison_result             res;
  Comparison_result             y_pos; 
  bool                          in_x_range;

  do
  {
    // Check if the query point is in the x-range of the curve associated
    // with the current halfedge.
    const X_monotone_curve_2& cv = (*curr).curve();

    if (inclusive)
    {
      // Check if p is lexicographically in the x-range of the curve.
      res1 = compare_xy (p, (*curr).source().point());
      res2 = compare_xy (p, (*curr).target().point());

      in_x_range = (res1 != res2 || res1 == EQUAL || res2 == EQUAL);
    }
    else
    {
      // Use a simple x-range query:
      in_x_range = is_in_x_range (cv, p);
    }

    // Skip the current edge if p is not in its x-range. 
    if (! in_x_range)
    {
      ++curr;
      continue;
    }
    
    // Get the relative position of p with respect to the current curve.
    res = compare_y_at_x (p, cv);

    if (res == point_above_under)
    {
      // p is in the x-range of the current halfedge, which lies above it
      // (if we shoot up) or below it (if we shoot the ray down).
      // Flip the inside_cc flag, as we have found a new edge on the boundary
      // that our ray intersects - unless the curve is vertical, in which case
      // the vertical ray overlaps it.
      if (! is_vertical (cv))
	inside_cc = !inside_cc;

      // Check if the current curve lies closer to p than the closest curve
      // so far (associated with closest_he).
      y_pos = EQUAL;
      if (! closest_he_valid ||
	  ((y_pos = compare_y_position (closest_he.curve(), cv)) ==
	   curve_above_under))
      {
	closest_he = *curr;
	closest_he_valid = true;
      }
      
      // Examine the case where the two curves are equal.
      if (closest_he_valid && y_pos == EQUAL &&
	  closest_he != *curr && closest_he != (*curr).twin())
      {
	// In this case closest_he and the current edge share an end-vertex
	// lying right above the query point p. This means we have the
	// following scenario:
        //                   
        //                   
        //  (.)----(.)---(.) 
        //          ^        
        //          |        
	//          |        
        //          p        
        //
	// In this case, we take closest_he to be the edge "left" edge in case
	// we shoot up, and the "right" edge, in case we shoot down.
	bool    curr_to_right;

	if ((*curr).source() == closest_he.source() || 
	    (*curr).source() == closest_he.target())
	{
	  // The common endpoint is the current source point. Check if the
	  // current halfedge extends to its left or to its right.
	  curr_to_right = (compare_xy ((*curr).source().point(),
				       (*curr).target().point()) == SMALLER);
	}
	else
	{
	  CGAL_assertion ((*curr).target() == closest_he.source() || 
			  (*curr).target() == closest_he.target());

	  // The common endpoint is the current target point. Check if the
	  // current halfedge extends to its left or to its right.
	  curr_to_right = (compare_xy ((*curr).source().point(),
				       (*curr).target().point()) == LARGER);
	}

	if ((shoot_up && !curr_to_right) ||
	    (!shoot_up && curr_to_right))
	{
	  closest_he = *curr;
	}
      }
    }
    else if (res == EQUAL)
    {
      // In this case p lies on the curve.
      if (!inclusive)
      {
	// If we do not include p itself (in vertical ray-shooting queries),
	// the only case we have to check is when p is located in the interior
	// of a vertical segment. Otherwise, we ignore this edge.
	if (is_vertical (cv) &&
	    (compare_xy (p, traits->construct_max_vertex_2_object() (cv)) 
	     == SMALLER))
	{
	  closest_he = *curr;
	  locate_type = WPL_EDGE;
	  return (true);
	}
      }
      else
      {
	// Check if p is one of the edge endpoints. If so, return the halfedge
	// whose target is the vertex that contains p.
	if (compare_xy (p, (*curr).source().point()) == EQUAL)
	{
	  // p lies on the source vertex:
	  closest_he = (*curr).twin();
	  locate_type = WPL_VERTEX;
	}
	else if (compare_xy (p, (*curr).target().point()) == EQUAL)
	{
	  // p lies on the target vertex:
	  closest_he = *curr;
	  locate_type = WPL_VERTEX;
	}
	else
	{
	  // p lies in the interior of the current edge:
	  closest_he = *curr;
	  locate_type = WPL_EDGE;
	}

	return (true);
      }
    }

    // Proceed to the next halfedge along the boundary.
    ++curr;

  } while (curr != circ);

  // Check if we found no intersections:
  if (! closest_he_valid)
  {
    locate_type = WPL_UFACE;
    return (false);
  }

  if (inclusive)
  {
    // We are in point-location mode.
    if (inside_cc)
    {
      // Make sure that the closest halfedge we return is directed from right
      // to left (if we shoot up) - or from left to right (if we shoot down).
      // This way, p is located in the incident face of closest_he.
      const bool    closest_he_directed_right =
	(compare_xy (closest_he.source().point(),
		     closest_he.target().point()) == SMALLER);

      if ((shoot_up && closest_he_directed_right) ||
	  (!shoot_up && !closest_he_directed_right))
      {
	closest_he = closest_he.twin();
      }

      locate_type = WPL_FACE;
    }
    else
    {
      locate_type = WPL_UFACE;
    }
  }
  else
  {
    // We are in vertical ray-shooting mode.
    if (is_vertical (closest_he.curve()))
    {
      // p is below (or above, if we shoot down) a vertical segment.
      // Direct the closest halfegde so that its target is the endpoint
      // closer to p.
      res = compare_xy (closest_he.source().point(), 
			closest_he.target().point());

      if ((shoot_up && res == SMALLER) ||
	  (!shoot_up && res == LARGER))
      {
	closest_he = closest_he.twin();
      }
      locate_type = WPL_VERTEX;
    }
    else
    {
      if (traits->compare_x_2_object() (closest_he.source().point(),
					p) == EQUAL)
      {
	// When we shoot a ray from p we ancounter the source vertex:
	closest_he = closest_he.twin();
	locate_type = WPL_VERTEX;
      }
      else if (traits->compare_x_2_object() (closest_he.target().point(),
					     p) == EQUAL)
      {
	// When we shoot a ray from p we ancounter the target vertex:
	locate_type = WPL_VERTEX;
      }
      else
      {
	locate_type = WPL_EDGE;
      }
    }
  }

  return (true);
}


CGAL_END_NAMESPACE

#endif
