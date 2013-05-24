// Copyright (c) 2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s)     : Ron Wein   <wein@post.tau.ac.il>
//                 (based on old version by Oren Nechushtan
//                                      and Iddo Hanniel)

#ifndef CGAL_ARR_WALK_ALONG_LINE_POINT_LOCATION_IMPL_H
#define CGAL_ARR_WALK_ALONG_LINE_POINT_LOCATION_IMPL_H

/*! \file
 * Member-function definitions for the
 * Arr_walk_along_line_point_location<Arrangement> class.
 */

namespace CGAL {

//-----------------------------------------------------------------------------
// Locate the arrangement feature containing the given point.
//
template <class Arrangement>
typename Arr_walk_along_line_point_location<Arrangement>::result_type
Arr_walk_along_line_point_location<Arrangement>::locate(const Point_2& p) const
{
  // Start from the initial face (as determined by the topology-traits class),
  // and an invalid halfedge representing the closest edge to p from above it
  // so far.
  typename Traits_adaptor_2::Equal_2  equal = geom_traits->equal_2_object();
  Inner_ccb_const_iterator            ic_it;
  Face_const_handle                   face;
  Halfedge_const_handle               closest_he;
  bool                                is_in_face;
  bool                                is_on_edge;
  bool                                closest_to_target;
  bool                                found_containing_hole;

  face = Face_const_handle (top_traits->initial_face());
  do {
    // Go over the holes in the current face.
    found_containing_hole = false;
    for (ic_it = face->inner_ccbs_begin();
         !found_containing_hole && ic_it != face->inner_ccbs_end();
         ++ic_it)
    {
      // Check if the point is inside the current hole.
      is_in_face = _is_in_connected_component (p, *ic_it,
                                               true,        // Shoot up.
                                               true,        // Including p.
                                               closest_he,
                                               is_on_edge,
                                               closest_to_target);

      // Check if the query point is located on the returned edge.
      if (is_on_edge) {
        // Check if the point is located on one of the edge endpoints.
        if (! closest_he->source()->is_at_open_boundary() &&
            equal (p, closest_he->source()->point()))
        {
          // The query point is located on the source vertex:
          Vertex_const_handle  vh = closest_he->source();
          return make_result(vh);
        }
        else if (! closest_he->target()->is_at_open_boundary() &&
                 equal (p, closest_he->target()->point()))
        {
          // The query point is located on the target vertex:
          Vertex_const_handle  vh = closest_he->target();
          return make_result(vh);
        }

        // The query point is located in the edge interior:
        return make_result(closest_he);
      }

      // Check if the point is contained in the interior of the current hole.
      if (is_in_face) {
        if (closest_to_target) {
          // Deal with the following scenario, where the definition of the
          // closest halfedge is not unique (all halfedges around x are
          // closest to p):
          //
          //          +--+--+
          //           \ | /
          //            \|/
          //     +-------x------+
          //     |              |
          //     |      (.)p    |
          //     |              |
          //     +--------------+
          //
          // In this case, we find the first halfegde whose target is x
          // in a clockwise direction from "6 o'clock" around x and take
          // its incident face.

          // Shoot up.
          closest_he = _first_around_vertex(closest_he->target(), true);
        }
        
        // Move inside the faces that constitute the hole, the first one
        // being incident face of the twin of closest halfedge found so far.
        CGAL_assertion (face != closest_he->twin()->face());

        face = closest_he->twin()->face();

        // Perform a vertical walk along the faces of the hole until locating
        // a face that contains the qeury point.
        do {
          CGAL_assertion_code
            (Halfedge_const_handle  old_closest_he = closest_he);

          is_in_face = _is_in_connected_component (p, face->outer_ccb(),
                                                   true,    // Shoot up.
                                                   true,    // Including p.
                                                   closest_he,
                                                   is_on_edge,
                                                   closest_to_target);

          // Check if the query point was located on an edge.
          if (is_on_edge) {
            // Check if the point is located on one of the edge endpoints.
            if (! closest_he->source()->is_at_open_boundary() &&
                equal (p, closest_he->source()->point()))
            {
              // The query point is located on the source vertex:
              Vertex_const_handle  vh = closest_he->source();
              return make_result(vh);
            }
            else if (! closest_he->target()->is_at_open_boundary() &&
                     equal (p, closest_he->target()->point()))
            {
              // The query point is located on the target vertex:
              Vertex_const_handle  vh = closest_he->target();
              return make_result(vh);
            }

            // The query point is located in the edge iterior:
            return make_result(closest_he);
          }

          // If the point is not contained in the face, move to the neighboring
          // face from below, using the closest halfedge located so far.
          if (! is_in_face) {
            if (closest_to_target) {
              // The query point lies below the target vertex of the closest
              // halfedge: In this case we have to locate the first halfedge
              // we encounter when going around this target vertex from
              // "6 o'clock" (see above).

              // Shoot up.
              closest_he = _first_around_vertex(closest_he->target(), true);
            }

            CGAL_assertion (old_closest_he != closest_he);
            face = closest_he->twin()->face();
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
  Isolated_vertex_const_iterator   iso_verts_it;

  for (iso_verts_it = face->isolated_vertices_begin();
       iso_verts_it != face->isolated_vertices_end(); ++iso_verts_it)
  {
    if (equal(p, iso_verts_it->point())) {
      Vertex_const_handle  vh = iso_verts_it;
      return make_result(vh);
    }
  }

  // The query point is contained in the face interior:
  return make_result(face);
}

//-----------------------------------------------------------------------------
// Locate the arrangement feature which a vertical ray emanating from the
// given point hits.
//
template <class Arrangement>
typename Arr_walk_along_line_point_location<Arrangement>::result_type
Arr_walk_along_line_point_location<Arrangement>::
_vertical_ray_shoot(const Point_2& p, bool shoot_up) const
{
  // Start from the initial face (as determined by the topology-traits class),
  // and an invalid halfedge representing the closest edge to p from above it
  // so far.
  typename Traits_adaptor_2::Is_vertical_2        is_vertical =
                                         geom_traits->is_vertical_2_object();
  Inner_ccb_const_iterator            ic_it;
  Face_const_handle                   face;
  Halfedge_const_handle               closest_he;
  bool                                is_in_face;
  bool                                is_on_edge;
  bool                                closest_to_target;
  bool                                found_containing_hole;

  face = Face_const_handle (top_traits->initial_face());
  do {
    // Go over the holes in the current face.
    found_containing_hole = false;
    for (ic_it = face->inner_ccbs_begin();
         !found_containing_hole && ic_it != face->inner_ccbs_end();
         ++ic_it)
    {
      // Check if the point is inside the current hole.
      is_in_face = _is_in_connected_component (p, *ic_it,
                                               shoot_up,
                                               false,     // Not including p.
                                               closest_he,
                                               is_on_edge,
                                               closest_to_target);

      // Check if the query point is located on the returned edge.
      // This can happen only if the edge is vertical.
      if (is_on_edge) {
        CGAL_assertion (is_vertical (closest_he->curve()));
        return make_result(closest_he);
      }

      // Check if the point is contained in the interior of the current hole.
      if (is_in_face) {
        if (closest_to_target) {
          // The query point lies below the target vertex of the closest
          // halfedge: In this case we have to locate the first halfedge
          // we encounter when going around this target vertex from
          // "6 o'clock" (when shooting up) or from "12 o'clock" (when
          // shooting down).
          closest_he = _first_around_vertex(closest_he->target(), shoot_up);
        }

        // Move inside the faces that constitute the hole, the first one
        // being incident face of the twin of closest halfedge found so far.
        CGAL_assertion (face != closest_he->twin()->face());

        face = closest_he->twin()->face();

        // Perform a vertical walk along the faces of the hole until locating
        // a face that contains the qeury point.
        do {
          CGAL_assertion_code (
            Halfedge_const_handle  old_closest_he = closest_he;
          );

          is_in_face = _is_in_connected_component (p, face->outer_ccb(),
                                                   shoot_up,
                                                   false,   // Not including p.
                                                   closest_he,
                                                   is_on_edge,
                                                   closest_to_target);

          // Check if the query point was located on an edge.
          // This can happen only if the edge is vertical.
          if (is_on_edge) {
            CGAL_assertion (is_vertical (closest_he->curve()));
            return make_result(closest_he);
          }

          // If the point is not contained in the face, move to the neighboring
          // face from below (or above, if we shoot downward), using the
          // closest halfedge located so far.
          if (! is_in_face) {
            if (closest_to_target) {
              // The query point lies below the target vertex of the closest
              // halfedge: In this case we have to locate the first halfedge
              // we encounter when going around this target vertex from
              // "6 o'clock" (when shooting up) or from "12 o'clock" (when
              // shooting down).
              closest_he = _first_around_vertex(closest_he->target(), shoot_up);
            }

            CGAL_assertion (old_closest_he != closest_he);
            face = closest_he->twin()->face();
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
  typename Traits_adaptor_2::Compare_x_2          compare_x =
    geom_traits->compare_x_2_object();
  typename Traits_adaptor_2::Compare_xy_2         compare_xy =
    geom_traits->compare_xy_2_object();
  typename Traits_adaptor_2::Compare_y_at_x_2     compare_y_at_x =
    geom_traits->compare_y_at_x_2_object();

  const Comparison_result point_above_under = (shoot_up ? SMALLER : LARGER);

  Isolated_vertex_const_iterator     iso_verts_it;
  Vertex_const_handle                closest_iso_v;
  const Vertex_const_handle          invalid_v;
  const Halfedge_const_handle        invalid_he;

  for (iso_verts_it = face->isolated_vertices_begin();
       iso_verts_it != face->isolated_vertices_end(); ++iso_verts_it)
  {
    // The current isolated vertex should have the same x-coordinate as the
    // query point in order to be below or above it.
    if (compare_x (p, iso_verts_it->point()) != EQUAL)
      continue;

    // Make sure the isolated vertex is above the query point (if we shoot up)
    // or below it (if we shoot down).
    if (compare_xy (p, iso_verts_it->point()) != point_above_under)
      continue;

    // Check if the current isolated vertex lies closer to the query point than
    // the closest feature so far.
    if (closest_iso_v == invalid_v) {
      // Compare the current isolated vertex with the closest halfedge.
      if (closest_he == invalid_he || closest_he->is_fictitious() ||
          compare_y_at_x(iso_verts_it->point(), closest_he->curve()) ==
          point_above_under)
      {
        // The isolated vertex is closer to the query point than thr closest
        // edge found so far (note that a fictitious edge means the edge is
        // at infinity).
        closest_iso_v = iso_verts_it;
      }
    }
    else if (compare_xy(iso_verts_it->point(), closest_iso_v->point()) ==
             point_above_under)
    {
      // The current isolated vertex is closer to the query point than the
      // closest isolated vertex found so far.
      closest_iso_v = iso_verts_it;
    }
  }

  if (closest_iso_v != invalid_v) {
    // The first object we encounter when we shoot a vertical ray from p is
    // an isolated vertex:
    return make_result(closest_iso_v);
  }

  // If we have no halfedge above, return the initial face.
  if (closest_he == invalid_he) {
    Face_const_handle  uf = Face_const_handle (top_traits->initial_face());

    return make_result(uf);
  }

  // If we reached here, closest_he is the closest edge from above (below)
  // the query point.
  if (closest_he->is_fictitious()) {
    // If we have a fictitious edge, the ray we shoot is completely contained
    // in an unbounded face. This face is incident to closest_he:
    if ((shoot_up && closest_he->direction() == ARR_LEFT_TO_RIGHT) ||
        (!shoot_up && closest_he->direction() == ARR_RIGHT_TO_LEFT))
      closest_he = closest_he->twin();

    Face_const_handle  uf = closest_he->face();

    CGAL_assertion (! uf->is_fictitious());
    return make_result(uf);
  }

  // Check if one of closest_he's end vertices lies directly above (below) the
  // query point, and if so, return this vertex.
  if (! is_vertical(closest_he->curve())) {
    if (! closest_he->source()->is_at_open_boundary() &&
        compare_x (p, closest_he->source()->point()) == EQUAL)
    {
      Vertex_const_handle  vh = closest_he->source();
      return make_result(vh);
    }

    if (! closest_he->target()->is_at_open_boundary() &&
        compare_x (p, closest_he->target()->point()) == EQUAL)
    {
      Vertex_const_handle  vh = closest_he->target();
      return make_result(vh);
    }
  }
  else {
    // The entire vertical segment is above (below) the query point: Return the
    // endpoint closest to it.
    const bool is_directed_up = (closest_he->direction() == ARR_LEFT_TO_RIGHT);

    if ((shoot_up && is_directed_up) || (! shoot_up && ! is_directed_up))
    {
      Vertex_const_handle  vh = closest_he->source();
      return make_result(vh);
    }
    else {
      Vertex_const_handle  vh = closest_he->target();
      return make_result(vh);
    }
  }

  // The interior of the edge is closest to the query point:
  return make_result(closest_he);
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
                            bool& is_on_edge,
                            bool& closest_to_target) const
{
  // As far as we know, we are not on an edge.
  is_on_edge = false;
  closest_to_target = false;
  
  // Set the results for comparison acording to the ray direction.
  const Comparison_result point_above_under = (shoot_up ? SMALLER : LARGER);
  const Comparison_result curve_above_under = (shoot_up ? LARGER : SMALLER);

  // Keep a counter of the number of halfedges of the connected component's
  // boundary that intersect an upward (or downward, if shoot_up is false)
  // vertical ray emanating from the query point p (except for some degenerate
  // cases that are explained below).
  unsigned int              n_ray_intersections = 0;

  typename Traits_adaptor_2::Equal_2              equal =
    geom_traits->equal_2_object();
  typename Traits_adaptor_2::Compare_x_2          compare_x =
    geom_traits->compare_x_2_object();
  typename Traits_adaptor_2::Is_vertical_2        is_vertical =
    geom_traits->is_vertical_2_object();
  typename Traits_adaptor_2::Compare_y_position_2 compare_y_position =
    geom_traits->compare_y_position_2_object();
  typename Traits_adaptor_2::Compare_y_at_x_right_2  compare_y_at_x_right =
    geom_traits->compare_y_at_x_right_2_object();
  typename Traits_adaptor_2::Compare_y_at_x_left_2   compare_y_at_x_left =
    geom_traits->compare_y_at_x_left_2_object();

  // Start from the first non-vertical segment in the connected component.
  const Halfedge_const_handle        invalid_he;
  Ccb_halfedge_const_circulator      first = circ;
  bool                               found_non_vertical = false;

  do {
    // In case of a fictitious edge, check whether it is horizontal; otherwise
    // skip it.
    if (first->is_fictitious()) {
      if ((first->source()->parameter_space_in_y() != ARR_INTERIOR) &&
          (first->target()->parameter_space_in_y() != ARR_INTERIOR))
      {
        found_non_vertical = true;
        break;
      }
      else {
        ++first;
        continue;
      }
    }

    // Stop if we found a non-vertical curve.
    if (! is_vertical(first->curve())) {
      found_non_vertical = true;
      break;
    }

    if (inclusive) {
      // Check if the current vertical curve contains the query point.
      if (! first->source()->is_at_open_boundary() &&
          compare_x(p, first->source()->point()) == EQUAL &&
          top_traits->compare_y_at_x(p, &(*first)) == EQUAL)
      {
        closest_he = first;
        is_on_edge = true;
        return true;
      }
    }
    else {
      // Check if the current vertical curve contains the query point in its
      // x-range.
      if (! first->source()->is_at_open_boundary() &&
          (compare_x(p, first->source()->point()) == EQUAL))
      {
        // Check if the current vertical curve contains the query point in its
        // iterior.
        const Comparison_result  res1 = 
          top_traits->compare_xy(p, &(*(first->source())));
        const Comparison_result  res2 =
          top_traits->compare_xy(p, &(*(first->target())));
        
        if (res1 != res2) {
          if (! ((res1 == EQUAL && res2 == curve_above_under) ||
                 (res1 == curve_above_under && res2 == EQUAL)))
          {
            closest_he = first;
            is_on_edge = true;
            return true;
          }
        }
        else if (res1 == point_above_under) {
          // Check if the vertical segment is the closest to the query point so
          // far.
          if (closest_he == invalid_he ||
              (closest_he != first->twin() &&
               ((! first->source()->is_at_open_boundary() &&
                 (top_traits->compare_y_at_x(first->source()->point(),
                                            &(*closest_he)) ==
                  point_above_under)) ||
                (! first->target()->is_at_open_boundary() &&
                 (top_traits->compare_y_at_x(first->target()->point(),
                                            &(*closest_he)) ==
                  point_above_under)))))
          {
            closest_he = first;
            closest_to_target =
              ((shoot_up && first->direction() == ARR_RIGHT_TO_LEFT) ||
               (! shoot_up && first->direction() == ARR_LEFT_TO_RIGHT));
          }
        }
      }
    }

    // Move to the next curve.
    ++first;

  } while (first != circ);

  if (! found_non_vertical)
    // In this case the entire component is comprised of vertical segments,
    // so it has an empty interior and p cannot lie inside it.
    return (false);

  // Go over all curves of the boundary, starting from the non-vertical curve
  // we have located, and count those which are above p.
  Ccb_halfedge_const_circulator  curr = first;
  Comparison_result   source_res, target_res;
  Comparison_result   res;
  Comparison_result   y_res;
  bool closest_in_ccb =
    ((closest_he != invalid_he) && (closest_he->face() == circ->face()));

  source_res = top_traits->compare_x(p, &(*(curr->source())));
  do {
    // Ignore the current edge if p is not in its x-range.
    target_res = top_traits->compare_x(p, &(*(curr->target())));

    if ((source_res == target_res) && (source_res != EQUAL)) {
      ++curr;
      source_res = target_res;
      continue;
    }

    // Check whether p lies above or below the current edge.
    res = top_traits->compare_y_at_x(p, &(*curr));

    if (res == EQUAL) {
      // The current edge contains the query point. If the seach is inclusive
      // we return the edge. Otherwise, we return it only if it is vertical,
      // and contains p in its interior.
      if (inclusive) {
        closest_he = curr;
        is_on_edge = true;
        return true;
      }
      else {
        if (! curr->is_fictitious() &&
            is_vertical(curr->curve()) &&
            (curr->source()->is_at_open_boundary() ||
             ! equal(curr->source()->point(), p)) &&
            (curr->target()->is_at_open_boundary() ||
             ! equal(curr->target()->point(), p)))
        {
          closest_he = curr;
          is_on_edge = true;
          return true;
        }
      }
    }

    // If the point is above the current edge (or below it, if we shoot down),
    // move to the next edge.
    if (res != point_above_under) {
      ++curr;
      source_res = target_res;
      continue;
    }

    // Note that we do not have count intersections (actually these are
    // overlaps) of the vertical ray we shoot with vertical segments along
    // the boundary.
    if (curr->is_fictitious() || ! is_vertical(curr->curve())) {
      // The current curve is not vertical. Check the query point is in the
      // semi-open x-range (source, target] of this curve and lies below it.
      if (source_res != EQUAL) {
        if (closest_he == invalid_he || (closest_he->twin() == curr)) {
          // 1. If we have no closests halfedge, we have just found one.
          // 2. If the closest halfedge is the twin of our current halfedge,
          // we can take our halfedge to be the closest one. This covers the
          // case where our closest halfedge is not in our CCB.
          closest_he = curr;
          closest_in_ccb = true;
          closest_to_target = (target_res == EQUAL);
        }
        else {
          // Compare with the vertically closest curve so far and detemine the
          // curve closest to p. We first check the case that the two curves
          // have a common endpoint (note that the two curves do not intersect
          // in their interiors). Observe that if such a common vertex exists,
          // it is certainly not a vertex at infinity, therefore it is
          // associated with a valid point.
          if (! closest_he->is_fictitious() && ! curr->is_fictitious() &&
              ((closest_he->source() == curr->source() &&
                closest_he->direction() == curr->direction()) ||
               (closest_he->source() == curr->target() &&
                closest_he->direction() != curr->direction())))
          {
            if (closest_he->direction() == ARR_LEFT_TO_RIGHT) {
              // Both curves extend to the right from a common point.
              y_res = compare_y_at_x_right(closest_he->curve(), curr->curve(), 
                                           closest_he->source()->point());
            }
            else {
              // Both curves extend to the left from a common point.
              y_res = compare_y_at_x_left(closest_he->curve(), curr->curve(), 
                                          closest_he->source()->point());
            }
          }
          else if (! closest_he->is_fictitious() && ! curr->is_fictitious() &&
                   ((closest_he->target() == curr->source() &&
                     closest_he->direction() != curr->direction()) ||
                    (closest_he->target() == curr->target() &&
                     closest_he->direction() == curr->direction())))
          {
            if (closest_he->direction() == ARR_LEFT_TO_RIGHT) {
              // Both curves extend to the left from a common point.
              y_res = compare_y_at_x_left(closest_he->curve(), curr->curve(), 
                                          closest_he->target()->point());
            }
            else {
              // Both curves extend to the right from a common point.
              y_res = compare_y_at_x_right(closest_he->curve(), curr->curve(), 
                                           closest_he->target()->point());
            }
          }
          else {
            // In case the two curves do not have a common endpoint, but
            // overlap in their x-range (both contain p), just compare their
            // positions. Note that in this case one of the edges may be
            // fictitious, so we preform the comparsion symbolically in this
            // case.
            if (closest_he->is_fictitious())
              y_res = curve_above_under;
            else if (curr->is_fictitious())
              y_res = point_above_under;
            else
              y_res = compare_y_position(closest_he->curve(), curr->curve());
          }

          // If the current curve is closer to p than the closest curve found
          // so far, assign its halfedge as the closest one.
          if (y_res == curve_above_under) {
            closest_he = curr;
            closest_in_ccb = true;
            closest_to_target = (target_res == EQUAL);
          }
        }

        // In the degenerate case where p lies below the target vertex of
        // the current halfedge, we have to be a bit careful:
        if (target_res == EQUAL) {
          // Locate the next halfedge along the boundary that does not
          // contain a vertical segment.
          Halfedge_const_handle  next_non_vert = curr;

          do {
            next_non_vert = next_non_vert->next();

            CGAL_assertion (next_non_vert != curr);
          }
          while ((! next_non_vert->is_fictitious() &&
                  is_vertical(next_non_vert->curve())) ||
                 (next_non_vert->is_fictitious() &&
                  next_non_vert->source()->parameter_space_in_x() != 
                  next_non_vert->target()->parameter_space_in_x()));

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
          target_res = top_traits->compare_x(p, &(*(next_non_vert->target())));

          CGAL_assertion (source_res != EQUAL && target_res != EQUAL);

          if (source_res != target_res)
            n_ray_intersections++;
        }
        else {
          // In this case p lies under the interior of the current x-montone
          // curve, so the vertical ray we shoot intersects it exactly once.
          n_ray_intersections++;
        }
      }
    }
    else {
      // Assign an edge associated with a vertical curve as the closest edge
      // only if the vertical curve has an endpoint closer to p than the
      // closest edge so far.
      if (closest_he == invalid_he ||
          (! closest_in_ccb && closest_he->twin() == curr) ||
          (! curr->source()->is_at_open_boundary() &&
           top_traits->compare_y_at_x(curr->source()->point(),
                &(*closest_he)) == point_above_under) ||
          (! curr->target()->is_at_open_boundary() &&
           top_traits->compare_y_at_x(curr->target()->point(),
                &(*closest_he)) == point_above_under))
      {
        closest_he = curr;
        closest_in_ccb = true;
        closest_to_target = 
              ((shoot_up && curr->direction() == ARR_RIGHT_TO_LEFT) ||
               (! shoot_up && curr->direction() == ARR_LEFT_TO_RIGHT));
      }
    }

    // Proceed to the next halfedge along the component boundary.
    ++curr;
    source_res = target_res;

  } while (curr != first);

  // The query point lies inside the connected components if and only if the
  // ray we shoot from it intersects the boundary an odd number of times.
  return ((n_ray_intersections % 2) != 0);
}

//-----------------------------------------------------------------------------
// Find the first halfedge around a given target vertex, when going clockwise
// from "6 o'clock" around this vertex (when shooting up) or starting from
// "12 o'clock (when shooting down).
//
template <class Arrangement>
typename Arr_walk_along_line_point_location<Arrangement>::Halfedge_const_handle
Arr_walk_along_line_point_location<Arrangement>::
_first_around_vertex (Vertex_const_handle v, bool shoot_up) const
{
  // Travrse the incident halfedges of the current vertex and locate the
  // lowest one to its left and the topmost to its right.
  typename Traits_adaptor_2::Compare_y_at_x_right_2 compare_y_at_x_right =
                                geom_traits->compare_y_at_x_right_2_object();
  typename Traits_adaptor_2::Compare_y_at_x_left_2  compare_y_at_x_left =
                                geom_traits->compare_y_at_x_left_2_object();

  Halfedge_const_handle   invalid_handle;
  Halfedge_const_handle   lowest_left;
  Halfedge_const_handle   top_right;

  typename Arrangement_2::Halfedge_around_vertex_const_circulator first = 
    v->incident_halfedges();
  typename Arrangement_2::Halfedge_around_vertex_const_circulator curr = first;

  do {
    // Check whether the current halfedge is defined to the left or to the
    // right of the given vertex.
    if (curr->direction() == ARR_LEFT_TO_RIGHT) {
      // The curve associated with the current halfedge is defined to the left
      // of v.
      if (lowest_left == invalid_handle ||
          (! curr->is_fictitious() &&
           (lowest_left->is_fictitious() ||
            compare_y_at_x_left (curr->curve(),
                                 lowest_left->curve(), 
                                 v->point()) == SMALLER)))
      {
        lowest_left = curr;
      }
    }
    else {
      // The curve associated with the current halfedge is defined to the right
      // of v.
      if (top_right == invalid_handle ||
          (! curr->is_fictitious() &&
           (top_right->is_fictitious() ||
            compare_y_at_x_right (curr->curve(),
                                  top_right->curve(), 
                                  v->point()) == LARGER)))
      {
        top_right = curr;
      }
    }

    curr++;
  } while (curr != first);

  if (shoot_up) {
    // As we start from "6 o'clock" in a clockwsie direction, the first
    // halfedge we encounter is the lowest to the left. However, if there
    // is no edge to the left, we first encounter the topmost halfedge to the
    // right.
    return (lowest_left != invalid_handle) ? lowest_left : top_right;
  }
  else {
    // As we start from "12 o'clock" in a clockwsie direction, the first
    // halfedge we encounter is the topmost to the right. However, if there
    // is no edge to the right, we first encounter the lowest halfedge to the
    // left.
    return (top_right != invalid_handle) ? (top_right) : (lowest_left);
  }
}

} //namespace CGAL

#endif
