// Copyright (c) 2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Ron Wein          <wein@post.tau.ac.il>
//                 Efi Fogel         <efif@post.tau.ac.il>

#ifndef CGAL_ARR_PLANAR_INC_INSERTION_ZONE_VISITOR_H
#define CGAL_ARR_PLANAR_INC_INSERTION_ZONE_VISITOR_H

#include <CGAL/license/Arrangement_on_surface_2.h>


/*! \file
 * Definition of the arr_planar_inc_insertion_zone_visitor_2 class.
 */

#include <CGAL/Arr_accessor.h>
#include <CGAL/Arrangement_2/Arr_traits_adaptor_2.h>

namespace CGAL {

/*! \class
 * A visitor class for Arrangement_zone_2, which performs incremental
 * insertion of an x-monotone curve into an arrangement.
 * The class should be templated by an Arrangement_2 class (that is, an
 * Arrangement_on_surface_2<GeomTraits, TopTraits> class, with the TopTraits
 * being a planar topology-traits class.
 */
template <class Arrangement_>
class Arr_inc_insertion_zone_visitor
{
public:

  typedef Arrangement_                                Arrangement_2;
  typedef typename Arrangement_2::Geometry_traits_2   Geometry_traits_2;

  typedef typename Arrangement_2::Vertex_handle       Vertex_handle;
  typedef typename Arrangement_2::Halfedge_handle     Halfedge_handle;
  typedef typename Arrangement_2::Face_handle         Face_handle;

  typedef typename Arrangement_2::Point_2             Point_2;
  typedef typename Arrangement_2::X_monotone_curve_2  X_monotone_curve_2;

  typedef std::pair<Halfedge_handle, bool>            Result;

protected:

  typedef Arr_traits_basic_adaptor_2<Geometry_traits_2>  Traits_adaptor_2;

private:

  Arrangement_2         *p_arr;         // The arrangement into which we
                                        // insert the curves.
  Traits_adaptor_2      *geom_traits;   // The arrangement geometry-traits.

  const Vertex_handle    invalid_v;     // An invalid vertex handle.
  const Halfedge_handle  invalid_he;    // An invalid halfedge handle.

  X_monotone_curve_2     sub_cv1;       // Auxiliary varibale (for splitting).
  X_monotone_curve_2     sub_cv2;       // Auxiliary varibale (for splitting).

public:

  /*! Constructor. */
  Arr_inc_insertion_zone_visitor () :
    p_arr (nullptr),
    geom_traits (nullptr),
    invalid_v (),
    invalid_he ()
  {}

  /*! Initialize the visitor with an arrangement object. */
  void init (Arrangement_2 *arr)
  {
    p_arr = arr;
    geom_traits = const_cast<Traits_adaptor_2*>
      (static_cast<const Traits_adaptor_2*> (p_arr->geometry_traits()));
  }

  /*!
   * Handle the a subcurve located in the interior of a given face.
   * \param cv The subcurve.
   * \param face The face containing cv's interior.
   * \param left_v The vertex that corresponds to the left endpoint of cv
   *               (or an invalid handle if no such arrangement vertex exists).
   * \param left_he The halfedge that contains the left endpoint of cv
   *               (or an invalid handle if no such halfedge exists).
   * \param right_v The vertex that corresponds to the right endpoint of cv
   *               (or an invalid handle if no such arrangement vertex exists).
   * \param right_he The halfedge that contains the right endpoint of cv
   *                 (or an invalid handle if no such halfedge exists).
   * \return A handle to the halfedge obtained from the insertion of the
   *         subcurve into the arrangement.
   */
  Result found_subcurve (const X_monotone_curve_2& cv, Face_handle face,
                         Vertex_handle left_v, Halfedge_handle left_he,
                         Vertex_handle right_v, Halfedge_handle right_he);

  /*!
   * Handle the a subcurve that overlaps a given edge.
   * \param cv The overlapping subcurve.
   * \param he The overlapped halfedge (directed from left to right).
   * \param left_v The vertex that corresponds to the left endpoint of cv
   *               (or an invalid handle if no such arrangement vertex exists).
   * \param right_v The vertex that corresponds to the right endpoint of cv
   *               (or an invalid handle if no such arrangement vertex exists).
   * \return A handle to the halfedge obtained from the insertion of the
   *         overlapping subcurve into the arrangement.
   */
  Result found_overlap (const X_monotone_curve_2& cv, Halfedge_handle he,
                        Vertex_handle left_v, Vertex_handle right_v);

private:

  /*!
   * Split an arrangement edge.
   * \param he The edge to split (one of the twin halfedges).
   * \param p The split point.
   * \param arr_access An arrangement accessor.
   */
  void _split_edge (Halfedge_handle he, const Point_2& p,
                    Arr_accessor<Arrangement_2>& arr_access);

};

//-----------------------------------------------------------------------------
// Memeber-function definitions:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Handle the a subcurve located in the interior of a given face.
//
template <class Arrangement>
typename Arr_inc_insertion_zone_visitor<Arrangement>::Result
Arr_inc_insertion_zone_visitor<Arrangement>::
found_subcurve (const X_monotone_curve_2& cv, Face_handle face,
                Vertex_handle left_v, Halfedge_handle left_he,
                Vertex_handle right_v, Halfedge_handle right_he)
{
  typename Traits_adaptor_2::Construct_min_vertex_2 min_vertex =
    geom_traits->construct_min_vertex_2_object();
  typename Traits_adaptor_2::Construct_max_vertex_2 max_vertex =
    geom_traits->construct_max_vertex_2_object();

  typename Traits_adaptor_2::Parameter_space_in_x_2 ps_in_x =
    geom_traits->parameter_space_in_x_2_object();
  typename Traits_adaptor_2::Parameter_space_in_y_2 ps_in_y =
    geom_traits->parameter_space_in_y_2_object();

  // Create an arrangement accessor.
  Arr_accessor<Arrangement_2>    arr_access (*p_arr);

  // Get the boundary conditions of the curve ends.
  const Arr_parameter_space bx_l = ps_in_x (cv, ARR_MIN_END);
  const Arr_parameter_space by_l = ps_in_y (cv, ARR_MIN_END);

  const Arr_parameter_space bx_r = ps_in_x (cv, ARR_MAX_END);
  const Arr_parameter_space by_r = ps_in_y (cv, ARR_MAX_END);

  // Check if the left and the right endpoints of cv should be associated
  // with arrangement vertices.
  const bool vertex_for_left =
    (left_v != invalid_v) || (left_he != invalid_he);
  const bool vertex_for_right =
    (right_v != invalid_v) || (right_he != invalid_he);

  // Find the previous halfedges for the left and right endpoints (if any).
  Halfedge_handle  prev_he_left;
  Halfedge_handle  prev_he_right;

  if (vertex_for_left) {
    // If we are given the previous halfedge, use it. Otherwise, we are given
    // the vertex and we should locate cv around it.
    if (left_he != invalid_he)
      prev_he_left = left_he;
    else if (! left_v->is_isolated())
      prev_he_left = arr_access.locate_around_vertex (left_v, cv);

    // In case the vertex does not exist, split left_he at cv's left endpoint
    // and create the vertex.
    if (left_v == invalid_v) {
      _split_edge (left_he, min_vertex (cv), arr_access);

      // Check whether we have just split the halfedge that right_he refers
      // to, and if this halfedge is directed from left to right.
      // If so, right_he's target is now the new vertex, and we have to
      // proceed to the next halfedge (which is going to be split).
      if ((right_he == left_he) && (left_he->direction() == ARR_LEFT_TO_RIGHT))
        right_he = right_he->next();
    }
  }
  else {
    // Check whether the left end of cv is bounded or not.
    if ((bx_l == ARR_LEFT_BOUNDARY) || (by_l != ARR_INTERIOR)) {
      // Use the arrangement accessor and obtain a vertex associated with
      // the unbounded left end (possibly with a predecessor halfedge).
      std::pair<Vertex_handle, Halfedge_handle>  pos =
        arr_access.place_and_set_curve_end (face, cv, ARR_MIN_END, bx_l, by_l);

      if (pos.second != invalid_he) {
        // Use the predecessor halfedge, if possible:
        prev_he_left = pos.second;

        /* The following cannot happen, because a vertex at infinity has a
         * single incident edge.
         *
         * Check whether we have just split the halfedge that right_he refers
         * to, and if this halfedge is directed from left to right.
         * If so, right_he's target is now the new vertex, and we have to
         * proceed to the next halfedge (which is going to be split).
         * if ((right_he == prev_he_left) &&
         *     (prev_he_left->direction() == ARR_LEFT_TO_RIGHT))
         *   right_he = right_he->next();
         */
      } else
        // Use just the isolated vertex:
        left_v = pos.first;
    }
  }

  if (vertex_for_right) {
    // If we are given the previous halfedge, use it. Otherwise, we are given
    // the vertex and we should locate cv around it.
    if (right_he != invalid_he)
      prev_he_right = right_he;
    else if (! right_v->is_isolated())
      prev_he_right = arr_access.locate_around_vertex (right_v, cv);

    // In case the vertex does not exist, split right_he at cv's right
    // endpoint and create the vertex.
    if (right_v == invalid_v) {
      _split_edge (right_he, max_vertex (cv), arr_access);

      // Check whether we have just split the halfedge that left_he refers to.
      // If so, prev_he_right's target is now the new vertex, and we have to
      // proceed to the next halfedge (whose target is right_v).
      if (right_he == prev_he_left)
        prev_he_left = prev_he_left->next();
    }
  }
  else {
    // Check whether the right end of cv is bounded or not.
    if ((bx_r == ARR_RIGHT_BOUNDARY) || (by_r != ARR_INTERIOR)) {
      // Use the arrangement accessor and obtain a vertex associated with
      // the unbounded right end (possibly with a predecessor halfedge).
      std::pair<Vertex_handle, Halfedge_handle>  pos =
        arr_access.place_and_set_curve_end (face, cv, ARR_MAX_END, bx_r, by_r);

      if (pos.second != invalid_he) {
        // Use the predecessor halfedge, if possible:
        prev_he_right = pos.second;

        // Check whether we have just split the halfedge that left_he refers to.
        // If so, prev_he_right's target is now the new vertex, and we have to
        // proceed to the next halfedge (whose target is right_v).
        if (prev_he_right == prev_he_left)
          prev_he_left = prev_he_left->next();
      } else
        // Use just the isolated vertex:
        right_v = pos.first;
    }
  }

  // Insert the curve into the arrangement.
  Halfedge_handle   inserted_he;

  if (prev_he_left == invalid_he) {
    // The left endpoint is associated with an isolated vertex, or is not
    // associated with any vertex. In the latter case, we create such a vertex
    // now.
    if (left_v == invalid_v) {
      left_v = ((bx_l == ARR_INTERIOR) && (by_l == ARR_INTERIOR)) ?
        arr_access.create_vertex (min_vertex (cv)) :
        arr_access.create_boundary_vertex (cv, ARR_MIN_END, bx_l, by_l);
    }

    if (prev_he_right == invalid_he) {
      // The right endpoint is associated with an isolated vertex, or is not
      // associated with any vertex. In the latter case, we create such a
      // vertex now.
      if (right_v == invalid_v) {
        right_v = ((bx_r == ARR_INTERIOR) && (by_r == ARR_INTERIOR)) ?
          arr_access.create_vertex (max_vertex (cv)) :
          arr_access.create_boundary_vertex (cv, ARR_MAX_END, bx_r, by_r);
      }

      // We should insert the curve in the interior of the face.
      inserted_he = arr_access.insert_in_face_interior_ex (face, cv, ARR_LEFT_TO_RIGHT,
                                                           left_v, right_v);
    }
    else {
      // The right endpoint is associated with an arrangement vertex, and
      // we have the predecessor halfedge for the insertion.
      inserted_he = arr_access.insert_from_vertex_ex (prev_he_right, cv, ARR_RIGHT_TO_LEFT,
                                                      left_v);

      // The returned halfedge is directed to the newly created vertex
      // (the left one), so we take its twin.
      inserted_he = inserted_he->twin();
    }
  }
  else {
    // We have a vertex associated with the left end of the curve, along
    // with a predecessor halfedge.
    if (prev_he_right == invalid_he) {
      // The right endpoint is associated with an isolated vertex, or is not
      // associated with any vertex. In the latter case, we create such a
      // vertex now.
      if (right_v == invalid_v) {
        right_v = ((bx_r == ARR_INTERIOR) && (by_r == ARR_INTERIOR)) ?
          arr_access.create_vertex (max_vertex (cv)) :
          arr_access.create_boundary_vertex (cv, ARR_MAX_END, bx_r, by_r);
      }

      // Use the left predecessor for the insertion.
      inserted_he = arr_access.insert_from_vertex_ex (prev_he_left, cv, ARR_LEFT_TO_RIGHT,
                                                      right_v);
    }
    else {
      // The right endpoint is associated with an arrangement vertex, and
      // we have the predecessor halfedge for the insertion.
      CGAL_assertion (prev_he_right != invalid_he);

      // Perform the insertion using the predecessor halfedges.
      inserted_he = p_arr->insert_at_vertices (cv, prev_he_left, prev_he_right);
    }
  }

  // Return the inserted halfedge, and indicate we should not halt the
  // zone-computation process.
  return Result (inserted_he, false);
}

//-----------------------------------------------------------------------------
// Handle the a subcurve located in the interior of a given face.
//
template <class Arrangement>
typename Arr_inc_insertion_zone_visitor<Arrangement>::Result
Arr_inc_insertion_zone_visitor<Arrangement>::
found_overlap (const X_monotone_curve_2& cv, Halfedge_handle he,
               Vertex_handle left_v, Vertex_handle right_v)
{
  // Get the boundary conditions of the curve ends.
  const Arr_parameter_space bx_l =
    geom_traits->parameter_space_in_x_2_object()(cv, ARR_MIN_END);
  const Arr_parameter_space by_l =
    geom_traits->parameter_space_in_y_2_object()(cv, ARR_MIN_END);

  const Arr_parameter_space bx_r =
    geom_traits->parameter_space_in_x_2_object()(cv, ARR_MAX_END);
  const Arr_parameter_space by_r =
    geom_traits->parameter_space_in_y_2_object()(cv, ARR_MAX_END);

  // Modify (perhaps split) the overlapping arrangement edge.
  Halfedge_handle   updated_he;

  if (left_v == invalid_v &&
      ! ((bx_l == ARR_LEFT_BOUNDARY) || (by_l != ARR_INTERIOR)))
  {
    // Split the curve associated with he at the left endpoint of cv.
    geom_traits->split_2_object()
      (he->curve(),
       geom_traits->construct_min_vertex_2_object() (cv), sub_cv1, sub_cv2);

    if (right_v == invalid_v &&
        ! ((bx_r == ARR_RIGHT_BOUNDARY) || (by_r != ARR_INTERIOR)))
    {
      // The overlapping curve is contained strictly in the interior of he:
      // Split he as an intermediate step.
      updated_he = p_arr->split_edge (he, sub_cv1, sub_cv2);
      updated_he = updated_he->next();

      // Split the left subcurve at the right endpoint of cv.
      geom_traits->split_2_object()
        (updated_he->curve(),
         geom_traits->construct_max_vertex_2_object() (cv), sub_cv1, sub_cv2);

      // Split updated_he once again, so that the left portion corresponds
      // to the overlapping curve and the right portion corresponds to
      // sub_cv2.
      updated_he = p_arr->split_edge (updated_he, cv, sub_cv2);
    }
    else {
      // Split he, such that the left portion corresponds to sub_cv1 and the
      // right portion corresponds to the overlapping curve.
      updated_he = p_arr->split_edge (he, sub_cv1, cv);
      updated_he = updated_he->next();
    }
  }
  else {
    if (right_v == invalid_v &&
        ! ((bx_r == ARR_RIGHT_BOUNDARY) || (by_r != ARR_INTERIOR)))
    {
      // Split the curve associated with he at the right endpoint of cv.
      geom_traits->split_2_object()
        (he->curve(),
         geom_traits->construct_max_vertex_2_object() (cv), sub_cv1, sub_cv2);

      // Split he, such that the left portion corresponds to the overlapping
      // curve and the right portion corresponds to sub_cv2.
      updated_he = p_arr->split_edge (he, cv, sub_cv2);
    }
    else {
      // The entire edge is overlapped: Modify the curve associated with cv
      // to be the overlapping curve.
      updated_he = p_arr->modify_edge (he, cv);
    }
  }

  // Return the updated halfedge, and indicate we should not halt the
  // zone-computation process.
  return Result (updated_he, false);
}

//-----------------------------------------------------------------------------
// Split an arrangement edge.
//
template <class Arrangement>
void Arr_inc_insertion_zone_visitor<Arrangement>::
_split_edge (Halfedge_handle he, const Point_2& p,
             Arr_accessor<Arrangement_2>& arr_access)
{
  // Split the curve at the split point.
  geom_traits->split_2_object() (he->curve(), p, sub_cv1, sub_cv2);

  // sub_cv1 is always to the left of the split point p and sub_cv2 lies to
  // its right. Thus, if the split edge is directed from left to right then
  // left end of sub_cv1 equals he's source, and if the edge is directed from
  // right to left, we have to reverse the subcurve order.
  if (he->direction() == ARR_LEFT_TO_RIGHT)
    arr_access.split_edge_ex (he, p, sub_cv1, sub_cv2);
  else
    arr_access.split_edge_ex (he, p, sub_cv2, sub_cv1);
}

} //namespace CGAL

#endif
