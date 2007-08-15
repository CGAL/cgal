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
// Author(s)     : Efi Fogel         <efif@post.tau.ac.il>

#ifndef CGAL_ARR_SPHERICAL_INC_INSERTION_ZONE_VISITOR_H
#define CGAL_ARR_SPHERICAL_INC_INSERTION_ZONE_VISITOR_H

/*! \file
 * Definition of the arr_spherical_inc_insertion_zone_visitor_2 class.
 */

#include <CGAL/Arr_accessor.h>
#include <CGAL/Arrangement_2/Arr_traits_adaptor_2.h>

CGAL_BEGIN_NAMESPACE

/*! A visitor class for Arrangement_zone_2, which performs incremental 
 * insertion of an x-monotone curve into an arrangement.
 * The class should be templated by an Arrangement_2 class (that is, an
 * Arrangement_on_surface_2<GeomTraits, TopTraits> class, with the TopTraits
 * being a spherical topology-traits class.
 */
template <typename T_Arrangement>
class Arr_spherical_inc_insertion_zone_visitor {
public:
  typedef T_Arrangement                                 Arrangement_2;
  typedef typename Arrangement_2::Geometry_traits_2     Geometry_traits_2;

  typedef typename Arrangement_2::Vertex_handle         Vertex_handle;
  typedef typename Arrangement_2::Halfedge_handle       Halfedge_handle;
  typedef typename Arrangement_2::Face_handle           Face_handle;

  typedef typename Arrangement_2::Point_2               Point_2;
  typedef typename Arrangement_2::X_monotone_curve_2    X_monotone_curve_2;

  typedef std::pair<Halfedge_handle, bool>              Result;

protected:
  typedef Arr_traits_basic_adaptor_2<Geometry_traits_2> Traits_adaptor_2;

private:
  //! The arrangement to insert curves into
  Arrangement_2 * m_arr;

  //! The arrangement geometry-traits.
  Traits_adaptor_2 * m_traits;

  //! An invalid vertex handle.
  const Vertex_handle m_invalid_v;

  //! An invalid halfedge handle.
  const Halfedge_handle m_invalid_he;

  //! Auxiliary varibale (for splitting).
  X_monotone_curve_2 m_sub_xc1;

  //! Auxiliary varibale (for splitting).
  X_monotone_curve_2 m_sub_xc2;

public:

  /*! Constructor */
  Arr_spherical_inc_insertion_zone_visitor() :
    m_arr(NULL),
    m_traits(NULL),
    m_invalid_v(),
    m_invalid_he()
  {}

  /*! Initialize the visitor with an arrangement object.
   * \param arr the arrangement
   */
  void init(Arrangement_2 * arr)
  {
    m_arr = arr;
    m_traits = const_cast<Traits_adaptor_2*> 
      (static_cast<const Traits_adaptor_2*>(m_arr->geometry_traits()));
  }

  /*! Handle the a subcurve located in the interior of a given face.
   * \param xc The subcurve.
   * \param face The face containing xc's interior.
   * \param left_v The vertex that corresponds to the left endpoint of xc
   *               (or an invalid handle if no such arrangement vertex exists).
   * \param left_he The halfedge that contains the left endpoint of xc
   *               (or an invalid handle if no such halfedge exists).
   * \param right_v The vertex that corresponds to the right endpoint of xc
   *               (or an invalid handle if no such arrangement vertex exists).
   * \param right_he The halfedge that contains the right endpoint of xc
   *                 (or an invalid handle if no such halfedge exists).
   * \return A handle to the halfedge obtained from the insertion of the
   *         subcurve into the arrangement.
   */
  Result found_subcurve(const X_monotone_curve_2 & xc,
                        Face_handle face,
                        Vertex_handle left_v, Halfedge_handle left_he,
                        Vertex_handle right_v, Halfedge_handle right_he);

  /*! Handle the a subcurve that overlaps a given edge.
   * \param xc The overlapping subcurve.
   * \param he The overlapped halfedge (directed from left to right).
   * \param left_v The vertex that corresponds to the left endpoint of xc
   *               (or an invalid handle if no such arrangement vertex exists).
   * \param right_v The vertex that corresponds to the right endpoint of xc
   *               (or an invalid handle if no such arrangement vertex exists).
   * \return A handle to the halfedge obtained from the insertion of the
   *         overlapping subcurve into the arrangement.
   */
  Result found_overlap(const X_monotone_curve_2 & xc, Halfedge_handle he,
                       Vertex_handle left_v, Vertex_handle right_v);

private:

  /*! Split an arrangement edge.
   * \param he The edge to split (one of the twin halfedges).
   * \param p The split point.
   * \param arr_access An arrangement accessor.
   */
  void _split_edge(Halfedge_handle he, const Point_2 & p,
                   Arr_accessor<Arrangement_2> & arr_access);

};

//! \brief handles the a subcurve located in the interior of a given face.
template <class Arrangement>
typename Arr_spherical_inc_insertion_zone_visitor<Arrangement>::Result
Arr_spherical_inc_insertion_zone_visitor<Arrangement>::
found_subcurve(const X_monotone_curve_2 & xc, Face_handle face,
               Vertex_handle left_v, Halfedge_handle left_he,
               Vertex_handle right_v, Halfedge_handle right_he)
{
  // Create an arrangement accessor.
  Arr_accessor<Arrangement_2> arr_access(*m_arr);
  
  // Get the boundary conditions of the curve ends.
  const Boundary_type bx_l = m_traits->boundary_in_x_2_object()(xc, MIN_END);
  const Boundary_type by_l = m_traits->boundary_in_y_2_object()(xc, MIN_END);

  const Boundary_type bx_r = m_traits->boundary_in_x_2_object()(xc, MAX_END);
  const Boundary_type by_r = m_traits->boundary_in_y_2_object()(xc, MAX_END);

  // Endpoints of xc should always be associated with arrangement vertices.
  const bool vertex_for_left = (left_v != m_invalid_v) ||
    (left_he != m_invalid_he);
  const bool vertex_for_right = (right_v != m_invalid_v) ||
    (right_he != m_invalid_he);

  CGAL_assertion(vertex_for_left);
  CGAL_assertion(vertex_for_right);
  
  // Find the previous halfedges for the left and right endpoints (if any).
  Halfedge_handle prev_he_left;
  Halfedge_handle prev_he_right;

  if (vertex_for_left) {
    // If we are given the previous halfedge, use it. Otherwise, we are given
    // the vertex and we should locate xc around it.
    if (left_he != m_invalid_he) prev_he_left = left_he;
    else if (!left_v->is_isolated())
      prev_he_left = arr_access.locate_around_vertex (left_v, xc);

    // In case the vertex does not exist, split left_he at xc's left endpoint
    // and create the vertex.
    if (left_v == m_invalid_v) {
      _split_edge(left_he, m_traits->construct_min_vertex_2_object()(xc),
                  arr_access);

      // Check if we have just split the halfedge that right_he refers to,
      // and if this halfedge is directed from left to right.
      // If so, right_he's target is now the new vertex, and we have to
      // proceed to the next halfedge (which is going to be split).
      if (right_he == left_he && left_he->direction() == LEFT_TO_RIGHT)
        right_he = right_he->next();
    }
  } else {
    // Check if the left end of xc is bounded of not.    
    if (bx_l == MINUS_INFINITY ||
        by_l == MINUS_INFINITY || by_l == PLUS_INFINITY)
    {
      // Use the arrangement accessor and obtain a vertex associated with
      // the unbounded left end (possibly with a predecessor halfedge).
      std::pair<Vertex_handle, Halfedge_handle> pos =
        arr_access.place_and_set_curve_end(face, xc, MIN_END, bx_l, by_l);

      // Use the predecessor halfedge, if possible, just the isolated vertex:
      if (pos.second != m_invalid_he) prev_he_left = pos.second;
      else left_v = pos.first;
    }
  }

  if (vertex_for_right) {
    // If we are given the previous halfedge, use it. Otherwise, we are given
    // the vertex and we should locate xc around it.
    if (right_he != m_invalid_he) prev_he_right = right_he;
    else if (!right_v->is_isolated())
      prev_he_right = arr_access.locate_around_vertex(right_v, xc);
    
    // In case the vertex does not exist, split right_he at xc's right
    // endpoint and create the vertex.
    if (right_v == m_invalid_v) {
      _split_edge(right_he, m_traits->construct_max_vertex_2_object()(xc),
                  arr_access);
      
      // Check if we have just split the halfedge that left_he refers to.
      // If so, prev_he_right's target is now the new vertex, and we have to
      // proceed to the next halfedge (whose target is right_v).
      if (right_he == prev_he_left) {
        prev_he_left = prev_he_left->next();
      }
    }
  }
  else {
    // Check if the right end of xc is bounded of not.    
    if (bx_r == PLUS_INFINITY ||
        by_r == MINUS_INFINITY || by_r == PLUS_INFINITY)
    {
      // Use the arrangement accessor and obtain a vertex associated with
      // the unbounded right end (possibly with a predecessor halfedge).
      std::pair<Vertex_handle, Halfedge_handle> pos =
        arr_access.place_and_set_curve_end(face, xc, MAX_END, bx_r, by_r);

      // Use the predecessor halfedge, if possible, or just the isolated vertex:
      if (pos.second != m_invalid_he) prev_he_right = pos.second;
      else right_v = pos.first;
    }
  }

  // Insert the curve into the arrangement.
  Halfedge_handle inserted_he;

  if (prev_he_left == m_invalid_he) {
    // The left endpoint is associated with an isolated vertex, or is not
    // associated with any vertex. In the latter case, we create such a vertex
    // now.
    if (left_v == m_invalid_v) {
      left_v =
        arr_access.create_vertex(m_traits->construct_min_vertex_2_object()(xc),
                                 bx_l, by_l);
    }

    if (prev_he_right == m_invalid_he) {
      // The right endpoint is associated with an isolated vertex, or is not
      // associated with any vertex. In the latter case, we create such a
      // vertex now.
      if (right_v == m_invalid_v) {
        right_v =
          arr_access.create_vertex
          (m_traits->construct_max_vertex_2_object()(xc), bx_r, by_r);
      }
     
      // We should insert the curve in the interior of the face.
      inserted_he = arr_access.insert_in_face_interior_ex(xc, face, left_v,
                                                          right_v, SMALLER);
    }
    else {
      // The right endpoint is associated with an arrangement vertex, and
      // we have the predecessor halfedge for the insertion.
      inserted_he = arr_access.insert_from_vertex_ex(xc, prev_he_right,
                                                     left_v, LARGER);

      // The returned halfedge is directed to the newly created vertex
      // (the left one), so we take its twin.
      inserted_he = inserted_he->twin();
    }
  } else {
    // We have a vertex associated with the left end of the curve, along
    // with a predecessor halfedge.
    if (prev_he_right == m_invalid_he) {
      // The right endpoint is associated with an isolated vertex, or is not
      // associated with any vertex. In the latter case, we create such a
      // vertex now.
      if (right_v == m_invalid_v) {
        right_v =
          arr_access.create_vertex
          (m_traits->construct_max_vertex_2_object()(xc), bx_r, by_r);
      }
     
      // Use the left predecessor for the insertion.
      inserted_he = arr_access.insert_from_vertex_ex(xc, prev_he_left, right_v,
                                                     SMALLER);
    } else {
      // The right endpoint is associated with an arrangement vertex, and
      // we have the predecessor halfedge for the insertion.
      CGAL_assertion(prev_he_right != m_invalid_he);

      // Perform the insertion using the predecessor halfedges.
      inserted_he = m_arr->insert_at_vertices(xc, prev_he_left, prev_he_right);
    }
  }

  // Return the inserted halfedge, and indicate we should not halt the
  // zone-computation process.
  return Result(inserted_he, false);
}

//-----------------------------------------------------------------------------
// Handle the a subcurve located in the interior of a given face.
//
template <class Arrangement>
typename Arr_spherical_inc_insertion_zone_visitor<Arrangement>::Result
Arr_spherical_inc_insertion_zone_visitor<Arrangement>::
found_overlap(const X_monotone_curve_2 & xc, Halfedge_handle he,
              Vertex_handle left_v, Vertex_handle right_v)
{
  // Modify (perhaps split) the overlapping arrangement edge.
  Halfedge_handle updated_he;

  if (left_v == m_invalid_v) {
    // Split the curve associated with he at the left endpoint of xc.
    m_traits->split_2_object()(he->curve(),
                               m_traits->construct_min_vertex_2_object()(xc),
                               m_sub_xc1, m_sub_xc2);

    if (right_v == m_invalid_v) {
      // The overlapping curve is contained strictly in the interior of he:
      // Split he as an intermediate step.
      updated_he = m_arr->split_edge(he, m_sub_xc1, m_sub_xc2);
      updated_he = updated_he->next();

      // Split the left subcurve at the right endpoint of xc.
      m_traits->split_2_object()(updated_he->curve(),
                                 m_traits->construct_max_vertex_2_object()(xc),
                                 m_sub_xc1, m_sub_xc2);

      // Split updated_he once again, so that the left portion corresponds
      // to the overlapping curve and the right portion corresponds to
      // m_sub_xc2.
      updated_he = m_arr->split_edge(updated_he, xc, m_sub_xc2);
    } else {
      // Split he, such that the left portion corresponds to m_sub_xc1 and the
      // right portion corresponds to the overlapping curve.
      updated_he = m_arr->split_edge(he, m_sub_xc1, xc);
      updated_he = updated_he->next();
    }
  } else {
    if (right_v == m_invalid_v) {
      // Split the curve associated with he at the right endpoint of xc.
      m_traits->split_2_object()(he->curve(),
                                 m_traits->construct_max_vertex_2_object()(xc),
                                 m_sub_xc1, m_sub_xc2);

      // Split he, such that the left portion corresponds to the overlapping
      // curve and the right portion corresponds to m_sub_xc2.
      updated_he = m_arr->split_edge(he, xc, m_sub_xc2);
    } else {
      // The entire edge is overlapped: Modify the curve associated with xc
      // to be the overlapping curve.
      updated_he = m_arr->modify_edge(he, xc);
    }
  }

  // Return the updated halfedge, and indicate we should not halt the
  // zone-computation process.
  return Result(updated_he, false);
}

//-----------------------------------------------------------------------------
// Split an arrangement edge.
//
template <class Arrangement>
void Arr_spherical_inc_insertion_zone_visitor<Arrangement>::
_split_edge(Halfedge_handle he, const Point_2 & p,
            Arr_accessor<Arrangement_2> & arr_access)
{
  // Split the curve at the split point.
  m_traits->split_2_object()(he->curve(), p, m_sub_xc1, m_sub_xc2);

  // m_sub_xc1 is always to the left of the split point p and m_sub_xc2 lies to
  // its right. Thus, if the split edge is directed from left to right then
  // left end of m_sub_xc1 equals he's source, and if the edge is directed from
  // right to left, we have to reverse the subcurve order.
  if (he->direction() == LEFT_TO_RIGHT)
    arr_access.split_edge_ex(he, p, m_sub_xc1, m_sub_xc2);
  else arr_access.split_edge_ex(he, p, m_sub_xc2, m_sub_xc1);
}

CGAL_END_NAMESPACE

#endif
