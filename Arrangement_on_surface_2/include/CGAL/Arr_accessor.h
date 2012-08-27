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
// Author(s)     : Ron Wein          <wein@post.tau.ac.il>
//                 Efi Fogel         <efif@post.tau.ac.il>

#ifndef CGAL_ARR_ACCESSOR_H
#define CGAL_ARR_ACCESSOR_H

/*! \file
 * Definition of the Arr_accessor<Arrangement> class.
 */

#include <CGAL/Arrangement_2/Arr_traits_adaptor_2.h>

namespace CGAL {

/*! \class
 * A class that provides access to some of the internal arrangement operations.
 * Used mostly by the global insertion functions and by the sweep-line visitors
 * for utilizing topological and geometrical information available during the
 * algorithms they perform.
 * The Arrangement parameter corresponds to an arrangement instantiation
 * (of the template Arrangement_on_surface_2).
 */
template <class Arrangement_>
class Arr_accessor
{
public:

  typedef Arrangement_                                  Arrangement_2;
  typedef Arr_accessor<Arrangement_2>                   Self;

  typedef typename Arrangement_2::Size                  Size;
  typedef typename Arrangement_2::Point_2               Point_2;
  typedef typename Arrangement_2::X_monotone_curve_2    X_monotone_curve_2;

  typedef typename Arrangement_2::Vertex_handle         Vertex_handle;
  typedef typename Arrangement_2::Vertex_const_handle   Vertex_const_handle;
  typedef typename Arrangement_2::Halfedge_handle       Halfedge_handle;
  typedef typename Arrangement_2::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Arrangement_2::Face_handle           Face_handle;
  typedef typename Arrangement_2::Face_const_handle     Face_const_handle;
  typedef typename Arrangement_2::Ccb_halfedge_circulator
                                                        Ccb_halfedge_circulator;

private:

  typedef typename Arrangement_2::DVertex               DVertex;
  typedef typename Arrangement_2::DHalfedge             DHalfedge;
  typedef typename Arrangement_2::DFace                 DFace;
  typedef typename Arrangement_2::DOuter_ccb            DOuter_ccb;
  typedef typename Arrangement_2::DInner_ccb            DInner_ccb;
  typedef typename Arrangement_2::DIso_vertex           DIso_vertex;

private:

  Arrangement_2  *p_arr;           // The associated arrangement.

public:

  /*! Constructor with an associated arrangement. */
  Arr_accessor (Arrangement_2& arr) :
    p_arr (&arr)
  {}

  /* Get the arrangement. */
  Arrangement_2& arrangement ()
  {
    return (*p_arr);
  }

  /* Get the arrangement (const version). */
  const Arrangement_2& arrangement() const
  {
    return (*p_arr);
  }

  /// \name Accessing the notification functions (for the global functions).
  //@{

  /*! Notify that a global operation is about to take place. */
  void notify_before_global_change ()
  {
    p_arr->_notify_before_global_change();
  }

  /*! Notify that a global operation was completed. */
  void notify_after_global_change ()
  {
    p_arr->_notify_after_global_change();
  }
  //@}

  /// \name Local operations and predicates for the arrangement.
  //@{

  /*!
   * Locate the arrangement feature that contains the given curve-end.
   * \param cv The curve.
   * \param ind ARR_MIN_END if we refer to cv's minimal end;
   *            ARR_MAX_END if we refer to its maximal end.
   * \param ps_x The boundary condition in x.
   * \param ps_y The boundary condition in y.
   * \pre The relevant end of cv has boundary conditions in x or in y.
   * \return An object that contains the curve end.
   *         This object may wrap a Face_const_handle (the general case),
   *         or a Halfedge_const_handle (in case of an overlap).
   */
  CGAL::Object locate_curve_end (const X_monotone_curve_2& cv,
                                 Arr_curve_end ind,
                                 Arr_parameter_space ps_x,
                                 Arr_parameter_space ps_y) const
  {
    CGAL_precondition (ps_x != ARR_INTERIOR || ps_y != ARR_INTERIOR);

    // Use the topology traits to locate the unbounded curve end.
    CGAL::Object  obj =
      p_arr->topology_traits()->locate_curve_end (cv, ind,
                                                  ps_x, ps_y);

    // Return a handle to the DCEL feature.
    DFace        *f;

    if (CGAL::assign (f, obj))
      return (CGAL::make_object (p_arr->_const_handle_for (f)));

    DHalfedge    *he;

    if (CGAL::assign (he, obj))
      return (CGAL::make_object (p_arr->_const_handle_for (he)));

    DVertex      *v;

    if (CGAL::assign (v, obj))
      return (CGAL::make_object (p_arr->_const_handle_for (v)));

    // We should never reach here:
    CGAL_error();
    return Object();
  }

  /*!
   * Locate the place for the given curve around the given vertex.
   * \param vh A handle for the arrangement vertex.
   * \param cv The given x-monotone curve.
   * \pre v is one of cv's endpoints.
   * \return A handle for a halfedge whose target is v, where cv should be
   *         inserted between this halfedge and the next halfedge around this
   *         vertex (in a clockwise order).
   */
  Halfedge_handle locate_around_vertex (Vertex_handle vh,
                                        const X_monotone_curve_2& cv) const
  {
    typedef
      Arr_traits_basic_adaptor_2<typename Arrangement_2::Geometry_traits_2>
      Traits_adaptor_2;

    const Traits_adaptor_2  *m_traits = 
      static_cast<const Traits_adaptor_2*> (p_arr->geometry_traits());

    Arr_curve_end                ind = ARR_MIN_END;

    if (m_traits->is_closed_2_object() (cv, ARR_MAX_END) &&
        m_traits->equal_2_object() (vh->point(),
                                    m_traits->construct_max_vertex_2_object()(cv)))
    {
      ind = ARR_MAX_END;
    }

    DHalfedge * he = p_arr->_locate_around_vertex(p_arr->_vertex (vh), cv, ind);

    CGAL_assertion (he != NULL);
    return (p_arr->_handle_for (he));
  }

  /*!
   * Locate the place for the given curve-end around the given vertex,
   * which lies on the boundary.
   * \param vh A handle for the arrangement vertex.
   * \param cv The curve.
   * \param ind ARR_MIN_END if we refer to cv's minimal end;
   *            ARR_MAX_END if we refer to its maximal end.
   * \param ps_x The boundary condition in x.
   * \param ps_y The boundary condition in y.
   * \pre The relevant end of cv has boundary conditions in x or in y.
   * \return A handle for a halfedge whose target is v, where cv should be
   *         inserted between this halfedge and the next halfedge around this
   *         vertex (in a clockwise order).
   */
  Halfedge_handle
      locate_around_boundary_vertex (Vertex_handle vh,
                                     const X_monotone_curve_2& cv,
                                     Arr_curve_end ind,
                                     Arr_parameter_space ps_x,
                                     Arr_parameter_space ps_y) const
  {
    CGAL_precondition (ps_x != ARR_INTERIOR || ps_y != ARR_INTERIOR);

    // Use the topology traits to locate the unbounded curve end.
    DHalfedge*  he = p_arr->topology_traits()->
                         locate_around_boundary_vertex (p_arr->_vertex (vh),
                                                        cv, ind, ps_x, ps_y);

    CGAL_assertion (he != NULL);
    return (p_arr->_handle_for (he));
  }

  /*!
   * Compute the distance (in halfedges) between two halfedges.
   * \param e1 A handle for the source halfedge.
   * \param e2 A handle for the destination halfedge.
   * \return In case e1 and e2 belong to the same connected component, the 
   *         function returns number of boundary halfedges between the two 
   *         halfedges. Otherwise, it returns (-1).
   */
  int halfedge_distance (Halfedge_const_handle e1,
                         Halfedge_const_handle e2) const
  {
    // If the two halfedges do not belong to the same component, return (-1).
    const DHalfedge  *he1 = p_arr->_halfedge (e1);
    const DHalfedge  *he2 = p_arr->_halfedge (e2);
    
    if (he1 == he2)
      return (0);

    const DInner_ccb *ic1 = (he1->is_on_inner_ccb()) ? he1->inner_ccb() : NULL;
    const DOuter_ccb *oc1 = (ic1 == NULL) ? he1->outer_ccb() : NULL;
    const DInner_ccb *ic2 = (he2->is_on_inner_ccb()) ? he2->inner_ccb() : NULL;
    const DOuter_ccb *oc2 = (ic2 == NULL) ? he2->outer_ccb() : NULL;

    if (oc1 != oc2 || ic1 != ic2)
      return (-1);

    // Compute the distance between the two halfedges.
    unsigned int         dist = p_arr->_halfedge_distance (he1, he2);
    return (static_cast<int> (dist));
  }

  /*!
   * Determine whether a given query halfedge lies in the interior of a new
   * face we are about to create, by connecting it with another halfedge
   * using a given x-monotone curve.
   * \param prev1 A handle for the query halfedge.
   * \param prev2 The other halfedge we are about to connect with prev1.
   * \param cv The x-monotone curve we use to connect prev1 and prev2.
   * \pre prev1 and prev2 belong to the same connected component, and by
   *      connecting them using cv we form a new face.
   * \return (true) if prev1 lies in the interior of the face we are about
   *         to create, (false) otherwise - in which case prev2 must lie
   *         inside this new face.
   */
  bool is_inside_new_face (Halfedge_handle prev1,
                           Halfedge_handle prev2,
                           const X_monotone_curve_2& cv) const
  {
    return (p_arr->_is_inside_new_face (p_arr->_halfedge (prev1),
                                        p_arr->_halfedge (prev2),
                                        cv));
  }

  /*!
   * Check if the given vertex represents one of the ends of a given curve.
   * \param v The vertex.
   * \param cv The curve.
   * \param ind ARR_MIN_END if we refer to cv's minimal end;
   *            ARR_MAX_END if we refer to its maximal end.
   * \param ps_x The boundary condition of the curve end in x.
   * \param ps_y The boundary condition of the curve end in y.
   * \return Whether v represents the left (or right) end of cv.
   */
  bool are_equal (Vertex_const_handle v,
                  const X_monotone_curve_2& cv, Arr_curve_end ind,
                  Arr_parameter_space ps_x, Arr_parameter_space ps_y) const
  {
    return (p_arr->topology_traits()->are_equal (p_arr->_vertex (v),
                                                 cv, ind, ps_x, ps_y));
  }

  /*!
   * Check whether the given halfedge lies on the outer boundary of its
   * incident face.
   * \param he The given halfedge.
   * \return (true) in case he lies on the outer boundary of its incident face;
   *         (false) if he lies on a hole inside this face.
   */
  bool is_on_outer_boundary (Halfedge_const_handle he) const
  {
    const DHalfedge    *p_he = p_arr->_halfedge (he);

    return (! p_he->is_on_inner_ccb());
  }

  /*!
   * Check whether the given halfedge lies on the inner boundary of its
   * incident face.
   * \param he The given halfedge.
   * \return (true) in case he lies on a hole inside its incident face;
   *         (false) if he lies on the outer boundary of this face.
   */
  bool is_on_inner_boundary (Halfedge_const_handle he) const
  {
    const DHalfedge    *p_he = p_arr->_halfedge (he);

    return (p_he->is_on_inner_ccb());
  }

  /*!
   * Create a new vertex and associate it with the given point.
   * \param p The point.
   * \return A handle for the newly created vertex.
   */
  Vertex_handle create_vertex (const Point_2& p)
  {
    DVertex* v = p_arr->_create_vertex (p);
    
    CGAL_assertion (v != NULL);
    return (p_arr->_handle_for (v));
  }
  
  /*!
   * Create a new boundary vertex.
   * \param cv The curve incident to the boundary.
   * \param ind The relevant curve-end.
   * \param ps_x The boundary condition in x.
   * \param by The boundary condition in y.
   * \param notify Should we send a notification to the topology traits
   *               on the creation of the vertex (true by default).
   * \pre Either ps_x or by does not equal ARR_INTERIOR.
   * \return A handle for the newly created vertex.
   */
  Vertex_handle create_boundary_vertex (const X_monotone_curve_2& cv,
                                        Arr_curve_end ind,
                                        Arr_parameter_space ps_x,
                                        Arr_parameter_space ps_y,
                                        bool notify = true)
  {
    DVertex   *v = p_arr->_create_boundary_vertex (cv, ind, ps_x, ps_y);

    CGAL_assertion (v != NULL);

    // Notify the topology traits on the creation of the boundary vertex.
    if (notify)
    {
      p_arr->topology_traits()->notify_on_boundary_vertex_creation(v, cv, ind,
                                                                   ps_x, ps_y);
    }

    return (p_arr->_handle_for (v));
  }

  /*!
   * Locate the arrangement features that will be used for inserting the
   * given curve end, which has a boundary condition, and set a proper vertex
   * there.
   * \param f The face that contains the curve end.
   * \param cv The x-monotone curve.
   * \param ind The curve end.
   * \param ps_x The boundary condition at the x-coordinate.
   * \param ps_y The boundary condition at the y-coordinate.
   * \return A pair of <Vertex_handle, Halfedge_handle>:
   *         The first element is the vertex that corresponds to the curve end.
   *         The second is its predecessor halfedge (if valid).
   */
  std::pair<Vertex_handle, Halfedge_handle>
  place_and_set_curve_end (Face_handle f,
                           const X_monotone_curve_2& cv, Arr_curve_end ind,
                           Arr_parameter_space ps_x, Arr_parameter_space ps_y)
  {
    DHalfedge  *pred;
    DVertex    *v = p_arr->_place_and_set_curve_end (p_arr->_face (f), cv, ind,
                                                     ps_x, ps_y, &pred);

    if (pred == NULL)
      // No predecessor halfedge, return just the vertex:
      return (std::make_pair (p_arr->_handle_for(v), Halfedge_handle()));

    // Return a pair of the vertex and predecessor halfedge:
    return (std::make_pair (p_arr->_handle_for(v), p_arr->_handle_for(pred)));
  }

  /*!
   * Insert an x-monotone curve into the arrangement, where the end vertices
   * are given by the target points of two given halfedges.
   * The two halfedges should be given such that in case a new face is formed,
   * it will be the incident face of the halfedge directed from the first
   * vertex to the second vertex.
   * \param cv the given curve.
   * \param prev1 The reference halfedge for the first vertex.
   * \param prev2 The reference halfedge for the second vertex.
   * \param res The comparsion result between the points associated with the
   *            target vertex of prev and the target vertex of prev2.
   * \param new_face Output - whether a new face has been created.
   * \return A handle for one of the halfedges corresponding to the inserted
   *         curve directed from prev1's target to prev2's target.
   *         In case a new face has been created, it is given as the incident
   *         face of this halfedge.
   */
  Halfedge_handle insert_at_vertices_ex (const X_monotone_curve_2& cv,
                                         Halfedge_handle prev1, 
                                         Halfedge_handle prev2,
                                         Comparison_result res,
                                         bool& new_face)
  {
    DHalfedge*  he = p_arr->_insert_at_vertices (cv,
                                                 p_arr->_halfedge (prev1),
                                                 p_arr->_halfedge (prev2),
                                                 res, new_face);

    CGAL_assertion (he != NULL);
    return (p_arr->_handle_for (he));
  }

  /*!
   * Insert an x-monotone curve into the arrangement, such that one of its
   * endpoints corresponds to a given arrangement vertex, given the exact
   * place for the curve in the circular list around this vertex. The other
   * endpoint corrsponds to a free vertex (a newly created vertex or an
   * isolated vertex).
   * \param cv The given x-monotone curve.
   * \param prev The reference halfedge. We should represent cv as a pair
   *             of edges, one of them should become prev's successor.
   * \param v The free vertex that corresponds to the other endpoint.
   * \param res The comparsion result between the points associated with
   *            the target vertex of prev and the vertex v.
   * \return A handle to one of the halfedges corresponding to the inserted
   *         curve, whose target is the vertex v.
   */
  Halfedge_handle insert_from_vertex_ex (const X_monotone_curve_2& cv,
                                         Halfedge_handle prev,
                                         Vertex_handle v,
                                         Comparison_result res)
  {
    DVertex    *p_v = p_arr->_vertex (v);

    if (p_v->is_isolated())
    {
      // Remove the isolated vertex record, which will not be isolated any
      // more.
      DIso_vertex  *iv = p_v->isolated_vertex();
      DFace        *f = iv->face();

      f->erase_isolated_vertex (iv);
      p_arr->_dcel().delete_isolated_vertex (iv);
    }

    DHalfedge*  he =
      p_arr->_insert_from_vertex (cv, p_arr->_halfedge (prev), p_v, res);

    CGAL_assertion (he != NULL);
    return (p_arr->_handle_for (he));
  }

  /*!
   * Insert an x-monotone curve into the arrangement, such that both its
   * endpoints correspond to free arrangement vertices (newly created vertices
   * or existing isolated vertices), so a new hole is formed in the face
   * that contains the two vertices.
   * \param cv The given x-monotone curve.
   * \param f The face containing the two end vertices.
   * \param v1 The free vertex that corresponds to the left endpoint of cv.
   * \param v2 The free vertex that corresponds to the right endpoint of cv.
   * \param res The comparsion result between the points associated with the
   *            vertices v1 and v2.
   * \return A handle to one of the halfedges corresponding to the inserted
   *         curve, directed from v1 to v2.
   */
  Halfedge_handle insert_in_face_interior_ex (const X_monotone_curve_2& cv,
                                              Face_handle f,
                                              Vertex_handle v1,
                                              Vertex_handle v2,
                                              Comparison_result res)
  {
    DVertex    *p_v1 = p_arr->_vertex (v1);
    DVertex    *p_v2 = p_arr->_vertex (v2);

    if (p_v1->is_isolated())
    {
      // Remove the isolated vertex record, which will not be isolated any
      // more.
      DIso_vertex  *iv1 = p_v1->isolated_vertex();
      DFace        *f1 = iv1->face();

      f1->erase_isolated_vertex (iv1);
      p_arr->_dcel().delete_isolated_vertex (iv1);
    }

    if (p_v2->is_isolated())
    {
      // Remove the isolated vertex record, which will not be isolated any
      // more.
      DIso_vertex  *iv2 = p_v2->isolated_vertex();
      DFace        *f2 = iv2->face();

      f2->erase_isolated_vertex (iv2);
      p_arr->_dcel().delete_isolated_vertex (iv2);
    }

    DHalfedge*  he = p_arr->_insert_in_face_interior (cv,
                                                      p_arr->_face (f),
                                                      p_v1,
                                                      p_v2,
                                                      res);

    CGAL_assertion (he != NULL);
    return (p_arr->_handle_for (he));
  
  }

  /*!
   * Insert the given vertex as an isolated vertex inside the given face.
   * \param f The face that should contain the isolated vertex.
   * \param v The isolated vertex.
   */
  void insert_isolated_vertex (Face_handle f, Vertex_handle v)
  {
    p_arr->_insert_isolated_vertex (p_arr->_face (f), p_arr->_vertex(v));
  }
  
  /*!
   * Relocate all holes and isolated vertices to their proper position,
   * immediately after a face has split due to the insertion of a new halfedge.
   * In case insert_at_vertices_ex() was invoked and indicated that a new face
   * has been created, this function should be called with the halfedge
   * returned by insert_at_vertices_ex().
   * \param new_he The new halfedge that caused the split, such that the new
   *               face lies to its left and the old face to its right.
   */
  void relocate_in_new_face (Halfedge_handle new_he)
  {
    p_arr->_relocate_in_new_face (p_arr->_halfedge (new_he));
    return;
  }

  void relocate_isolated_vertices_in_new_face (Halfedge_handle new_he)
  {
    p_arr->_relocate_isolated_vertices_in_new_face (p_arr->_halfedge(new_he));
    return;
  }

  void relocate_holes_in_new_face (Halfedge_handle new_he)
  {
    p_arr->_relocate_holes_in_new_face (p_arr->_halfedge(new_he));
    return;
  }

  /*!
   * Move an outer CCB from one face to another.
   * \param from_face The source face.
   * \param to_face The destination face.
   * \param ccb A CCB circulator that corresponds to component to move.
   */
  void move_outer_ccb (Face_handle from_face, Face_handle to_face,
                       Ccb_halfedge_circulator ccb)
  {
    p_arr->_move_outer_ccb (p_arr->_face (from_face),
                            p_arr->_face (to_face),
                            p_arr->_halfedge (ccb));
    return;
  }

  /*!
   * Move an inner CCB from one face to another.
   * \param from_face The source face.
   * \param to_face The destination face.
   * \param ccb A CCB circulator that corresponds to component to move.
   */
  void move_inner_ccb (Face_handle from_face, Face_handle to_face,
                       Ccb_halfedge_circulator ccb)
  {
    p_arr->_move_inner_ccb (p_arr->_face (from_face),
                            p_arr->_face (to_face),
                            p_arr->_halfedge (ccb));
    return;
  }
  
  /*!
   * Move an isolated vertex from one face to another.
   * \param from_face The source face.
   * \param to_face The destination face.
   * \param v The isolated vertex to move.
   */
  void move_isolated_vertex (Face_handle from_face, Face_handle to_face,
                             Vertex_handle v)
  {
    p_arr->_move_isolated_vertex (p_arr->_face (from_face),
                                  p_arr->_face (to_face),
                                  p_arr->_vertex (v));
    return;
  }

  /*!
   * Remove an isolated vertex from its face.
   * \param v The isolated vertex to remove.
   */
  void remove_isolated_vertex_ex (Vertex_handle v)
  {
    CGAL_precondition (v->is_isolated());
    DVertex *iso_v = p_arr->_vertex (v);

    p_arr->_remove_isolated_vertex (iso_v);
    return;
  }

  /*!
   * Modify the point associated with a given vertex. The point may be
   * geometrically different than the one currently associated with the vertex.
   * \param v The vertex to modify.
   * \param p The new point to associate with v.
   * \return A handle for the modified vertex (same as v).
   */
  Vertex_handle modify_vertex_ex (Vertex_handle v,
                                  const Point_2& p)
  {
    p_arr->_modify_vertex (p_arr->_vertex (v),
                           p);

    return (v);
  }
        
  /*!
   * Modify the x-monotone curve associated with a given edge. The curve may be
   * geometrically different than the one currently associated with the edge.
   * \param e The edge to modify.
   * \param cv The new x-monotone curve to associate with e.
   * \return A handle for the modified edge (same as e).
   */
  Halfedge_handle modify_edge_ex (Halfedge_handle e,
                                  const X_monotone_curve_2& cv)
  {
    p_arr->_modify_edge (p_arr->_halfedge (e), cv);

    return (e);
  }
          
  /*!
   * Split a given edge into two at a given point, and associate the given
   * x-monotone curves with the split edges.
   * \param e The edge to split (one of the pair of twin halfegdes).
   * \param p The split point.
   * \param cv1 The curve that should be associated with the first split edge,
   *            whose source equals e's source and its target is p.
   * \param cv2 The curve that should be associated with the second split edge,
   *            whose source is p and its target equals e's target.
   * \return A handle for the first split halfedge, whose source equals the
   *         source of e, and whose target is the split point.
   */
  Halfedge_handle split_edge_ex (Halfedge_handle e,
                                 const Point_2& p,
                                 const X_monotone_curve_2& cv1, 
                                 const X_monotone_curve_2& cv2)
  {
    DHalfedge*  he = p_arr->_split_edge (p_arr->_halfedge (e), p, cv1, cv2);

    CGAL_assertion (he != NULL);
    return (p_arr->_handle_for (he));
  }

  /*!
   * Split a given edge into two at the given vertex, and associate the given
   * x-monotone curves with the split edges.
   * \param e The edge to split (one of the pair of twin halfegdes).
   * \param v The split vertex.
   * \param cv1 The curve that should be associated with the first split edge,
   *            whose source equals e's source and its target is v's point.
   * \param cv2 The curve that should be associated with the second split edge,
   *            whose source is v's point and its target equals e's target.
   * \return A handle for the first split halfedge, whose source equals the
   *         source of e, and whose target is the split vertex v.
   */
  Halfedge_handle split_edge_ex (Halfedge_handle e,
                                 Vertex_handle v,
                                 const X_monotone_curve_2& cv1,
                                 const X_monotone_curve_2& cv2)
  {
    DHalfedge*  he = p_arr->_split_edge (p_arr->_halfedge (e),
                                         p_arr->_vertex (v),
                                         cv1, cv2);

    CGAL_assertion (he != NULL);
    return (p_arr->_handle_for (he));
  }

  /*!
   * Split a fictitious edge at the given vertex.
   * \param e The edge to split (one of the pair of twin halfegdes).
   * \param v The split vertex.
   * \return A handle for the first split halfedge, whose source equals the
   *         source of e, and whose target is the split vertex v.
   */
  Halfedge_handle split_fictitious_edge (Halfedge_handle e, Vertex_handle v)
  {
    CGAL_precondition (e->is_fictitious());

    DHalfedge  *he =  
      p_arr->topology_traits()->split_fictitious_edge (p_arr->_halfedge (e),
                                                       p_arr->_vertex (v));

    return (p_arr->_handle_for (he));
  }

  /*!
   * Remove a pair of twin halfedges from the arrangement.
   * \param e A handle for one of the halfedges to be removed.
   * \param remove_source Should the source vertex of e be removed if it
   *                      becomes isolated (true by default).
   * \param remove_target Should the target vertex of e be removed if it
   *                      becomes isolated (true by default).
   * \pre In case the removal causes the creation of a new hole, e should 
   *      point at this hole.
   * \return A handle for the remaining face.
   */
  Face_handle remove_edge_ex (Halfedge_handle e,
                              bool remove_source = true,
                              bool remove_target = true)
  {
    DFace*      f = p_arr->_remove_edge (p_arr->_halfedge (e),
                                         remove_source, remove_target);
    
    CGAL_assertion (f != NULL);
    return (p_arr->_handle_for (f));
  }

  /*!
   * Check if the two given halfedges lie on the same inner component.
   * \param e1 A handle for the first halfedge.
   * \param e2 A handle for the second halfedge.
   * \return Whether e1 and e2 lie on the same inner component.
   */
  bool are_on_same_inner_component (Halfedge_handle e1, Halfedge_handle e2)
  {
     DHalfedge        *he1 = p_arr->_halfedge (e1);
     DHalfedge        *he2 = p_arr->_halfedge (e2);

    const DInner_ccb  *ic1 = (he1->is_on_inner_ccb()) ? he1->inner_ccb() : NULL;

    if (ic1 == NULL)
      return (false);

    const DInner_ccb  *ic2 = (he2->is_on_inner_ccb()) ? he2->inner_ccb() : NULL;

    return (ic1 == ic2);
  }

  /*!
   * Check if the two given halfedges lie on the same outer component.
   * \param e1 A handle for the first halfedge.
   * \param e2 A handle for the second halfedge.
   * \return Whether e1 and e2 lie on the same outer component.
   */
  bool are_on_same_outer_component (Halfedge_handle e1, Halfedge_handle e2)
  {
     DHalfedge        *he1 = p_arr->_halfedge (e1);
     DHalfedge        *he2 = p_arr->_halfedge (e2);

    const DOuter_ccb  *oc1 = (he1->is_on_outer_ccb()) ? he1->outer_ccb() : NULL;

    if (oc1 == NULL)
      return (false);

    const DOuter_ccb  *oc2 = (he2->is_on_outer_ccb()) ? he2->outer_ccb() : NULL;

    return (oc1 == oc2);
  }
  //@}

  /// \name Traversal methods for the BOOST graph traits.
  //@{

  /*! \class
   * An iterator for traversing all arrangement vertices, including vertices
   * at infinity (not including fictitious vertices).
   */
  typedef typename Arrangement_2::_Is_valid_vertex       Is_valid_vertex;
  typedef typename Arrangement_2::_Valid_vertex_iterator Valid_vertex_iterator;

  /*! Get an iterator for the first valid arrangement vertex. */
  Valid_vertex_iterator valid_vertices_begin()
  { 
    return (Valid_vertex_iterator
            (p_arr->topology_traits()->dcel().vertices_begin(),
             p_arr->topology_traits()->dcel().vertices_end(),
             Is_valid_vertex (p_arr->topology_traits()))); 
  }

  /*! Get a past-the-end iterator for the valid arrangement vertices. */
  Valid_vertex_iterator valid_vertices_end()
  { 
    return (Valid_vertex_iterator
            (p_arr->topology_traits()->dcel().vertices_end(),
             p_arr->topology_traits()->dcel().vertices_end(),
             Is_valid_vertex (p_arr->topology_traits()))); 
  }

  /*! Get the number of valid arrangement vertices. */
  Size number_of_valid_vertices () const
  {
    return (p_arr->topology_traits()->number_of_valid_vertices());
  }
  //@}

  /// \name Functions used by the arrangement reader and writer.
  //@{
  typedef typename Arrangement_2::Dcel                Dcel;
  typedef typename Arrangement_2::DVertex_const_iter  Dcel_vertex_iterator;
  typedef typename Arrangement_2::DEdge_const_iter    Dcel_edge_iterator;
  typedef typename Arrangement_2::DFace_const_iter    Dcel_face_iterator;
  typedef typename Arrangement_2::DOuter_ccb_const_iter
                                                      Dcel_outer_ccb_iterator;
  typedef typename Arrangement_2::DInner_ccb_const_iter
                                                      Dcel_inner_ccb_iterator;
  typedef typename Arrangement_2::DIso_vertex_const_iter
                                                      Dcel_iso_vertex_iterator;

  typedef DVertex                               Dcel_vertex;
  typedef DHalfedge                             Dcel_halfedge;
  typedef DFace                                 Dcel_face;
  typedef DOuter_ccb                            Dcel_outer_ccb;
  typedef DInner_ccb                            Dcel_inner_ccb;
  typedef DIso_vertex                           Dcel_isolated_vertex;

  /*!
   * Get the arrangement DCEL.
   */
  const Dcel& dcel () const
  {
    return (p_arr->_dcel());
  }

  /*!
   * Clear the entire arrangment.
   */
  void clear_all ()
  {
    p_arr->clear();
    p_arr->_dcel().delete_all();
    return;
  }

   /*!
   * Set the boundary of a vertex
   * \param p A vertex
   * \param ps_x The boundary condition at x.
   * \param ps_y The boundary condition at y.
   * \return A pointer to the created DCEL vertex.
   */
  Dcel_vertex* set_vertex_boundary (const Vertex_handle v,
                           Arr_parameter_space ps_x, Arr_parameter_space ps_y)
  {
    Dcel_vertex     *v_to_set = p_arr->_vertex (v);

    v_to_set->set_boundary (ps_x, ps_y);

    return (v_to_set);
  }

  /*!
   * Create a new vertex.
   * \param p A pointer to the point (may be NULL in case of a vertex at
   *          infinity).
   * \param ps_x The boundary condition at x.
   * \param ps_y The boundary condition at y.
   * \return A pointer to the created DCEL vertex.
   */
  Dcel_vertex* new_vertex (const Point_2 *p,
                           Arr_parameter_space ps_x, Arr_parameter_space ps_y)
  {
    Dcel_vertex     *new_v = p_arr->_dcel().new_vertex();

    if (p != NULL)
    {
      typename Dcel::Vertex::Point  *p_pt = p_arr->_new_point(*p);
      new_v->set_point (p_pt);
    }
    else
    {
      CGAL_precondition (p_arr->is_open(ps_x, ps_y));
      new_v->set_point (NULL);
    }

    new_v->set_boundary (ps_x, ps_y);
    return (new_v);
  }

  /*!
   * Create a new edge (halfedge pair), associated with the given curve.
   * \param cv A pointer to the x-monotone curve (may be NULL in case of
   *           a fictitious edge).
   * \return A pointer to one of the created DCEL halfedge.
   */
  Dcel_halfedge* new_edge (const X_monotone_curve_2 *cv)
  {
    Dcel_halfedge       *new_he = p_arr->_dcel().new_edge();

    if (cv != NULL)
    {
      typename Dcel::Halfedge::X_monotone_curve  *p_cv = p_arr->_new_curve(*cv);
      new_he->set_curve (p_cv);
    }
    else
    {
      new_he->set_curve (NULL);
    }

    return (new_he);
  }

  /*!
   * Create a new face.
   * \return A pointer to the created DCEL face.
   */
  Dcel_face* new_face ()
  {
    return (p_arr->_dcel().new_face());
  }

  /*!
   * Create a new outer CCB.
   * \return A pointer to the created DCEL outer CCB.
   */
  Dcel_outer_ccb* new_outer_ccb ()
  {
    return (p_arr->_dcel().new_outer_ccb());
  }

  /*!
   * Create a new inner CCB.
   * \return A pointer to the created DCEL inner CCB.
   */
  Dcel_inner_ccb* new_inner_ccb ()
  {
    return (p_arr->_dcel().new_inner_ccb());
  }

  /*!
   * Create a new isolated vertex.
   * \return A pointer to the created DCEL isolated vertex.
   */
  Dcel_isolated_vertex* new_isolated_vertex ()
  {
    return (p_arr->_dcel().new_isolated_vertex());
  }

  /*!
   * Update the topology traits after the DCEL has been updated.
   */
  void dcel_updated()
  {
    p_arr->topology_traits()->dcel_updated();
    return;
  }
  //@}

};

} //namespace CGAL

#endif
