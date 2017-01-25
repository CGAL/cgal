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

#ifndef CGAL_ARR_OBSERVER_H
#define CGAL_ARR_OBSERVER_H

#include <CGAL/license/Arrangement_on_surface_2.h>


#include <CGAL/Arr_enums.h>

/*! \file
 * Definition of the Arr_observer<Arrangement> base class.
 */

namespace CGAL {

/*! \class
 * A base class for arrangement observers.
 * The Arrangement parameter corresponds to an arrangement instantiation.
 */
template <class Arrangement_>
class Arr_observer
{
public:

  typedef Arrangement_                                     Arrangement_2;
  typedef Arr_observer<Arrangement_2>                      Self;

  typedef typename Arrangement_2::Point_2                  Point_2;
  typedef typename Arrangement_2::X_monotone_curve_2       X_monotone_curve_2;

  typedef typename Arrangement_2::Vertex_handle            Vertex_handle;
  typedef typename Arrangement_2::Halfedge_handle          Halfedge_handle;
  typedef typename Arrangement_2::Face_handle              Face_handle;
  typedef typename Arrangement_2::Ccb_halfedge_circulator  
                                                      Ccb_halfedge_circulator;

private:

  Arrangement_2  *p_arr;           // The associated arrangement.

  /*! Copy constructor - not supported. */
  Arr_observer (const Self& );

  /*! Assignment operator - not supported. */
  Self& operator= (const Self& );

public:

  /// \name Construction and destruction functions.
  //@{

  /*! Default constructor. */
  Arr_observer () :
    p_arr (NULL)
  {}

  /*! Constructor with an associated arrangement. */
  Arr_observer (Arrangement_2& arr) :
    p_arr (&arr)
  {
    // Register the observer object in the arrangement.
    p_arr->_register_observer (this);
  }

  /*! Destructor. */
  virtual ~Arr_observer ()
  {
    // Unregister the observer object from the arrangement.
    if (p_arr != NULL)
      p_arr->_unregister_observer (this);
  }
  //@}

  /// \name Modifying the associated arrangement.
  //@{

  /*! Get the associated arrangement (non-const version). */
  const Arrangement_2* arrangement () const
  {
    return (p_arr);
  }

  /*! Get the associated arrangement (non-const version). */
  Arrangement_2* arrangement ()
  {
    return (p_arr);
  }

  /*!
   * Attach the observer to an arrangement. 
   * \pre The observer is not already attached to an arrangement.
   */
  void attach (Arrangement_2& arr)
  {
    // Do nothing if the associated arrangement is not changed.
    if (p_arr == &arr)
      return;

    // The observer is not already attached to an arrangement.
    CGAL_precondition (p_arr == NULL);

    if (p_arr != NULL)
      return;

    // Notify the concrete oberver (the sub-class) about the attachment.
    before_attach (arr);

    // Register the observer object in the new arrangement.
    p_arr = &arr;
    p_arr->_register_observer (this);

    // Notify the concrete oberver that the attachment took place.
    after_attach();

    return;
  }

  /*! Detach the observer from the arrangement. */
  void detach ()
  {
    if (p_arr == NULL)
      return;

    // Notify the concrete oberver (the sub-class) about the detachment.
    before_detach ();

    // Unregister the observer object from the current arrangement, and mark
    // that the oberver is not attached to an arrangement.
    p_arr->_unregister_observer (this);
    p_arr = NULL;
   
    // Notify the concrete oberver that the detachment took place.
    after_detach();

    return;
  }
  //@}

  /// \name Notification functions on global arrangement operations.
  //@{

  /*! 
   * Notification before the arrangement is assigned with another
   * arrangement.
   * \param arr The arrangement to be copied.
   */
  virtual void before_assign (const Arrangement_2& /* arr */)
  {}

  /*!
   * Notification after the arrangement has been assigned with another
   * arrangement.
   */
  virtual void after_assign ()
  {}

  /*! Notification before the arrangement is cleared. */
  virtual void before_clear ()
  {}

  /*!
   * Notification after the arrangement is cleared.
   */
  virtual void after_clear ()
  {}

  /*! Notification before a global operation modifies the arrangement. */
  virtual void before_global_change ()
  {}

  /*! Notification after a global operation is completed. */
  virtual void after_global_change ()
  {}
  //@}

  /// \name Notification functions on observer attachment or detachment.
  //@{

  /*! 
   * Notification before the observer is attached to an arrangement.
   * \param arr The arrangement we are about to attach the observer to.
   */
  virtual void before_attach (const Arrangement_2& /* arr */)
  {}

  /*!
   * Notification after the observer has been attached to an arrangement.
   */
  virtual void after_attach ()
  {}

  /*! 
   * Notification before the observer is detached from the arrangement.
   */
  virtual void before_detach ()
  {}

  /*!
   * Notification after the observer has been detached to the arrangement.
   */
  virtual void after_detach ()

  {}
  //@}

  /// \name Notification functions on local changes in the arrangement.
  //@{

  /*!
   * Notification before the creation of a new vertex.
   * \param p The point to be associated with the vertex.
   *          This point cannot lies on the surface boundaries.
   */
  virtual void before_create_vertex (const Point_2& /* p */)
  {}

  /*!
   * Notification after the creation of a new vertex.
   * \param v A handle to the created vertex.
   */
  virtual void after_create_vertex (Vertex_handle /* v */)
  {}

  /*!
   * Notification before the creation of a new boundary vertex.
   * \param cv The curve incident to the surface boundary.
   * \param ind The relevant curve-end.
   * \param ps_x The boundary condition of the vertex in x.
   * \param ps_y The boundary condition of the vertex in y.
   */
  virtual void before_create_boundary_vertex (const X_monotone_curve_2& /*cv*/,
                                              Arr_curve_end /* ind */,
                                              Arr_parameter_space /* ps_x */,
                                              Arr_parameter_space /* ps_y */)
  {}

  /*!
   * Notification after the creation of a new vertex at infinity.
   * \param v A handle to the created vertex.
   */
  virtual void after_create_boundary_vertex (Vertex_handle /* v */)
  {}

  /*!
   * Notification before the creation of a new edge.
   * \param c The x-monotone curve to be associated with the edge.
   * \param v1 A handle to the first end-vertex of the edge.
   * \param v2 A handle to the second end-vertex of the edge.
   */
  virtual void before_create_edge (const X_monotone_curve_2& /* c */,
                                   Vertex_handle /* v1 */,
                                   Vertex_handle /* v2 */)
  {}

  /*!
   * Notification after the creation of a new edge.
   * \param e A handle to one of the twin halfedges that were created.
   */
  virtual void after_create_edge (Halfedge_handle /* e */)
  {}

  /*!
   * Notification before the modification of an existing vertex.
   * \param v A handle to the vertex to be updated.
   * \param p The point to be associated with the vertex.
   */
  virtual void before_modify_vertex (Vertex_handle /* v */,
                                     const Point_2& /* p */)
  {}

  /*!
   * Notification after a vertex was modified.
   * \param v A handle to the updated vertex.
   */
  virtual void after_modify_vertex (Vertex_handle /* v */)
  {}

  /*!
   * Notification before the modification of an existing edge.
   * \param e A handle to one of the twin halfedges to be updated.
   * \param c The x-monotone curve to be associated with the edge.
   */
  virtual void before_modify_edge (Halfedge_handle /* e */,
                                   const X_monotone_curve_2& /* c */)
  {}

  /*!
   * Notification after an edge was modified.
   * \param e A handle to one of the twin halfedges that were updated.
   */
  virtual void after_modify_edge (Halfedge_handle /* e */)
  {}

  /*!
   * Notification before the splitting of an edge into two.
   * \param e A handle to one of the existing halfedges.
   * \param v A vertex representing the split point.
   * \param c1 The x-monotone curve to be associated with the first edge.
   * \param c2 The x-monotone curve to be associated with the second edge.
   */
  virtual void before_split_edge (Halfedge_handle /* e */,
                                  Vertex_handle /* v */,
                                  const X_monotone_curve_2& /* c1 */,
                                  const X_monotone_curve_2& /* c2 */)
  {}

  /*!
   * Notification after an edge was split.
   * \param e1 A handle to one of the twin halfedges forming the first edge.
   * \param e2 A handle to one of the twin halfedges forming the second edge.
   */
  virtual void after_split_edge (Halfedge_handle /* e1 */,
                                 Halfedge_handle /* e2 */)
  {}

  /*!
   * Notification before the splitting of a fictitious edge into two.
   * \param e A handle to one of the existing halfedges.
   * \param v A vertex representing the unbounded split point.
   */
  virtual void before_split_fictitious_edge (Halfedge_handle /* e */,
                                             Vertex_handle /* v */)
  {}

  /*!
   * Notification after a fictitious edge was split.
   * \param e1 A handle to one of the twin halfedges forming the first edge.
   * \param e2 A handle to one of the twin halfedges forming the second edge.
   */
  virtual void after_split_fictitious_edge (Halfedge_handle /* e1 */,
                                            Halfedge_handle /* e2 */)
  {}

  /*!
   * Notification before the splitting of a face into two.
   * \param f A handle to the existing face.
   * \param e The new edge whose insertion causes the face to split.
   */
  virtual void before_split_face (Face_handle /* f */,
                                  Halfedge_handle /* e */)
  {}

  /*!
   * Notification after a face was split.
   * \param f A handle to the face we have just split.
   * \param new_f A handle to the new face that has been created.
   * \param is_hole Whether the new face forms a hole inside f.
   */
  virtual void after_split_face (Face_handle /* f */,
                                 Face_handle /* new_f */,
                                 bool /* is_hole */)
  {}

  /*!
   * Notification before the splitting of an outer CCB into two.
   * \param f A handle to the face that owns the outer CCB.
   * \param h A circulator representing the component boundary.
   * \param e The new edge whose removal causes the outer CCB to split.
   */
  virtual void before_split_outer_ccb (Face_handle /* f */,
                                       Ccb_halfedge_circulator /* h */,
                                       Halfedge_handle /* e */)
  {}

  /*!
   * Notification after an outer CCB was split.
   * \param f A handle to the face that owns the outer CCBs.
   * \param h1 A circulator representing the boundary of the first component.
   * \param h2 A circulator representing the boundary of the second component.
   */
  virtual void after_split_outer_ccb (Face_handle /* f */,
                                      Ccb_halfedge_circulator /* h1 */,
                                      Ccb_halfedge_circulator /* h2 */)
  {}

  /*!
   * Notification before the splitting of an inner CCB into two.
   * \param f A handle to the face containing the inner CCB.
   * \param h A circulator representing the component boundary.
   * \param e The new edge whose removal causes the inner CCB to split.
   */
  virtual void before_split_inner_ccb (Face_handle /* f */,
                                       Ccb_halfedge_circulator /* h */,
                                       Halfedge_handle /* e */)
  {}

  /*!
   * Notification after an inner CCB was split.
   * \param f A handle to the face containing the inner CCBs.
   * \param h1 A circulator representing the boundary of the first component.
   * \param h2 A circulator representing the boundary of the second component.
   */
  virtual void after_split_inner_ccb (Face_handle /* f */,
                                      Ccb_halfedge_circulator /* h1 */,
                                      Ccb_halfedge_circulator /* h2 */)
  {}

  /*!
   * Notification before the creation of a new outer CCB of a face.
   * \param f A handle to the face that owns the outer CCB.
   * \param e A halfedge along the new outer CCB.
   */
  virtual void before_add_outer_ccb (Face_handle /* f */,
                                     Halfedge_handle /* e */)
  {}

  /*!
   * Notification after an outer CCB was added to a face.
   * \param h A circulator representing the boundary of the new outer CCB.
   */
  virtual void after_add_outer_ccb (Ccb_halfedge_circulator /* h */)
  {}

  /*!
   * Notification before the creation of a new inner CCB inside a face.
   * \param f A handle to the face containing the inner CCB.
   * \param e The new halfedge that forms the new inner CCB.
   */
  virtual void before_add_inner_ccb (Face_handle /* f */,
                                     Halfedge_handle /* e */)
  {}

  /*!
   * Notification after an inner CCB was created inside a face.
   * \param h A circulator representing the boundary of the new inner CCB.
   */
  virtual void after_add_inner_ccb (Ccb_halfedge_circulator /* h */)
  {}

  /*!
   * Notification before the creation of a new isolated vertex inside a face.
   * \param f A handle to the face containing the isolated vertex.
   * \param v The isolated vertex.
   */
  virtual void before_add_isolated_vertex (Face_handle /* f */,
                                           Vertex_handle /* v */)
  {}

  /*!
   * Notification after an isolated vertex was created inside a face.
   * \param v The isolated vertex.
   */
  virtual void after_add_isolated_vertex (Vertex_handle /* v */)
  {}

  /*!
   * Notification before the merging of two edges.
   * \param e1 A handle to one of the halfedges forming the first edge.
   * \param e2 A handle to one of the halfedges forming the second edge.
   * \param c The x-monotone curve to be associated with the merged edge.
   */
  virtual void before_merge_edge (Halfedge_handle /* e1 */,
                                  Halfedge_handle /* e2 */,
                                  const X_monotone_curve_2& /* c */)
  {}

  /*!
   * Notification after an edge was merged.
   * \param e A handle to one of the twin halfedges forming the merged edge.
   */
  virtual void after_merge_edge (Halfedge_handle /* e */)
  {}

  /*!
   * Notification before the merging of two fictitious edges.
   * \param e1 A handle to one of the halfedges forming the first edge.
   * \param e2 A handle to one of the halfedges forming the second edge.
   */
  virtual void before_merge_fictitious_edge (Halfedge_handle /* e1 */,
                                             Halfedge_handle /* e2 */)
  {}

  /*!
   * Notification after a fictitious edge was merged.
   * \param e A handle to one of the twin halfedges forming the merged edge.
   */
  virtual void after_merge_fictitious_edge (Halfedge_handle /* e */)
  {}

  /*!
   * Notification before the merging of two faces.
   * \param f1 A handle to the first face.
   * \param f2 A handle to the second face.
   * \param e The edge whose removal causes the faces to merge.
   */
  virtual void before_merge_face (Face_handle /* f1 */,
                                  Face_handle /* f2 */,
                                  Halfedge_handle /* e */)
  {}

  /*!
   * Notification after a face was merged.
   * \param f A handle to the merged face.
   */
  virtual void after_merge_face (Face_handle /* f */)
  {}

  /*!
   * Notification before the merging of two outer CCBs.
   * \param f A handle to the face that owns the outer CCBs.
   * \param h1 A circulator representing the boundary of the first component.
   * \param h2 A circulator representing the boundary of the second component.
   * \param e The edge whose insertion or removal causes the CCBs to merge.
   */
  virtual void before_merge_outer_ccb (Face_handle /* f */,
                                       Ccb_halfedge_circulator /* h1 */,
                                       Ccb_halfedge_circulator /* h2 */,
                                       Halfedge_handle /* e */)
  {}

  /*!
   * Notification after an outer CCB was merged.
   * \param f A handle to the face that owns the outer CCBs.
   * \param h A circulator representing the boundary of the merged component.
   */
  virtual void after_merge_outer_ccb (Face_handle /* f */,
                                      Ccb_halfedge_circulator /* h */)
  {}

  /*!
   * Notification before the merging of two inner CCBs (holes).
   * \param f A handle to the face that contains the inner CCBs.
   * \param h1 A circulator representing the boundary of the first component.
   * \param h2 A circulator representing the boundary of the second component.
   * \param e The edge whose insertion causes the inner CCBs to merge.
   */
  virtual void before_merge_inner_ccb (Face_handle /* f */,
                                       Ccb_halfedge_circulator /* h1 */,
                                       Ccb_halfedge_circulator /* h2 */,
                                       Halfedge_handle /* e */)
  {}

  /*!
   * Notification after an inner CCB was merged.
   * \param f A handle to the face that contains the inner CCBs.
   * \param h A circulator representing the boundary of the merged component.
   */
  virtual void after_merge_inner_ccb (Face_handle /* f */,
                                      Ccb_halfedge_circulator /* h */)
  {}

  /*!
   * Notification before an outer CCB is moved from one face to another.
   * \param from_f A handle to the face that currently owns the outer CCB.
   * \param to_f A handle to the face that should own the outer CCB.
   * \param h A circulator representing the boundary of the component.
   */
  virtual void before_move_outer_ccb (Face_handle /* from_f */,
                                      Face_handle /* to_f */,
                                      Ccb_halfedge_circulator /* h */)
  {}

  /*!
   * Notification after an outer CCB is moved from one face to another.
   * \param h A circulator representing the boundary of the component.
   */
  virtual void after_move_outer_ccb (Ccb_halfedge_circulator /* h */)
  {}


  /*!
   * Notification before an inner CCB is moved from one face to another.
   * \param from_f A handle to the face currently containing the inner CCB.
   * \param to_f A handle to the face that should contain the inner CCB.
   * \param h A circulator representing the boundary of the component.
   */
  virtual void before_move_inner_ccb (Face_handle /* from_f */,
                                      Face_handle /* to_f */,
                                      Ccb_halfedge_circulator /* h */)
  {}

  /*!
   * Notification after an inner CCB is moved from one face to another.
   * \param h A circulator representing the boundary of the component.
   */
  virtual void after_move_inner_ccb (Ccb_halfedge_circulator /* h */)
  {}

  /*!
   * Notification before an isolated vertex is moved from one face to another.
   * \param from_f A handle to the face currently containing the vertex.
   * \param to_f A handle to the face that should contain the vertex.
   * \param v The isolated vertex.
   */
  virtual void before_move_isolated_vertex (Face_handle /* from_f */,
                                            Face_handle /* to_f */,
                                            Vertex_handle /* v */)
  {}

  /*!
   * Notification after an isolated vertex is moved from one face to another.
   * \param v The isolated vertex.
   */
  virtual void after_move_isolated_vertex (Vertex_handle /* v */)
  {}

  /*!
   * Notificaion before the removal of a vertex.
   * \param v A handle to the vertex to be deleted.
   */
  virtual void before_remove_vertex (Vertex_handle /* v */)
  {}

  /*!
   * Notificaion after the removal of a vertex.
   */
  virtual void after_remove_vertex ()
  {}

  /*!
   * Notification before the removal of an edge.
   * \param e A handle to one of the twin halfedges to be deleted.
   */
  virtual void before_remove_edge (Halfedge_handle /* e */)
  {}

  /*!
   * Notificaion after the removal of an edge.
   */
  virtual void after_remove_edge ()
  {}

  /*!
   * Notification before the removal of an outer CCB.
   * \param f The face that owns the outer CCB.
   * \param h A circulator representing the boundary of the component.
   */
  virtual void before_remove_outer_ccb (Face_handle /* f */,
                                        Ccb_halfedge_circulator /* h */)
  {}

  /*!
   * Notificaion after the removal of an outer CCB.
   * \param f The face that used to own the outer CCB.
   */
  virtual void after_remove_outer_ccb (Face_handle /* f */)
  {}

  /*!
   * Notification before the removal of an inner CCB.
   * \param f The face containing the inner CCB.
   * \param h A circulator representing the boundary of the component.
   */
  virtual void before_remove_inner_ccb (Face_handle /* f */,
                                        Ccb_halfedge_circulator /* h */)
  {}

  /*!
   * Notificaion after the removal of an inner CCB.
   * \param f The face that used to contain the inner CCB.
   */
  virtual void after_remove_inner_ccb (Face_handle /* f */)
  {}

  //@}

};

} //namespace CGAL

#endif
