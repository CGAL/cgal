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
#ifndef CGAL_ARR_OBSERVER_H
#define CGAL_ARR_OBSERVER_H

/*! \file
 * Definition of the Arr_observer<Arrangement> base class.
 */

CGAL_BEGIN_NAMESPACE

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

protected:

  Arrangement_2  *p_arr;           // The associated arrangement.

private:

  /*! Copy constructor - not supported. */
  Arr_observer (const Self& );

  /*! Assignment operator - not supported. */
  Self& operator= (const Self& );

public:

  /// \name Construction and destruction functions.
  //@{

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
    p_arr->_unregister_observer (this);
  }
  //@}

  /// \name Notification functions on global arrangement operations.
  //@{

  /*! 
   * Notification before the arrangement is assigned with another
   * arrangement.
   * \param arr The arrangement to be copied.
   */
  virtual void before_assign (const Arrangement_2& arr)
  {}

  /*!
   * Notification after the arrangement has been assigned with another
   * arrangement.
   * \param u A handle to the unbounded face.
   */
  virtual void after_assign ()
  {}

  /*! Notification before the arrangement is cleared. */
  virtual void before_clear ()
  {}

  /*!
   * Notification after the arrangement is cleared.
   * \param u A handle to the unbounded face.
   */
  virtual void after_clear (Face_handle /* u */)
  {}

  /*! Notification before a global operation modifies the arrangement. */
  virtual void before_global_change ()
  {}

  /*! Notification after a global operation is completed. */
  virtual void after_global_change ()
  {}
  //@}

  /// \name Notification functions on local changes in the arrangement.
  //@{

  /*!
   * Notification before the creation of a new edge.
   * \param c The x-monotone curve to be associated with the edge.
   */
  virtual void before_create_edge (const X_monotone_curve_2& /* c */)
  {}

  /*!
   * Notification after the creation of a new edge.
   * \param e A handle to one of the twin halfedges that were created.
   */
  virtual void after_create_edge (Halfedge_handle /* e */)
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
   * \param c1 The x-monotone curve to be associated with the first edge.
   * \param c2 The x-monotone curve to be associated with the second edge.
   */
  virtual void before_split_edge (Halfedge_handle /* e */,
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
   * Notification before the creation of a new hole inside a face.
   * \param f A handle to the face containing the hole.
   * \param e The new halfedge that forms the new hole.
   */
  virtual void before_add_hole (Face_handle /* f */,
                                Halfedge_handle /* e */)
  {}

  /*!
   * Notification after a hole was created inside a face.
   * \param h A circulator representing the boundary of the new hole.
   */
  virtual void after_add_hole (Ccb_halfedge_circulator /* h */)
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
   * Notification before a hole is moved from one face to another.
   * \param from_f A handle to the face currently containing the hole.
   * \param to_f A handle to the face that should contain the hole.
   * \param h A circulator representing the boundary of the hole.
   */
  virtual void before_move_hole (Face_handle /* from_f */,
                                 Face_handle /* to_f */,
                                 Ccb_halfedge_circulator /* h */)
  {}

  /*!
   * Notification after a hole is moved from one face to another.
   * \param h A circulator representing the boundary of the hole.
   */
  virtual void after_move_hole (Ccb_halfedge_circulator /* h */)
  {}

  /*!
   * Notification before the removal of an edge.
   * \param e A handle to one of the twin halfedges to be deleted.
   */
  virtual void before_remove_edge (Halfedge_handle /* e */)
  {}

  /*!
   * Notification before the removal of a hole.
   * \param h A circulator representing the boundary of the hole.
   */
  virtual void before_remove_hole (Ccb_halfedge_circulator /* h */)
  {}
  //@}

};

CGAL_END_NAMESPACE

#endif
