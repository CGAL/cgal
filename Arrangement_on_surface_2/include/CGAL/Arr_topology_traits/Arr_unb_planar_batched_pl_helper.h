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
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>
//                 Ron Wein        <wein@post.tau.ac.il>
//                 Efi Fogel       <efif@post.tau.ac.il>

#ifndef CGAL_ARR_UNB_PLANAR_BATCHED_PL_HELPER_H
#define CGAL_ARR_UNB_PLANAR_BATCHED_PL_HELPER_H

/*!
 * Definition of the Arr_unb_planar_batched_pl_helper class-template.
 */

namespace CGAL {

#include <CGAL/Sweep_line_empty_visitor.h>

/*! \class Arr_unb_planar_batched_pl_helper
 * A helper class for the batched point-location sweep-line visitor, suitable
 * for an Arrangement_on_surface_2 instantiated with a topology-traits class
 * for unbounded curves in the plane.
 */
template <class Traits_, class Arrangement_>
class Arr_unb_planar_batched_pl_helper
{
public:

  typedef Traits_                                      Traits_2;
  typedef Arrangement_                                 Arrangement_2;

  typedef typename Arrangement_2::Face_const_handle    Face_const_handle;

  typedef Sweep_line_empty_visitor<Traits_2>           Base_visitor;
  typedef typename Base_visitor::Event                 Event;
  typedef typename Base_visitor::Subcurve              Subcurve;

protected:

  typedef typename Arrangement_2::Topology_traits       Topology_traits;
  typedef typename Arrangement_2::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Arrangement_2::Vertex_const_handle   Vertex_const_handle;

  // Data members:
  const Topology_traits  *m_top_traits; // The topology-traits class.
  Halfedge_const_handle   m_top_fict;   // The current top fictitious halfedge.

public:

  /*!
   * Constructor.
   * \param arr The arrangement.
   */
  Arr_unb_planar_batched_pl_helper (const Arrangement_2 *arr) :
    m_top_traits (arr->topology_traits())
  {}

  /// \name Notification functions.
  //@{

  /* A notification issued before the sweep process starts. */
  void before_sweep ();

  /*!
   * A notification invoked after the sweep-line finishes handling the given
   * event.
   */
  void after_handle_event (Event* event);
  //@}

  /*! Get the current top face. */
  Face_const_handle top_face () const
  {
    return (m_top_fict->face());
  }
};

//-----------------------------------------------------------------------------
// Memeber-function definitions:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// A notification issued before the sweep process starts.
//
template <class Tr, class Arr>
void Arr_unb_planar_batched_pl_helper<Tr, Arr>::before_sweep ()
{
  // Initialize the fictitious halfedge lying on the top edge of the
  // fictitious face. We start from the leftmost halfedge, which is
  // incident to the top-left vertex and directed from right to left.
  Vertex_const_handle  v_tl = 
    Vertex_const_handle (m_top_traits->top_left_vertex());

  m_top_fict = v_tl->incident_halfedges();
  if (m_top_fict->direction() == ARR_LEFT_TO_RIGHT)
    m_top_fict = m_top_fict->next()->twin();
  
  CGAL_assertion_code (
    Vertex_const_handle  v_tr = 
      Vertex_const_handle (m_top_traits->top_right_vertex());
  );
  CGAL_assertion
    ((m_top_fict->source() == v_tr) ||
     (m_top_fict->source()->parameter_space_in_x() == ARR_INTERIOR &&
      m_top_fict->source()->parameter_space_in_y() == ARR_TOP_BOUNDARY));

  return;
}

//-----------------------------------------------------------------------------
// A notification invoked after the sweep-line finishes handling the given
// event.
//
template <class Tr, class Arr>
void Arr_unb_planar_batched_pl_helper<Tr, Arr>::
after_handle_event (Event* event)
{
  // If the event is at infinity and occurs on the top edge of the fictitious
  // face (namely x is finite and y = +oo), we have to update the fictitious
  // edge we keep.
  if (event->is_closed())
    return;

  if (event->parameter_space_in_x() != ARR_INTERIOR)
    return;
    
  if (event->parameter_space_in_y() == ARR_TOP_BOUNDARY)
    m_top_fict = m_top_fict->twin()->next()->twin();
    
  return;
}

} //namespace CGAL

#endif
