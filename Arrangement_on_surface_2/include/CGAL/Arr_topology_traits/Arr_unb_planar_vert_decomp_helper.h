// Copyright (c) 2007,2009,2010,2011 Tel-Aviv University (Israel).
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
//                 Efi Fogel  <efif@post.tau.ac.il>

#ifndef CGAL_ARR_UNB_PLANAR_VERT_DECOMP_HELPER_H
#define CGAL_ARR_UNB_PLANAR_VERT_DECOMP_HELPER_H

/*! \file
 * Definition of the Arr_unb_planar_vert_decomp_helper class-template.
 */

namespace CGAL {

#include <CGAL/Sweep_line_empty_visitor.h>

/*! \class Arr_unb_planar_vert_decomp_helper
 * A helper class for the vertical decomposition sweep-line visitor, suitable
 * for an Arrangement_on_surface_2 instantiated with a topology-traits class
 * for unbounded curves in the plane.
 */
template <class Traits_, class Arrangement_>
class Arr_unb_planar_vert_decomp_helper
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
  const Topology_traits *m_top_traits;  // The topology-traits class.
  Halfedge_const_handle  m_top_fict;    // The current top fictitious halfedge.
  Halfedge_const_handle  m_bottom_fict; // The current bottom fictitious
                                        // halfedge.

public:

  /*!
   * Constructor.
   * \param arr The arrangement.
   */
  Arr_unb_planar_vert_decomp_helper (const Arrangement_2 *arr) :
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

  /*! Get the current top object. */
  CGAL::Object top_object () const
  {
    // Wrap the top fictitious halfedge by a CGAL object.
    return (CGAL::make_object (m_top_fict));
  }

  /*! Get the current bottom object. */
  CGAL::Object bottom_object () const
  {
    // Wrap the bottom fictitious halfedge by a CGAL object.
    return (CGAL::make_object (m_bottom_fict));
  }
};

//-----------------------------------------------------------------------------
// Memeber-function definitions:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// A notification issued before the sweep process starts.
//
template <class Tr, class Arr>
void Arr_unb_planar_vert_decomp_helper<Tr, Arr>::before_sweep ()
{
  // Initialize the fictitious halfedge lying on the top edge of the
  // bounding rectangle. We start from the leftmost halfedge, which is
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

  // Initialize the fictitious halfedge lying on the bottom edge of the
  // bounding rectangle. We start from the leftmost halfedge, which is
  // incident to the bottom-left vertex and whose source is not at -oo.
  Vertex_const_handle  v_bl = 
    Vertex_const_handle (m_top_traits->bottom_left_vertex());

  m_bottom_fict = v_bl->incident_halfedges();
  if (m_bottom_fict->source()->parameter_space_in_x() == ARR_LEFT_BOUNDARY)
      m_bottom_fict = m_bottom_fict->next()->twin();

  CGAL_assertion_code (
    Vertex_const_handle  v_br = 
      Vertex_const_handle (m_top_traits->bottom_right_vertex());
  );
  CGAL_assertion
    ((m_bottom_fict->source() == v_br) ||
     (m_bottom_fict->source()->parameter_space_in_x() == ARR_INTERIOR &&
      m_bottom_fict->source()->parameter_space_in_y() == ARR_BOTTOM_BOUNDARY));

  return;
}

//-----------------------------------------------------------------------------
// A notification invoked after the sweep-line finishes handling the given
// event.
//
template <class Tr, class Arr>
void Arr_unb_planar_vert_decomp_helper<Tr, Arr>::
after_handle_event (Event* event)
{
  // If the event is at infinity and occurs on the top edge of the fictitious
  // face (namely x is finite and y = +oo), we have to update the fictitious
  // halfedges we keep.
  if (event->is_closed())
    return;

  if (event->parameter_space_in_x() != ARR_INTERIOR)
    return;

  Arr_parameter_space     ps_y = event->parameter_space_in_y();

  if (ps_y == ARR_TOP_BOUNDARY)
  {
    // Update the fictitious top halfedge.
    m_top_fict = m_top_fict->twin()->next()->twin();
  }
  else
  {
    CGAL_assertion (ps_y == ARR_BOTTOM_BOUNDARY);

    // Update the fictitious bottom halfedge.
    m_bottom_fict = m_bottom_fict->prev();
  }

  return;
}

} //namespace CGAL

#endif
