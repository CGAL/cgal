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

#ifndef CGAL_ARR_UNB_PLANAR_INSERTION_HELPER_H
#define CGAL_ARR_UNB_PLANAR_INSERTION_HELPER_H

#include <CGAL/license/Arrangement_on_surface_2.h>


/*!
 * Definition of the Arr_unb_planar_insertion_helper class-template.
 */

#include <CGAL/Sweep_line_2/Arr_construction_sl_visitor.h>
#include <CGAL/Arr_topology_traits/Arr_unb_planar_construction_helper.h>

namespace CGAL {

/*! \class Arr_unb_planar_insertion_helper
 * A helper class for the insertion sweep-line visitors, suitable
 * for an Arrangement_on_surface_2 instantiated with a topology-traits class
 * for unbounded curves in the plane.
 */
template <class Traits_, class Arrangement_, class Event_, class Subcurve_> 
class Arr_unb_planar_insertion_helper :
  public Arr_unb_planar_construction_helper<Traits_, Arrangement_,
                                            Event_, Subcurve_>
{
public:

  typedef Traits_                                      Traits_2;
  typedef Arrangement_                                 Arrangement_2;
  typedef Event_                                       Event;
  typedef Subcurve_                                    Subcurve;

  typedef typename Traits_2::X_monotone_curve_2        X_monotone_curve_2;
  typedef typename Traits_2::Point_2                   Point_2;


  typedef Arr_unb_planar_construction_helper<Traits_2, Arrangement_2, Event,
                                             Subcurve> Base;

  typedef Sweep_line_empty_visitor<Traits_2, Subcurve, Event>
                                                       Base_visitor;

  typedef Arr_unb_planar_insertion_helper<Traits_2, Arrangement_2, Event,
                                          Subcurve>    Self;

  typedef Arr_construction_sl_visitor<Self>            Parent_visitor;

  typedef typename Arrangement_2::Face_handle          Face_handle;

  typedef typename Base::Indices_list                  Indices_list;
  typedef typename Base::Halfedge_indices_map          Halfedge_indices_map;

protected:

  typedef typename Base::Topology_traits               Topology_traits;
  typedef typename Base::Vertex_handle                 Vertex_handle;
  typedef typename Base::Halfedge_handle               Halfedge_handle;
  
public:
 
  /*! Constructor. */
  Arr_unb_planar_insertion_helper (Arrangement_2 *arr) :
    Base (arr)
  {}

  /*! Destructor. */
  virtual ~Arr_unb_planar_insertion_helper(){}

  /// \name Notification functions.
  //@{

  /* A notification issued before the sweep process starts. */
  virtual void before_sweep ();

  /*!
   * A notification invoked before the sweep-line starts handling the given
   * event.
   */
  virtual void before_handle_event (Event* event);
  //@}
};

//-----------------------------------------------------------------------------
// Memeber-function definitions:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// A notification issued before the sweep process starts.
//
template <class Tr, class Arr, class Evnt, class Sbcv> 
void Arr_unb_planar_insertion_helper<Tr,Arr,Evnt,Sbcv>::before_sweep ()
{
  // Obtain the four fictitious vertices that form the "corners" of the
  // fictitious face in the DCEL.
  Vertex_handle   v_bl =
    Vertex_handle (this->m_top_traits->bottom_left_vertex());
  Vertex_handle   v_tl = 
    Vertex_handle (this->m_top_traits->top_left_vertex());
  Vertex_handle   v_br = 
    Vertex_handle (this->m_top_traits->bottom_right_vertex());

  // Get the fictitous halfedges incident to these vertices, and lying on
  // the left, right, top and bottom edges of the fictitious face.
  //
  //            m_th
  //  v_tl (.)<---x      (.) v_tr
  //                      ^
  //        x             | m_rh
  //   m_lh |             x
  //        v              
  //  v_bl (.)----->x    (.) v_br
  //              m_bh
  //
  this->m_lh = v_bl->incident_halfedges();

  if (this->m_lh->source()->parameter_space_in_x() != ARR_LEFT_BOUNDARY)
    this->m_lh = this->m_lh->next()->twin();
   
  this->m_bh = this->m_lh->next();
  
  this->m_th = v_tl->incident_halfedges();
  if (this->m_th->source()->parameter_space_in_x() == ARR_LEFT_BOUNDARY)
    this->m_th = this->m_th->next()->twin();

  this->m_rh = v_br->incident_halfedges();
  if (this->m_rh->source()->parameter_space_in_x() == ARR_RIGHT_BOUNDARY)
    this->m_rh = this->m_rh->twin();
  else
    this->m_rh = this->m_rh->next();

  CGAL_assertion_code (
    Face_handle fict_face = Face_handle (this->m_top_traits->fictitious_face());
  );
  CGAL_assertion (this->m_lh->direction() == ARR_RIGHT_TO_LEFT);
  CGAL_assertion (this->m_lh->face() != fict_face);
  CGAL_assertion (this->m_lh->target() == v_bl);

  CGAL_assertion (this->m_bh->direction() == ARR_LEFT_TO_RIGHT);
  CGAL_assertion (this->m_bh->face() != fict_face);
  CGAL_assertion (this->m_bh->source() == v_bl);

  CGAL_assertion (this->m_rh->direction() == ARR_LEFT_TO_RIGHT);
  CGAL_assertion (this->m_rh->face() != fict_face);
  CGAL_assertion (this->m_rh->source() == v_br);

  CGAL_assertion (this->m_th->direction() == ARR_RIGHT_TO_LEFT);
  CGAL_assertion (this->m_th->face() != fict_face);
  CGAL_assertion (this->m_th->target() == v_tl);
}

//-----------------------------------------------------------------------------
// A notification invoked before the sweep-line starts handling the given
// event.
//
template <class Tr, class Arr, class Evnt, class Sbcv> 
void Arr_unb_planar_insertion_helper<Tr,Arr,Evnt,Sbcv>::
before_handle_event (Event* event)
{
  if (event->is_closed())
    return;

  // In case the event lies at inifinity, check whether its incident curve
  // is already in the arrangement.
  if (event->curve().halfedge_handle() == Halfedge_handle())
  {
    // The curve is not in the arrangement, use the base construction helper
    // to handle the event:
    Base::before_handle_event (event);
    return;
  }

  // The curve is already in the arrangement, but has an infinite end,
  // so we have to update the fictitious halfedges.
  const Arr_parameter_space ps_x = event->parameter_space_in_x();
      
  if (ps_x == ARR_LEFT_BOUNDARY)
  {
    // The event lies on the left fictitious halfedge.
    this->m_lh = this->m_lh->twin()->next()->twin();
    this->m_prev_minus_inf_x_event = NULL;
  }
  else if (ps_x == ARR_RIGHT_BOUNDARY)
  {
    // The event lies on the right fictitious halfedge.
    this->m_rh = this->m_rh->twin()->prev()->twin();
  }
  else
  {
    const Arr_parameter_space ps_y = event->parameter_space_in_y();
    
    if (ps_y == ARR_BOTTOM_BOUNDARY)
    {
      // The event lies on the bottom fictitious halfedge.
      this->m_bh = this->m_bh->twin()->prev()->twin();
    }
    else
    {
      // The event lies on the top fictitious halfedge.
      CGAL_assertion (ps_y == ARR_TOP_BOUNDARY);
      this->m_th = this->m_th->twin()->next()->twin();
      this->m_prev_plus_inf_y_event = NULL;
    }
  }
 
  return;
}


} //namespace CGAL

#endif
