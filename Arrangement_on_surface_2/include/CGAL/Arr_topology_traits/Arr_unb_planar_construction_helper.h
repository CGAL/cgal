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

#ifndef CGAL_ARR_UNB_PLANAR_CONSTRUCTION_HELPER_H
#define CGAL_ARR_UNB_PLANAR_CONSTRUCTION_HELPER_H

/*!
 * Definition of the Arr_unb_planar_construction_helper class-template.
 */

#include <CGAL/Sweep_line_empty_visitor.h>
#include <CGAL/Unique_hash_map.h>

namespace CGAL {

/*! \class Arr_unb_planar_construction_helper
 * A helper class for the construction sweep-line visitor, suitable
 * for an Arrangement_on_surface_2 instantiated with a topology-traits class
 * for unbounded curves in the plane.
 */
template <class Traits_, class Arrangement_, class Event_, class Subcurve_> 
class Arr_unb_planar_construction_helper
{
public:

  typedef Traits_                                      Traits_2;
  typedef Arrangement_                                 Arrangement_2;
  typedef Event_                                       Event;
  typedef Subcurve_                                    Subcurve;

  typedef typename Traits_2::X_monotone_curve_2        X_monotone_curve_2;
  typedef typename Traits_2::Point_2                   Point_2;

  typedef Sweep_line_empty_visitor<Traits_2, Subcurve, Event>
                                                       Base_visitor;

  typedef typename Arrangement_2::Halfedge_handle      Halfedge_handle;
  typedef typename Arrangement_2::Face_handle          Face_handle;
  
  typedef typename Subcurve::Halfedge_indices_list     Indices_list;
  typedef Unique_hash_map<Halfedge_handle, 
                          Indices_list>                Halfedge_indices_map;

protected:

  typedef typename Arrangement_2::Topology_traits      Topology_traits;
  typedef typename Arrangement_2::Vertex_handle        Vertex_handle;

  // Data members:
  Topology_traits*         m_top_traits;  // The topology-traits class.
  Arr_accessor<Arrangement_2>
                           m_arr_access;  // An arrangement accessor.

  Halfedge_handle          m_lh;          // The current left fictitious
                                          // halfedge (on x = -oo).
  Halfedge_handle          m_th;          // The current top fictitious
                                          // halfedge (on y = +oo).
  Halfedge_handle          m_bh;          // The current bottom fictitious
                                          // halfedge (on y = -oo).
  Halfedge_handle          m_rh;          // The current right fictitious
                                          // halfedge (on x = +oo).

  Indices_list      m_subcurves_at_ubf;       // Indices of the curves that
                                              // "see" m_th from below.
  Event*            m_prev_minus_inf_x_event; // The previous event at x = -oo.
  Event*            m_prev_plus_inf_y_event;  // The previous event at y = +oo.

  Halfedge_indices_map*    m_he_ind_map_p;    // A pointer to a map of
                                              // halfedges to indices lists
                                              // (stored in the visitor class).
  
public:
 
  /*! Constructor. */
  Arr_unb_planar_construction_helper(Arrangement_2* arr) :
    m_top_traits(arr->topology_traits()),
    m_arr_access(*arr),
    m_prev_minus_inf_x_event(NULL),
    m_prev_plus_inf_y_event(NULL),
    m_he_ind_map_p(NULL)
  {}

  /*! Destructor. */
  virtual ~Arr_unb_planar_construction_helper(){}

  /// \name Notification functions.
  //@{

  /* A notification issued before the sweep process starts. */
  virtual void before_sweep();

  /*!
   * A notification invoked before the sweep-line starts handling the given
   * event.
   */
  virtual void before_handle_event(Event* event);

  /*! A notification invoked when a new subcurve is created. */
  virtual void add_subcurve(Halfedge_handle /* he */,
                            Subcurve* /* sc */)
  {}

  /*! Collect a subcurve index that does not see any status-line from below. */
  void add_subcurve_in_top_face(unsigned int index)
  {
    m_subcurves_at_ubf.push_back(index);
  }

  /*! Get the indices of the halfedges below the subcurve. */
  Indices_list& halfedge_indices_list()
  {
    return (m_subcurves_at_ubf);
  }


  /*! A notification invoked before the given event it deallocated. */
  void before_deallocate_event(Event* event)
  {
    // The last event at y = +oo may be deallocated if it has no incident
    // right subcurves, so we should not keep a pointer to it.
    if (event == m_prev_plus_inf_y_event)
      m_prev_plus_inf_y_event = NULL;
  }
  //@} 
  
  /*!
   * Set the map that maps each halfedge to the list of subcurve indices
   * that "see" the halfedge from below.
   */
  void set_halfedge_indices_map(Halfedge_indices_map& table)
  {
    m_he_ind_map_p = &table;
  }

  /*!
   * Determine if we should swap the order of predecessor halfedges when
   * calling insert_at_vertices_ex() .
   */
  bool swap_predecessors(Event* event) const
  {
    // If we insert an edge whose right end lies on the top edge of the
    // ficititous bounding rectangle, we have to flip the order of predecessor
    // halfegdes.
    return ((event->parameter_space_in_x() == ARR_INTERIOR) &&
            (event->parameter_space_in_y() == ARR_TOP_BOUNDARY));
  }

  /*! Get the current top face. */
  Face_handle top_face() const
  {
    return (m_th->face());
  }
};

//-----------------------------------------------------------------------------
// Memeber-function definitions:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// A notification issued before the sweep process starts.
//
template <class Tr, class Arr, class Evnt, class Sbcv> 
void Arr_unb_planar_construction_helper<Tr,Arr,Evnt,Sbcv>::before_sweep()
{
  // Obtain the four fictitious vertices that form the "corners" of the
  // fictitious face in the DCEL.
  Vertex_handle v_bl = Vertex_handle(m_top_traits->bottom_left_vertex());
  Vertex_handle v_tl = Vertex_handle(m_top_traits->top_left_vertex());
  CGAL_assertion_code
    (Vertex_handle v_br = Vertex_handle(m_top_traits->bottom_right_vertex());
     Vertex_handle v_tr = Vertex_handle(m_top_traits->top_right_vertex()));
  
  // Get the fictitous halfedges incident to these vertices, representing
  // the left, right, top and bottom edges of the fictitious face.
  //
  //              m_th
  //  v_tl (.)<----------(.) v_tr
  //        |             ^
  //   m_lh |             | m_rh
  //        v             |
  //  v_bl (.)---------->(.) v_br
  //              m_bh
  //
  m_lh = v_tl->incident_halfedges();
  m_lh = (m_lh->source() != v_bl) ? m_lh->next() : m_lh->twin();
  
  m_bh = m_lh->next();
  m_rh = m_bh->next();
  m_th = m_rh->next();

  CGAL_assertion_code
    (Face_handle  fict_face = Face_handle(m_top_traits->fictitious_face()));
  CGAL_assertion(m_lh->direction() == ARR_RIGHT_TO_LEFT);
  CGAL_assertion(m_lh->face() != fict_face);
  CGAL_assertion(m_lh->source() == v_tl && m_lh->target() == v_bl);
  
  CGAL_assertion(m_bh->direction() == ARR_LEFT_TO_RIGHT);
  CGAL_assertion(m_bh->face() != fict_face);
  CGAL_assertion(m_bh->source() == v_bl && m_bh->target() == v_br);
  
  CGAL_assertion(m_rh->direction() == ARR_LEFT_TO_RIGHT);
  CGAL_assertion(m_rh->face() != fict_face);
  CGAL_assertion(m_rh->source() == v_br && m_rh->target() == v_tr);
  
  CGAL_assertion(m_th->direction() == ARR_RIGHT_TO_LEFT);
  CGAL_assertion(m_th->face() != fict_face);
  CGAL_assertion(m_th->source() == v_tr && m_th->target() == v_tl);
}

//-----------------------------------------------------------------------------
// A notification invoked before the sweep-line starts handling the given
// event.
//
template <class Tr, class Arr, class Evnt, class Sbcv> 
void Arr_unb_planar_construction_helper<Tr, Arr, Evnt, Sbcv>::
before_handle_event(Event* event)
{
  if (event->is_closed())
    return;

  // As the event lieas at infinity, it must have only one (right or left)
  // incident curve.
  CGAL_assertion(((event->number_of_left_curves() == 0) &&
                  (event->number_of_right_curves() == 1)) ||
                 ((event->number_of_left_curves() == 1) &&
                  (event->number_of_right_curves() == 0)));
  Arr_curve_end  ind = (event->number_of_left_curves() == 0 &&
                        event->number_of_right_curves() == 1) ?
    ARR_MIN_END : ARR_MAX_END;
  const X_monotone_curve_2&  xc = (ind == ARR_MIN_END) ?
    (*(event->right_curves_begin()))->last_curve() :
    (*(event->left_curves_begin()))->last_curve();

  const Arr_parameter_space  ps_x = event->parameter_space_in_x();
  const Arr_parameter_space  ps_y = event->parameter_space_in_y();

  // Create a vertex at infinity and split the corresponding fictitious edge.
  Vertex_handle v_at_inf =
    m_arr_access.create_boundary_vertex(xc, ind, ps_x, ps_y, false);

  switch (ps_x) {
   case ARR_LEFT_BOUNDARY:
    // The event lies on the left fictitious halfedge.
    m_arr_access.split_fictitious_edge(m_lh, v_at_inf);
    event->set_halfedge_handle(m_lh);

    // Update the incident halfedge of the previous vertex at x = -oo
    // (m_lh used to be incident to it, but now we have split it).
    if (m_prev_minus_inf_x_event != NULL)
      m_prev_minus_inf_x_event->set_halfedge_handle(m_lh->next());
    m_prev_minus_inf_x_event = event;
    return;

   case ARR_RIGHT_BOUNDARY:
    // The event lies on the right fictitious halfedge.
    m_arr_access.split_fictitious_edge(m_rh, v_at_inf);
    event->set_halfedge_handle(m_rh);
    m_rh = m_rh->next();
    return;

  case ARR_INTERIOR:
    break;

  default:
    CGAL_error();
  }

  switch (ps_y) {
   case ARR_BOTTOM_BOUNDARY:
    // The event lies on the bottom fictitious halfedge.
    m_arr_access.split_fictitious_edge(m_bh, v_at_inf);
    event->set_halfedge_handle(m_bh);
    m_bh = m_bh->next();
    return;

   case ARR_TOP_BOUNDARY:
    {
      // The event lies on the top fictitious halfedge.
      m_arr_access.split_fictitious_edge(m_th, v_at_inf);
      event->set_halfedge_handle(m_th);

      // Update the incident halfedge of the previous vertex at y = +oo
      // (m_th used to be incident to it, but now we have split it).
      if(m_prev_plus_inf_y_event != NULL)
        m_prev_plus_inf_y_event->set_halfedge_handle(m_th->next());
      m_prev_plus_inf_y_event = event;

      // Associate all curve indices of subcurves that "see" m_th from
      // below with the left portion of the split halfedge (m_th->next()).
      if (m_he_ind_map_p != NULL)
      {
        Indices_list& list_ref = (*m_he_ind_map_p)[m_th->next()];
        list_ref.clear();
        list_ref.splice(list_ref.end(), m_subcurves_at_ubf);
      }
      else
      {
        m_subcurves_at_ubf.clear();
      }
      CGAL_assertion(m_subcurves_at_ubf.empty());
    }
    return;

  case ARR_INTERIOR:
  default:
    // We are not supposed to reach here at all.
    CGAL_error();
  }
}

} //namespace CGAL

#endif
