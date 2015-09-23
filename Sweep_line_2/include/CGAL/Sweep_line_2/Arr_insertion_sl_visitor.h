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
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>
//                 Ron Wein <wein@post.tau.ac.il>
//                 Efi Fogel <efif@post.tau.ac.il>

#ifndef CGAL_ARR_INSERTION_SL_VISITOR_H
#define CGAL_ARR_INSERTION_SL_VISITOR_H

/*!
 * Definition of the Arr_insertion_sl_visitor class-template.
 */

#include <CGAL/Sweep_line_2/Arr_basic_insertion_sl_visitor.h>

namespace CGAL {

/*! \class Arr_insertion_sl_visitor
 * A sweep-line visitor for inserting new curves into an existing arrangement
 * embedded on a surface.
 */
template <typename Helper_> 
class Arr_insertion_sl_visitor : 
  public Arr_basic_insertion_sl_visitor<Helper_>
{
public:
  typedef Helper_                                       Helper;
 
  typedef Arr_basic_insertion_sl_visitor<Helper>        Base;

  typedef typename Base::Traits_2                       Traits_2;
  typedef typename Base::Arrangement_2                  Arrangement_2;
  typedef typename Base::Event                          Event;
  typedef typename Base::Subcurve                       Subcurve;
  
  typedef typename Base::Halfedge_handle                Halfedge_handle;
  typedef typename Base::X_monotone_curve_2             X_monotone_curve_2;
  typedef typename Base::Point_2                        Point_2;

private:
  X_monotone_curve_2 sub_cv1;         // Auxiliary variables
  X_monotone_curve_2 sub_cv2;         // (used for splitting curves).

public:
  /*! Constructor. */
  Arr_insertion_sl_visitor (Arrangement_2* arr) : Base(arr) {}

  /// \name Edge-split functions (to be overridden by the child visitor).
  //@{

  /*!
   * Check if the halfedge associated with the given subcurve will be split
   * at the given event.
   * \param sc The subcurve.
   * \param event The event.
   */
  virtual bool is_split_event(Subcurve* sc, Event* event);

  /*!
   * Split the given edge edge.
   * \param he The edge to split.
   * \param sc The associated subcurve.
   * \param The split point.
   * \return A handle to the split edge.
   */
  virtual Halfedge_handle split_edge(Halfedge_handle he, Subcurve* sc,
                                     const Point_2& pt);
  //@}   
};

//-----------------------------------------------------------------------------
// Memeber-function definitions:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Check if the halfedge associated with the given subcurve will be split
// at the given event.
//
template <typename Hlpr> 
bool Arr_insertion_sl_visitor<Hlpr>::is_split_event(Subcurve* sc, Event* event)
{
  if (sc->last_curve().halfedge_handle() == Halfedge_handle(NULL))
    return false;

  if (! sc->originating_subcurve1()) {
    return (reinterpret_cast<Event*>(sc->left_event()) != 
            this->current_event());
  }
  return
    (this->is_split_event
     (reinterpret_cast<Subcurve*>(sc->originating_subcurve1()), event) || 
     this->is_split_event
     (reinterpret_cast<Subcurve*>(sc->originating_subcurve2()), event));
  }

//-----------------------------------------------------------------------------
// Split an edge.
//
template <typename Hlpr>
typename Arr_insertion_sl_visitor<Hlpr>::Halfedge_handle
Arr_insertion_sl_visitor<Hlpr>::split_edge(Halfedge_handle he, Subcurve* sc,
                                           const Point_2& pt)
{
  // Make sure that the halfedge associated with sc is the directed from
  // right to left, since we always "look" above , and the incident face 
  // is on the left of the  halfedge
  CGAL_assertion (he->direction() == ARR_RIGHT_TO_LEFT);

  this->traits()->split_2_object()(he->curve(), pt, sub_cv2, sub_cv1);
  Halfedge_handle new_he =  
    this->m_arr_access.split_edge_ex(he, pt.base(),
                                     sub_cv1.base(), sub_cv2.base());
  Event* last_event_on_sc = reinterpret_cast<Event*>(sc->last_event());
  if (last_event_on_sc->halfedge_handle() == he)
    last_event_on_sc->set_halfedge_handle(new_he->next());

  return new_he;
}

} //namespace CGAL

#endif
